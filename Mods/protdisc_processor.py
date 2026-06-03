# Created by @thibaultcolombani with the help of Claude Code (Anthropic)
"""
protdisc_processor.py - ProteomeDiscoverer data processor for MassSleuth.

Takes the dict produced by load_folder() for engine='protdisc' and normalises
the PSMs DataFrame (with optional InputFiles merge) to the canonical column
schema used by all other processors and HTML tabs.

Works with any PD search node: Comet, Sequest HT, MSFragger-PD, etc.
The output file format (*PSMs*.txt, *InputFiles*.txt, ...) is the same
regardless of the underlying search algorithm.

PD file structure
-----------------
*PSMs*.txt       Required. PSM-level results, tab-separated.
*InputFiles*.txt Optional. Maps integer File ID → raw file path.

Key PSMs.txt columns used
--------------------------
Sequence                   bare peptide sequence
Annotated Sequence         modified sequence (PD bracket notation)
Modifications              semicolon-separated named mods, e.g.:
                             "mTRAQ D0 [N-Term]; mTRAQ D4 [K5]"
                             "N-Term(mTRAQ D0); K5(mTRAQ D4)"
Charge                     precursor charge
RT [min]                   retention time in minutes (also: "Retention Time [min]")
Intensity                  precursor MS1 intensity
Percolator q-Value         PSM-level FDR (also: "q-Value")
Percolator PEP             posterior error probability (also: "PEP")
Master Protein Accessions  semicolon-separated protein accessions
File ID                    integer join key to InputFiles
Confidence                 "High" / "Medium" / "Low" — used when no q-value column

InputFiles.txt columns used
----------------------------
File ID / Study File ID    integer join key
File Name                  raw file path → stem becomes 'Raw file'

Canonical columns produced
---------------------------
Raw file           ← File Name stem (from InputFiles) or Spectrum File stem
Sequence           ← Sequence
Modified sequence  ← Annotated Sequence
Retention time     ← RT [min] (already in minutes — no conversion)
Charge             ← Charge (Int32)
Intensity          ← Intensity (Float64)
Q.Value            ← Percolator q-Value (or q-Value)
PEP                ← Percolator PEP (or PEP)
Protein.Group      ← Master Protein Accessions

Derived columns
---------------
seqcharge          "<Sequence>_<Charge>"
Precursor.Id       "<Modified sequence>_<Charge>"  (or seqcharge if absent)
C_term_aa          last residue of Sequence
N_term_aa          first residue of Sequence
fully_labeled      1 if N-term tagged AND all K in sequence are tagged
partly_labeled     1 if N-term tagged AND at least one K is untagged

Tag name detection
------------------
Auto-detected from Modifications: modification names that appear on BOTH
N-terminal AND K side-chain positions (amine-reactive = labels both sites).
Multiple distinct N-terminal variant names of the same base tag (D0/D4/D8)
→ multiplexed; Channel = N-terminal modification name in plex_report.

Quality filter (priority order — first match wins)
---------------------------------------------------
1. Percolator q-Value ≤ 0.01
2. q-Value ≤ 0.01
3. Confidence == 'High'
4. Percolator PEP ≤ 0.01
"""

from __future__ import annotations

import re
from collections import Counter
from typing import Dict, List, Optional, Tuple

import polars as pl


# ---------------------------------------------------------------------------
# Module-level helpers
# ---------------------------------------------------------------------------

def _find_col(df: pl.DataFrame, candidates: List[str]) -> Optional[str]:
    """Return the first candidate column name present in df, or None."""
    for c in candidates:
        if c in df.columns:
            return c
    return None


def _merge_input_files(psms: pl.DataFrame,
                       input_files: Optional[pl.DataFrame]) -> pl.DataFrame:
    """
    Derive 'Raw file' by joining PSMs (File ID) with InputFiles (File Name stem).

    Falls back to extracting stem from 'Spectrum File' column; last resort:
    cast 'File ID' to string.
    """
    if input_files is not None and 'File ID' in psms.columns:
        id_col = _find_col(input_files, ['File ID', 'Study File ID'])

        if id_col is not None and 'File Name' in input_files.columns:
            file_map = (
                input_files
                .select([
                    pl.col(id_col).cast(pl.String).alias('__fid'),
                    pl.col('File Name').alias('__fname'),
                ])
                .unique(subset=['__fid'])
                .with_columns(
                    pl.col('__fname')
                      .str.replace_all(r'\', '/', literal=False)
                      .str.extract(r'([^/]+)$', group_index=1)
                      .str.replace(r'\.[^.]+$', '', literal=False)
                      .alias('__raw')
                )
            )
            psms = (
                psms.with_columns(pl.col('File ID').cast(pl.String).alias('__fid'))
                    .join(file_map.select(['__fid', '__raw']),
                          on='__fid', how='left')
                    .rename({'__raw': 'Raw file'})
                    .drop('__fid')
            )
            return psms

    # Fallback: extract from Spectrum File column
    if 'Spectrum File' in psms.columns:
        return psms.with_columns(
            pl.col('Spectrum File')
              .str.replace_all(r'\', '/', literal=False)
              .str.extract(r'([^/]+)$', group_index=1)
              .str.replace(r'\.[^.]+$', '', literal=False)
              .alias('Raw file')
        )

    if 'File ID' in psms.columns:
        # When multiple PSM source files are concatenated their File IDs overlap.
        # Prefix with the PSM filename stem so runs from different files stay separate.
        if '_source_file' in psms.columns and psms['_source_file'].n_unique() > 1:
            stem_expr = (
                pl.col('_source_file')
                  .str.replace_all(r'\', '/', literal=False)
                  .str.extract(r'([^/]+)\.[^./]+$', group_index=1)
                  .str.replace(r'_PSMs.*$', '', literal=False)
            )
            return psms.with_columns(
                (stem_expr + '_' + pl.col('File ID').cast(pl.String)).alias('Raw file')
            )
        return psms.with_columns(
            pl.col('File ID').cast(pl.String).alias('Raw file')
        )

    return psms


def _detect_tag_names(mods_series: pl.Series) -> List[str]:
    """
    Auto-detect amine-reactive tag names from PD Modifications column.

    Supports two PD notation formats:
      "mTRAQ D0 [N-Term]; mTRAQ D4 [K5]"   — name then bracketed position
      "N-Term(mTRAQ D0); K5(mTRAQ D4)"      — position then name in parens

    Returns sorted list of N-terminal modification names whose base name
    (first word) also appears on K side-chains (= amine-reactive tag).
    Multiple distinct N-terminal variants of the same base (D0/D4/D8)
    → all returned so plex_report can be built with each as a Channel.
    """
    f1_nt_re = re.compile(r'^(.+?)\s*\[N-[Tt]erm\]')
    f1_k_re  = re.compile(r'^(.+?)\s*\[K\d+\]')
    f2_nt_re = re.compile(r'N-[Tt]erm\(([^)]+)\)')
    f2_k_re  = re.compile(r'K\d*\(([^)]+)\)')   # PD format: K8(TMTpro) or K(TMTpro)

    nterm_names: Counter = Counter()
    k_names:     Counter = Counter()

    for mods in mods_series:
        if not mods:
            continue
        for entry in re.split(r';\s*', str(mods)):
            entry = entry.strip()
            if not entry:
                continue
            m = f1_nt_re.match(entry)
            if m:
                nterm_names[m.group(1).strip()] += 1
                continue
            m = f1_k_re.match(entry)
            if m:
                k_names[m.group(1).strip()] += 1
                continue
            m = f2_nt_re.search(entry)
            if m:
                nterm_names[m.group(1).strip()] += 1
                continue
            m = f2_k_re.search(entry)
            if m:
                k_names[m.group(1).strip()] += 1

    # Base name = first word; must be > 2 chars to exclude single-letter designators
    k_bases = {n.split()[0] for n in k_names if n.split() and len(n.split()[0]) > 2}

    return sorted(
        nname for nname in nterm_names
        if nname.split() and nname.split()[0] in k_bases
    )


def _add_derived(df: pl.DataFrame,
                 tag_names: Optional[List[str]] = None) -> pl.DataFrame:
    """
    Add seqcharge, Precursor.Id, terminal AAs, and labeling columns.

    Labeling uses name-agnostic patterns from the 'Modifications' column,
    matching the approach of the Tag-Labeling reference script:
      - N-term labeled: ANY N-terminal modification present
        (N-Term(name) or name [N-Term])
      - K-position labeled: ANY K-position modification present
        (K8(name), K(name), or name [K8])
    Works for any amine-reactive tag (TMTpro, mTRAQ, PSMtag, iTRAQ, …).

    tag_names acts as a gate: labeling flags are only added when an
    amine-reactive tag was pre-detected on both N-term and K (non-empty
    tag_names from _detect_tag_names). This prevents the Labeling tab
    from appearing in label-free experiments where N-terminal acetylation
    would create false positives.

      fully_labeled  = N-term has any mod AND every K in sequence has a K-position mod
      partly_labeled = N-term has any mod AND at least one K has no K-position mod
    """
    exprs = []

    if 'Sequence' in df.columns and 'Charge' in df.columns:
        bare = pl.col('Sequence') + '_' + pl.col('Charge').cast(pl.String)
        exprs.append(bare.alias('seqcharge'))
        if 'Modified sequence' in df.columns:
            exprs.append(
                (pl.col('Modified sequence') + '_' + pl.col('Charge').cast(pl.String))
                .alias('Precursor.Id')
            )
        else:
            exprs.append(bare.alias('Precursor.Id'))

    if 'Sequence' in df.columns:
        exprs += [
            pl.col('Sequence').str.slice(-1).alias('C_term_aa'),
            pl.col('Sequence').str.slice(0, 1).alias('N_term_aa'),
        ]

    if exprs:
        df = df.with_columns(exprs)

    # Only add labeling flags when an amine-reactive tag was confirmed present
    if 'Modifications' not in df.columns or 'Sequence' not in df.columns or not tag_names:
        return df

    mods      = pl.col('Modifications').fill_null('')
    n_k_total = pl.col('Sequence').str.count_matches('K')

    # Name-agnostic N-term detection (both PD notation formats):
    #   Format 2: N-Term(name)   → N-[Tt]erm\s*\(
    #   Format 1: name [N-Term]  → \[N-[Tt]erm\]
    nterm_any = mods.str.contains(r'N-[Tt]erm\s*\(|\[N-[Tt]erm\]', literal=False)

    # Name-agnostic K-position count (both PD notation formats):
    #   Format 2: K8(name) or K(name)  → K\d*\s*\(
    #   Format 1: name [K8]            → \[K\d+\]
    n_k_tagged = (
        mods.str.count_matches(r'K\d*\s*\(', literal=False) +
        mods.str.count_matches(r'\[K\d+\]', literal=False)
    )

    df = df.with_columns([
        (nterm_any & (n_k_tagged >= n_k_total)).cast(pl.Int8).alias('fully_labeled'),
        (nterm_any & (n_k_tagged < n_k_total)).cast(pl.Int8).alias('partly_labeled'),
    ])
    return df


def _build_plex_report(df: pl.DataFrame,
                       tag_names: List[str]) -> Optional[pl.DataFrame]:
    """
    Build plex_report from multiplexed data.

    Extracts the N-terminal modification name from 'Modifications' as 'Channel'.
    Supports both PD notation formats.
    Filters to rows where Channel is in the detected tag_names list.
    """
    if 'Modifications' not in df.columns or not tag_names:
        return None

    plex = df.with_columns(
        pl.when(
            pl.col('Modifications').str.contains(r'\[N-[Tt]erm\]', literal=False)
        ).then(
            pl.col('Modifications')
              .str.extract(r'^([^;]+?)\s*\[N-[Tt]erm\]', group_index=1)
              .str.strip_chars()
        ).otherwise(
            pl.col('Modifications')
              .str.extract(r'N-[Tt]erm\(([^)]+)\)', group_index=1)
              .str.strip_chars()
        ).alias('Channel')
    )
    plex = plex.filter(pl.col('Channel').is_in(tag_names))

    if plex.is_empty():
        return None

    n_ch = plex['Channel'].n_unique()
    print(f"  [protdisc] plex_report: {len(plex):,} rows, "
          f"{n_ch} channels: {sorted(tag_names)}")
    return plex


_ABD_PAT = re.compile(r'^Abundance[:\s]+(.+)$')


def _build_plex_report_tmt(df: pl.DataFrame) -> Optional[pl.DataFrame]:
    """
    Build plex_report from PD TMT-style 'Abundance: NNN' wide columns.

    Unpivots to long format (one row per PSM × channel).
    Channel name = column suffix after 'Abundance: ' (e.g. '126', '127N', '134CD').
    Active channels: ≥ 5% of PSMs have non-zero signal (minimum 10 rows).
    """
    abd_cols = {m.group(1): c for c in df.columns if (m := _ABD_PAT.match(c))}
    if not abd_cols:
        return None

    n = len(df)
    threshold = max(10, int(n * 0.05))
    active = {
        ch: col for ch, col in abd_cols.items()
        if df.filter(pl.col(col).cast(pl.Float64, strict=False) > 0).height >= threshold
    }
    if not active:
        return None

    print(f"  [protdisc] TMT Abundance: {len(active)} active channel(s) detected")

    carry = [c for c in ('Raw file', 'Precursor.Id', 'seqcharge',
                         'Sequence', 'Modified sequence', 'Charge',
                         'Protein.Group', 'PEP', 'Q.Value', 'Retention time',
                         'Precursor.Mz', 'fully_labeled', 'partly_labeled')
             if c in df.columns]

    def _ch_key(ch: str):
        m = re.match(r'^(\d+)', ch)
        return (int(m.group(1)) if m else 0, ch)

    frames = []
    for ch, col in sorted(active.items(), key=lambda x: _ch_key(x[0])):
        ch_df = (
            df.select(carry + [col])
              .with_columns(
                  pl.col(col).cast(pl.Float64, strict=False).alias('Intensity')
              )
              .drop(col)
              # For TMT, reporter ions ARE the MS2 signal: expose as Precursor.Quantity
              # so plex_diagnostic has_ms2 = True and MS2 dist/CV sections render.
              .with_columns(pl.col('Intensity').alias('Precursor.Quantity'))
              .with_columns(pl.lit(ch).alias('Channel'))
              .filter(pl.col('Intensity') > 0)
        )
        frames.append(ch_df)

    if not frames:
        return None

    result = pl.concat(frames, how='diagonal')
    print(f"  [protdisc] plex_report: {len(result):,} rows  "
          f"({len(active)} channels × {n:,} PSMs)")
    return result


def _filter_quality(df: pl.DataFrame) -> Tuple[pl.DataFrame, List[str]]:
    """
    Apply quality filter using the best available column (priority order):
      1. Percolator q-Value ≤ 0.01
      2. q-Value ≤ 0.01
      3. Confidence == 'High'
      4. Percolator PEP ≤ 0.01
    """
    filters: List[str] = []

    # By this point normalize() has renamed q-Value → Q.Value and PEP stays PEP
    if 'Q.Value' in df.columns:
        before = len(df)
        df = df.filter(pl.col('Q.Value').cast(pl.Float64, strict=False) <= 0.01)
        print(f"  [protdisc] q-Value ≤ 0.01: {before:,} → {len(df):,}")
        filters.append('q-Value ≤ 0.01')
    elif 'Confidence' in df.columns:
        before = len(df)
        df = df.filter(pl.col('Confidence') == 'High')
        print(f"  [protdisc] Confidence == High: {before:,} → {len(df):,}")
        filters.append('Confidence = High')
    elif 'PEP' in df.columns:
        before = len(df)
        df = df.filter(pl.col('PEP').cast(pl.Float64, strict=False) <= 0.01)
        print(f"  [protdisc] PEP ≤ 0.01: {before:,} → {len(df):,}")
        filters.append('PEP ≤ 0.01')

    return df, filters


# ---------------------------------------------------------------------------
# Public processor class
# ---------------------------------------------------------------------------

class ProtdiscProcessor:
    """
    Normalises ProteomeDiscoverer output loaded by data_loader.load_folder().

    Parameters
    ----------
    data : dict[str, pl.DataFrame]
        Output of load_folder() when engine == 'protdisc'.
        Required key: 'psms'         (*PSMs*.txt content).
        Optional key: 'input_files'  (*InputFiles*.txt — maps File ID → filename).
    """

    def __init__(self, data: Dict[str, pl.DataFrame]):
        if 'psms' not in data or data['psms'] is None:
            raise ValueError(
                "ProtdiscProcessor requires a 'psms' DataFrame in data dict."
            )

        self.report:        pl.DataFrame           = data['psms']
        self._input_files:  Optional[pl.DataFrame] = data.get('input_files')
        self.features:      Optional[pl.DataFrame] = None
        self.fill_times:    Optional[pl.DataFrame] = None
        self.tic:           Optional[pl.DataFrame] = None
        self.sn:            Optional[pl.DataFrame] = None
        self.ms1_extracted: Optional[pl.DataFrame] = None
        self.plex_report:   Optional[pl.DataFrame] = None
        self.raw_report:    Optional[pl.DataFrame] = None
        self.enzymes:       Optional[List[str]]    = None
        self.filters:       List[str]              = []

    # ------------------------------------------------------------------
    # Main pipeline
    # ------------------------------------------------------------------

    def normalize(self) -> 'ProtdiscProcessor':
        """
        Run the full normalisation pipeline:
          1.  Merge with InputFiles → 'Raw file' from File Name stem
          2.  Map Retention time column (PD already outputs minutes)
          3.  Rename Q.Value / PEP / Protein.Group / Modified sequence / Precursor.Mz
          4.  Cast Charge (Int32) and Intensity (Float64)
          5.  Detect amine-reactive tag names from Modifications (pre-filter)
          6.  Add derived columns (seqcharge, Precursor.Id, terminal AAs, labels)
          7.  Save raw_report snapshot
          8.  Apply quality filter
          9.  Deduplicate by (Raw file, Precursor.Id) keeping best intensity
         10.  Build plex_report when > 1 distinct tag variant detected
        Returns self for chaining.
        """
        df = self.report

        # 1. Run names from InputFiles
        df = _merge_input_files(df, self._input_files)

        # 2. Retention time — PD writes minutes; just cast + rename
        rt_col = _find_col(df, ['RT [min]', 'Retention Time [min]', 'RTmin',
                                 'RT in min', 'Retention time'])
        if rt_col and rt_col != 'Retention time':
            df = (df.with_columns(
                    pl.col(rt_col).cast(pl.Float64, strict=False).alias('Retention time')
                  ).drop(rt_col))

        # 3. Column renames
        rmap: Dict[str, str] = {}

        modseq_col = _find_col(df, ['Annotated Sequence', 'AnnotatedSequence',
                                     'Annotated sequence'])
        if modseq_col and modseq_col != 'Modified sequence':
            rmap[modseq_col] = 'Modified sequence'

        qval_col = _find_col(df, ['Percolator q-Value', 'q-Value',
                                   'Percolator q Value'])
        if qval_col and qval_col != 'Q.Value':
            rmap[qval_col] = 'Q.Value'

        pep_col = _find_col(df, ['Percolator PEP', 'PEP'])
        if pep_col and pep_col != 'PEP':
            rmap[pep_col] = 'PEP'

        prot_col = _find_col(df, ['Master Protein Accessions',
                                   'MasterProteinAccessions',
                                   'Protein Accessions'])
        if prot_col and prot_col != 'Protein.Group':
            rmap[prot_col] = 'Protein.Group'

        mz_col = _find_col(df, ['m/z [Da]', 'mz in Da', 'Theo. MH+ [Da]', 'm/z'])
        if mz_col and mz_col != 'Precursor.Mz':
            rmap[mz_col] = 'Precursor.Mz'

        if rmap:
            df = df.rename(rmap)

        # Strip PD flanking residue notation: "[R].PEPTIDE.[K]" → "PEPTIDE"
        # PD includes the preceding/following AA in Annotated Sequence as context.
        # The same peptide can have different flanking residues across fractions
        # → without stripping, Precursor.Id differs between runs → cross-run
        # intersection collapses to nearly zero.
        if 'Modified sequence' in df.columns:
            df = df.with_columns(
                pl.col('Modified sequence')
                  .str.replace(r'^\[[A-Z\-]+\]\.', '', literal=False)
                  .str.replace(r'\.\[[A-Z\-]+\]$', '', literal=False)
                  .alias('Modified sequence')
            )

        # 4. Type casts
        if 'Charge' in df.columns:
            df = df.with_columns(pl.col('Charge').cast(pl.Int32, strict=False))
        if 'Intensity' in df.columns:
            df = df.with_columns(pl.col('Intensity').cast(pl.Float64, strict=False))

        # 5. Tag name detection (pre-filter for maximum PSM coverage)
        tag_names: List[str] = []
        if 'Modifications' in df.columns:
            tag_names = _detect_tag_names(df['Modifications'])
            if tag_names:
                print(f"  [protdisc] tag names detected: {tag_names}")

        # 6. Derived columns
        df = _add_derived(df, tag_names)

        # 7. Pre-filter snapshot
        self.raw_report = df

        # 8. Quality filter
        df, self.filters = _filter_quality(df)

        # 9. Deduplicate by (Raw file, Precursor.Id)
        if 'Precursor.Id' in df.columns and 'Raw file' in df.columns:
            before   = len(df)
            sort_col = 'Intensity' if 'Intensity' in df.columns else None
            if sort_col:
                df = (df.sort(sort_col, descending=True, nulls_last=True)
                        .unique(subset=['Raw file', 'Precursor.Id'], keep='first'))
            else:
                df = df.unique(subset=['Raw file', 'Precursor.Id'], keep='first')
            removed = before - len(df)
            if removed:
                print(f"  [protdisc] dedup: {before:,} → {len(df):,} "
                      f"(dropped {removed:,} duplicate PSMs)")

        self.report = df
        print(f"  [protdisc] normalized report: {len(self.report):,} rows "
              f"× {self.report.width} columns")

        # 10. Plex report
        if len(tag_names) > 1:
            # PSMtag / mTRAQ: distinct N-terminal isotope variants are channels
            self.plex_report = _build_plex_report(self.report, tag_names)
        else:
            # TMT / iTRAQ: all channels share one mod name; channels live in Abundance: NNN cols
            abd_cols = [c for c in self.report.columns if _ABD_PAT.match(c)]
            if abd_cols:
                # Add Precursor.Quantity = total reporter ion signal (sum across channels)
                # so MS2-based tabs (IDs missed cleavages, protein groups) can use it
                # as a quantification filter when MS1 Intensity is absent.
                self.report = self.report.with_columns(
                    pl.sum_horizontal(
                        [pl.col(c).cast(pl.Float64, strict=False).fill_null(0.0)
                         for c in abd_cols]
                    ).alias('Precursor.Quantity')
                )
            plex = _build_plex_report_tmt(self.report)
            if plex is not None:
                self.plex_report = plex
                if abd_cols:
                    self.report = self.report.drop(abd_cols)

        return self

    def filter_enzyme(self, enzymes: List[str]) -> 'ProtdiscProcessor':
        """Filter to fully specific peptides for the given enzyme(s)."""
        from Mods.enzyme import filter_peptides
        self.report  = filter_peptides(self.report, enzymes)
        self.enzymes = list(enzymes)
        return self

    # ------------------------------------------------------------------
    # Convenience accessors
    # ------------------------------------------------------------------

    def runs(self) -> List[str]:
        """Sorted list of unique run names."""
        return sorted(self.report['Raw file'].unique().to_list())

    def summary_stats(self) -> Dict:
        df = self.report
        return {
            'rows':       len(df),
            'runs':       df['Raw file'].n_unique()      if 'Raw file'      in df.columns else None,
            'precursors': df['Precursor.Id'].n_unique()  if 'Precursor.Id'  in df.columns else None,
            'seqcharges': df['seqcharge'].n_unique()     if 'seqcharge'     in df.columns else None,
            'proteins':   df['Protein.Group'].n_unique() if 'Protein.Group' in df.columns else None,
        }
