# Created by @thibaultcolombani with the help of Claude Code (Anthropic)
"""
fragpipe_processor.py - FragPipe data processor for MassSleuth.

Takes the dict produced by load_folder() for engine='fragpipe' and normalises
the 'psms' DataFrame (psm.tsv) to the canonical column schema used by all
other processors and HTML tabs.

Input psm.tsv columns used
---------------------------
Spectrum File     path to interact-<name>.pep.xml  → stem stripped of interact- prefix
Peptide           bare peptide sequence
Modified Peptide  MSFragger notation: n[NOM]PEPTIDEK[NOM] (nominal masses in brackets)
                  n  = N-terminal mod,  K/R/... = side-chain mod
                  e.g. n[309]AEGGGGGGRPGAPAAGDGK[436]
Charge            precursor charge state
Retention         retention time in seconds  → converted to minutes (÷60)
Intensity         MS1 peak intensity
Qvalue            PSM-level FDR from PeptideProphet
Probability       PeptideProphet probability of being correct (1 − PEP)
Assigned Modifications
                  accurate-mass mod list, e.g.:
                    "N-term(308.1161), 5K(308.1161), 3C(57.0214)"
                  N-terminal: "N-term(mass)" or "1X(mass)" (position 1)
                  K side-chain: "NNK(mass)" where NN = 1-based residue position
                  This column is used for tag-mass auto-detection and labeling flags.
Protein ID        semicolon-separated protein accessions
Is Decoy          bool; rows where True are removed

Canonical columns produced in report
--------------------------------------
Raw file          ← Spectrum File stem, interact- prefix stripped
Sequence          ← Peptide (bare sequence)
Modified sequence ← Modified Peptide (MSFragger notation, kept as-is)
Charge            ← Charge (Int32)
Retention time    ← Retention / 60  (minutes)
Intensity         ← Intensity
Q.Value           ← Qvalue
PEP               ← 1 − Probability  (clipped to [0, 1])
Protein.Group     ← Protein ID
Precursor.Mz      ← Calibrated Observed M/Z  (when present)
Hyperscore        ← Hyperscore  (when present)

Derived columns
---------------
seqcharge      "<Sequence>_<Charge>"  (charge-specific, tag-free)
Precursor.Id   same as seqcharge
C_term_aa      last residue of Sequence
N_term_aa      first residue of Sequence
fully_labeled  1 if N-term has auto-detected tag mass AND all K residues carry it
partly_labeled 1 if N-term tagged but at least one K is untagged
               Absent when no amine-reactive tag is detected (label-free data).

Tag mass detection
------------------
Auto-detected from "Assigned Modifications" by finding mass values that appear
on BOTH "N-term(mass)" positions AND "NNK(mass)" K side-chains.  Excludes
single-site mods: Cys CAM (+57.02, C only), Met oxidation (+15.99, M only),
N-terminal acetylation (+42, N-term only), etc.

N-terminal labeling is detected as:
  • "N-term(mass)"       — PeptideProphet standard notation
  • "1X(mass)"           — position-1 format some FragPipe versions emit

Quality filters applied during normalize()
------------------------------------------
Is Decoy == True rows removed
Qvalue ≤ 0.01
"""

from __future__ import annotations

import re
from collections import Counter
from typing import Dict, List, Optional, Tuple

import polars as pl


# ---------------------------------------------------------------------------
# Module-level helpers
# ---------------------------------------------------------------------------

def _extract_run_names(df: pl.DataFrame) -> pl.DataFrame:
    """
    Derive 'Raw file' from the 'Spectrum File' column.

    PeptideProphet writes full paths like:
        interact-2025-03-27_HY_T6d0-4ng_DDA1.pep.xml
        C:/path/to/interact-samplename.pep.xml

    DIA analyses add _rankN suffixes for chimeric spectrum ranks:
        interact-2025-04-28_NAME_rank1.pep.xml  → NAME
        interact-2025-04-28_NAME_rank2.pep.xml  → NAME  (same run)

    Steps:
      1. Normalise path separators  (backslash → forward slash)
      2. Extract the filename (last path component)
      3. Strip .pep.xml suffix (and any remaining single-dot extension)
      4. Strip leading interact- prefix
      5. Strip trailing _rankN suffix (DIA chimeric ranks)
    """
    if 'Spectrum File' not in df.columns:
        return df
    df = df.with_columns(
        pl.col('Spectrum File')
          .str.replace_all(r'\\', '/', literal=False)
          .str.extract(r'([^/]+)$', group_index=1)          # filename only
          .str.replace(r'\.pep\.xml$', '', literal=False)   # strip .pep.xml
          .str.replace(r'\.[^./]+$', '', literal=False)     # strip remaining ext
          .str.replace(r'^interact-', '', literal=False)    # strip interact- prefix
          .str.replace(r'_rank\d+$', '', literal=False)     # strip DIA _rankN suffix
          .alias('Raw file')
    ).drop('Spectrum File')
    return df


def _detect_tag_masses(mods_series: pl.Series) -> List[str]:
    """
    Auto-detect amine-reactive tag masses from the 'Assigned Modifications' column.

    FragPipe format: "N-term(308.1161), 5K(308.1161), 3C(57.0214)"
      N-terminal: "N-term(mass)"  or  "1X(mass)" (first-position variant)
      K side-chain: "<pos>K(mass)"

    Returns sorted list of mass strings that appear on BOTH N-terminal and K
    positions across all rows (amine-reactive = hits both sites).
    """
    nterm_re  = re.compile(r'N-term\(([0-9]+\.?[0-9]*)\)')
    first_re  = re.compile(r'(?:^|,\s*)1[A-Z]\(([0-9]+\.?[0-9]*)\)')
    k_re      = re.compile(r'\d+K\(([0-9]+\.?[0-9]*)\)')

    nterm_masses: Counter = Counter()
    k_masses:     Counter = Counter()

    for mods in mods_series:
        if not mods:
            continue
        s = str(mods)
        m = nterm_re.search(s)
        if m:
            nterm_masses[m.group(1)] += 1
        else:
            m2 = first_re.search(s)
            if m2:
                nterm_masses[m2.group(1)] += 1
        for km in k_re.finditer(s):
            k_masses[km.group(1)] += 1

    tag = set(nterm_masses) & set(k_masses)
    return sorted(tag)


def _add_derived(df: pl.DataFrame,
                 tag_masses: Optional[List[str]] = None) -> pl.DataFrame:
    """
    Add seqcharge, Precursor.Id, terminal AAs, and labeling columns.

    Labeling logic (uses 'Assigned Modifications'):
      fully_labeled  = N-term carries any detected tag mass
                       AND every K in the sequence carries any tag mass
      partly_labeled = N-term carries any detected tag mass
                       AND at least one K in the sequence does not carry a tag mass

    For R-terminated peptides (0 K residues), fully_labeled = N-term_tagged,
    partly_labeled = False — consistent with the convention in other processors.
    """
    exprs = []

    if 'Sequence' in df.columns and 'Charge' in df.columns:
        bare = pl.col('Sequence') + '_' + pl.col('Charge').cast(pl.String)
        exprs.append(bare.alias('seqcharge'))
        # Precursor.Id is modification-aware (Modified sequence + charge);
        # seqcharge stays tag-free (bare sequence + charge) for the Labeling tab.
        if 'Modified sequence' in df.columns:
            exprs.append((pl.col('Modified sequence') + '_' + pl.col('Charge').cast(pl.String)).alias('Precursor.Id'))
        else:
            exprs.append(bare.alias('Precursor.Id'))

    if 'Sequence' in df.columns:
        exprs += [
            pl.col('Sequence').str.slice(-1).alias('C_term_aa'),
            pl.col('Sequence').str.slice(0, 1).alias('N_term_aa'),
        ]

    if exprs:
        df = df.with_columns(exprs)

    # Labeling flags from Assigned Modifications
    if 'Assigned Modifications' in df.columns and 'Sequence' in df.columns and tag_masses:
        mods      = pl.col('Assigned Modifications').fill_null('')
        n_k_total = pl.col('Sequence').str.count_matches('K')

        nterm_pats:   list = []
        k_count_exprs: list = []

        for mass in tag_masses:
            e = mass.replace('.', r'\.')
            # N-terminal: "N-term(mass)" or "1X(mass)" (first-position variant)
            nterm_pats.append(mods.str.contains(rf'N-term\({e}\)', literal=False))
            nterm_pats.append(mods.str.contains(rf'(?:^|, )1[A-Z]\({e}\)', literal=False))
            # K side-chain at any position
            k_count_exprs.append(mods.str.count_matches(rf'\d+K\({e}\)', literal=False))

        nterm_any  = pl.any_horizontal(nterm_pats)
        n_k_tagged = pl.sum_horizontal(k_count_exprs)

        fully  = nterm_any & (n_k_tagged >= n_k_total)
        partly = nterm_any & (n_k_tagged < n_k_total)

        df = df.with_columns([
            fully.cast(pl.Int8).alias('fully_labeled'),
            partly.cast(pl.Int8).alias('partly_labeled'),
        ])

    return df


def _build_plex_report(df: pl.DataFrame,
                       tag_masses: List[str]) -> Optional[pl.DataFrame]:
    """
    Build plex_report for multiplexed data by extracting the N-terminal tag
    mass (from 'Assigned Modifications') as the 'Channel' column.

    Filters to rows whose N-terminal mod matches one of the detected tag masses.
    Returns None when no rows match or 'Assigned Modifications' is absent.
    """
    if 'Assigned Modifications' not in df.columns or not tag_masses:
        return None

    plex = df.with_columns(
        pl.col('Assigned Modifications')
          .str.extract(r'N-term\(([0-9]+\.?[0-9]*)\)', group_index=1)
          .alias('Channel')
    )
    plex = plex.filter(pl.col('Channel').is_in(tag_masses))

    if plex.is_empty():
        return None

    n_ch = plex['Channel'].n_unique()
    print(f"  [fragpipe] plex_report: {len(plex):,} rows, {n_ch} channels: {sorted(tag_masses)}")
    return plex


def _filter_quality(df: pl.DataFrame) -> Tuple[pl.DataFrame, List[str]]:
    """
    Remove decoy rows and apply Qvalue ≤ 0.01 filter.
    Returns (filtered_df, list_of_filter_descriptions).
    """
    filters: List[str] = []

    if 'Is Decoy' in df.columns:
        before = len(df)
        df = df.filter(~pl.col('Is Decoy'))
        removed = before - len(df)
        if removed:
            print(f"  [fragpipe] removed {removed:,} decoy rows → {len(df):,} remaining")

    if 'Q.Value' in df.columns:
        before = len(df)
        df = df.filter(pl.col('Q.Value') <= 0.01)
        print(f"  [fragpipe] Qvalue ≤ 0.01: {before:,} → {len(df):,}")
        filters.append('Qvalue ≤ 0.01')

    return df, filters


# ---------------------------------------------------------------------------
# Public processor class
# ---------------------------------------------------------------------------

class FragpipeProcessor:
    """
    Normalises and enriches FragPipe output loaded by data_loader.load_folder().

    Parameters
    ----------
    data : dict[str, pl.DataFrame]
        Output of load_folder() when engine == 'fragpipe'.
        Expected key: 'psms' (required) — psm.tsv content.
    """

    def __init__(self, data: Dict[str, pl.DataFrame]):
        if 'psms' not in data or data['psms'] is None:
            raise ValueError("FragpipeProcessor requires a 'psms' DataFrame in data dict.")

        self.report:        pl.DataFrame           = data['psms']
        # Attributes mirroring other processor interfaces (FragPipe DDA has none of these)
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

    def normalize(self) -> 'FragpipeProcessor':
        """
        Run the full normalization pipeline:
          1.  Extract 'Raw file' from 'Spectrum File' (strip interact- prefix)
          2.  Rename 'Peptide' → 'Sequence'
          3.  Rename 'Modified Peptide' → 'Modified sequence'
          4.  Convert Retention (seconds) → Retention time (minutes)
          5.  Derive PEP = 1 − Probability
          6.  Rename Qvalue → Q.Value, Protein ID → Protein.Group, etc.
          7.  Cast Charge to Int32
          8.  Detect amine-reactive tag masses from Assigned Modifications
          9.  Add derived columns (seqcharge, Precursor.Id, terminal AAs, labels)
         10.  Save raw_report (pre-quality-filter, for Summary tab)
         11.  Apply quality filters (Is Decoy removed, Qvalue ≤ 0.01)
         12.  Build plex_report when > 1 tag mass detected
        Returns self for chaining.
        """
        df = self.report

        # --- Run names ---
        df = _extract_run_names(df)

        # --- Sequence ---
        if 'Peptide' in df.columns:
            df = df.rename({'Peptide': 'Sequence'})

        # --- Modified sequence ---
        if 'Modified Peptide' in df.columns:
            df = df.rename({'Modified Peptide': 'Modified sequence'})

        # --- Retention time: seconds → minutes ---
        if 'Retention' in df.columns:
            df = df.with_columns(
                (pl.col('Retention').cast(pl.Float64, strict=False) / 60.0)
                  .alias('Retention time')
            ).drop('Retention')

        # --- PEP = 1 - Probability ---
        if 'Probability' in df.columns:
            df = df.with_columns(
                (1.0 - pl.col('Probability').cast(pl.Float64, strict=False))
                  .clip(0.0, 1.0)
                  .alias('PEP')
            ).drop('Probability')

        # --- Column renames ---
        rmap: Dict[str, str] = {}
        if 'Qvalue' in df.columns:
            rmap['Qvalue'] = 'Q.Value'
        if 'Protein ID' in df.columns:
            rmap['Protein ID'] = 'Protein.Group'
        if 'Calibrated Observed M/Z' in df.columns:
            rmap['Calibrated Observed M/Z'] = 'Precursor.Mz'
        if rmap:
            df = df.rename(rmap)

        # --- Charge: cast to Int32 ---
        if 'Charge' in df.columns:
            df = df.with_columns(pl.col('Charge').cast(pl.Int32, strict=False))

        # --- Tag mass detection (pre-filter = max coverage) ---
        tag_masses: List[str] = []
        if 'Assigned Modifications' in df.columns:
            tag_masses = _detect_tag_masses(df['Assigned Modifications'])
            if tag_masses:
                print(f"  [fragpipe] tag masses detected: {tag_masses}")

        # --- Derived columns ---
        df = _add_derived(df, tag_masses)

        self.raw_report = df  # pre-quality-filter snapshot for Summary tab

        # --- Quality filter ---
        df, self.filters = _filter_quality(df)

        # --- Deduplication: keep best-intensity PSM per (run, Precursor.Id) ---
        if 'Precursor.Id' in df.columns and 'Intensity' in df.columns:
            before = len(df)
            df = (df.sort('Intensity', descending=True, nulls_last=True)
                    .unique(subset=['Raw file', 'Precursor.Id'], keep='first'))
            removed = before - len(df)
            if removed:
                print(f"  [fragpipe] dedup: {before:,} → {len(df):,} (dropped {removed:,} lower-intensity duplicate PSMs)")

        self.report = df
        print(f"  [fragpipe] normalized report: {len(self.report):,} rows × {self.report.width} columns")

        # --- Plex report (only for multiplexed data with > 1 channel mass) ---
        if len(tag_masses) > 1:
            self.plex_report = _build_plex_report(self.report, tag_masses)

        return self

    def filter_enzyme(self, enzymes: List[str]) -> 'FragpipeProcessor':
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
