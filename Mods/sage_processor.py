# Created by @thibaultcolombani with the help of Claude Code (Anthropic)
"""
sage_processor.py - Sage data processor for MassSleuth.

Takes the dict produced by load_folder() for engine='sage' and normalises
the 'psms' DataFrame (concatenated *.sage.tsv files) to the canonical column
schema.  When lfq.tsv is also present it is merged to supply proper MS1
intensities; ms2_intensity is kept as Precursor.Quantity for the Ion Sampling
MS2 sections.

Input *.sage.tsv columns used
------------------------------
filename          full path to source mzML/raw file → stem as 'Raw file'
peptide           modified sequence, Sage/ProForma bracket notation:
                    N-terminal: [+mass]-PEPTIDE  or  [+mass]PEPTIDE
                    Side-chain: PEPTIDEK[+mass]
                  Retained as-is in 'Modified sequence'.
charge            precursor charge state
rt                retention time in minutes (per Sage documentation)
ms2_intensity     total MS2 spectrum intensity → Precursor.Quantity
spectrum_q        spectrum-level q-value → Q.Value  (primary quality filter)
peptide_q         peptide-level q-value
protein_q         protein-level q-value → Protein.Q.Value
proteins          semicolon-separated protein accessions → Protein.Group
precursor_mz      precursor m/z → Precursor.Mz (when present)
posterior_error   ln(PEP): natural-log posterior error probability.
                  Converted to linear with exp() → PEP column [0, 1].
label             1 = target, -1 = decoy

Input lfq.tsv columns used (when present)
------------------------------------------
peptide           ProForma modified sequence (matches results.sage.tsv peptide)
charge            always -1 in current Sage output; ignored for joining
                  (join is on peptide × Raw file only)
q_value           LFQ-specific q-value (typically 0.09–0.28; not used for
                  filtering — rely on PSM spectrum_q instead)
spectral_angle    normalized spectral contrast angle for the MS1 integration
<filename>.mzML   MS1 peak-area intensity for that file (wide format)

Canonical columns produced in report
--------------------------------------
Raw file          ← stem of filename
Sequence          ← peptide stripped of all [mod] tokens and leading '-'
Modified sequence ← peptide as-is (ProForma bracket notation)
Charge            ← charge (cast to Int32)
Retention time    ← rt (already in minutes per Sage docs; no conversion)
Intensity         ← lfq MS1 area when available, else ms2_intensity
Precursor.Quantity← ms2_intensity (MS2 total intensity per spectrum)
Q.Value           ← spectrum_q
Protein.Q.Value   ← protein_q
Protein.Group     ← proteins
Precursor.Mz      ← precursor_mz (when present)
PEP               ← exp(posterior_error)  (linear PEP, [0, 1])
Hyperscore        ← hyperscore (when present)

Derived columns
---------------
seqcharge      "<Sequence>_<Charge>"  (tag-free, charge-specific)
Precursor.Id   same as seqcharge
C_term_aa      last residue of Sequence
N_term_aa      first residue of Sequence
fully_labeled  1 if all expected amine sites (N-term + all K) carry the
               auto-detected tag mass (mass present on BOTH N-term and K positions)
partly_labeled 1 if at least one but not all amine sites carry the tag mass
               Columns absent when no tag mass is detected (unlabeled data)

Quality filters applied during normalize()
------------------------------------------
Decoy rows (label == -1) removed
spectrum_q ≤ 0.01
"""

from __future__ import annotations

import re
from collections import Counter
from typing import Dict, List, Optional, Tuple

import polars as pl


# ---------------------------------------------------------------------------
# Column rename map  (Sage raw name → canonical name)
# ---------------------------------------------------------------------------

_COL_MAP: Dict[str, str] = {
    'charge':       'Charge',
    'proteins':     'Protein.Group',
    'precursor_mz': 'Precursor.Mz',
    'spectrum_q':   'Q.Value',
    'protein_q':    'Protein.Q.Value',
    'hyperscore':   'Hyperscore',
    # rt → Retention time handled separately (direct rename; already in minutes)
    # ms2_intensity → Precursor.Quantity handled separately
    # posterior_error → PEP handled separately (requires exp() conversion)
}

# Regex patterns for Sage ProForma sequence parsing
_RE_BRACKET_MOD = r'\[[^\]]+\]'   # any [mod] token
_RE_LEADING_DASH = r'^-'           # optional dash between N-term mod and sequence


# LFQ meta-columns (everything else is a filename-intensity column)
_LFQ_META = {'peptide', 'charge', 'proteins', 'q_value', 'score', 'spectral_angle'}


# ---------------------------------------------------------------------------
# Module-level transformation helpers
# ---------------------------------------------------------------------------

def _rename(df: pl.DataFrame, col_map: Dict[str, str]) -> pl.DataFrame:
    rmap = {k: v for k, v in col_map.items() if k in df.columns}
    return df.rename(rmap) if rmap else df


def _extract_run_names(df: pl.DataFrame) -> pl.DataFrame:
    """Replace full-path 'filename' with the mzML stem as 'Raw file'."""
    if 'filename' not in df.columns:
        return df
    df = df.with_columns(
        pl.col('filename')
          .str.replace_all(r'\\', '/', literal=False)
          .str.extract(r'([^/]+)\.[^./]+$', group_index=1)
          .alias('Raw file')
    )
    return df.drop('filename')


def _extract_sequence(df: pl.DataFrame) -> pl.DataFrame:
    """
    Derive 'Sequence' (stripped) from the 'peptide' ProForma column.
    Removes all [mod] tokens and any leading dash left by N-terminal mod removal.
    Renames 'peptide' → 'Modified sequence'.
    """
    if 'peptide' not in df.columns:
        return df
    df = df.with_columns(
        pl.col('peptide')
          .str.replace_all(_RE_BRACKET_MOD, '', literal=False)
          .str.replace(_RE_LEADING_DASH, '', literal=False)
          .alias('Sequence')
    )
    return df.rename({'peptide': 'Modified sequence'})


def _cast_columns(df: pl.DataFrame) -> pl.DataFrame:
    """Cast charge to Int32."""
    col = 'charge' if 'charge' in df.columns else ('Charge' if 'Charge' in df.columns else None)
    if col:
        df = df.with_columns(pl.col(col).cast(pl.Int32))
    return df


def _convert_posterior_error(df: pl.DataFrame) -> pl.DataFrame:
    """Drop posterior_error (ln PEP). Sage's primary quality metric is spectrum_q → Q.Value."""
    if 'posterior_error' not in df.columns:
        return df
    return df.drop('posterior_error')


def _assign_ms2_as_precursor_quantity(df: pl.DataFrame) -> pl.DataFrame:
    """
    Rename ms2_intensity → Precursor.Quantity.
    This enables Ion Sampling MS2 sections (which expect Precursor.Quantity).
    The canonical Intensity column will later be filled from lfq MS1 data.
    """
    if 'ms2_intensity' in df.columns:
        return df.rename({'ms2_intensity': 'Precursor.Quantity'})
    return df


def _filter_decoys(df: pl.DataFrame) -> pl.DataFrame:
    """Remove rows where label == -1 (Sage decoy convention)."""
    if 'label' not in df.columns:
        return df
    before = len(df)
    df = df.filter(pl.col('label') != -1)
    removed = before - len(df)
    if removed:
        print(f"  [sage] removed {removed:,} decoy rows → {len(df):,} remaining")
    return df


def _filter_quality(df: pl.DataFrame) -> Tuple[pl.DataFrame, List[str]]:
    """Filter spectrum_q ≤ 0.01 (already renamed to Q.Value)."""
    filters: List[str] = []
    if 'Q.Value' in df.columns:
        before = len(df)
        df = df.filter(pl.col('Q.Value') <= 0.01)
        print(f"  [sage] spectrum_q ≤ 0.01: {before:,} → {len(df):,}")
        filters.append('spectrum_q ≤ 0.01')
    return df, filters


def _detect_tag_masses(mod_seqs: pl.Series) -> List[str]:
    """
    Auto-detect amine-reactive tag masses from ProForma-annotated sequences.

    A tag mass is one that appears on BOTH N-terminal positions (^[+X]) AND
    K side-chain positions (K[+X]).  This correctly excludes single-site mods
    such as N-terminal acetylation, Met oxidation, and Cys carbamidomethylation,
    which appear at only one residue type.

    Works for any tag family (PSMtag, TMT, mTRAQ, iTRAQ, …).
    Returns sorted list of detected mass strings (e.g. ['308.1161']).
    """
    _nterm_re = re.compile(r'^\[\+([0-9]+\.?[0-9]*)\]')
    _k_re     = re.compile(r'K\[\+([0-9]+\.?[0-9]*)\]')

    nterm_masses: Counter = Counter()
    k_masses:     Counter = Counter()

    for seq in mod_seqs:
        if not seq:
            continue
        m = _nterm_re.match(seq)
        if m:
            nterm_masses[m.group(1)] += 1
        for m in _k_re.finditer(seq):
            k_masses[m.group(1)] += 1

    tag = set(nterm_masses) & set(k_masses)
    return sorted(tag)


def _add_derived(df: pl.DataFrame, tag_masses: Optional[List[str]] = None) -> pl.DataFrame:
    """Add seqcharge, Precursor.Id, terminal amino acids, and labeling columns."""
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
        exprs.append(pl.col('Sequence').str.slice(-1).alias('C_term_aa'))
        exprs.append(pl.col('Sequence').str.slice(0, 1).alias('N_term_aa'))

    if exprs:
        df = df.with_columns(exprs)

    # Labeling detection: auto-detect tag mass as the modification mass that
    # appears on BOTH N-terminal and K side-chain positions. This works for any
    # amine-reactive tag (PSMtag, TMT, mTRAQ, …) and correctly excludes
    # single-site mods (acetylation on N-term only, oxidation on M, CAM on C).
    # If no such mass is found the data is unlabeled → columns not added.
    #
    # For multiplexed data (multiple tag masses detected), fully_labeled requires
    # the N-terminal mass and ALL K masses to be the SAME plex channel — a K
    # carrying a different channel mass is NOT fully labeled.
    if 'Modified sequence' in df.columns and 'Sequence' in df.columns:
        if tag_masses is None:
            tag_masses = _detect_tag_masses(df['Modified sequence'])
        if tag_masses:
            mod_seq   = pl.col('Modified sequence')
            n_k_total = pl.col('Sequence').str.count_matches('K')
            escaped   = [m.replace('.', r'\.') for m in tag_masses]
            alt       = '|'.join(escaped)
            nterm_any = mod_seq.str.contains(rf'^\[\+(?:{alt})\]', literal=False)

            # Count K residues carrying any plex tag (regardless of which channel)
            n_k_tagged = pl.sum_horizontal([
                mod_seq.str.count_matches(rf'K\[\+{e}\]') for e in escaped
            ])
            fully = nterm_any & (n_k_tagged >= n_k_total)
            partly = nterm_any & (n_k_tagged < n_k_total)

            df = df.with_columns([
                fully.cast(pl.Int8).alias('fully_labeled'),
                partly.cast(pl.Int8).alias('partly_labeled'),
            ])

    return df


def _build_plex_report(df: pl.DataFrame, tag_masses: List[str]) -> Optional[pl.DataFrame]:
    """
    Build plex_report for multiplexed Sage data by extracting the N-terminal tag
    mass as 'Channel'. Each PSM row already represents one channel observation.

    Filters to rows where the N-terminal mass is one of the detected tag masses
    (excludes PSMs with no N-terminal tag, acetylated N-terms, etc.).
    Returns None if no rows match.
    """
    if 'Modified sequence' not in df.columns or not tag_masses:
        return None

    plex = df.with_columns(
        pl.col('Modified sequence')
          .str.extract(r'^\[\+([0-9]+\.?[0-9]*)\]', group_index=1)
          .alias('Channel')
    )
    plex = plex.filter(pl.col('Channel').is_in(tag_masses))

    if plex.is_empty():
        return None

    # Use tag-free seqcharge as Precursor.Id for cross-channel comparisons.
    # Modified sequence includes the channel-specific N-terminal mass, so
    # modification-aware IDs are disjoint across channels → intersection = 0.
    if 'seqcharge' in plex.columns:
        plex = plex.with_columns(pl.col('seqcharge').alias('Precursor.Id'))

    n_ch = plex['Channel'].n_unique()
    print(f"  [sage] plex_report: {len(plex):,} rows, {n_ch} channels: {sorted(tag_masses)}")
    return plex


# ---------------------------------------------------------------------------
# LFQ helpers
# ---------------------------------------------------------------------------

def _melt_lfq(lfq: pl.DataFrame) -> Optional[pl.DataFrame]:
    """
    Melt wide-format lfq.tsv to long format: one row per (peptide, Raw file).

    The filename columns (everything not in _LFQ_META) are the actual mzML
    file names used as column headers.  Their stems are extracted as 'Raw file'.
    The 'charge' column in lfq.tsv is always -1 and is dropped.

    Returns None if no filename columns are found.
    """
    file_cols = [c for c in lfq.columns if c not in _LFQ_META]
    if not file_cols:
        return None

    keep_meta = [c for c in ('peptide', 'spectral_angle') if c in lfq.columns]

    melted = (
        lfq.select(keep_meta + file_cols)
           .unpivot(
               on=file_cols,
               index=keep_meta,
               variable_name='_filename',
               value_name='Intensity',
           )
    )

    # Extract run name from column header (mzML filename → stem)
    melted = melted.with_columns(
        pl.col('_filename')
          .str.replace_all(r'\\', '/', literal=False)
          .str.extract(r'([^/]+)\.[^./]+$', group_index=1)
          .alias('Raw file')
    ).drop('_filename')

    # Rename peptide → Modified sequence for joining
    melted = melted.rename({'peptide': 'Modified sequence'})

    # Ensure numeric type (unpivot may produce String when columns have mixed dtypes)
    melted = melted.with_columns(pl.col('Intensity').cast(pl.Float64, strict=False))
    # Keep only rows with a positive MS1 intensity
    melted = melted.filter(pl.col('Intensity') > 0)

    print(f"  [sage] lfq: {len(melted):,} peptide × run entries with MS1 intensity")
    return melted


def _merge_lfq(df: pl.DataFrame, lfq_long: pl.DataFrame) -> pl.DataFrame:
    """
    Left-join lfq MS1 intensities onto the PSM report by (Modified sequence,
    Raw file).  Creates 'Intensity' from lfq MS1 where available; falls back
    to Precursor.Quantity (ms2_intensity) for PSMs not covered by lfq.
    """
    lfq_select = ['Modified sequence', 'Raw file', 'Intensity']
    if 'spectral_angle' in lfq_long.columns:
        lfq_select.append('spectral_angle')

    merged = df.join(
        lfq_long.select(lfq_select),
        on=['Modified sequence', 'Raw file'],
        how='left',
    )

    # Where lfq has no MS1 → fall back to ms2_intensity (Precursor.Quantity)
    if 'Precursor.Quantity' in merged.columns:
        merged = merged.with_columns(
            pl.when(pl.col('Intensity').is_not_null() & (pl.col('Intensity') > 0))
              .then(pl.col('Intensity'))
              .otherwise(pl.col('Precursor.Quantity'))
              .alias('Intensity')
        )

    n_ms1 = merged.filter(pl.col('Intensity') == merged.join(
        lfq_long.select(['Modified sequence', 'Raw file', 'Intensity']),
        on=['Modified sequence', 'Raw file'], how='left')['Intensity']
    ).height if False else None  # skip expensive re-join; just report coverage

    ms1_covered = merged['Intensity'].is_not_null().sum()
    print(f"  [sage] MS1 intensity (lfq): {ms1_covered:,} / {len(merged):,} PSMs covered")
    return merged


# ---------------------------------------------------------------------------
# Public processor class
# ---------------------------------------------------------------------------

class SageProcessor:
    """
    Normalises and enriches Sage output loaded by data_loader.load_folder().

    Parameters
    ----------
    data : dict[str, pl.DataFrame]
        Output of load_folder() when engine == 'sage'.
        Expected key: 'psms' (required) — concatenated *.sage.tsv files.
        Optional key: 'lfq' — wide-format lfq.tsv with MS1 intensities.
    """

    def __init__(self, data: Dict[str, pl.DataFrame]):
        if 'psms' not in data or data['psms'] is None:
            raise ValueError("SageProcessor requires a 'psms' DataFrame in data dict.")

        self.report:        pl.DataFrame           = data['psms']
        self.lfq:           Optional[pl.DataFrame] = data.get('lfq')
        # Attributes mirroring other processor interfaces
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

    def normalize(self) -> 'SageProcessor':
        """
        Run the full normalization pipeline:
          1. Extract run names from filename paths
          2. Extract Sequence from peptide; rename peptide → Modified sequence
          3. Cast charge → Int32
          4. Convert posterior_error: exp() → linear PEP
          5. Rename ms2_intensity → Precursor.Quantity
          6. Remove decoy rows (label == -1)
          7. Rename columns to canonical names
          8. Add derived columns (seqcharge, Precursor.Id, terminal AAs, labels)
          9. Save raw_report (pre-quality-filter, for Summary tab)
         10. Apply quality filter (spectrum_q ≤ 0.01)
         11. Merge lfq MS1 intensities as Intensity (if lfq.tsv was loaded)
        Returns self for chaining.
        """
        df = self.report

        df = _extract_run_names(df)
        df = _extract_sequence(df)
        df = _cast_columns(df)
        df = _convert_posterior_error(df)
        df = _assign_ms2_as_precursor_quantity(df)
        df = _filter_decoys(df)
        df = _rename(df, _COL_MAP)
        # rt is already in minutes per Sage documentation
        if 'rt' in df.columns:
            df = df.rename({'rt': 'Retention time'})

        # Detect tag masses once from pre-filter data (maximum coverage),
        # then reuse for labeling flags and plex_report construction.
        tag_masses: List[str] = []
        if 'Modified sequence' in df.columns:
            tag_masses = _detect_tag_masses(df['Modified sequence'])
            if tag_masses:
                print(f"  [sage] tag masses detected: {tag_masses}")

        df = _add_derived(df, tag_masses)

        self.raw_report = df  # pre-quality-filter snapshot for Summary tab

        df, self.filters = _filter_quality(df)

        # Merge lfq MS1 intensities
        if self.lfq is not None and not self.lfq.is_empty():
            lfq_long = _melt_lfq(self.lfq)
            if lfq_long is not None:
                df = _merge_lfq(df, lfq_long)
        else:
            # No lfq: Intensity = Precursor.Quantity (ms2_intensity)
            if 'Precursor.Quantity' in df.columns and 'Intensity' not in df.columns:
                df = df.with_columns(pl.col('Precursor.Quantity').alias('Intensity'))

        self.report = df
        print(f"  [sage] normalized report: {len(self.report):,} rows × {self.report.width} columns")

        # Build plex_report when multiple tag masses are present (multiplexed)
        if len(tag_masses) > 1:
            self.plex_report = _build_plex_report(self.report, tag_masses)

        return self

    def filter_enzyme(self, enzymes: List[str]) -> 'SageProcessor':
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
