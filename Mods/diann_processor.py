# Created by @thibaultcolombani with the help of Claude Code (Anthropic)
"""
diann_processor.py - DIA-NN data processor for MassSleuth.

Takes the dict produced by load_folder() for engine='diann' and normalises
the report DataFrame to the canonical column schema, adds derived columns
for QC and labeling-efficiency analysis, and exposes the supplementary
DataFrames (features, fill_times, tic, sn, ms1_extracted).

Canonical columns produced in report
-------------------------------------
Raw file          ← Run  (+ channel suffix when multiplexed)
Sequence          ← Stripped.Sequence
Modified sequence ← Modified.Sequence
Charge            ← Precursor.Charge
Retention time    ← RT  (minutes)
Retention time start / stop  ← RT.Start / RT.Stop
Retention length  ← FWHM
Intensity         ← Ms1.Area
PEP               (unchanged)
Q.Value           (unchanged)
Protein.Group     (unchanged)
Protein.Ids       (unchanged)
Protein.Names     (unchanged)
Genes             (unchanged)
IM                (unchanged, ion mobility)
Channel           (unchanged)
Channel.Q.Value   (unchanged)
Translated.Q.Value (unchanged)
PG.Q.Value        (unchanged)
Protein.Q.Value   (unchanged)
GG.Q.Value        (unchanged)
Proteotypic       (unchanged)
Decoy             (unchanged, removed rows kept via filter)
Precursor.Id      (unchanged)
Precursor.Mz      (unchanged)

Derived columns
---------------
seqcharge      "ModifiedSequence_Charge"
C_term_aa      last residue of Sequence (A–Z)
N_term_aa      first residue of Sequence (A–Z)
fully_labeled  1 if all expected amine sites (N-term + all K side chains) carry a (tag) label
partly_labeled 1 if at least one amine site is labeled but not all

Quality filters applied during normalize()
------------------------------------------
Non-multiplexed : PEP ≤ 0.01  AND  Q.Value ≤ 0.01
Multiplexed     : PEP ≤ 0.01  AND  Channel.Q.Value ≤ 0.01
"""

from __future__ import annotations

from typing import Dict, List, Optional, Tuple

import polars as pl


# ---------------------------------------------------------------------------
# Column rename map  (DIA-NN raw name → canonical name)
# ---------------------------------------------------------------------------

_COL_MAP: Dict[str, str] = {
    'Run':               'Raw file',
    'Stripped.Sequence': 'Sequence',
    'Modified.Sequence': 'Modified sequence',
    'Precursor.Charge':  'Charge',
    'RT':                'Retention time',
    'RT.Start':          'Retention time start',
    'RT.Stop':           'Retention time stop',
    'FWHM':              'Retention length',
    'Ms1.Area':          'Intensity',
}

# Features.tsv column rename map
_FEATURES_COL_MAP: Dict[str, str] = {
    'Raw.file': 'Raw file',
    'rtStart':  'Retention time start',
    'rtApex':   'Retention time apex',
    'rtEnd':    'Retention time end',
    'mz':       'Precursor.Mz',
}


# ---------------------------------------------------------------------------
# Module-level transformation helpers
# ---------------------------------------------------------------------------

def _rename(df: pl.DataFrame, col_map: Dict[str, str]) -> pl.DataFrame:
    rename = {k: v for k, v in col_map.items() if k in df.columns}
    return df.rename(rename) if rename else df


def _filter_decoys(df: pl.DataFrame) -> pl.DataFrame:
    if 'Decoy' in df.columns:
        before = len(df)
        df = df.filter(pl.col('Decoy') == 0)
        print(f"  [diann] removed {before - len(df):,} decoy rows → {len(df):,} remaining")
    return df


def _normalize_channel_labels(df: pl.DataFrame) -> pl.DataFrame:
    """Normalize PSMtag channel labels: 'd0' → '0', 'd4' → '4', etc."""
    if 'Channel' not in df.columns:
        return df
    return df.with_columns(
        pl.col('Channel').str.replace(r'^d(\d+)$', r'$1').alias('Channel')
    )


def _is_multiplexed(df: pl.DataFrame) -> bool:
    if 'Channel' not in df.columns:
        return False
    return df.filter(pl.col('Channel').is_not_null() & (pl.col('Channel') != '')).height > 0


def _filter_quality(df: pl.DataFrame, is_multiplexed: bool) -> Tuple[pl.DataFrame, List[str]]:
    """
    Apply quality filters and return the filtered DataFrame plus a list of
    human-readable filter strings for display in the report header.

    Non-multiplexed: PEP ≤ 0.01  AND  Q.Value ≤ 0.01
    Multiplexed    : PEP ≤ 0.01  AND  Channel.Q.Value ≤ 0.01
    """
    filters: List[str] = []

    if 'PEP' in df.columns:
        before = len(df)
        df = df.filter(pl.col('PEP') <= 0.01)
        print(f"  [diann] PEP ≤ 0.01: {before:,} → {len(df):,}")
        filters.append('PEP ≤ 0.01')

    if is_multiplexed:
        if 'Channel.Q.Value' in df.columns:
            before = len(df)
            df = df.filter(pl.col('Channel.Q.Value') <= 0.01)
            print(f"  [diann] Channel.Q.Value ≤ 0.01: {before:,} → {len(df):,}")
            filters.append('Channel.Q.Value ≤ 0.01')
    else:
        if 'Q.Value' in df.columns:
            before = len(df)
            df = df.filter(pl.col('Q.Value') <= 0.01)
            print(f"  [diann] Q.Value ≤ 0.01: {before:,} → {len(df):,}")
            filters.append('Q.Value ≤ 0.01')

    return df, filters


def _collapse_channels(df: pl.DataFrame) -> pl.DataFrame:
    """
    For multiplexed runs, collapse per-channel rows to one representative row
    per (Raw file, Precursor.Id) by averaging numeric measurements.
    Non-multiplexed rows (Channel empty or absent) pass through unchanged.
    """
    if 'Channel' not in df.columns:
        return df
    has_ch = pl.col('Channel').is_not_null() & (pl.col('Channel') != '')
    mux    = df.filter(has_ch)
    non_mux = df.filter(~has_ch)
    if mux.is_empty():
        return df

    group_keys = [c for c in ['Raw file', 'Precursor.Id'] if c in mux.columns]
    avg_cols   = [c for c in ['Intensity', 'Retention time', 'Retention time start',
                               'Retention time stop', 'Retention length',
                               'Q.Value', 'PEP', 'Channel.Q.Value', 'IM']
                  if c in mux.columns]
    first_cols = [c for c in mux.columns if c not in group_keys + avg_cols]

    collapsed = (
        mux.group_by(group_keys)
           .agg(
               [pl.col(c).mean() for c in avg_cols] +
               [pl.col(c).first() for c in first_cols]
           )
           .select(df.columns)   # restore original column order
    )
    before = mux.height
    print(f"  [diann] channel collapse: {before:,} rows → {collapsed.height:,} "
          f"(÷{before // max(collapsed.height, 1)} channels avg)")
    return pl.concat([non_mux, collapsed], how='diagonal')


def _add_derived(df: pl.DataFrame) -> pl.DataFrame:
    """Add seqcharge, C_term_aa, N_term_aa, and labeling columns."""
    exprs = []

    # seqcharge
    if 'Modified sequence' in df.columns and 'Charge' in df.columns:
        exprs.append(
            (pl.col('Modified sequence') + '_' + pl.col('Charge').cast(pl.String))
            .alias('seqcharge')
        )

    # Terminal amino acids (from stripped Sequence)
    if 'Sequence' in df.columns:
        exprs.append(pl.col('Sequence').str.slice(-1).alias('C_term_aa'))
        exprs.append(pl.col('Sequence').str.slice(0, 1).alias('N_term_aa'))

    if exprs:
        df = df.with_columns(exprs)

    # Labeling columns — DIA-NN plex tag notation: (tag) on N-term, K(tag) on K side-chain
    if 'Modified sequence' in df.columns and 'Sequence' in df.columns:
        mod_seq = pl.col('Modified sequence')

        # N-terminal alpha-amine: sequence starts with optional UniMod mods then (tag)
        has_nterm_label = mod_seq.str.contains(r'^(?:\(UniMod:\d+\))*\(tag\)', literal=False)

        # K side-chain epsilon-amines: count labeled K residues
        n_k_labeled = mod_seq.str.count_matches(r'K\(tag\)')

        # Total K residues in stripped sequence
        n_k_total = pl.col('Sequence').str.count_matches('K')

        n_expected = 1 + n_k_total
        n_actual   = has_nterm_label.cast(pl.Int32) + n_k_labeled

        df = df.with_columns([
            (n_actual >= n_expected).cast(pl.Int8).alias('fully_labeled'),
            ((n_actual > 0) & (n_actual < n_expected)).cast(pl.Int8).alias('partly_labeled'),
        ])

    return df


# ---------------------------------------------------------------------------
# Public processor class
# ---------------------------------------------------------------------------

class DiannProcessor:
    """
    Normalises and enriches DIA-NN output loaded by data_loader.load_folder().

    Parameters
    ----------
    data : dict[str, pl.DataFrame]
        Output of load_folder() when engine == 'diann'.
        Expected keys: 'report' (required), optionally 'features',
        'fill_times', 'tic', 'sn', 'ms1_extracted'.
    """

    def __init__(self, data: Dict[str, pl.DataFrame]):
        if 'report' not in data or data['report'] is None:
            raise ValueError("DiannProcessor requires a 'report' DataFrame in data dict.")

        self.report:        pl.DataFrame           = data['report']
        self.features:      Optional[pl.DataFrame] = data.get('features')
        self.fill_times:    Optional[pl.DataFrame] = data.get('fill_times')
        self.tic:           Optional[pl.DataFrame] = data.get('tic')
        self.sn:            Optional[pl.DataFrame] = data.get('sn')
        self.ms1_extracted: Optional[pl.DataFrame] = data.get('ms1_extracted')
        self.enzymes:       Optional[List[str]]    = None
        self.filters:       List[str]              = []
        self.raw_report:    Optional[pl.DataFrame] = None  # pre-quality-filter, for Summary
        self.plex_report:   Optional[pl.DataFrame] = None  # per-channel rows, for Plex DIA tab

    # ------------------------------------------------------------------
    # Main pipeline
    # ------------------------------------------------------------------

    def normalize(self) -> 'DiannProcessor':
        """
        Run the full normalization pipeline in place:
          1. Rename columns to canonical names
          2. Remove decoy rows
          3. Apply quality filters (PEP + Q.Value or Channel.Q.Value)
          4. Handle multiplexed channels (append suffix to Raw file)
          5. Add derived columns (seqcharge, C_term_aa, labeling flags)
        Also renames features.tsv columns if present.
        Returns self for chaining.
        """
        df = self.report
        df = _rename(df, _COL_MAP)
        df = _normalize_channel_labels(df)
        df = _filter_decoys(df)

        is_mux = _is_multiplexed(df)
        if is_mux:
            n_ch = df.filter(pl.col('Channel').is_not_null() & (pl.col('Channel') != ''))['Channel'].n_unique()
            print(f"  [diann] multiplexed data — {n_ch} channel(s) detected")

        df = _add_derived(df)
        self.raw_report = _collapse_channels(df)  # pre-quality-filter, for Summary

        df, self.filters = _filter_quality(df, is_mux)

        if is_mux:
            self.plex_report = df              # per-channel rows, for Plex DIA tab
            df = _collapse_channels(df)        # merge channels for all other tabs

        self.report = df

        if self.features is not None:
            self.features = _rename(self.features, _FEATURES_COL_MAP)

        print(f"  [diann] normalized report: {len(self.report):,} rows × {self.report.width} columns")
        return self

    def filter_enzyme(self, enzymes: List[str]) -> 'DiannProcessor':
        """
        Filter self.report to fully specific peptides for the given enzyme(s).
        Must be called after normalize(). Returns self for chaining.
        """
        from Mods.enzyme import filter_peptides
        self.report  = filter_peptides(self.report, enzymes)
        self.enzymes = list(enzymes)
        return self

    # ------------------------------------------------------------------
    # Convenience accessors
    # ------------------------------------------------------------------

    def runs(self) -> list[str]:
        """Sorted list of unique run names."""
        return sorted(self.report['Raw file'].unique().to_list())

    def summary(self) -> Dict:
        df = self.report
        out: Dict = {
            'rows':        len(df),
            'runs':        df['Raw file'].n_unique() if 'Raw file' in df.columns else None,
            'precursors':  df['Precursor.Id'].n_unique() if 'Precursor.Id' in df.columns else None,
            'seqcharges':  df['seqcharge'].n_unique() if 'seqcharge' in df.columns else None,
            'proteins':    df['Protein.Group'].n_unique() if 'Protein.Group' in df.columns else None,
        }
        if self.features is not None:
            out['features_rows'] = len(self.features)
        return out
