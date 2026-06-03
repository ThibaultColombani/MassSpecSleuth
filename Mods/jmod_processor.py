# Created by @thibaultcolombani with the help of Claude Code (Anthropic)
"""
jmod_processor.py - Jmod data processor for MassSleuth.

Takes the dict produced by load_folder() for engine='jmod' and normalises
the 'ids' DataFrame to the canonical column schema, adds derived columns
for QC and labeling-efficiency analysis.

Input parquet / CSV columns
---------------------------
stripped_seq       raw peptide sequence
z                  charge state (float → Int32)
untag_prec         tag-free precursor id  "<stripped_seq>_<z>"
file_name          full path to source mzML → stem extracted as run name
channel            integer channel index (0 = LF or untagged; multiple = plex)
is_decoy           bool
Qvalue             channel-level FDR Q-value  (≈ DIA-NN Channel.Q.Value)
Protein_Qvalue     protein-level Q-value
PredVal            PSM confidence score (~0–1, NOT intensity)
protein            protein group
BestChannel_Qvalue best-channel Q-value (≈ DIA-NN Q.Value; used for non-mux filter)
plex_Area          channel-level MS intensity (actual counts)
seq                modified sequence with tag tokens — UNIQUE PER CHANNEL
                   e.g. "A(PSMtag_5plex-0)VAFR" vs "A(PSMtag_5plex-4)VAFR"
                   Do NOT use for cross-channel comparisons; use untag_prec instead.
silac_channel      SILAC channel index (0.0 = no SILAC); kept separate from plex channels
untag_seq          modified sequence without tag tokens
rt                 retention time (minutes)
mz                 precursor m/z
coeff              MS2 quantification coefficient (≈ DIA-NN Precursor.Quantity; CSV only)

Canonical columns produced in report
--------------------------------------
Raw file           ← stem of file_name
Sequence           ← stripped_seq
Modified sequence  ← seq  (channel-unique; use seqcharge for cross-channel work)
Charge             ← z  (cast to Int32)
Retention time     ← rt
Precursor.Mz       ← mz
Precursor.Quantity ← coeff  (CSV only; MS2-based quantification)
Intensity          ← plex_Area  (channel-level MS intensity)
Q.Value            ← BestChannel_Qvalue  (≈ DIA-NN Q.Value)
Channel.Q.Value    ← Qvalue              (≈ DIA-NN Channel.Q.Value)
Protein.Q.Value    ← Protein_Qvalue
Protein.Group      ← protein
Precursor.Id       ← untag_prec
Channel            ← channel  (cast to String, e.g. "0", "4", "8" …)

Derived columns
---------------
seqcharge      = Precursor.Id  (tag-free, channel-independent)
C_term_aa      last residue of Sequence
N_term_aa      first residue of Sequence
R_labeled      1 if C-term R and N-terminal PSMtag label present
K_full_labeled 1 if C-term K and both N-term and K-side PSMtag labels present
K_half_labeled 1 if C-term K and only one PSMtag label site present

Quality filters applied during normalize()
------------------------------------------
Non-multiplexed : BestChannel_Qvalue (→ Q.Value) ≤ 0.01
Multiplexed     : Qvalue (→ Channel.Q.Value) ≤ 0.01
"""

from __future__ import annotations

from pathlib import PurePosixPath, PureWindowsPath
from typing import Dict, List, Optional

import polars as pl


# ---------------------------------------------------------------------------
# Column rename map  (Jmod raw name → canonical name)
# ---------------------------------------------------------------------------

_COL_MAP: Dict[str, str] = {
    'stripped_seq':       'Sequence',
    'z':                  'Charge',
    'seq':                'Modified sequence',
    'Qvalue':             'Channel.Q.Value',  # channel-level q-value (≈ DIA-NN Channel.Q.Value)
    'BestChannel_Qvalue': 'Q.Value',           # best-channel q-value  (≈ DIA-NN Q.Value)
    'coeff':              'Precursor.Quantity', # MS2 quantification   (≈ DIA-NN Precursor.Quantity)
    'Protein_Qvalue':     'Protein.Q.Value',
    'protein':            'Protein.Group',
    'untag_prec':         'Precursor.Id',
    'rt':                 'Retention time',
    'mz':                 'Precursor.Mz',
}
# Intensity is assigned separately: MS1_Area (all_IDs.csv) takes priority over plex_Area (parquet)


# ---------------------------------------------------------------------------
# Module-level transformation helpers
# ---------------------------------------------------------------------------

def _rename(df: pl.DataFrame, col_map: Dict[str, str]) -> pl.DataFrame:
    rmap = {k: v for k, v in col_map.items() if k in df.columns}
    return df.rename(rmap) if rmap else df


def _extract_run_names(df: pl.DataFrame) -> pl.DataFrame:
    """Replace full-path 'file_name' with the mzML stem as 'Raw file'."""
    if 'file_name' not in df.columns:
        return df
    # Works for both Unix (/a/b/name.mzML) and Windows paths (C:\a\b\name.mzML)
    df = df.with_columns(
        pl.col('file_name')
          .str.replace_all(r'\', '/', literal=False)
          .str.extract(r'([^/]+)\.[^./]+$', group_index=1)
          .alias('Raw file')
    )
    return df.drop('file_name')


def _cast_columns(df: pl.DataFrame) -> pl.DataFrame:
    """Cast z → Int32 and channel integer → String Channel."""
    exprs = []
    if 'z' in df.columns:
        exprs.append(pl.col('z').cast(pl.Int32))
    if 'channel' in df.columns:
        exprs.append(pl.col('channel').cast(pl.String).alias('Channel'))
    if exprs:
        df = df.with_columns(exprs)
    if 'channel' in df.columns and 'Channel' in df.columns:
        df = df.drop('channel')
    return df


def _filter_decoys(df: pl.DataFrame) -> pl.DataFrame:
    col = 'is_decoy' if 'is_decoy' in df.columns else ('decoy' if 'decoy' in df.columns else None)
    if col is None:
        return df
    before = len(df)
    df = df.filter(pl.col(col) == False)  # noqa: E712
    print(f"  [jmod] removed {before - len(df):,} decoy rows → {len(df):,} remaining")
    return df


def _is_multiplexed(df: pl.DataFrame) -> bool:
    """True if more than one distinct channel value exists in the dataset."""
    col = 'Channel' if 'Channel' in df.columns else None
    if col is None:
        return False
    return df[col].n_unique() > 1


def _assign_intensity(df: pl.DataFrame):
    """
    Assign the canonical Intensity column and return (df, has_full_data).

    Priority:
      MS1_Area (all_IDs.csv)  → Intensity  [has_full_data = True]
      plex_Area (parquet)     → Intensity  [has_full_data = False]

    When MS1_Area is used, plex_Area is kept separately for the plex_report,
    where channel-level area is needed instead of the precursor-level MS1 area.
    """
    if 'MS1_Area' in df.columns:
        df = df.rename({'MS1_Area': 'Intensity'})
        return df, True
    if 'plex_Area' in df.columns:
        df = df.rename({'plex_Area': 'Intensity'})
        return df, False
    return df, False


def _filter_quality(df: pl.DataFrame, is_mux: bool):
    """Apply quality filter; return (filtered_df, filter_strings)."""
    filters: List[str] = []

    if is_mux:
        if 'Channel.Q.Value' in df.columns:
            before = len(df)
            df = df.filter(pl.col('Channel.Q.Value') <= 0.01)
            print(f"  [jmod] Channel.Q.Value (Qvalue) ≤ 0.01: {before:,} → {len(df):,}")
            filters.append('Channel Q ≤ 0.01')
    else:
        if 'Q.Value' in df.columns:
            before = len(df)
            df = df.filter(pl.col('Q.Value') <= 0.01)
            print(f"  [jmod] Q.Value (BestChannel_Qvalue) ≤ 0.01: {before:,} → {len(df):,}")
            filters.append('Q.Value ≤ 0.01')

    return df, filters


def _add_derived(df: pl.DataFrame) -> pl.DataFrame:
    """Add seqcharge, C_term_aa, N_term_aa, and labeling columns."""
    exprs = []

    # seqcharge = tag-free Precursor.Id (channel-independent, stable across runs)
    if 'Precursor.Id' in df.columns:
        exprs.append(pl.col('Precursor.Id').alias('seqcharge'))

    if 'Sequence' in df.columns:
        exprs.append(pl.col('Sequence').str.slice(-1).alias('C_term_aa'))
        exprs.append(pl.col('Sequence').str.slice(0, 1).alias('N_term_aa'))

    if exprs:
        df = df.with_columns(exprs)

    # Labeling columns — Jmod tag notation: X(PSMtag_<label>-<channel>)PEPTIDE
    if 'Modified sequence' in df.columns and 'Sequence' in df.columns:
        mod_seq = pl.col('Modified sequence')

        # N-terminal alpha-amine: first residue carries a (PSMtag...) modification
        has_nterm_label = mod_seq.str.contains(r'^[A-Z]\(PSMtag', literal=False)

        # K side-chain epsilon-amines: count K residues with PSMtag modifications
        n_k_labeled = mod_seq.str.count_matches(r'K\(PSMtag')

        # Total K residues in stripped sequence (each is a potential labeling site)
        n_k_total = pl.col('Sequence').str.count_matches('K')

        # Expected labels = 1 (N-term) + all K side chains; actual = N-term + labeled K
        n_expected = 1 + n_k_total
        n_actual   = has_nterm_label.cast(pl.Int32) + n_k_labeled

        df = df.with_columns([
            (n_actual >= n_expected).cast(pl.Int8).alias('fully_labeled'),
            ((n_actual > 0) & (n_actual < n_expected)).cast(pl.Int8).alias('partly_labeled'),
        ])

    return df


def _collapse_channels(df: pl.DataFrame) -> pl.DataFrame:
    """
    Collapse per-channel rows to one representative row per (Raw file, Precursor.Id).
    Only applied to runs that have more than one distinct channel.
    """
    if 'Channel' not in df.columns or 'Precursor.Id' not in df.columns:
        return df

    ch_counts = (df.group_by('Raw file')
                   .agg(pl.col('Channel').n_unique().alias('n_ch')))
    multi_ch_runs = ch_counts.filter(pl.col('n_ch') > 1)['Raw file'].to_list()

    if not multi_ch_runs:
        return df

    is_mux = pl.col('Raw file').is_in(multi_ch_runs)
    mux     = df.filter(is_mux)
    non_mux = df.filter(~is_mux)

    group_keys = ['Raw file', 'Precursor.Id']
    avg_cols   = [c for c in ['Intensity', 'Q.Value', 'Channel.Q.Value', 'plex_Area',
                               'Precursor.Quantity']
                  if c in mux.columns]
    first_cols = [c for c in mux.columns if c not in group_keys + avg_cols]

    collapsed = (
        mux.group_by(group_keys)
           .agg(
               [pl.col(c).mean() for c in avg_cols] +
               [pl.col(c).first() for c in first_cols]
           )
           .select(df.columns)
    )
    before = mux.height
    print(f"  [jmod] channel collapse: {before:,} rows → {collapsed.height:,} "
          f"(÷{before // max(collapsed.height, 1)} channels avg)")
    return pl.concat([non_mux, collapsed], how='diagonal')


# ---------------------------------------------------------------------------
# Public processor class
# ---------------------------------------------------------------------------

class JmodProcessor:
    """
    Normalises and enriches Jmod output loaded by data_loader.load_folder().

    Parameters
    ----------
    data : dict[str, pl.DataFrame]
        Output of load_folder() when engine == 'jmod'.
        Expected key: 'ids' (required) — all_IDs_filtered.parquet or *_IDs.parquet.
    """

    def __init__(self, data: Dict[str, pl.DataFrame]):
        if 'ids' not in data or data['ids'] is None:
            raise ValueError("JmodProcessor requires an 'ids' DataFrame in data dict.")

        self.report:        pl.DataFrame           = data['ids']
        self.features:      Optional[pl.DataFrame] = None
        self.fill_times:    Optional[pl.DataFrame] = None
        self.tic:           Optional[pl.DataFrame] = None
        self.sn:            Optional[pl.DataFrame] = None
        self.ms1_extracted: Optional[pl.DataFrame] = None
        self.enzymes:       Optional[List[str]]    = None
        self.filters:       List[str]              = []
        self.raw_report:    Optional[pl.DataFrame] = None
        self.plex_report:   Optional[pl.DataFrame] = None
        self.has_full_data: bool                   = False  # True when loaded from all_IDs.csv

    # ------------------------------------------------------------------
    # Main pipeline
    # ------------------------------------------------------------------

    def normalize(self) -> 'JmodProcessor':
        """
        Run the full normalization pipeline in place:
          1. Extract run names from file_name paths
          2. Cast channel → String 'Channel', z → Int32
          3. Remove decoy rows
          4. Rename columns to canonical names
          5. Add derived columns (seqcharge, terminal AAs, labeling flags)
          6. Save raw_report (pre-quality-filter, channel-collapsed)
          7. Apply quality filter (Qvalue ≤ 0.01 or BestChannel_Qvalue ≤ 0.01)
          8. For multiplexed: build plex_report and collapse channels in report
        Returns self for chaining.
        """
        df = self.report

        df = _extract_run_names(df)
        df = _cast_columns(df)
        df = _filter_decoys(df)
        df = _rename(df, _COL_MAP)
        df, self.has_full_data = _assign_intensity(df)
        if self.has_full_data:
            print(f"  [jmod] full data (all_IDs.csv) — MS1_Area, rt_error, ms1_cor, tic available")
        df = _add_derived(df)

        is_mux = _is_multiplexed(df)
        if is_mux:
            n_ch = df['Channel'].n_unique()
            print(f"  [jmod] multiplexed data — {n_ch} channel(s) detected")

        # raw_report: decoy-filtered, channel-collapsed, pre-quality-filter (for Summary)
        self.raw_report = _collapse_channels(df)

        df, self.filters = _filter_quality(df, is_mux)

        if is_mux:
            ch_counts = (df.group_by('Raw file')
                           .agg(pl.col('Channel').n_unique().alias('n_ch')))
            plex_runs = ch_counts.filter(pl.col('n_ch') > 1)['Raw file'].to_list()
            if plex_runs:
                plex_df = df.filter(pl.col('Raw file').is_in(plex_runs))
                # plex_report uses channel-level plex_Area as Intensity (not precursor-level MS1_Area)
                if 'plex_Area' in plex_df.columns:
                    plex_df = plex_df.with_columns(pl.col('plex_Area').alias('Intensity'))
                self.plex_report = plex_df
            df = _collapse_channels(df)
            # plex_Area is meaningless after channel collapse — drop it
            if 'plex_Area' in df.columns:
                df = df.drop('plex_Area')

        # plex_Area is only useful in plex_report (channel-level); drop from main report
        if 'plex_Area' in df.columns:
            df = df.drop('plex_Area')

        self.report = df
        print(f"  [jmod] normalized report: {len(self.report):,} rows × {self.report.width} columns")
        return self

    def filter_enzyme(self, enzymes: List[str]) -> 'JmodProcessor':
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

    def summary(self) -> Dict:
        df = self.report
        return {
            'rows':       len(df),
            'runs':       df['Raw file'].n_unique()      if 'Raw file'      in df.columns else None,
            'precursors': df['Precursor.Id'].n_unique()  if 'Precursor.Id'  in df.columns else None,
            'seqcharges': df['seqcharge'].n_unique()     if 'seqcharge'     in df.columns else None,
            'proteins':   df['Protein.Group'].n_unique() if 'Protein.Group' in df.columns else None,
        }
