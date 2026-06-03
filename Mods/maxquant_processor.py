# Created by @thibaultcolombani with the help of Claude Code (Anthropic)
"""
maxquant_processor.py - MaxQuant data processor for MassSleuth.

Takes the dict produced by load_folder() for engine='maxquant' and normalises
the evidence DataFrame to the canonical column schema, adds derived columns
for QC and labeling-efficiency analysis.

Input evidence.txt columns used
---------------------------------
Raw file          run name (already canonical)
Sequence          stripped peptide sequence (already canonical)
Modified sequence modified sequence with MaxQuant notation, e.g. _PEPTIDEK(tag)R_
Charge            precursor charge state (already canonical)
Retention time    apex RT in minutes (already canonical)
Retention length  peak FWHM in minutes (already canonical)
Intensity         MS1 feature intensity (already canonical)
PEP               posterior error probability (already canonical)
Proteins          protein group accessions → Protein.Group
m/z               precursor m/z → Precursor.Mz
Reverse           '+' for decoy entries (used for filtering, then dropped)
Potential contaminant '+' for contaminant entries (used for filtering, then dropped)

Canonical columns produced in report
--------------------------------------
Raw file          (unchanged)
Sequence          (unchanged)
Modified sequence (unchanged)
Charge            (unchanged)
Retention time    (unchanged, minutes)
Retention length  (unchanged, minutes)
Intensity         (unchanged)
PEP               (unchanged)
Protein.Group     ← Proteins
Precursor.Mz      ← m/z

Derived columns
---------------
seqcharge      "ModifiedSequence_Charge"
Precursor.Id   same as seqcharge (no native precursor ID in MaxQuant)
C_term_aa      last residue of Sequence
N_term_aa      first residue of Sequence
fully_labeled  1 if all expected amine sites carry a modification
partly_labeled 1 if at least one amine site is labeled but not all

Quality filters applied during normalize()
------------------------------------------
Reverse == '+' rows removed (decoys)
Potential contaminant == '+' rows removed
PEP ≤ 0.01
"""

from __future__ import annotations

import re
from typing import Dict, List, Optional, Tuple

import polars as pl


# ---------------------------------------------------------------------------
# Column rename map  (MaxQuant raw name → canonical name)
# ---------------------------------------------------------------------------

_COL_MAP: Dict[str, str] = {
    'Proteins': 'Protein.Group',
    'm/z':      'Precursor.Mz',
}


# ---------------------------------------------------------------------------
# Module-level transformation helpers
# ---------------------------------------------------------------------------

def _rename(df: pl.DataFrame, col_map: Dict[str, str]) -> pl.DataFrame:
    rmap = {k: v for k, v in col_map.items() if k in df.columns}
    return df.rename(rmap) if rmap else df


def _filter_decoys(df: pl.DataFrame) -> pl.DataFrame:
    if 'Reverse' not in df.columns:
        return df
    before = len(df)
    df = df.filter(pl.col('Reverse').is_null() | (pl.col('Reverse') != '+'))
    removed = before - len(df)
    if removed:
        print(f"  [maxquant] removed {removed:,} decoy rows → {len(df):,} remaining")
    return df.drop('Reverse') if 'Reverse' in df.columns else df


def _filter_contaminants(df: pl.DataFrame) -> pl.DataFrame:
    col = 'Potential contaminant'
    if col not in df.columns:
        return df
    before = len(df)
    df = df.filter(pl.col(col).is_null() | (pl.col(col) != '+'))
    removed = before - len(df)
    if removed:
        print(f"  [maxquant] removed {removed:,} contaminant rows → {len(df):,} remaining")
    return df.drop(col) if col in df.columns else df


def _filter_quality(df: pl.DataFrame) -> Tuple[pl.DataFrame, List[str]]:
    filters: List[str] = []
    if 'PEP' in df.columns:
        before = len(df)
        df = df.filter(pl.col('PEP') <= 0.01)
        print(f"  [maxquant] PEP ≤ 0.01: {before:,} → {len(df):,}")
        filters.append('PEP ≤ 0.01')
    return df, filters


def _add_derived(df: pl.DataFrame) -> pl.DataFrame:
    """Add seqcharge, Precursor.Id, C_term_aa, N_term_aa, and labeling columns."""
    exprs = []

    # seqcharge and Precursor.Id from Modified sequence + Charge
    if 'Modified sequence' in df.columns and 'Charge' in df.columns:
        sc = (pl.col('Modified sequence') + '_' + pl.col('Charge').cast(pl.String))
        exprs.append(sc.alias('seqcharge'))
        exprs.append(sc.alias('Precursor.Id'))

    # Terminal amino acids from stripped Sequence
    if 'Sequence' in df.columns:
        exprs.append(pl.col('Sequence').str.slice(-1).alias('C_term_aa'))
        exprs.append(pl.col('Sequence').str.slice(0, 1).alias('N_term_aa'))

    if exprs:
        df = df.with_columns(exprs)

    # Labeling columns — use MaxQuant's pre-computed Variable_ count columns when present.
    # These are tag-name agnostic: any variable amine-reactive mod (PSMtag, mTRAQ, iTRAQ…)
    # produces columns like "Variable_<any name> Lys" / "Variable_<any name> N-term".
    # Fixed mods (TMT bulk) are not visible here; handled post-hoc in normalize() via plex_report.
    if 'Sequence' in df.columns:
        lys_col   = next((c for c in df.columns
                          if c.startswith('Variable_') and c.endswith(' Lys')), None)
        nterm_col = next((c for c in df.columns
                          if c.startswith('Variable_') and c.endswith(' N-term')), None)

        if lys_col is not None or nterm_col is not None:
            n_k_labeled     = pl.col(lys_col).fill_null(0)   if lys_col   else pl.lit(0).cast(pl.Int32)
            n_nterm_labeled = pl.col(nterm_col).fill_null(0) if nterm_col else pl.lit(0).cast(pl.Int32)
            n_k_total       = pl.col('Sequence').str.count_matches('K')
            acetyl = (pl.col('Acetyl (Protein N-term)').fill_null(0)
                      if 'Acetyl (Protein N-term)' in df.columns
                      else pl.lit(0).cast(pl.Int32))
            n_expected = (1 - acetyl) + n_k_total
            n_actual   = n_nterm_labeled + n_k_labeled
            df = df.with_columns([
                (n_actual >= n_expected).cast(pl.Int8).alias('fully_labeled'),
                ((n_actual > 0) & (n_actual < n_expected)).cast(pl.Int8).alias('partly_labeled'),
            ])

    return df


# ---------------------------------------------------------------------------
# TMT plex report builder
# ---------------------------------------------------------------------------

_REP_PAT = re.compile(r'^Reporter intensity corrected (\d+)$')


def _build_plex_report(df: pl.DataFrame) -> Optional[pl.DataFrame]:
    """
    Pivot 'Reporter intensity corrected N' wide columns to long format.
    Returns None if no reporter columns are present.
    Active channels: ≥ 5% of PSMs have non-zero signal (minimum 10 rows).
    """
    rep_cols = {m.group(1): c for c in df.columns if (m := _REP_PAT.match(c))}
    if not rep_cols:
        return None

    n = len(df)
    threshold = max(10, int(n * 0.05))
    active = {idx: col for idx, col in rep_cols.items()
              if df.filter(pl.col(col) > 0).height >= threshold}
    if not active:
        return None

    print(f"  [maxquant] TMT: {len(active)} active reporter channel(s) detected")

    carry = [c for c in ('Raw file', 'Precursor.Id', 'seqcharge',
                          'Sequence', 'Modified sequence', 'Charge',
                          'Protein.Group', 'PEP', 'Retention time', 'Precursor.Mz',
                          'fully_labeled', 'partly_labeled')
             if c in df.columns]

    frames = []
    for idx, col in sorted(active.items(), key=lambda x: int(x[0])):
        ch_df = (df.select(carry + [col])
                   .rename({col: 'Intensity'})
                   .with_columns(pl.col('Intensity').alias('Precursor.Quantity'))
                   .with_columns(pl.lit(idx).alias('Channel'))
                   .filter(pl.col('Intensity') > 0))
        frames.append(ch_df)

    if not frames:
        return None

    result = pl.concat(frames, how='diagonal')
    print(f"  [maxquant] plex_report: {len(result):,} rows  "
          f"({len(active)} channels × {n:,} PSMs)")
    return result


# ---------------------------------------------------------------------------
# Public processor class
# ---------------------------------------------------------------------------

class MaxquantProcessor:
    """
    Normalises and enriches MaxQuant output loaded by data_loader.load_folder().

    Parameters
    ----------
    data : dict[str, pl.DataFrame]
        Output of load_folder() when engine == 'maxquant'.
        Expected key: 'evidence' (required). Optional: 'msms', 'msmsScans',
        'allPeptides', 'summary', 'parameters', 'msScans'.
    """

    def __init__(self, data: Dict[str, pl.DataFrame]):
        if 'evidence' not in data or data['evidence'] is None:
            raise ValueError("MaxquantProcessor requires an 'evidence' DataFrame in data dict.")

        self.report:        pl.DataFrame           = data['evidence']
        self.msms:          Optional[pl.DataFrame] = data.get('msms')
        self.msms_scans:    Optional[pl.DataFrame] = data.get('msmsScans')
        self.all_peptides:  Optional[pl.DataFrame] = data.get('allPeptides')
        self.summary:       Optional[pl.DataFrame] = data.get('summary')
        self.parameters:    Optional[pl.DataFrame] = data.get('parameters')
        self.ms_scans:      Optional[pl.DataFrame] = data.get('msScans')

        # Attributes mirroring DiannProcessor / JmodProcessor interface
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

    def normalize(self) -> 'MaxquantProcessor':
        """
        Run the full normalization pipeline in place:
          1. Rename Proteins → Protein.Group, m/z → Precursor.Mz
          2. Remove decoy rows (Reverse == '+')
          3. Remove contaminant rows (Potential contaminant == '+')
          4. Add derived columns (seqcharge, Precursor.Id, terminal AAs, labeling flags)
          5. Save raw_report (pre-quality-filter, for Summary)
          6. Apply quality filter (PEP ≤ 0.01)
        Returns self for chaining.
        """
        df = self.report
        df = _rename(df, _COL_MAP)
        df = _filter_decoys(df)
        df = _filter_contaminants(df)
        df = _add_derived(df)

        self.raw_report = df  # pre-quality-filter snapshot for Summary tab

        df, self.filters = _filter_quality(df)

        # Mean of all reporter channels → Precursor.Quantity (MS2 proxy) for main report
        rep_cols_all = [c for c in df.columns if _REP_PAT.match(c)]
        if rep_cols_all:
            n_rep = len(rep_cols_all)
            df = df.with_columns(
                (pl.sum_horizontal(
                    [pl.col(c).cast(pl.Float64, strict=False).fill_null(0.0) for c in rep_cols_all]
                ) / n_rep).alias('Precursor.Quantity')
            )

        # Build long-format plex_report from TMT reporter intensity columns
        self.plex_report = _build_plex_report(df)

        # Fixed-TMT fallback: reporter ions detected but no Variable_ labeling columns
        # exist → all passing PSMs are fully labeled (MaxQuant required the fixed mod).
        if self.plex_report is not None:
            has_var_lys = any(c.startswith('Variable_') and c.endswith(' Lys')
                              for c in df.columns)
            if not has_var_lys:
                df = df.with_columns([
                    pl.lit(1).cast(pl.Int8).alias('fully_labeled'),
                    pl.lit(0).cast(pl.Int8).alias('partly_labeled'),
                ])

        # Drop the wide reporter intensity columns from the main report
        rep_drop = [c for c in df.columns if c.startswith('Reporter intensity')]
        if rep_drop:
            df = df.drop(rep_drop)

        self.report = df
        print(f"  [maxquant] normalized report: {len(self.report):,} rows × {self.report.width} columns")
        return self

    def filter_enzyme(self, enzymes: List[str]) -> 'MaxquantProcessor':
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
