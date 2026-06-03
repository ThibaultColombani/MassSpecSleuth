# Created by @thibaultcolombani with the help of Claude Code (Anthropic)
"""
combined_processor.py — Thin adapter that merges normalized processors into one.

Use when loading data from multiple folders or multiple search engines.
All input processors must already have .normalize() called.
"""
from __future__ import annotations

from typing import Dict, List, Optional

import polars as pl

_INT_TYPES  = {pl.Int8, pl.Int16, pl.Int32, pl.Int64}
_FLOAT_TYPES = {pl.Float32, pl.Float64}

_ENGINE_LABELS = {
    'DiannProcessor':     'DIA-NN',
    'JmodProcessor':      'JMOD',
    'MaxquantProcessor':  'MQ',
    'SageProcessor':      'SAGE',
    'FragpipeProcessor':  'FP',
    'ProtdiscProcessor':  'PD',
}


def _engine_label(p) -> str:
    return _ENGINE_LABELS.get(type(p).__name__, type(p).__name__.replace('Processor', '').upper())


def _tag_runs(df: Optional[pl.DataFrame], suffix: str) -> Optional[pl.DataFrame]:
    """Append ' [suffix]' to Raw file values in df."""
    if df is None or 'Raw file' not in df.columns:
        return df
    return df.with_columns(
        (pl.col('Raw file') + f' [{suffix}]').alias('Raw file')
    )


def _harmonize(frames: list) -> list:
    """Cast shared columns with mismatched integer/float types to the widest type."""
    if len(frames) <= 1:
        return frames

    # Collect all types per column across all frames
    col_types: dict = {}
    for df in frames:
        for col, dtype in zip(df.columns, df.dtypes):
            col_types.setdefault(col, set()).add(dtype)

    # Determine casts needed
    result = []
    for df in frames:
        casts = []
        for col in df.columns:
            types = col_types[col]
            if len(types) == 1:
                continue
            if types <= _INT_TYPES:
                casts.append(pl.col(col).cast(pl.Int64))
            elif types <= (_INT_TYPES | _FLOAT_TYPES):
                casts.append(pl.col(col).cast(pl.Float64))
        result.append(df.with_columns(casts) if casts else df)
    return result


class CombinedProcessor:
    """
    Merges one or more normalized processor objects (DiannProcessor, JmodProcessor, …)
    into a single object with the same interface expected by the HTML exporters.

    Parameters
    ----------
    processors : list
        Normalized processor objects (any mix of engines).
    """

    def __init__(self, processors: list):
        if not processors:
            raise ValueError("CombinedProcessor requires at least one processor.")

        # Tag runs with engine labels whenever multiple different engines are combined
        unique_labels = list(dict.fromkeys(_engine_label(p) for p in processors))
        should_tag = len(unique_labels) > 1

        if should_tag:
            labels = [_engine_label(p) for p in processors]
            print(f"  [combined] multi-engine report — tagging runs with engine labels: "
                  f"{', '.join(unique_labels)}")
            reports   = [_tag_runs(p.report,     lb) for p, lb in zip(processors, labels)]
            raws      = [_tag_runs(getattr(p, 'raw_report',  None), lb) for p, lb in zip(processors, labels)]
            plexes    = [_tag_runs(getattr(p, 'plex_report', None), lb) for p, lb in zip(processors, labels)]
        else:
            reports   = [p.report for p in processors]
            raws      = [getattr(p, 'raw_report',  None) for p in processors]
            plexes    = [getattr(p, 'plex_report', None) for p in processors]

        self.report = pl.concat(_harmonize(reports), how='diagonal')

        raw_frames = [df for df in raws if df is not None]
        self.raw_report: Optional[pl.DataFrame] = (
            pl.concat(_harmonize(raw_frames), how='diagonal') if raw_frames else None
        )

        plex_frames = [df for df in plexes if df is not None]
        self.plex_report: Optional[pl.DataFrame] = (
            pl.concat(_harmonize(plex_frames), how='diagonal') if plex_frames else None
        )

        # Merge filter strings, deduplicated and in order
        seen: dict = {}
        for p in processors:
            for f in getattr(p, 'filters', []):
                seen[f] = None
        self.filters: List[str] = list(seen)

        # Merge enzyme lists from all processors (deduplicated, order-preserving)
        seen_enz: dict = {}
        for p in processors:
            for e in (getattr(p, 'enzymes', None) or []):
                seen_enz[e] = None
        self.enzymes: Optional[List[str]] = list(seen_enz) or None
        self.has_full_data: bool = any(
            getattr(p, 'has_full_data', False) for p in processors
        )

        # Optional DIA-NN supplementary frames (pass through if present)
        self.features      = next((p.features      for p in processors if getattr(p, 'features',      None) is not None), None)
        self.fill_times    = next((p.fill_times     for p in processors if getattr(p, 'fill_times',    None) is not None), None)
        self.tic           = next((p.tic            for p in processors if getattr(p, 'tic',           None) is not None), None)
        self.sn            = next((p.sn             for p in processors if getattr(p, 'sn',            None) is not None), None)
        self.ms1_extracted = next((p.ms1_extracted  for p in processors if getattr(p, 'ms1_extracted', None) is not None), None)

    # ------------------------------------------------------------------
    # Interface required by HTML exporters
    # ------------------------------------------------------------------

    def runs(self) -> List[str]:
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
