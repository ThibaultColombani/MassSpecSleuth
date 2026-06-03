# Created by @thibaultcolombani with the help of Claude Code (Anthropic)
"""
data_loader.py - Multi-engine file discovery, detection, and loading for MassSleuth.

Given a folder path, scans for search engine output files, detects the engine,
and loads ALL relevant files for that engine into a dict of named Polars DataFrames.


Supported engines: DIA-NN, MaxQuant, Sage, FragPipe, ProteomeDiscoverer, Jmod

Usage
-----
    from Mods.data_loader import load_folder

    engine, data = load_folder("/path/to/results")
    # engine  → 'diann' | 'maxquant' | 'sage' | 'fragpipe' | 'protdisc' | 'jmod' | 'unknown'
    # data    → dict[str, pl.DataFrame]  — only keys for files that were found

    # DIA-NN keys:       'report', 'ms1_extracted', 'features', 'fill_times', 'tic', 'sn'
    # MaxQuant keys:     'evidence', 'msms', 'msmsScans', 'allPeptides', 'summary',
    #                    'parameters', 'msScans'
    # Sage keys:         'psms', 'lfq'
    # FragPipe keys:     'psms'
    # ProteomeDisc keys: 'psms', 'protein_groups', 'proteins', 'input_files'
    # Jmod keys:         'ids'
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import polars as pl


# ---------------------------------------------------------------------------
# File catalog definitions (one entry per logical file slot per engine)
# ---------------------------------------------------------------------------

@dataclass
class FileSpec:
    key: str               # logical name used as dict key in the output
    filenames: List[str]   # candidate filenames or glob patterns (checked in order)
    required: bool = False # if True, missing file triggers a warning


_ENGINE_CATALOG: Dict[str, List[FileSpec]] = {

    'diann': [
        # Primary report — load all *report*.parquet files (concatenated), fall back to tsv
        FileSpec('report',        ['*report*.parquet', 'report.tsv'],                      required=True),
        # Supplementary instrument files (present when DIA-NN writes them)
        FileSpec('ms1_extracted', ['report.pr_matrix_channels_ms1_extracted.tsv'],         required=False),
        FileSpec('features',      ['features.tsv'],                                        required=False),
        FileSpec('fill_times',    ['fill_times.tsv'],                                      required=False),
        FileSpec('tic',           ['tic.tsv'],                                             required=False),
        FileSpec('sn',            ['sn.tsv'],                                              required=False),
    ],

    'maxquant': [
        FileSpec('evidence',    ['evidence.txt'],    required=True),
        FileSpec('msms',        ['msms.txt'],        required=False),
        FileSpec('msmsScans',   ['msmsScans.txt'],   required=False),
        FileSpec('allPeptides', ['allPeptides.txt'], required=False),
        FileSpec('summary',     ['summary.txt'],     required=False),
        FileSpec('parameters',  ['parameters.txt'],  required=False),
        FileSpec('msScans',     ['msScans.txt'],     required=False),
    ],

    'sage': [
        # Each raw file produces its own *.sage.tsv; they are concatenated
        FileSpec('psms', ['*.sage.tsv'], required=True),
        # LFQ table: wide format, one MS1-intensity column per file
        FileSpec('lfq',  ['lfq.tsv'],   required=False),
    ],

    'fragpipe': [
        FileSpec('psms', ['psm.tsv'], required=True),
    ],

    'protdisc': [
        FileSpec('psms',           ['*PSMs*.txt'],          required=True),
        FileSpec('protein_groups', ['*ProteinGroups*.txt'], required=False),
        FileSpec('proteins',       ['*Proteins*.txt'],      required=False),
        FileSpec('input_files',    ['*InputFiles*.txt'],    required=False),
    ],

    'jmod': [
        # Preferred: all_IDs.csv (full 25-col subset, rich QC data)
        # Fallback:  all_IDs_filtered.parquet (17 cols, minimal)
        # Legacy:    *_IDs.parquet (old single-mode format)
        FileSpec('ids', ['all_IDs.csv', 'all_IDs_filtered.parquet', '*_IDs.parquet'], required=True),
    ],
}


# ---------------------------------------------------------------------------
# Engine detection — priority order: most specific signature first
# ---------------------------------------------------------------------------

_DETECTION_SIGNATURES = [
    ('sage',     lambda name: name.endswith('.sage.tsv')),
    ('fragpipe', lambda name: name == 'psm.tsv'),
    ('jmod',     lambda name: name.endswith('_IDs.parquet') or name in ('all_IDs_filtered.parquet', 'all_IDs.csv')),
    ('diann',    lambda name: ('report' in name and name.endswith('.parquet')) or name == 'report.tsv'),
    ('maxquant', lambda name: name == 'evidence.txt'),
    ('protdisc', lambda name: name.endswith('.txt') and
                              any(kw in name for kw in ('PSMs', 'InputFiles', 'ProteinGroups', 'Proteins'))),
]


def _list_files(root: Path) -> List[Path]:
    """All non-hidden files under root (recursive)."""
    if root.is_file():
        return [root]
    return [p for p in root.rglob('*') if p.is_file() and not p.name.startswith('.')]


def _detect_engine(files: List[Path]) -> Tuple[str, List[Path]]:
    """
    Identify the search engine from a list of candidate paths.

    Returns (engine_name, matched_primary_files).
    Returns ('unknown', []) if nothing matches.
    If multiple engines are detected, warns and picks the highest-priority one.
    """
    all_matched = [(engine, [f for f in files if sig(f.name)])
                   for engine, sig in _DETECTION_SIGNATURES]
    all_matched = [(e, m) for e, m in all_matched if m]

    if not all_matched:
        return 'unknown', []

    if len(all_matched) > 1:
        detected = ', '.join(e for e, _ in all_matched)
        print(f"  Warning: multiple search engines detected ({detected})"
              f" — using {all_matched[0][0].upper()}")

    return all_matched[0]


# ---------------------------------------------------------------------------
# File loading helpers (Polars)
# ---------------------------------------------------------------------------

# Columns to select when loading Jmod all_IDs.csv (641 MB, 139 cols → load only what's needed)
_JMOD_ALL_IDS_COLS = [
    # Core — same as all_IDs_filtered.parquet
    'stripped_seq', 'z', 'untag_prec', 'file_name', 'channel', 'is_decoy',
    'Qvalue', 'Protein_Qvalue', 'PredVal', 'protein', 'BestChannel_Qvalue',
    'plex_Area', 'seq', 'silac_channel', 'untag_seq', 'rt', 'mz',
    # Extended — only in all_IDs.csv
    'coeff',                                                          # MS2 quant (→ Precursor.Quantity)
    'MS1_Area', 'MS1_Int', 'tic', 'ms1_cor', 'rt_error', 'mz_error', 'iso_cor', 'cosine',
]


def _load_parquet(path: Path) -> Optional[pl.DataFrame]:
    try:
        return pl.read_parquet(path)
    except Exception as e:
        print(f"  Warning: could not load {path.name}: {e}")
        return None


def _load_tsv(path: Path, ignore_errors: bool = True) -> Optional[pl.DataFrame]:
    try:
        return pl.read_csv(path, separator='	', infer_schema_length=10_000,
                           ignore_errors=ignore_errors)
    except Exception as e:
        print(f"  Warning: could not load {path.name}: {e}")
        return None


def _load_jmod_csv(path: Path) -> Optional[pl.DataFrame]:
    """Load Jmod all_IDs.csv selecting only the QC-relevant columns."""
    try:
        # Determine which of the desired columns are actually present
        header = pl.read_csv(path, n_rows=0).columns
        cols = [c for c in _JMOD_ALL_IDS_COLS if c in header]
        missing = [c for c in _JMOD_ALL_IDS_COLS if c not in header]
        if missing:
            print(f"    Note: {path.name} missing expected columns: {missing}")
        return pl.read_csv(path, columns=cols, infer_schema_length=10_000,
                           ignore_errors=True)
    except Exception as e:
        print(f"  Warning: could not load {path.name}: {e}")
        return None


def _load_file(path: Path) -> Optional[pl.DataFrame]:
    if path.suffix == '.parquet':
        return _load_parquet(path)
    if path.name == 'all_IDs.csv':
        return _load_jmod_csv(path)
    return _load_tsv(path)


def _concat(frames: List[pl.DataFrame], key: str) -> pl.DataFrame:
    if len(frames) == 1:
        return frames[0]
    combined = pl.concat(frames, how='diagonal_relaxed')
    print(f"    [{key}] Concatenated {len(frames)} files → {len(combined):,} rows")
    return combined


# ---------------------------------------------------------------------------
# Catalog-based file loading
# ---------------------------------------------------------------------------

def _resolve_spec(spec: FileSpec, root: Path) -> List[Path]:
    """
    Find files matching any candidate in spec.filenames anywhere under root
    (recursive).  Supports exact names and glob patterns (e.g. '*.sage.tsv').
    Returns the first pattern that produces at least one match; later patterns
    in the list are not tried (allows 'preferred' → 'fallback' ordering).
    """
    found: List[Path] = []
    for pattern in spec.filenames:
        found.extend(sorted(root.rglob(pattern)))
        if found:
            break   # first matching pattern wins
    return found


def _load_catalog(engine: str, root: Path) -> Dict[str, pl.DataFrame]:
    """
    Load all catalog files for `engine` by recursively searching under root.
    Multiple matching files with the same key are concatenated.

    Returns a dict where absent/empty files are simply not included.
    """
    specs = _ENGINE_CATALOG.get(engine, [])
    per_key: Dict[str, List[pl.DataFrame]] = {s.key: [] for s in specs}

    for spec in specs:
        paths = _resolve_spec(spec, root)

        if not paths:
            if spec.required:
                print(f"  Warning: required file '{spec.filenames[0]}' not found under {root.name}")
            continue

        for path in paths:
            df = _load_file(path)
            if df is not None and len(df) > 0:
                df = df.with_columns(pl.lit(str(path)).alias('_source_file'))
                per_key[spec.key].append(df)
                print(f"    {spec.key}: {path.name}  ({len(df):,} rows)")

    result: Dict[str, pl.DataFrame] = {}
    for spec in specs:
        frames = per_key[spec.key]
        if frames:
            result[spec.key] = _concat(frames, spec.key)

    return result


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def load_folder(path: str) -> Tuple[str, Dict[str, pl.DataFrame]]:
    """
    Scan a folder, detect the search engine, and load all relevant files.

    Parameters
    ----------
    path : str
        Path to a folder (or single file) containing search engine output.

    Returns
    -------
    engine : str
        One of 'diann', 'maxquant', 'sage', 'fragpipe', 'protdisc', 'jmod', 'unknown'.
    data : dict[str, pl.DataFrame]
        Named DataFrames — keys present only for files that were found and loaded.
        Each DataFrame has an extra '_source_file' column.
    """
    root = Path(path)
    if not root.exists():
        raise FileNotFoundError(f"Path does not exist: {path}")

    print(f"\nScanning: {root}")
    all_files = _list_files(root)
    print(f"  {len(all_files)} file(s) found")

    engine, primary_files = _detect_engine(all_files)

    if engine == 'unknown':
        print("  No recognised search engine output found.")
        return 'unknown', {}

    print(f"  Engine: {engine.upper()}")

    search_root = root.parent.resolve() if root.is_file() else root.resolve()
    data = _load_catalog(engine, search_root)

    if not data:
        print("  No data could be loaded.")
    else:
        loaded_keys = ', '.join(
            f"{k} ({len(df):,} rows)" for k, df in data.items()
        )
        print(f"  Loaded: {loaded_keys}")

    return engine, data
