# Created by @thibaultcolombani with the help of Claude Code (Anthropic)
"""
enzyme.py — Protease cleavage rules and peptide specificity filtering.

Usage in any processor:
    from Mods.enzyme import filter_peptides, ENZYMES
    df = filter_peptides(df, ['trypsin'])          # single enzyme
    df = filter_peptides(df, ['trypsin', 'glu-c']) # combination
"""
from __future__ import annotations

from typing import Dict, List, Set

import polars as pl

# ---------------------------------------------------------------------------
# Enzyme catalogue
# ---------------------------------------------------------------------------
# 'cterm': expected C-terminal residues of fully specific peptides
# 'nterm': expected N-terminal residues (for N-terminal cutting enzymes)
# When multiple enzymes are combined, the union of each side is used and a
# peptide is kept if it matches either side's residue set.
ENZYMES: Dict[str, Dict[str, Set[str]]] = {
    # ── C-terminal cutters ────────────────────────────────────────────────
    'trypsin':        {'cterm': {'R', 'K'}},
    'trypsin/p':      {'cterm': {'R', 'K'}},     # allows cleavage before P
    'lys-c':          {'cterm': {'K'}},
    'lys-c/p':        {'cterm': {'K'}},
    'arg-c':          {'cterm': {'R'}},
    'glu-c':          {'cterm': {'E', 'D'}},     # phosphate buffer (E + D)
    'glu-c-bic':      {'cterm': {'E'}},          # ammonium bicarbonate (E only)
    'v8':             {'cterm': {'E', 'D'}},     # alias for glu-c phosphate
    'chymotrypsin':   {'cterm': {'F', 'Y', 'W', 'L', 'M'}},
    'chymotrypsin/p': {'cterm': {'F', 'Y', 'W', 'L', 'M'}},
    'clostripain':    {'cterm': {'R'}},           # highly Arg-specific, like Arg-C
    'cnbr':           {'cterm': {'M'}},           # chemical, cleaves after Met
    'bnps-skatole':   {'cterm': {'W'}},           # chemical, cleaves after Trp
    'ntcb':           {'cterm': {'C'}},           # 2-nitro-5-thiocyanatobenzoic acid
    'thermolysin':    {'cterm': {'A', 'F', 'I', 'L', 'M', 'V'}},  # hydrophobic bulky
    'pepsin-ph2':     {'cterm': {'F', 'L'}},
    'pepsin-ph13':    {'cterm': {'F', 'L', 'W', 'Y', 'A', 'E', 'Q'}},
    'elastase':       {'cterm': {'A', 'G', 'S', 'V'}},
    'proteinase-k':   {'cterm': set('ACDEFGHIKLMNPQRSTVWY')},  # near-non-specific
    # ── N-terminal cutters ────────────────────────────────────────────────
    'asp-n':          {'nterm': {'D', 'C'}},     # cleaves N-terminal to Asp (and Cys)
    'asp-n/d':        {'nterm': {'D'}},           # strict Asp-N, Asp only
    'lys-n':          {'nterm': {'K'}},           # cleaves N-terminal to Lys
    'formic-acid':    {'nterm': {'D'}},           # chemical, cleaves Asp-Pro bonds
}


def enzyme_names() -> List[str]:
    """Sorted list of all known enzyme names."""
    return sorted(ENZYMES)


def filter_peptides(
    df: pl.DataFrame,
    enzymes: List[str],
    sequence_col: str = 'Sequence',
) -> pl.DataFrame:
    """
    Keep only fully specific peptides for the given enzyme(s).

    C-terminal and N-terminal residue constraints from all specified enzymes
    are unioned, then a row is kept if it satisfies either side.

    Parameters
    ----------
    df : pl.DataFrame
        Must contain `sequence_col` (or pre-derived 'C_term_aa' / 'N_term_aa').
    enzymes : list[str]
        One or more enzyme names from ENZYMES.
    sequence_col : str
        Column holding the bare (unmodified) peptide sequence.

    Returns
    -------
    pl.DataFrame
        Filtered frame with a print summary of removed rows.
    """
    cterm: Set[str] = set()
    nterm: Set[str] = set()
    unknown = [e for e in enzymes if e.lower() not in ENZYMES]
    if unknown:
        raise ValueError(
            f"Unknown enzyme(s): {unknown}. "
            f"Available: {enzyme_names()}"
        )
    for name in enzymes:
        enz = ENZYMES[name.lower()]
        cterm |= enz.get('cterm', set())
        nterm |= enz.get('nterm', set())

    conditions = []
    if cterm:
        c_col = (
            pl.col('C_term_aa') if 'C_term_aa' in df.columns
            else pl.col(sequence_col).str.slice(-1)
        )
        conditions.append(c_col.is_in(sorted(cterm)))
    if nterm:
        n_col = (
            pl.col('N_term_aa') if 'N_term_aa' in df.columns
            else pl.col(sequence_col).str.slice(0, 1)
        )
        conditions.append(n_col.is_in(sorted(nterm)))

    if not conditions:
        return df

    mask = conditions[0]
    for c in conditions[1:]:
        mask = mask | c

    before = df.height
    filtered = df.filter(mask)
    removed = before - filtered.height
    label = '+'.join(enzymes)
    print(
        f"  [enzyme] {label}: removed {removed:,} non-specific rows "
        f"({100*removed/before:.1f}%) → {filtered.height:,} remaining"
    )
    return filtered
