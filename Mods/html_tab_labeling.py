# Created by @thibaultcolombani with the help of Claude Code (Anthropic)
"""
html_tab_labeling.py — Labeling efficiency tab.

Adapts breakdown to the enzyme filter:
  - With --enzyme trypsin → breakdown by R and K
  - With --enzyme lys-c   → breakdown by K only
  - No enzyme             → C-terminal AAs with ≥ 1% of all precursors, sorted by frequency desc.

"Fully labeled" = all expected amine sites have a PSMtag (N-term + all K side chains).
"Partly labeled" = at least one site is labeled but not all.
"""
from __future__ import annotations
import math
from typing import Dict, List
import polars as pl

TAB_ID    = 'labeling'
TAB_LABEL = 'Labeling'

_MIN_AA_FRAC = 0.01   # minimum fraction of total precursors for an AA to appear in the no-enzyme breakdown

def _get_cleavage_res(processor, breakdown_aas: List[str]) -> frozenset:
    """Return the set of C-terminal cleavage residues for the active enzyme(s).

    With --enzyme: uses the enzyme catalogue.
    No enzyme: falls back to the breakdown AAs themselves (AAs that are common
    C-terminal residues are likely cleavage sites, e.g. K and R for trypsin).
    """
    from Mods.enzyme import ENZYMES
    enzymes = getattr(processor, 'enzymes', None)
    if enzymes:
        cterm: set = set()
        for name in enzymes:
            cterm |= ENZYMES.get(name.lower(), {}).get('cterm', set())
        if cterm:
            return frozenset(cterm)
    return frozenset(breakdown_aas)


def _drop_miscleaved(adf: pl.DataFrame, aa: str,
                     cleavage_res: frozenset) -> pl.DataFrame:
    """Remove miscleaved peptides for a given C-terminal AA group.

    A fully specific peptide ending in `aa` should contain no other residue
    from `cleavage_res` inside it — any extra occurrence is an uncleaved
    internal cleavage site (missed cleavage).
    """
    if aa not in cleavage_res or 'Sequence' not in adf.columns:
        return adf
    pattern = '[' + ''.join(sorted(cleavage_res)) + ']'
    return adf.filter(pl.col('Sequence').str.count_matches(pattern) == 1)


def is_visible(processor) -> bool:
    df = processor.report
    for col in ('fully_labeled', 'partly_labeled'):
        if col in df.columns and int(df[col].sum()) > 0:
            return True
    return False


def _breakdown_aas(processor, df: pl.DataFrame) -> List[str]:
    """
    C-terminal AAs to show.
    - Enzyme specified: the enzyme's C-terminal residues (intersected with data).
    - No enzyme: C-term AAs with ≥ _MIN_AA_FRAC fraction of total precursors.
    """
    from Mods.enzyme import ENZYMES
    enzymes = getattr(processor, 'enzymes', None)
    if enzymes:
        cterm: set = set()
        for name in enzymes:
            cterm |= ENZYMES.get(name.lower(), {}).get('cterm', set())
        if cterm and 'C_term_aa' in df.columns:
            present = set(df['C_term_aa'].drop_nulls().unique().to_list())
            aas = sorted(cterm & present)
            # R first, K second, then remaining alphabetically
            return [aa for aa in ('R', 'K') if aa in aas] + [aa for aa in aas if aa not in ('R', 'K')]

    if 'C_term_aa' not in df.columns:
        return []
    min_count = max(1, int(len(df) * _MIN_AA_FRAC))
    vc = (df['C_term_aa'].drop_nulls()
            .value_counts()
            .sort('count', descending=True))
    aas = [str(row[0]) for row in vc.iter_rows() if row[1] >= min_count]
    # R first, K second, then remaining in frequency order
    ordered = []
    for aa in ('R', 'K'):
        if aa in aas:
            ordered.append(aa)
    ordered += [aa for aa in aas if aa not in ('R', 'K')]
    return ordered


def compute(processor) -> Dict:
    df   = processor.report
    runs = processor.runs()

    breakdown_aas = _breakdown_aas(processor, df)
    cleavage_res  = _get_cleavage_res(processor, breakdown_aas)

    rows: List[Dict] = []
    for run in runs:
        rdf = df.filter(pl.col('Raw file') == run)
        row: Dict = {'run': run, 'tot': rdf.height}
        for aa in breakdown_aas:
            adf      = (rdf.filter(pl.col('C_term_aa') == aa)
                        if 'C_term_aa' in rdf.columns else rdf.clear())
            adf      = _drop_miscleaved(adf, aa, cleavage_res)
            aa_tot   = adf.height
            aa_full  = int(adf['fully_labeled'].sum())   if 'fully_labeled'  in adf.columns else 0
            aa_partly = int(adf['partly_labeled'].sum()) if 'partly_labeled' in adf.columns else 0
            key = aa.lower()
            row[f'{key}_tot']        = aa_tot
            row[f'{key}_full']       = aa_full
            row[f'{key}_partly']     = aa_partly
            row[f'{key}_full_pct']   = round(100 * aa_full   / aa_tot, 1) if aa_tot else 0
            row[f'{key}_partly_pct'] = round(100 * aa_partly / aa_tot, 1) if aa_tot else 0
        rows.append(row)

    # Only show "Partly" in the legend for AAs where at least one run has partly > 0
    partly_aas = [aa for aa in breakdown_aas
                  if any(row.get(f'{aa.lower()}_partly', 0) > 0 for row in rows)]

    result: Dict = {'rows': rows, 'breakdown_aas': breakdown_aas, 'partly_aas': partly_aas}

    if breakdown_aas and 'Retention time' in df.columns:
        result['lab_chrom']      = _lab_chrom(df, runs, breakdown_aas, cleavage_res)
        result['lab_chrom_prec'] = _lab_chrom_prec(df, runs, breakdown_aas, cleavage_res)

    if breakdown_aas and 'seqcharge' in df.columns:
        result['lab_prec'] = _lab_prec(df, runs, _is_reporter_ion_plex(processor))

    return result


# ---------------------------------------------------------------------------
# Python data helpers
# ---------------------------------------------------------------------------

def _lab_chrom(df: pl.DataFrame, runs: List[str], breakdown_aas: List[str],
               cleavage_res: frozenset = frozenset()) -> Dict:
    """RT histogram counts per labeled AA per run (0.5-min bins)."""
    gmin   = float(df['Retention time'].min())
    gmax   = float(df['Retention time'].max())
    bsize  = 0.5
    n_bins = max(1, math.ceil((gmax - gmin) / bsize))
    centers = [round(gmin + (i + 0.5) * bsize, 3) for i in range(n_bins)]

    df2 = df.with_columns(
        ((pl.col('Retention time') - gmin) / bsize).floor().cast(pl.Int32).alias('_b')
    )

    def _bcounts(frame: pl.DataFrame) -> List[int]:
        if frame.is_empty():
            return [0] * n_bins
        grp = frame.group_by('_b').agg(pl.len().alias('c')).sort('_b')
        out = [0] * n_bins
        for b, c in grp.iter_rows():
            if 0 <= b < n_bins:
                out[b] = int(c)
        return out

    labeled_any = pl.lit(False)
    for col in ['fully_labeled', 'partly_labeled']:
        if col in df.columns:
            labeled_any = labeled_any | (pl.col(col) > 0)

    data = []
    for run in runs:
        rdf    = df2.filter(pl.col('Raw file') == run)
        aa_cnt = {}
        for aa in breakdown_aas:
            aa_labeled = (rdf.filter((pl.col('C_term_aa') == aa) & labeled_any)
                          if 'C_term_aa' in df.columns else rdf.clear())
            aa_labeled = _drop_miscleaved(aa_labeled, aa, cleavage_res)
            aa_cnt[aa.lower()] = _bcounts(aa_labeled)
        max_c = max([max(c) for c in aa_cnt.values() if c] + [1])
        data.append({'run': run, **aa_cnt, 'max_count': max_c})

    return {'data': data, 'centers': centers, 'rt_range': [gmin, gmax]}


def _lab_chrom_prec(df: pl.DataFrame, runs: List[str], breakdown_aas: List[str],
                    cleavage_res: frozenset = frozenset()) -> Dict:
    """Per-run per-AA per-bin precursor ID arrays for grouped chromatography view.

    Uses Sequence + '_' + Charge as the ID key so cross-engine groups count
    unique peptides correctly (DIA-NN and Jmod use different seqcharge formats).
    """
    if not all(c in df.columns for c in ['Retention time', 'Sequence', 'Charge', 'C_term_aa']):
        return {}

    gmin   = float(df['Retention time'].min())
    gmax   = float(df['Retention time'].max())
    bsize  = 0.5
    n_bins = max(1, math.ceil((gmax - gmin) / bsize))

    labeled_any = pl.lit(False)
    for col in ['fully_labeled', 'partly_labeled']:
        if col in df.columns:
            labeled_any = labeled_any | (pl.col(col) > 0)

    labeled = df.filter(labeled_any)

    # Remove miscleaved peptides before binning: for each C-terminal residue
    # that is itself a cleavage site, keep only peptides with exactly one
    # occurrence of any cleavage residue in the sequence (the terminal one).
    if cleavage_res and 'Sequence' in labeled.columns and 'C_term_aa' in labeled.columns:
        cr_pattern = '[' + ''.join(sorted(cleavage_res)) + ']'
        labeled = labeled.filter(
            pl.col('C_term_aa').is_in(list(cleavage_res)).not_() |
            (pl.col('Sequence').str.count_matches(cr_pattern) == 1)
        )

    labeled = (labeled
               .with_columns(
                   (pl.col('Sequence') + '_' + pl.col('Charge').cast(pl.String)).alias('_prec_key')
               )
               .select(['Raw file', '_prec_key', 'Retention time', 'C_term_aa'])
               .with_columns(
                   ((pl.col('Retention time') - gmin) / bsize)
                   .floor().cast(pl.Int32).alias('_b')
               ))

    if labeled.is_empty():
        return {'n_bins': n_bins}

    all_sc   = labeled['_prec_key'].drop_nulls().unique().sort().to_list()
    sc_to_id = {sc: i for i, sc in enumerate(all_sc)}

    result: Dict = {'n_bins': n_bins}
    for run in runs:
        rdf = labeled.filter(pl.col('Raw file') == run)
        run_data: Dict = {aa.lower(): [set() for _ in range(n_bins)] for aa in breakdown_aas}

        sc_list  = rdf['_prec_key'].to_list()
        b_list   = rdf['_b'].to_list()
        ct_list  = rdf['C_term_aa'].to_list()

        for sc, b, ct in zip(sc_list, b_list, ct_list):
            if sc is None or b is None or not (0 <= int(b) < n_bins):
                continue
            sc_id = sc_to_id.get(str(sc))
            if sc_id is None or ct not in breakdown_aas:
                continue
            run_data[ct.lower()][int(b)].add(sc_id)

        result[run] = {k: [sorted(s) for s in sets] for k, sets in run_data.items()}

    return result


def _is_reporter_ion_plex(processor) -> bool:
    """True when plex_report Intensity == Precursor.Quantity (TMT reporter ions, not PSMtag MS1)."""
    plex = getattr(processor, 'plex_report', None)
    if plex is None:
        return False
    if 'Intensity' not in plex.columns or 'Precursor.Quantity' not in plex.columns:
        return False
    sample = plex.filter(pl.col('Intensity') > 0).head(500)
    if sample.is_empty():
        return False
    n_diff = int(sample.select((pl.col('Intensity') != pl.col('Precursor.Quantity')).sum())[0, 0])
    return n_diff <= 5


def _lab_prec(df: pl.DataFrame, runs: List[str], force_ms2: bool = False) -> Dict:
    """Fully labeled precursors: {run: {seq_charge_key: [rt_length_s_or_null, log10_int_or_null]}}.

    Intensity source: MS1 (Intensity) when available and non-zero; falls back to
    MS2 (Precursor.Quantity) for TMT-like experiments where Intensity = 0.
    Result includes 'intensity_label': 'MS1' or 'MS2' to drive HTML titles.

    Key is Sequence + '_' + Charge (stripped, modification-free) so cross-engine
    intersection works regardless of tag-notation differences between DIA-NN and Jmod.
    """
    if 'fully_labeled' not in df.columns:
        return {}
    if 'Sequence' not in df.columns or 'Charge' not in df.columns:
        return {}
    labeled = df.filter(pl.col('fully_labeled') > 0)
    labeled = labeled.with_columns(
        (pl.col('Sequence') + '_' + pl.col('Charge').cast(pl.String)).alias('_prec_key')
    )

    has_rl = 'Retention length' in labeled.columns
    # For TMT/reporter-ion plex (force_ms2=True): use Precursor.Quantity (MS2) even when MS1 exists.
    # Otherwise: use MS1 (Intensity) when available and non-zero; fall back to MS2.
    has_ms1_data = ('Intensity' in labeled.columns and (labeled['Intensity'].max() or 0) > 0)
    has_ms2_data = ('Precursor.Quantity' in labeled.columns and
                    (labeled['Precursor.Quantity'].max() or 0) > 0)
    if force_ms2 and has_ms2_data:
        int_col, intensity_label = 'Precursor.Quantity', 'MS2'
    elif has_ms1_data:
        int_col, intensity_label = 'Intensity', 'MS1'
    elif has_ms2_data:
        int_col, intensity_label = 'Precursor.Quantity', 'MS2'
    else:
        int_col, intensity_label = None, 'MS1'

    sel = ['Raw file', '_prec_key']
    if has_rl:   sel.append('Retention length')
    if int_col:  sel.append(int_col)
    slim = labeled.select(sel)

    result: Dict = {'intensity_label': intensity_label}
    for run in runs:
        rdf = slim.filter(pl.col('Raw file') == run)
        prec_map: Dict = {}
        if rdf.is_empty():
            result[run] = prec_map
            continue
        sc_list  = rdf['_prec_key'].to_list()
        rl_list  = rdf['Retention length'].to_list() if has_rl  else [None] * len(sc_list)
        int_list = rdf[int_col].to_list()             if int_col else [None] * len(sc_list)
        for sc, rl, intensity in zip(sc_list, rl_list, int_list):
            if sc is None:
                continue
            rl_s     = round(float(rl) * 60, 3) if rl is not None and float(rl) > 0 else None
            log10_i  = round(math.log10(float(intensity)), 3) if intensity is not None and float(intensity) > 0 else None
            prec_map[str(sc)] = [rl_s, log10_i]
        result[run] = prec_map
    return result


# ---------------------------------------------------------------------------
# HTML
# ---------------------------------------------------------------------------

def html_content(data: Dict) -> str:
    bk = data.get('breakdown_aas', [])
    bk_str = ', '.join(bk) if bk else '—'
    n_runs      = len(data.get('rows', []))
    default_cols = min(max(1, n_runs), 6)
    cols_opts = ''.join(
        f'<option value="{n}"{" selected" if n == default_cols else ""}>{n}</option>'
        for n in range(1, 7)
    )

    html = f"""
<div class="run-selector" style="margin-bottom:0">
  <label style="font-weight:700">View mode:</label>
  <button id="lab-mode-ind" class="btn btn-primary"   onclick="setLabMode('individual')">Individual</button>
  <button id="lab-mode-grp" class="btn btn-secondary" onclick="setLabMode('grouped')">Grouped</button>
  <span id="lab-group-hint" style="display:none;font-size:.82em;color:#666;font-style:italic">
    Runs grouped by base name (split at first <code style="font-style:normal">-</code>). Rename runs in the Files tab to control groups.
  </span>
  <label style="margin-left:16px">Samples per row:</label>
  <select id="lab-chrom-cols-sel" onchange="updateLabChromLayout()">
    {cols_opts}
  </select>
</div>
<div class="section">
  <h3>Tag Labeling Efficiency</h3>
  <p class="help">
    Breakdown by C-terminal residue ({bk_str}).
    <em>Fully labeled</em>: all expected amine sites (N-terminus + all K side chains) carry a PSMtag.
    <em>Partly labeled</em>: at least one site labeled, but not all.
    Bars for each residue are stacked (Full + Partly).
    Apply <code>--enzyme</code> to restrict to specific residues.
    Without <code>--enzyme</code>, residues with ≥1% of precursors are shown.
  </p>
  <div class="chart-grid chart-grid-1">
    <div class="chart-box">
      <h4>Labeling Efficiency (%)</h4>
      <button class="export-btn" onclick="exportChart('ch_lab_pct')">PNG</button>
      <canvas id="ch_lab_pct"></canvas>
    </div>
  </div>
</div>
<div class="section">
  <h3>Precursor Counts</h3>
  <p class="help">Total precursors and per-residue breakdown: total, fully labeled, partly labeled.</p>
  <div class="chart-grid chart-grid-1">
    <div class="chart-box">
      <h4>Precursor Counts</h4>
      <button class="export-btn" onclick="exportChart('ch_lab_cnt')">PNG</button>
      <canvas id="ch_lab_cnt"></canvas>
    </div>
  </div>
</div>"""

    if 'lab_chrom' in data:
        html += f"""
<div class="section">
  <h3>Chromatography — Labeled Peptides</h3>
  <p class="help">Retention time distribution of labeled precursors per C-terminal residue per run (0.5-min bins, labeled = fully or partly).</p>
  <div id="lab-chrom-grid" class="chart-grid" style="grid-template-columns:repeat({default_cols},1fr)"></div>
</div>"""

    if 'lab_prec' in data:
        int_lbl = data['lab_prec'].get('intensity_label', 'MS1')
        html += f"""
<div class="section">
  <h3>Retention Length — Fully Labeled Peptides</h3>
  <p class="help">Peak width FWHM (seconds) for fully labeled precursors. Mean ± SD per run. Right panel: precursors present in all active runs only.</p>
  <div class="chart-grid chart-grid-2">
    <div class="chart-box">
      <h4>All Fully Labeled Peptides</h4>
      <button class="export-btn" onclick="exportChart('ch_lab_rl')">PNG</button>
      <canvas id="ch_lab_rl"></canvas>
    </div>
    <div class="chart-box">
      <h4>Intersected Fully Labeled Peptides</h4>
      <button class="export-btn" onclick="exportChart('ch_lab_rl_int')">PNG</button>
      <canvas id="ch_lab_rl_int"></canvas>
    </div>
  </div>
</div>
<div class="section">
  <h3>{int_lbl} Intensity — Fully Labeled Peptides</h3>
  <p class="help">log₁₀({int_lbl} intensity) for fully labeled precursors. Mean ± SD per run. Right panel: precursors present in all active runs only.</p>
  <div class="chart-grid chart-grid-2">
    <div class="chart-box">
      <h4>All Fully Labeled Peptides</h4>
      <button class="export-btn" onclick="exportChart('ch_lab_ms1')">PNG</button>
      <canvas id="ch_lab_ms1"></canvas>
    </div>
    <div class="chart-box">
      <h4>Intersected Fully Labeled Peptides</h4>
      <button class="export-btn" onclick="exportChart('ch_lab_ms1_int')">PNG</button>
      <canvas id="ch_lab_ms1_int"></canvas>
    </div>
  </div>
</div>
<div class="section">
  <h3>{int_lbl} Ratio — Intersected Fully Labeled Peptides</h3>
  <p class="help">log₂(sample / control) of {int_lbl} intensity for fully labeled precursors present in all active runs. Requires ≥2 active runs.</p>
  <div class="run-selector">
    <label id="lab-ratio-ctrl-lbl">Control run:</label>
    <select id="lab-ratio-ctrl" onchange="updateLabRatio()"></select>
  </div>
  <div class="chart-grid chart-grid-1">
    <div class="chart-box">
      <h4 id="lab-ratio-title">{int_lbl} Ratio (sample / control)</h4>
      <button class="export-btn" onclick="exportChart('ch_lab_ratio')">PNG</button>
      <canvas id="ch_lab_ratio"></canvas>
    </div>
  </div>
</div>"""

    return html


# ---------------------------------------------------------------------------
# JavaScript
# ---------------------------------------------------------------------------

def javascript() -> str:
    return """
// ── Labeling tab ─────────────────────────────────────────────────────────────
let _labMode = 'individual';

// Color pairs per C-terminal amino acid: [fullColor, partlyColor].
// Full colors match _AA_COL for consistency with ion sampling / chromatography tabs.
// Partly colors use a lighter, distinct hue so "Partly" is never confused with "Total"
// (which is derived from the full color at 0.40 opacity).
const _LAB_AA_COLS = {
  R: ['rgba(231,76,60,0.8)',   'rgba(255,138,128,0.70)'],  // dark red  / salmon
  K: ['rgba(39,174,96,0.8)',   'rgba(102,204,102,0.70)'],  // dark green / lime
  E: ['rgba(255,159,64,0.8)',  'rgba(255,210,160,0.70)'],  // orange    / pale orange
  D: ['rgba(54,162,235,0.8)',  'rgba(144,202,249,0.70)'],  // blue      / sky blue
  F: ['rgba(121,85,72,0.8)',   'rgba(188,143,143,0.70)'],  // brown     / rose brown
  Y: ['rgba(96,125,139,0.8)',  'rgba(176,190,197,0.70)'],  // slate     / light slate
  W: ['rgba(241,196,15,0.8)',  'rgba(255,236,153,0.70)'],  // gold      / pale gold
  L: ['rgba(52,152,219,0.8)',  'rgba(144,202,249,0.70)'],  // blue      / sky blue
  M: ['rgba(155,89,182,0.8)',  'rgba(206,161,215,0.70)'],  // purple    / lavender
};
const _LAB_DEFAULT_FULL   = 'rgba(149,165,166,0.8)';
const _LAB_DEFAULT_PARTLY = 'rgba(149,165,166,0.35)';
function _labAACol(aa, idx) { return (_LAB_AA_COLS[aa] || [_LAB_DEFAULT_FULL, _LAB_DEFAULT_PARTLY])[idx]; }

function setLabMode(mode) {
  _labMode = mode;
  document.getElementById('lab-mode-ind').className = mode === 'individual' ? 'btn btn-primary' : 'btn btn-secondary';
  document.getElementById('lab-mode-grp').className = mode === 'grouped'    ? 'btn btn-primary' : 'btn btn-secondary';
  const hint = document.getElementById('lab-group-hint');
  if (hint) hint.style.display = mode === 'grouped' ? 'inline' : 'none';
  const ctrlLbl = document.getElementById('lab-ratio-ctrl-lbl');
  if (ctrlLbl) ctrlLbl.textContent = mode === 'grouped' ? 'Control group:' : 'Control run:';
  _redrawLabelingCharts();
}

function _parseSampleBaseName(r) {
  const n = label(r); const i = n.indexOf('-');
  return i < 0 ? n : n.substring(0, i);
}
function _getLabGroups() {
  const ar = getActiveRuns(), order = [], map = {};
  ar.forEach(r => {
    const bn = _parseSampleBaseName(r);
    if (!map[bn]) { map[bn] = []; order.push(bn); }
    map[bn].push(r);
  });
  return order.map(bn => ({ baseName: bn, runs: map[bn] }));
}

function initLabeling() { _redrawLabelingCharts(); }

function _redrawLabelingCharts() {
  const ar  = getActiveRuns();
  const lab = ar.map(r => DATA.labeling.rows.find(d => d.run === r)).filter(Boolean);
  if (_labMode === 'individual') {
    _drawLabBarsInd(ar, lab);
    if (DATA.labeling.lab_chrom) _initLabChrom();
    if (DATA.labeling.lab_prec)  { _drawLabScatter(); _initLabRatioDropdownInd(); }
  } else {
    const groups = _getLabGroups();
    _drawLabBarsGrp(groups);
    if (DATA.labeling.lab_chrom) _initLabChromGrp(groups);
    if (DATA.labeling.lab_prec)  { _drawLabScatterGrp(groups); _initLabRatioDropdownGrp(groups); }
  }
}

// ── Build labeling datasets dynamically from breakdown_aas ────────────────────
function _labMakeDatasets(lab, keyFull, keyPartly, aaList) {
  const ds = [];
  const partlyAAs = new Set(DATA.labeling.partly_aas || []);
  aaList.forEach(aa => {
    const aL = aa.toLowerCase();
    ds.push({
      label: `${aa} Fully Labeled %`, stack: aa,
      data: lab.map(d => d[`${aL}_${keyFull}`] || 0),
      backgroundColor: _labAACol(aa, 0),
      borderColor: _labAACol(aa, 0).replace('0.8', '1').replace('0.35','1'), borderWidth: 1,
    });
    if (partlyAAs.has(aa)) {
      ds.push({
        label: `${aa} Partly Labeled %`, stack: aa,
        data: lab.map(d => d[`${aL}_${keyPartly}`] || 0),
        backgroundColor: _labAACol(aa, 1),
        borderColor: _labAACol(aa, 1).replace('0.70', '1'), borderWidth: 1,
      });
    }
  });
  return ds;
}

function _labMakeCountDatasets(lab, aaList) {
  const ds = [
    { label:'Total', data:lab.map(d=>d.tot), backgroundColor:'rgba(128,128,128,0.6)', borderColor:'rgba(128,128,128,1)', borderWidth:1 },
  ];
  const partlyAAs = new Set(DATA.labeling.partly_aas || []);
  aaList.forEach(aa => {
    const aL = aa.toLowerCase();
    ds.push({ label:`Total ${aa}`,  data:lab.map(d=>d[`${aL}_tot`] ||0), backgroundColor:_labAACol(aa,0).replace('0.8','0.4'), borderColor:_labAACol(aa,0).replace('0.8','0.7'), borderWidth:1 });
    ds.push({ label:`${aa} Fully`,  data:lab.map(d=>d[`${aL}_full`]||0), backgroundColor:_labAACol(aa,0),                      borderColor:_labAACol(aa,0).replace('0.8','1'),   borderWidth:1 });
    if (partlyAAs.has(aa)) {
      ds.push({ label:`${aa} Partly`, data:lab.map(d=>d[`${aL}_partly`]||0), backgroundColor:_labAACol(aa,1), borderColor:_labAACol(aa,1).replace('0.35','0.7'), borderWidth:1 });
    }
  });
  return ds;
}

// ── Bar charts — Individual ───────────────────────────────────────────────────
function _drawLabBarsInd(ar, lab) {
  const bkAAs = DATA.labeling.breakdown_aas || [];
  const names = lab.map(d => label(d.run));
  barChart('ch_lab_pct', names, _labMakeDatasets(lab, 'full_pct', 'partly_pct', bkAAs), {
    scales: {
      x: { stacked: true },
      y: { stacked: true, max: 100, title: { display: true, text: 'Labeling Efficiency (%)' } },
    },
    aspectRatio: 2.5,
    plugins: {
      tooltip: {
        mode: 'index',
        callbacks: {
          label: (ctx) => {
            const d   = lab[ctx.dataIndex];
            const val = ctx.parsed.y.toFixed(1);
            const lbl = ctx.dataset.label;
            const aa  = lbl.split(' ')[0];
            const aL  = aa.toLowerCase();
            if (lbl.includes('Fully'))  return `${aa} Full: ${val}% (${d[aL+'_full']}/${d[aL+'_tot']})`;
            if (lbl.includes('Partly')) return `${aa} Partly: ${val}% (${d[aL+'_partly']}/${d[aL+'_tot']})`;
            return `${lbl}: ${val}%`;
          }
        }
      }
    }
  });
  barChart('ch_lab_cnt', names, _labMakeCountDatasets(lab, bkAAs),
    { aspectRatio: 2.5, y: { title: { display: true, text: 'Precursor Count' } } });
}

// ── SD whisker plugin ─────────────────────────────────────────────────────────
const _labSdPlugin = {
  id: '_labSd',
  afterDatasetsDraw(chart) {
    const {ctx, scales: {y}} = chart;
    if (!y) return;
    const zeroY = y.getPixelForValue(0);
    chart.data.datasets.forEach((ds, di) => {
      if (!ds._sd) return;
      const meta = chart.getDatasetMeta(di);
      if (meta.hidden) return;
      ctx.save();
      ctx.strokeStyle = 'rgba(30,30,30,0.75)'; ctx.lineWidth = 1.5; ctx.lineCap = 'round';
      ds._sd.forEach((sd, bi) => {
        if (!sd || !isFinite(sd) || sd <= 0) return;
        const bar = meta.data[bi];
        if (!bar || bar.hidden) return;
        const sdPx = Math.abs(zeroY - y.getPixelForValue(sd));
        const topY = bar.y, capW = Math.min((bar.width || 12) / 2, 6);
        ctx.beginPath();
        ctx.moveTo(bar.x, topY); ctx.lineTo(bar.x, topY - sdPx);
        ctx.moveTo(bar.x - capW, topY - sdPx); ctx.lineTo(bar.x + capW, topY - sdPx);
        ctx.stroke();
      });
      ctx.restore();
    });
  }
};

// ── Bar charts — Grouped ──────────────────────────────────────────────────────
function _drawLabBarsGrp(groups) {
  const bkAAs  = DATA.labeling.breakdown_aas || [];
  const names  = groups.map(g => g.baseName);
  const allRows = DATA.labeling.rows;

  const gStats = groups.map(g => {
    const rows = g.runs.map(r => allRows.find(d => d.run === r)).filter(Boolean);
    if (!rows.length) return null;
    const mean = k => rows.reduce((s, d) => s + (d[k] || 0), 0) / rows.length;
    const sd   = k => { const m = mean(k); return Math.sqrt(rows.reduce((s, d) => s + ((d[k]||0)-m)**2, 0) / rows.length); };
    const stat = { tot: mean('tot'), tot_sd: sd('tot') };
    bkAAs.forEach(aa => {
      const aL = aa.toLowerCase();
      ['tot','full','partly'].forEach(k => {
        stat[`${aL}_${k}`]    = mean(`${aL}_${k}`);
        stat[`${aL}_${k}_sd`] = sd(`${aL}_${k}`);
      });
      ['full_pct','partly_pct'].forEach(k => {
        stat[`${aL}_${k}`]    = mean(`${aL}_${k}`);
      });
      const totPcts = rows.map(d => (d[`${aL}_full_pct`]||0) + (d[`${aL}_partly_pct`]||0));
      const totM    = totPcts.reduce((a,b)=>a+b,0) / totPcts.length;
      stat[`${aL}_tot_pct_sd`] = Math.sqrt(totPcts.reduce((a,b)=>a+(b-totM)**2,0) / totPcts.length);
    });
    return stat;
  });

  // Efficiency chart (grouped, with SD on top of stack per AA)
  destroyChart('ch_lab_pct');
  const ctx1 = document.getElementById('ch_lab_pct');
  if (ctx1) {
    const efDs = [];
    const partlyAAs = new Set(DATA.labeling.partly_aas || []);
    const sdVals = aa => gStats.map(d => d?.[`${aa.toLowerCase()}_tot_pct_sd`] ?? 0);
    bkAAs.forEach(aa => {
      const aL = aa.toLowerCase();
      const hasPartly = partlyAAs.has(aa);
      efDs.push({ label:`${aa} Fully Labeled %`, stack:aa,
        data: gStats.map(d => +(d?.[`${aL}_full_pct`] ??0).toFixed(1)),
        ...(!hasPartly ? {_sd: sdVals(aa)} : {}),
        backgroundColor: _labAACol(aa,0), borderColor: _labAACol(aa,0).replace('0.8','1'), borderWidth:1 });
      if (hasPartly) {
        efDs.push({ label:`${aa} Partly Labeled %`, stack:aa,
          data: gStats.map(d => +(d?.[`${aL}_partly_pct`] ??0).toFixed(1)),
          _sd: sdVals(aa),
          backgroundColor: _labAACol(aa,1), borderColor: _labAACol(aa,1).replace('0.35','0.7'), borderWidth:1 });
      }
    });
    CHARTS['ch_lab_pct'] = new Chart(ctx1, {
      type:'bar', data:{labels:names, datasets:efDs},
      options:{
        responsive:true, maintainAspectRatio:true, aspectRatio:2.5,
        scales:{x:{stacked:true,ticks:{minRotation:90,maxRotation:90,font:{size:10}}},
                y:{stacked:true,beginAtZero:true,max:100,title:{display:true,text:'Labeling Efficiency — Group Mean (%)'}}},
        plugins:{legend:{display:true},
          tooltip:{mode:'index',callbacks:{label:(ctx)=>{
            const d=gStats[ctx.dataIndex]; if(!d) return '';
            const v=ctx.parsed.y.toFixed(1), n=groups[ctx.dataIndex]?.runs.length||1;
            const sd=(ctx.dataset._sd?.[ctx.dataIndex]??0).toFixed(1);
            return `${ctx.dataset.label}: ${v}%${ctx.dataset._sd?` ±${sd}%`:''} (n=${n})`;
          }}}
        }
      },
      plugins:[_labSdPlugin],
    });
  }

  // Counts chart
  destroyChart('ch_lab_cnt');
  const ctx2 = document.getElementById('ch_lab_cnt');
  if (ctx2) {
    const cntDs = [
      { label:'Total', data:gStats.map(d=>Math.round(d?.tot??0)), _sd:gStats.map(d=>d?.tot_sd??0),
        backgroundColor:'rgba(128,128,128,0.6)', borderColor:'rgba(128,128,128,1)', borderWidth:1 },
    ];
    bkAAs.forEach(aa => {
      const aL = aa.toLowerCase();
      cntDs.push({ label:`Total ${aa}`,  data:gStats.map(d=>Math.round(d?.[`${aL}_tot`]   ??0)), _sd:gStats.map(d=>d?.[`${aL}_tot_sd`]    ??0), backgroundColor:_labAACol(aa,0).replace('0.8','0.4'), borderColor:_labAACol(aa,0).replace('0.8','0.7'), borderWidth:1 });
      cntDs.push({ label:`${aa} Fully`,  data:gStats.map(d=>Math.round(d?.[`${aL}_full`]  ??0)), _sd:gStats.map(d=>d?.[`${aL}_full_sd`]   ??0), backgroundColor:_labAACol(aa,0),                       borderColor:_labAACol(aa,0).replace('0.8','1'),   borderWidth:1 });
      cntDs.push({ label:`${aa} Partly`, data:gStats.map(d=>Math.round(d?.[`${aL}_partly`]??0)), _sd:gStats.map(d=>d?.[`${aL}_partly_sd`] ??0), backgroundColor:_labAACol(aa,1),                       borderColor:_labAACol(aa,1).replace('0.35','0.7'),borderWidth:1 });
    });
    CHARTS['ch_lab_cnt'] = new Chart(ctx2, {
      type:'bar', data:{labels:names, datasets:cntDs},
      options:{
        responsive:true, maintainAspectRatio:true, aspectRatio:2.5,
        scales:{x:{ticks:{minRotation:90,maxRotation:90,font:{size:10}}},
                y:{beginAtZero:true,title:{display:true,text:'Precursor Count — Group Mean'}}},
        plugins:{legend:{display:true},
          tooltip:{mode:'index',callbacks:{label:(ctx)=>{
            const n=groups[ctx.dataIndex]?.runs.length||1;
            const v=Math.round(ctx.parsed.y).toLocaleString();
            const sd=Math.round(ctx.dataset._sd?.[ctx.dataIndex]??0).toLocaleString();
            return `${ctx.dataset.label}: ${v} ±${sd} (n=${n})`;
          }}}
        }
      },
      plugins:[_labSdPlugin],
    });
  }
}

// ── Chromatography helpers ────────────────────────────────────────────────────
function _labNiceTicks(maxVal) {
  if (!maxVal || maxVal <= 0) return { step: 1, axMax: 4 };
  const rawStep = maxVal / 4;
  const mag     = Math.pow(10, Math.floor(Math.log10(rawStep)));
  const step    = Math.ceil(rawStep / mag) * mag;
  return { step, axMax: step * 4 };
}

function _labRtChart(canvasId, datasets, rtMin, rtMax, xMax, xStep) {
  destroyChart(canvasId);
  const ctx = document.getElementById(canvasId);
  if (!ctx) return;
  const yMin = Math.floor(rtMin / 5) * 5;
  const yMax = Math.ceil(rtMax  / 5) * 5;
  CHARTS[canvasId] = new Chart(ctx, {
    type: 'line', data: { datasets },
    options: {
      responsive:true, maintainAspectRatio:false, parsing:false,
      layout:{padding:{bottom:8,right:12}},
      plugins:{legend:{display:true,position:'top',labels:{boxWidth:12,font:{size:10}}}},
      elements:{point:{radius:2.5,hoverRadius:4},line:{borderWidth:1.5}},
      scales:{
        x:{type:'linear',min:0,max:xMax,title:{display:true,text:'Count'},
           ticks:{stepSize:xStep,precision:0,maxRotation:0,minRotation:0,font:{size:9}}},
        y:{type:'linear',reverse:true,min:yMin,max:yMax,
           title:{display:true,text:'RT (min)'},
           ticks:{stepSize:0.5,callback:(v)=>Math.abs(Math.round(v*10)%50)<1?v:null},
           grid:{color:(ctx)=>Math.abs(Math.round(ctx.tick.value*10)%50)<1?'#ddd':'transparent'}},
      }
    }
  });
}

function _labChromBox(canvasId, title) {
  const box = document.createElement('div');
  box.className = 'chart-box';
  box.innerHTML =
    `<h4>${title}</h4>` +
    `<button class="export-btn" onclick="exportChart('${canvasId}')">PNG</button>` +
    `<div style="position:relative;height:375px"><canvas id="${canvasId}"></canvas></div>`;
  return box;
}

// Per-AA line color for chromatography (same as bar chart full color)
function _labChromColor(aa) { return (_LAB_AA_COLS[aa] || [_LAB_DEFAULT_FULL])[0].replace('0.8','1'); }

// ── Chromatography — Individual ───────────────────────────────────────────────
function _initLabChrom() {
  const d      = DATA.labeling.lab_chrom;
  if (!d || !d.data) return;
  const bkAAs  = DATA.labeling.breakdown_aas || [];
  const ar     = getActiveRuns();
  const defCols = Math.min(ar.length || 1, 6);
  const colsSel = document.getElementById('lab-chrom-cols-sel');
  if (colsSel) colsSel.value = String(defCols);
  const grid = document.getElementById('lab-chrom-grid');
  if (!grid) return;
  grid.style.gridTemplateColumns = `repeat(${defCols},1fr)`;
  grid.innerHTML = '';
  const [rtMin, rtMax] = d.rt_range || [0, 120];
  const centers        = d.centers  || [];
  const maxCount = Math.max(1, ...ar.map(r => (d.data.find(x=>x.run===r)||{}).max_count||0));
  const {step:xStep, axMax:xMax} = _labNiceTicks(maxCount);
  ar.forEach((run,idx) => { if (d.data.find(x=>x.run===run)) grid.appendChild(_labChromBox(`ch_lab_rt_${idx}`,label(run))); });
  ar.forEach((run,idx) => {
    const rd = d.data.find(x=>x.run===run);
    if (!rd) return;
    const ds = bkAAs.map(aa => ({
      label: `${aa} Labeled`,
      data:  centers.map((c,j)=>({x:rd[aa.toLowerCase()]?.[j]||0, y:c})),
      borderColor: _labChromColor(aa), backgroundColor:'transparent', borderWidth:1.5,
    }));
    _labRtChart(`ch_lab_rt_${idx}`, ds, rtMin, rtMax, xMax, xStep);
  });
}

// ── Chromatography — Grouped ──────────────────────────────────────────────────
function _initLabChromGrp(groups) {
  const d    = DATA.labeling.lab_chrom;
  const prec = DATA.labeling.lab_chrom_prec;
  if (!d || !d.centers) return;
  const bkAAs   = DATA.labeling.breakdown_aas || [];
  const defCols = Math.min(groups.length||1, 6);
  const colsSel = document.getElementById('lab-chrom-cols-sel');
  if (colsSel) colsSel.value = String(defCols);
  const grid = document.getElementById('lab-chrom-grid');
  if (!grid) return;
  grid.style.gridTemplateColumns = `repeat(${defCols},1fr)`;
  grid.innerHTML = '';
  const [rtMin, rtMax] = d.rt_range || [0, 120];
  const centers = d.centers || [];
  const nBins   = centers.length;

  let groupChrom;
  if (prec && prec.n_bins) {
    groupChrom = groups.map(g => {
      const aaSets = Object.fromEntries(bkAAs.map(aa => [aa.toLowerCase(), Array.from({length:nBins},()=>new Set())]));
      g.runs.forEach(r => {
        const rd = prec[r];
        if (!rd) return;
        bkAAs.forEach(aa => {
          const aL = aa.toLowerCase();
          (rd[aL]||[]).forEach((ids,bi) => { if(bi<nBins) ids.forEach(id=>aaSets[aL][bi].add(id)); });
        });
      });
      const aaCounts = Object.fromEntries(bkAAs.map(aa=>[aa.toLowerCase(), aaSets[aa.toLowerCase()].map(s=>s.size)]));
      const maxC = Math.max(...Object.values(aaCounts).flatMap(v=>v), 1);
      return { baseName:g.baseName, ...aaCounts, max_count:maxC };
    });
  } else {
    groupChrom = groups.map(g => {
      const aaCounts = Object.fromEntries(bkAAs.map(aa=>[aa.toLowerCase(), Array(nBins).fill(0)]));
      g.runs.forEach(r => {
        const rd = (d.data||[]).find(x=>x.run===r);
        if (!rd) return;
        bkAAs.forEach(aa => {
          const aL = aa.toLowerCase();
          (rd[aL]||[]).forEach((v,i)=>{ aaCounts[aL][i] += v||0; });
        });
      });
      const maxC = Math.max(...Object.values(aaCounts).flatMap(v=>v), 1);
      return { baseName:g.baseName, ...aaCounts, max_count:maxC };
    });
  }

  const maxCount = Math.max(1, ...groupChrom.map(g=>g.max_count));
  const {step:xStep, axMax:xMax} = _labNiceTicks(maxCount);
  groupChrom.forEach((g,idx) => grid.appendChild(_labChromBox(`ch_lab_rt_${idx}`,g.baseName)));
  groupChrom.forEach((g,idx) => {
    const ds = bkAAs.map(aa => ({
      label: `${aa} Labeled (unique)`,
      data:  centers.map((c,j)=>({x:g[aa.toLowerCase()]?.[j]||0, y:c})),
      borderColor: _labChromColor(aa), backgroundColor:'transparent', borderWidth:1.5,
    }));
    _labRtChart(`ch_lab_rt_${idx}`, ds, rtMin, rtMax, xMax, xStep);
  });
}

// ── Layout update ─────────────────────────────────────────────────────────────
function updateLabChromLayout() {
  const n = parseInt(document.getElementById('lab-chrom-cols-sel').value);
  const grid = document.getElementById('lab-chrom-grid');
  if (grid) grid.style.gridTemplateColumns = `repeat(${n},1fr)`;
  for (let i=0;i<50;i++) destroyChart(`ch_lab_rt_${i}`);
  if (!DATA.labeling.lab_chrom?.data) return;
  requestAnimationFrame(() => {
    if (_labMode === 'individual') _initLabChrom();
    else _initLabChromGrp(_getLabGroups());
  });
}

// ── Dot / scatter mean±SD ─────────────────────────────────────────────────────
// runColors: optional array of hex colors, one per entry in dotData, for consistent per-run coloring.
function _labDotChart(canvasId, names, dotData, yTitle, aspectRatio=2.0, runColors=null) {
  destroyChart(canvasId);
  const ctx = document.getElementById(canvasId);
  if (!ctx) return;
  const valid = dotData.map((d,i)=>d?{...d,i}:null).filter(Boolean);
  if (!valid.length) return;
  const yCandidates = valid.flatMap(d=>[d.mean-d.std,d.mean+d.std]).filter(isFinite);
  const rawMin = Math.min(...yCandidates), rawMax = Math.max(...yCandidates);
  const pad = (rawMax-rawMin)*0.12||0.5;
  const _DOTCOLS = ['#3498db','#e74c3c','#27ae60','#f39c12','#9b59b6','#1abc9c','#e67e22','#2980b9','#c0392b','#16a085'];
  const datasets = [];
  valid.forEach((d,vi) => {
    const col = (runColors && runColors[d.i]) ? runColors[d.i] : _DOTCOLS[vi%_DOTCOLS.length];
    datasets.push({label:names[d.i],type:'scatter',data:[{x:d.i,y:d.mean}],
      backgroundColor:col+'bb',borderColor:col,borderWidth:2,pointRadius:8,pointHoverRadius:10,showLine:false});
    datasets.push({label:'_err',type:'line',data:[{x:d.i,y:d.mean-d.std},{x:d.i,y:d.mean+d.std}],
      borderColor:col,backgroundColor:'transparent',pointRadius:0,showLine:true,borderWidth:2,tension:0});
  });
  CHARTS[canvasId] = new Chart(ctx, {
    data:{datasets},
    options:{
      responsive:true,maintainAspectRatio:true,aspectRatio,parsing:false,
      plugins:{
        legend:{display:false},
        tooltip:{filter:item=>item.dataset.label!=='_err',callbacks:{
          title:(ctxs)=>names[Math.round(ctxs[0].parsed.x)]||'',
          label:(ctx)=>{ const d=valid.find(d=>d.i===Math.round(ctx.parsed.x));
            return d?`${yTitle}: ${d.mean.toFixed(3)} ± ${d.std.toFixed(3)} (n=${d.n})`:''; }
        }}
      },
      scales:{
        x:{type:'linear',min:-0.5,max:dotData.length-0.5,
           afterBuildTicks(scale){scale.ticks=names.map((_,i)=>({value:i}));},
           ticks:{callback:(v)=>{const i=Math.round(v);return(i>=0&&i<names.length)?names[i]:'';},
                  maxRotation:90,minRotation:90,font:{size:9}},grid:{display:false}},
        y:{min:rawMin-pad,max:rawMax+pad,title:{display:true,text:yTitle}},
      },
    },
  });
}

function _labRunStats(entries, idx) {
  const vals = entries.map(e=>Array.isArray(e)?e[idx]:null).filter(v=>v!=null&&isFinite(v));
  if (!vals.length) return null;
  const mean = vals.reduce((a,b)=>a+b,0)/vals.length;
  const std  = Math.sqrt(vals.reduce((a,b)=>a+(b-mean)**2,0)/vals.length);
  return {mean:+mean.toFixed(3),std:+std.toFixed(3),n:vals.length};
}

function _drawLabScatter() {
  const ar=getActiveRuns(), names=ar.map(r=>label(r)), prec=DATA.labeling.lab_prec||{};
  const runColors = ar.map(r => color(r));
  const allRlData  = ar.map(run=>_labRunStats(Object.values(prec[run]||{}),0));
  const allMs1Data = ar.map(run=>_labRunStats(Object.values(prec[run]||{}),1));
  let common = null;
  if (ar.length>=2) {
    // Filter to array entries only (excludes 'intensity_label' string key)
    let s=new Set(Object.keys(prec[ar[0]]||{}).filter(k=>Array.isArray((prec[ar[0]]||{})[k])));
    for(let i=1;i<ar.length;i++) s=new Set([...s].filter(k=>Array.isArray((prec[ar[i]]||{})[k])));
    if(s.size>0) common=[...s];
  }
  const intRlData  = common?ar.map(run=>_labRunStats(common.map(sc=>prec[run][sc]),0)):ar.map(()=>null);
  const intMs1Data = common?ar.map(run=>_labRunStats(common.map(sc=>prec[run][sc]),1)):ar.map(()=>null);
  const intLbl = (DATA.labeling.lab_prec?.intensity_label || 'MS1');
  const intAxisLbl = `log₁₀(${intLbl})`;
  _labDotChart('ch_lab_rl',     names,allRlData,  'FWHM (s)',  1.33, runColors);
  _labDotChart('ch_lab_rl_int', names,intRlData,  'FWHM (s)',  1.33, runColors);
  _labDotChart('ch_lab_ms1',    names,allMs1Data, intAxisLbl,  1.33, runColors);
  _labDotChart('ch_lab_ms1_int',names,intMs1Data, intAxisLbl,  1.33, runColors);
}

function _drawLabScatterGrp(groups) {
  const prec=DATA.labeling.lab_prec||{}, names=groups.map(g=>g.baseName);
  const grpColors = groups.map(g => color(g.runs[0]));
  function _grpStats(g,idx) {
    const vals=g.runs.flatMap(r=>Object.values(prec[r]||{}).map(v=>Array.isArray(v)?v[idx]:null).filter(v=>v!=null&&isFinite(v)));
    if(!vals.length) return null;
    const mean=vals.reduce((a,b)=>a+b,0)/vals.length;
    return {mean:+mean.toFixed(3),std:+Math.sqrt(vals.reduce((a,b)=>a+(b-mean)**2,0)/vals.length).toFixed(3),n:vals.length};
  }
  const ar=getActiveRuns(); let common=null;
  if(ar.length>=2){let s=new Set(Object.keys(prec[ar[0]]||{}).filter(k=>Array.isArray((prec[ar[0]]||{})[k])));for(let i=1;i<ar.length;i++)s=new Set([...s].filter(k=>Array.isArray((prec[ar[i]]||{})[k])));if(s.size>0)common=[...s];}
  function _grpStatsInt(g,idx){
    if(!common) return null;
    const vals=g.runs.flatMap(r=>common.map(sc=>{const v=prec[r]?.[sc];return Array.isArray(v)?v[idx]:null;}).filter(v=>v!=null&&isFinite(v)));
    if(!vals.length) return null;
    const mean=vals.reduce((a,b)=>a+b,0)/vals.length;
    return {mean:+mean.toFixed(3),std:+Math.sqrt(vals.reduce((a,b)=>a+(b-mean)**2,0)/vals.length).toFixed(3),n:vals.length};
  }
  const intLbl2 = (DATA.labeling.lab_prec?.intensity_label || 'MS1');
  const intAxisLbl2 = `log₁₀(${intLbl2})`;
  _labDotChart('ch_lab_rl',     names,groups.map(g=>_grpStats(g,0)),    'FWHM (s)',  1.33, grpColors);
  _labDotChart('ch_lab_rl_int', names,groups.map(g=>_grpStatsInt(g,0)), 'FWHM (s)',  1.33, grpColors);
  _labDotChart('ch_lab_ms1',    names,groups.map(g=>_grpStats(g,1)),    intAxisLbl2, 1.33, grpColors);
  _labDotChart('ch_lab_ms1_int',names,groups.map(g=>_grpStatsInt(g,1)), intAxisLbl2, 1.33, grpColors);
}

// ── Ratio ─────────────────────────────────────────────────────────────────────
function _initLabRatioDropdownInd() {
  const ar=getActiveRuns(); if(ar.length<2) return;
  const sel=document.getElementById('lab-ratio-ctrl'); if(!sel) return;
  sel.innerHTML=ar.map((r,i)=>`<option value="${r}"${i===0?' selected':''}>${label(r)}</option>`).join('');
  updateLabRatio();
}
function _updateLabRatioInd() {
  const sel=document.getElementById('lab-ratio-ctrl'), ctrl=sel?.value; if(!ctrl) return;
  const ar=getActiveRuns(), prec=DATA.labeling.lab_prec||{}, names=ar.map(r=>label(r));
  const titleEl=document.getElementById('lab-ratio-title');
  const _intLblR = DATA.labeling.lab_prec?.intensity_label || 'MS1';
  if(titleEl) titleEl.textContent=`${_intLblR} Ratio — log₂(sample / ${label(ctrl)})`;
  let common=new Set(Object.keys(prec[ar[0]]||{}).filter(k=>Array.isArray((prec[ar[0]]||{})[k])));
  for(let i=1;i<ar.length;i++) common=new Set([...common].filter(k=>Array.isArray((prec[ar[i]]||{})[k])));
  if(!common.size){destroyChart('ch_lab_ratio');return;}
  const LG2=Math.log10(2);
  const dotData=ar.map(run=>{
    const log2FCs=[...common].map(sc=>{
      const cLog=(prec[ctrl]?.[sc]||[])[1], rLog=(prec[run]?.[sc]||[])[1];
      return(cLog==null||rLog==null||!isFinite(cLog)||!isFinite(rLog))?null:(rLog-cLog)/LG2;
    }).filter(v=>v!=null&&isFinite(v));
    if(!log2FCs.length) return null;
    const mean=log2FCs.reduce((a,b)=>a+b,0)/log2FCs.length;
    return {mean:+mean.toFixed(3),std:+Math.sqrt(log2FCs.reduce((a,b)=>a+(b-mean)**2,0)/log2FCs.length).toFixed(3),n:log2FCs.length};
  });
  _labDotChart('ch_lab_ratio',names,dotData,'log₂(sample / control)',4.0, ar.map(r=>color(r)));
}
function _initLabRatioDropdownGrp(groups) {
  if(groups.length<2) return;
  const sel=document.getElementById('lab-ratio-ctrl'); if(!sel) return;
  sel.innerHTML=groups.map((g,i)=>`<option value="${g.baseName}"${i===0?' selected':''}>${g.baseName}</option>`).join('');
  updateLabRatio();
}
function _updateLabRatioGrp() {
  const sel=document.getElementById('lab-ratio-ctrl'), ctrlBn=sel?.value; if(!ctrlBn) return;
  const groups=_getLabGroups(), prec=DATA.labeling.lab_prec||{}, ar=getActiveRuns();
  const names=groups.map(g=>g.baseName);
  const titleEl=document.getElementById('lab-ratio-title');
  const _intLblRG = DATA.labeling.lab_prec?.intensity_label || 'MS1';
  if(titleEl) titleEl.textContent=`${_intLblRG} Ratio — log₂(group / ${ctrlBn})`;
  let common=new Set(Object.keys(prec[ar[0]]||{}).filter(k=>Array.isArray((prec[ar[0]]||{})[k])));
  for(let i=1;i<ar.length;i++) common=new Set([...common].filter(k=>Array.isArray((prec[ar[i]]||{})[k])));
  if(!common.size){destroyChart('ch_lab_ratio');return;}
  const ctrlGrp=groups.find(g=>g.baseName===ctrlBn); if(!ctrlGrp) return;
  const LG2=Math.log10(2);
  const dotData=groups.map(g=>{
    const log2FCs=[...common].map(sc=>{
      const ctrlVals=ctrlGrp.runs.map(r=>(prec[r]?.[sc]||[])[1]).filter(v=>v!=null&&isFinite(v));
      const smpVals=g.runs.map(r=>(prec[r]?.[sc]||[])[1]).filter(v=>v!=null&&isFinite(v));
      if(!ctrlVals.length||!smpVals.length) return null;
      return(smpVals.reduce((a,b)=>a+b,0)/smpVals.length - ctrlVals.reduce((a,b)=>a+b,0)/ctrlVals.length)/LG2;
    }).filter(v=>v!=null&&isFinite(v));
    if(!log2FCs.length) return null;
    const mean=log2FCs.reduce((a,b)=>a+b,0)/log2FCs.length;
    return {mean:+mean.toFixed(3),std:+Math.sqrt(log2FCs.reduce((a,b)=>a+(b-mean)**2,0)/log2FCs.length).toFixed(3),n:log2FCs.length};
  });
  _labDotChart('ch_lab_ratio',names,dotData,'log₂(group / control)',4.0, groups.map(g=>color(g.runs[0])));
}
function updateLabRatio() {
  if(_labMode==='individual') _updateLabRatioInd(); else _updateLabRatioGrp();
}
"""
