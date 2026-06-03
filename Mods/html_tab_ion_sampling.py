# Created by @thibaultcolombani with the help of Claude Code (Anthropic)
"""html_tab_ion_sampling.py — Ion Sampling tab: MS1/MS2 intensity distributions, fragment m/z, ratios."""
from __future__ import annotations
import math
from typing import Dict, List, Optional
import polars as pl

TAB_ID    = 'ion_sampling'
TAB_LABEL = 'Ion Sampling'


def is_visible(processor) -> bool:
    return ('Precursor.Quantity' in processor.report.columns or
            getattr(processor, 'fill_times', None) is not None or
            'Intensity' in processor.report.columns)


def compute(processor) -> Dict:
    df      = processor.report
    runs    = processor.runs()
    aa_keep = _enzyme_cterm_aas(getattr(processor, 'enzymes', None))
    result  = {
        'ms1_intensity':  _ms1_intensity(df, runs, aa_keep),
        'prec_data':      _prec_data(df, runs, 'Intensity'),
        'prec_aa_map':    _prec_aa_map(df),
        'ms2_intensity':  _ms2_intensity(df, runs, aa_keep),
        'ms2_prec_data':  _prec_data(df, runs, 'Precursor.Quantity'),
        'ms2_vs_mz':      _ms2_vs_mz(df, runs),
        'frag_mz_map':    _frag_mz_map(df),
        'ms1_fill_dist':  None,
        'ms1_fill_grad':  None,
        'ms2_fill_dist':  None,
        'ms2_fill_grad':  None,
    }
    ft = getattr(processor, 'fill_times', None)
    if isinstance(ft, pl.DataFrame) and len(ft) > 0:
        result['ms1_fill_dist'] = _fill_dist(ft, runs, level=1)
        result['ms1_fill_grad'] = _fill_grad(ft, runs, level=1)
        result['ms2_fill_dist'] = _fill_dist(ft, runs, level=2)
        result['ms2_fill_grad'] = _fill_grad(ft, runs, level=2)
    return result


def _enzyme_cterm_aas(enzymes):
    if not enzymes:
        return None
    from Mods.enzyme import ENZYMES as _E
    result = set()
    for e in enzymes:
        result.update(_E.get(e, {}).get('cterm', set()))
    return frozenset(result) if result else None


# ─── Shared low-level helpers ─────────────────────────────────────────────────

def _median(sv: list):
    n = len(sv)
    if n == 0:
        return None
    mid = n // 2
    return round(sv[mid] if n % 2 else (sv[mid - 1] + sv[mid]) / 2, 4)


def _hist_params(vals_sorted: list, n_bins: int):
    gmin = vals_sorted[0]
    gmax = vals_sorted[-1]
    bsize = (gmax - gmin) / n_bins if gmax > gmin else 1.0
    centers = [round(gmin + (i + 0.5) * bsize, 4) for i in range(n_bins)]
    return gmin, gmax, bsize, centers


def _bin(vals: list, gmin: float, gmax: float, bsize: float, n_bins: int) -> list:
    counts = [0] * n_bins
    for v in vals:
        if v < gmin or v > gmax:
            continue
        counts[min(n_bins - 1, int((v - gmin) / bsize))] += 1
    return counts


def _aa_intensity_section(df_run, col, ordered_aas, gmin, gmax, bsize, n_bins):
    """Build enzyme-aware per-run dict for one run's intensity column."""
    aa_series, median_per_aa = {}, {}
    for aa in ordered_aas:
        lv = [math.log10(float(v)) for v in
              df_run.filter(pl.col('C_term_aa') == aa)[col].to_list() if v and float(v) > 0]
        aa_series[aa]     = _bin(lv, gmin, gmax, bsize, n_bins)
        median_per_aa[aa] = _median(sorted(lv))
    max_count = max((max(v) for v in aa_series.values() if v), default=1)
    return {'aa_series': aa_series, 'median_per_aa': median_per_aa, 'max_count': max_count}


# ─── Section 1: MS1 Intensity Distribution ────────────────────────────────────

def _ms1_intensity(df: pl.DataFrame, runs: List[str], aa_keep=None) -> Dict:
    if 'Intensity' not in df.columns:
        return {'data': [], 'aas': None, 'centers': [], 'global_range': [0, 10]}

    all_log = []
    for run in runs:
        vals = (df.filter((pl.col('Raw file') == run) & (pl.col('Intensity') > 0))
                  .select('Intensity').drop_nulls()['Intensity'].to_list())
        all_log.extend(math.log10(float(v)) for v in vals if v > 0)

    if not all_log:
        return {'data': [], 'aas': None, 'centers': [], 'global_range': [0, 10]}

    n_bins = 50
    gmin, gmax, bsize, centers = _hist_params(sorted(all_log), n_bins)

    if aa_keep and 'C_term_aa' in df.columns:
        all_aas = sorted(df['C_term_aa'].drop_nulls().unique().to_list())
        all_aas = [a for a in all_aas if a in aa_keep]
        priority = ['R', 'K']
        ordered_aas = [a for a in priority if a in all_aas] + sorted(a for a in all_aas if a not in priority)
    else:
        ordered_aas = None

    data = []
    for run in runs:
        rdf = (df.filter((pl.col('Raw file') == run) & (pl.col('Intensity') > 0))
                 .drop_nulls('Intensity'))
        if ordered_aas:
            entry = _aa_intensity_section(rdf, 'Intensity', ordered_aas, gmin, gmax, bsize, n_bins)
            data.append({'run': run, **entry})
        else:
            lv     = [math.log10(float(v)) for v in rdf['Intensity'].to_list() if v > 0]
            counts = _bin(lv, gmin, gmax, bsize, n_bins)
            data.append({'run': run, 'counts': counts, 'max_count': max(counts) if counts else 1,
                         'median': _median(sorted(lv))})

    return {'data': data, 'aas': ordered_aas, 'centers': centers, 'global_range': [gmin, gmax]}


# ─── Section 6: MS2 Intensity Distribution ────────────────────────────────────

def _ms2_intensity(df: pl.DataFrame, runs: List[str], aa_keep=None) -> Dict:
    col = 'Precursor.Quantity'
    if col not in df.columns:
        return {'data': [], 'aas': None, 'centers': [], 'global_range': [0, 10]}

    all_log = []
    for run in runs:
        vals = (df.filter((pl.col('Raw file') == run) & (pl.col(col) > 0))
                  .select(col).drop_nulls()[col].to_list())
        all_log.extend(math.log10(float(v)) for v in vals if v > 0)

    if not all_log:
        return {'data': [], 'aas': None, 'centers': [], 'global_range': [0, 10]}

    n_bins = 50
    gmin, gmax, bsize, centers = _hist_params(sorted(all_log), n_bins)

    if aa_keep and 'C_term_aa' in df.columns:
        all_aas = sorted(df['C_term_aa'].drop_nulls().unique().to_list())
        all_aas = [a for a in all_aas if a in aa_keep]
        priority = ['R', 'K']
        ordered_aas = [a for a in priority if a in all_aas] + sorted(a for a in all_aas if a not in priority)
    else:
        ordered_aas = None

    data = []
    for run in runs:
        rdf = (df.filter((pl.col('Raw file') == run) & (pl.col(col) > 0))
                 .drop_nulls(col))
        if rdf.is_empty():
            continue  # no MS2 data for this run (e.g., Jmod) — omit so JS skips it
        if ordered_aas and 'C_term_aa' in rdf.columns:
            entry = _aa_intensity_section(rdf, col, ordered_aas, gmin, gmax, bsize, n_bins)
            data.append({'run': run, **entry})
        else:
            lv = [math.log10(float(v)) for v in rdf[col].to_list() if v and float(v) > 0]
            if not lv:
                continue
            counts = _bin(lv, gmin, gmax, bsize, n_bins)
            data.append({'run': run, 'counts': counts,
                         'max_count': max(counts) if counts else 1,
                         'median': _median(sorted(lv))})

    return {'data': data, 'aas': ordered_aas, 'centers': centers, 'global_range': [gmin, gmax]}


# ─── Precursor data for live JS intersection ──────────────────────────────────

def _prec_data(df: pl.DataFrame, runs: List[str], col: str = 'Intensity') -> Dict:
    """Per-run {seq_charge_key: log10_value} for live JavaScript intersection.

    Uses Sequence + '_' + Charge as the key so cross-engine intersection works
    (DIA-NN and Jmod use different Precursor.Id formats).
    """
    if col not in df.columns:
        return {}
    if 'Sequence' not in df.columns or 'Charge' not in df.columns:
        return {}
    df2 = df.with_columns(
        (pl.col('Sequence') + '_' + pl.col('Charge').cast(pl.String)).alias('_prec_key')
    )
    result = {}
    for run in runs:
        rdf = (df2.filter((pl.col('Raw file') == run) & (pl.col(col) > 0))
                   .select(['_prec_key', col]).drop_nulls()
                   .group_by('_prec_key').agg(pl.col(col).mean()))
        prec_map = {}
        for pid, intensity in zip(rdf['_prec_key'].to_list(), rdf[col].to_list()):
            if intensity and float(intensity) > 0:
                prec_map[pid] = round(math.log10(float(intensity)), 3)
        result[run] = prec_map
    return result


def _prec_aa_map(df: pl.DataFrame) -> Dict:
    """Global {seq_charge_key: C_term_aa} for enzyme-aware intersected charts."""
    if 'C_term_aa' not in df.columns or 'Sequence' not in df.columns or 'Charge' not in df.columns:
        return {}
    result = (df.with_columns(
                  (pl.col('Sequence') + '_' + pl.col('Charge').cast(pl.String)).alias('_prec_key')
              )
               .select(['_prec_key', 'C_term_aa'])
               .drop_nulls()
               .group_by('_prec_key')
               .agg(pl.col('C_term_aa').first().alias('aa')))
    return {pid: aa
            for pid, aa in zip(result['_prec_key'].to_list(), result['aa'].to_list())
            if aa is not None}


# ─── Section 12: MS2 Intensity vs Precursor M/Z (per-run grid) ───────────────

def _ms2_vs_mz(df: pl.DataFrame, runs: List[str]) -> Dict:
    """Per-run, bin Precursor.Mz in 25 Da bins, median log10(Precursor.Quantity)."""
    if 'Precursor.Mz' not in df.columns or 'Precursor.Quantity' not in df.columns:
        return {}

    BIN_SIZE = 25.0
    valid = df.filter((pl.col('Precursor.Quantity') > 0) & pl.col('Precursor.Mz').is_not_null())
    if valid.is_empty():
        return {}

    mz_min = float(valid['Precursor.Mz'].min())
    mz_max = float(valid['Precursor.Mz'].max())
    n_bins = max(1, int((mz_max - mz_min) / BIN_SIZE) + 1)
    centers = [round(mz_min + (i + 0.5) * BIN_SIZE, 1) for i in range(n_bins)]

    all_meds: List[float] = []
    data = []
    for run in runs:
        rdf = valid.filter(pl.col('Raw file') == run).select(['Precursor.Mz', 'Precursor.Quantity'])
        if rdf.is_empty():
            continue  # no MS2 data for this run (e.g. Jmod) — omit so JS creates no chart box

        mzs = rdf['Precursor.Mz'].to_list()
        pqs = rdf['Precursor.Quantity'].to_list()
        bins: List[List[float]] = [[] for _ in range(n_bins)]
        for mz, pq in zip(mzs, pqs):
            if pq and float(pq) > 0:
                bi = min(n_bins - 1, int((float(mz) - mz_min) / BIN_SIZE))
                if bi >= 0:
                    bins[bi].append(math.log10(float(pq)))

        vals = []
        for bv in bins:
            med = _median(sorted(bv)) if bv else None
            vals.append(med)
            if med is not None:
                all_meds.append(med)
        data.append({'run': run, 'values': vals})

    if not all_meds:
        return {}
    return {
        'data':         data,
        'centers':      centers,
        'mz_range':     [round(mz_min, 1), round(mz_max, 1)],
        'global_range': [round(min(all_meds), 3), round(max(all_meds), 3)],
    }


# ─── Section 13+14: Fragment M/Z map (Best.Fr.Mz per precursor) ──────────────

def _frag_mz_map(df: pl.DataFrame) -> Dict:
    """Global median Best.Fr.Mz per Precursor.Id for live JS fragment-m/z charts."""
    if 'Best.Fr.Mz' not in df.columns or 'Precursor.Id' not in df.columns:
        return {}
    result = (df.filter(pl.col('Best.Fr.Mz').is_not_null() & (pl.col('Best.Fr.Mz') > 0))
               .select(['Precursor.Id', 'Best.Fr.Mz'])
               .group_by('Precursor.Id')
               .agg(pl.col('Best.Fr.Mz').median().alias('mz')))
    return {pid: round(float(mz), 1)
            for pid, mz in zip(result['Precursor.Id'].to_list(), result['mz'].to_list())
            if mz is not None}


# ─── Fill time sections (MS1 level=1, MS2 level=2) ───────────────────────────

def _fill_dist(ft: pl.DataFrame, runs: List[str], level: int = 1) -> Optional[Dict]:
    run_col = next((c for c in ['Raw.file', 'Raw file', 'Filename', 'Run'] if c in ft.columns), None)
    if run_col is None or 'Fill.Time' not in ft.columns:
        return None

    ft_lvl = ft.filter(pl.col('Ms.Level') == level) if 'Ms.Level' in ft.columns else ft
    if ft_lvl.is_empty():
        return None

    all_vals = sorted(ft_lvl.select('Fill.Time').drop_nulls()['Fill.Time'].cast(pl.Float64).to_list())
    if not all_vals:
        return None

    n_bins = 15
    gmin, gmax, bsize, centers = _hist_params(all_vals, n_bins)

    data = []
    for run in runs:
        vals = (ft_lvl.filter(pl.col(run_col) == run)
                      .select('Fill.Time').drop_nulls()['Fill.Time'].cast(pl.Float64).to_list())
        if not vals:
            continue
        counts = _bin(vals, gmin, gmax, bsize, n_bins)
        data.append({'run': run, 'counts': counts, 'max_count': max(counts) if counts else 1})

    return ({'data': data, 'centers': centers, 'global_range': [gmin, gmax]} if data else None)


def _fill_grad(ft: pl.DataFrame, runs: List[str], level: int = 1) -> Optional[Dict]:
    run_col = next((c for c in ['Raw.file', 'Raw file', 'Filename', 'Run'] if c in ft.columns), None)
    if run_col is None or 'Fill.Time' not in ft.columns or 'RT.Start' not in ft.columns:
        return None

    ft_lvl = (ft.filter(pl.col('Ms.Level') == level) if 'Ms.Level' in ft.columns else ft).with_columns(
        pl.col('Fill.Time').cast(pl.Float64),
        pl.col('RT.Start').cast(pl.Float64),
    )
    if ft_lvl.is_empty():
        return None

    rt_all = ft_lvl['RT.Start'].drop_nulls().to_list()
    ft_all = ft_lvl['Fill.Time'].drop_nulls().to_list()
    if not rt_all:
        return None

    rt_min, rt_max = min(rt_all), max(rt_all)
    ft_min, ft_max = min(ft_all), max(ft_all)
    bin_size = (rt_max - rt_min) / 15 if rt_max > rt_min else 1.0

    data = []
    for run in runs:
        rft = ft_lvl.filter(pl.col(run_col) == run)
        if rft.is_empty():
            continue
        bins = []
        for i in range(15):
            rt_lo = rt_min + i * bin_size
            vals = (rft.filter((pl.col('RT.Start') >= rt_lo) &
                                (pl.col('RT.Start') < rt_lo + bin_size))
                       ['Fill.Time'].drop_nulls().to_list())
            if not vals:
                continue
            mean_v = sum(vals) / len(vals)
            sd_v   = (sum((v - mean_v) ** 2 for v in vals) / len(vals)) ** 0.5
            bins.append({'rt_center': round(rt_lo + bin_size / 2, 3),
                         'mean': round(mean_v, 3), 'sd': round(sd_v, 3)})
        if bins:
            data.append({'run': run, 'bins': bins})

    if not data:
        return None
    return {
        'data':     data,
        'rt_range': [rt_min, rt_max],
        'ft_range': [max(0.0, ft_min - 1.0), ft_max + 1.0],
    }


# ─── HTML ─────────────────────────────────────────────────────────────────────

def html_content(data: Dict) -> str:
    n_runs = len((data.get('ms1_intensity') or {}).get('data') or [])
    default_cols = min(max(1, n_runs), 6)
    cols_opts = ''.join(
        f'<option value="{n}"{" selected" if n == default_cols else ""}>{n}</option>'
        for n in range(1, 7)
    )

    def _grid_section(h3, badge_id, help_txt, grid_id, msg_id, default_cols=default_cols):
        return f"""
<div class="section">
  <h3>{h3} <span id="{badge_id}" style="font-weight:normal;font-size:.85em;color:#555"></span></h3>
  <p class="help">{help_txt}</p>
  <div id="{grid_id}" class="chart-grid" style="grid-template-columns:repeat({default_cols},1fr)"></div>
  <p id="{msg_id}" style="display:none;color:#888;font-style:italic;text-align:center;padding:20px 0"></p>
</div>"""

    def _fill_section(h3, help_txt, grid_id):
        return f"""
<div class="section">
  <h3>{h3}</h3>
  <p class="help">{help_txt}</p>
  <div id="{grid_id}" class="chart-grid" style="grid-template-columns:repeat({default_cols},1fr)"></div>
</div>"""

    # MS1 fill-time sections (conditional)
    sec_ms1_fdt = ''
    if (data.get('ms1_fill_dist') or {}).get('data'):
        sec_ms1_fdt = _fill_section(
            'MS1 Fill Time Distribution',
            'Distribution of MS1 ion accumulation times per run. '
            'Fill times close to the instrument maximum indicate ion-limited sampling.',
            'ion-fdt-grid')
    sec_ms1_fgr = ''
    if (data.get('ms1_fill_grad') or {}).get('data'):
        sec_ms1_fgr = _fill_section(
            'MS1 Fill Times along Gradient',
            'Mean MS1 fill time across 15 retention-time bins (band = ±½ SD). '
            'A systematic rise late in the gradient may indicate increasing sample complexity.',
            'ion-fgr-grid')

    # MS2 fill-time sections (conditional)
    sec_ms2_fdt = ''
    if (data.get('ms2_fill_dist') or {}).get('data'):
        sec_ms2_fdt = _fill_section(
            'MS2 Fill Time Distribution',
            'Distribution of MS2 ion accumulation times per run. '
            'Fill times near the instrument maximum indicate DIA windows are ion-limited.',
            'ion-fd2-grid')
    sec_ms2_fgr = ''
    if (data.get('ms2_fill_grad') or {}).get('data'):
        sec_ms2_fgr = _fill_section(
            'MS2 Fill Times along Gradient',
            'Mean MS2 fill time across 15 retention-time bins (band = ±½ SD). '
            'A rise mid-gradient reflects the apex of peptide elution.',
            'ion-fg2-grid')

    return f"""
<div class="section">
  <div class="run-selector">
    <label>Samples per row:</label>
    <select id="ion-cols-sel" onchange="updateIonLayout()">
      {cols_opts}
    </select>
  </div>
</div>

<div class="section">
  <h3>MS1 Intensity Distribution</h3>
  <p class="help">Log₁₀ MS1 precursor intensity per run (full range, 50 bins). Peaks should overlap across runs for consistent ionisation efficiency.</p>
  <div id="ion-ms1-grid" class="chart-grid" style="grid-template-columns:repeat({default_cols},1fr)"></div>
</div>

{_grid_section(
    'MS1 Intensity – Intersected Precursors', 'ion-int-badge',
    'Log₁₀ MS1 intensity for precursors detected in every active run, coloured by C-terminal residue when an enzyme is selected. '
    'Updates when runs are hidden/shown in the Files tab. Requires ≥2 active runs with ≥20 shared precursors.',
    'ion-int-grid', 'ion-int-msg')}

{_grid_section(
    'Normalized MS1 Intensity – Intersected Precursors', 'ion-nrm-badge',
    'Log₂ fold-change of MS1 intensity vs. the first active run, for intersected precursors. '
    'Axis fixed [−3, +3]. A narrow peak at 0 indicates stable ionisation. Updates live with run selection.',
    'ion-nrm-grid', 'ion-nrm-msg')}

{sec_ms1_fdt}{sec_ms1_fgr}

<div class="section">
  <h3>MS2 Intensity Distribution</h3>
  <p class="help">Log₁₀ MS2 precursor quantity per run (full range, 50 bins), coloured by C-terminal residue when an enzyme is selected.</p>
  <div id="ion-ms2-grid" class="chart-grid" style="grid-template-columns:repeat({default_cols},1fr)"></div>
</div>

{_grid_section(
    'MS2 Intensity – Intersected Precursors', 'ion-i2t-badge',
    'Log₁₀ MS2 precursor quantity for precursors detected in every active run, coloured by C-terminal residue. '
    'Updates live with run selection. Requires ≥2 active runs with ≥20 shared precursors.',
    'ion-i2t-grid', 'ion-i2t-msg')}

{_grid_section(
    'Normalized MS2 Intensity – Intersected Precursors', 'ion-n2m-badge',
    'Log₂ fold-change of MS2 precursor quantity vs. the first active run, for intersected precursors. '
    'Axis fixed [−3, +3]. Updates live with run selection.',
    'ion-n2m-grid', 'ion-n2m-msg')}

{sec_ms2_fdt}{sec_ms2_fgr}

{_grid_section(
    'MS2 Intensity vs Precursor M/Z', 'ion-mzc-badge',
    'Median log₁₀ MS2 intensity (Precursor.Quantity) in 25-Da precursor m/z bins per run. '
    'Y-axis = precursor m/z; X-axis = log₁₀ intensity. Systematic shifts across the m/z range '
    'may reflect DIA window boundaries or instrument bias.',
    'ion-mz2-grid', 'ion-mzc-msg')}

{_grid_section(
    'Best Fragment M/Z vs MS2 Intensity – Intersected Precursors', 'ion-fmz-badge',
    'Median log₁₀ MS2 precursor quantity in 5-Da bins of best-fragment m/z for precursors detected in every active run. '
    '"Best fragment" (Best.Fr.Mz) = the single most abundant fragment ion per precursor as reported by DIA-NN. '
    'Note: the DIA-NN parquet does not contain Fragment.Quant.Raw / Fragment.Info, so only the dominant fragment '
    'is available here — all other fragment ions are not used. '
    'Y-axis = fragment m/z; X-axis = log₁₀ precursor quantity. '
    'Requires ≥2 active runs with ≥20 shared precursors.',
    'ion-fmz-grid', 'ion-fmz-msg')}

{_grid_section(
    'Fragment M/Z Normalized to First Run – Intersected Precursors', 'ion-fnm-badge',
    'Log₂ fold-change of MS2 precursor quantity vs. the first active run per 5-Da best-fragment m/z bin. '
    'Uses Best.Fr.Mz (dominant fragment only — full fragment-ion data not available in the parquet). '
    'Y-axis = fragment m/z; X-axis = log₂ FC. A flat line at 0 indicates uniform intensity across the m/z range. '
    'Requires ≥2 active runs with ≥20 shared precursors.',
    'ion-fnm-grid', 'ion-fnm-msg')}

{_grid_section(
    'MS2/MS1 Intensity Ratio – Intersected Precursors', 'ion-rto-badge',
    'Distribution of log₂(Precursor.Quantity / MS1 Intensity) for precursors detected in every active run '
    'with both MS1 and MS2 quantification. A narrow, stable distribution indicates consistent fragmentation '
    'efficiency. Updates live with run selection.',
    'ion-rto-grid', 'ion-rto-msg')}"""


# ─── JavaScript ───────────────────────────────────────────────────────────────

def javascript() -> str:
    return r"""
// _ionChart: distribution chart (x=count, y=value) reused for all ion sections.
// yTitle: y-axis label. yStepSize: override y tick step (use 1 for log2 FC).
function _ionChart(canvasId, datasets, gMin, gMax, maxCount, yTitle, yStepSize) {
  if (yTitle === undefined) yTitle = 'log₁₀(Intensity)';
  destroyChart(canvasId);
  const ctx = document.getElementById(canvasId);
  if (!ctx) return;
  const {step: xStep, axMax: xMax} = _niceTicks(maxCount);
  const yMin  = Math.floor(gMin);
  const yMax  = Math.ceil(gMax);
  const yStep = yStepSize !== undefined ? yStepSize : Math.max(1, Math.ceil((yMax - yMin) / 5));
  const realDs = datasets.filter(d => !d.label.startsWith('_med'));
  const medLabelPlugin = {
    id: 'medLabel',
    afterDraw(chart) {
      const c = chart.ctx;
      const yScale = chart.scales.y;
      chart.data.datasets.forEach((ds, i) => {
        if (!ds.label.startsWith('_med')) return;
        if (chart.getDatasetMeta(i).hidden) return;
        const val = ds.data[0].y;
        if (val < yScale.min || val > yScale.max) return;
        const y = yScale.getPixelForValue(val);
        c.save();
        c.fillStyle = ds.borderColor;
        c.font = 'bold 9px sans-serif';
        c.textAlign = 'left';
        c.fillText(val.toFixed(2), chart.chartArea.left + 3, y - 3);
        c.restore();
      });
    }
  };
  CHARTS[canvasId] = new Chart(ctx, {
    type: 'line',
    plugins: [medLabelPlugin],
    data: { datasets },
    options: {
      responsive: true,
      maintainAspectRatio: false,
      parsing: false,
      layout: { padding: { bottom: 8, right: 12 } },
      plugins: {
        legend: {
          display: realDs.length > 1, position: 'top',
          labels: { boxWidth: 12, font: { size: 10 },
                    filter: (item) => !item.text.startsWith('_med') },
          onClick: (evt, item, legend) => {
            Chart.defaults.plugins.legend.onClick(evt, item, legend);
            const chart = legend.chart;
            const ds = chart.data.datasets[item.datasetIndex];
            const medIdx = chart.data.datasets.findIndex(d => d.label === '_med_' + ds.label);
            if (medIdx >= 0) {
              chart.getDatasetMeta(medIdx).hidden = chart.getDatasetMeta(item.datasetIndex).hidden;
              chart.update();
            }
          }
        },
        tooltip: { filter: (item) => !item.dataset.label.startsWith('_med') },
      },
      elements: { line: { borderWidth: 1.5 } },
      scales: {
        x: { type: 'linear', min: 0, max: xMax,
             title: { display: true, text: 'Count' },
             ticks: { stepSize: xStep, precision: 0, maxRotation: 0, minRotation: 0,
                      font: { size: 9 } } },
        y: { type: 'linear', min: yMin, max: yMax,
             title: { display: true, text: yTitle },
             ticks: { stepSize: yStep, precision: 1 } },
      }
    }
  });
}

// Per-run M/Z profile chart: Y = m/z bins, X = log intensity (or log2 FC).
// rd.values[i] = intensity value for m/z bin i (null = no data in that bin).
// opts.showDots: show small dots at each bin point (like _ionChart style).
// opts.noMedian: suppress the vertical median dashed line.
// opts.medLabelRight: place the median value label to the right of the dashed bar.
function _mzProfileChart(canvasId, rd, mzCenters, xRange, yRange, xTitle, opts) {
  opts = opts || {};
  destroyChart(canvasId);
  const ctx = document.getElementById(canvasId);
  if (!ctx) return;
  const [xMin, xMax] = xRange;
  const [yMin, yMax] = yRange;

  const pts = mzCenters.map((c, i) =>
    rd.values[i] != null ? {x: rd.values[i], y: c} : null
  ).filter(p => p !== null);
  if (!pts.length) return;

  const xStep  = Math.max(0.2, Math.ceil((xMax - xMin) / 5 * 5) / 5);
  const mzStep = Math.ceil((yMax - yMin) / 4 / 50) * 50 || 100;

  const plugins = [];
  if (!opts.noMedian) {
    const sortedX = pts.map(p => p.x).sort((a, b) => a - b);
    const n = sortedX.length;
    const med = n % 2 ? sortedX[Math.floor(n / 2)] : (sortedX[n / 2 - 1] + sortedX[n / 2]) / 2;
    plugins.push({
      id: 'mzMed',
      afterDraw(chart) {
        const xScale = chart.scales.x;
        if (med < xScale.min || med > xScale.max) return;
        const xPx = xScale.getPixelForValue(med);
        const c = chart.ctx;
        c.save();
        c.strokeStyle = color(rd.run);
        c.setLineDash([5, 5]);
        c.lineWidth = 1;
        c.beginPath();
        c.moveTo(xPx, chart.chartArea.top);
        c.lineTo(xPx, chart.chartArea.bottom);
        c.stroke();
        c.fillStyle = color(rd.run);
        c.font = 'bold 9px sans-serif';
        c.textAlign = opts.medLabelRight ? 'left' : 'center';
        c.fillText(med.toFixed(2), opts.medLabelRight ? xPx + 4 : xPx, chart.chartArea.top + 11);
        c.restore();
      }
    });
  }

  const showDots = !!opts.showDots;
  CHARTS[canvasId] = new Chart(ctx, {
    type: 'line',
    plugins,
    data: { datasets: [{
      label: label(rd.run),
      data: pts,
      borderColor: color(rd.run),
      backgroundColor: 'transparent',
      borderWidth: 1.5,
      pointRadius: showDots ? 2.5 : 0,
      pointHoverRadius: 4,
      showLine: true,
    }]},
    options: {
      responsive: true,
      maintainAspectRatio: false,
      parsing: false,
      layout: { padding: { bottom: 8, right: 16, top: 4 } },
      plugins: {
        legend: { display: false },
        tooltip: { callbacks: {
          title: (items) => `M/Z: ${items[0].parsed.y.toFixed(0)}`,
          label: (item)  => `${xTitle || 'Value'}: ${item.parsed.x.toFixed(3)}`
        }}
      },
      elements: { line: { borderWidth: 1.5, tension: 0.15 } },
      scales: {
        x: { type: 'linear', min: xMin, max: xMax,
             title: { display: true, text: xTitle || 'Value' },
             ticks: { stepSize: xStep, precision: 2, maxRotation: 0, minRotation: 0,
                      font: { size: 9 } } },
        y: { type: 'linear', min: Math.floor(yMin / mzStep) * mzStep,
             max: Math.ceil(yMax / mzStep) * mzStep,
             title: { display: true, text: 'Precursor M/Z' },
             ticks: { stepSize: mzStep, precision: 0 } },
      }
    }
  });
}

// Fill times along gradient: mean line + ±½SD band
function _fillGradChart(canvasId, bins, rtRange, ftRange) {
  destroyChart(canvasId);
  const ctx = document.getElementById(canvasId);
  if (!ctx || !bins || !bins.length) return;
  const [rtMin, rtMax] = rtRange;
  const [ftMin, ftMax] = ftRange;
  CHARTS[canvasId] = new Chart(ctx, {
    type: 'line',
    data: { datasets: [
      { label: '_band_hi',
        data: bins.map(b => ({x: b.rt_center, y: b.mean + b.sd / 2})),
        fill: '+1', backgroundColor: 'rgba(52,152,219,0.2)',
        borderColor: 'transparent', borderWidth: 0, pointRadius: 0, pointHoverRadius: 0 },
      { label: '_band_lo',
        data: bins.map(b => ({x: b.rt_center, y: Math.max(0, b.mean - b.sd / 2)})),
        fill: false, borderColor: 'transparent', borderWidth: 0, pointRadius: 0, pointHoverRadius: 0 },
      { label: 'Mean',
        data: bins.map(b => ({x: b.rt_center, y: b.mean})),
        fill: false, borderColor: '#3498db', backgroundColor: 'transparent',
        borderWidth: 1.5, pointRadius: 2.5, pointHoverRadius: 4 }
    ]},
    options: {
      responsive: true,
      maintainAspectRatio: false,
      parsing: false,
      layout: { padding: { bottom: 8, right: 12 } },
      plugins: {
        legend: { display: false },
        tooltip: {
          filter: (item) => item.dataset.label === 'Mean',
          callbacks: {
            title: (items) => `RT: ${items[0].parsed.x.toFixed(1)} min`,
            label:  (item) => `Fill time: ${item.parsed.y.toFixed(1)} ms`
          }
        }
      },
      elements: { line: { borderWidth: 1.5, tension: 0.3 }, point: { radius: 0 } },
      scales: {
        x: { type: 'linear', min: rtMin, max: rtMax,
             title: { display: true, text: 'Retention Time (min)' },
             ticks: { maxRotation: 0, minRotation: 0, font: { size: 9 } } },
        y: { type: 'linear', min: Math.max(0, ftMin), max: ftMax,
             title: { display: true, text: 'Fill Time (ms)' },
             ticks: { precision: 0 } }
      }
    }
  });
}

function _ionDatasets(rd, aas, centers, run) {
  if (aas && aas.length > 0) {
    return aas.flatMap((aa, i) => {
      const sets = [{
        label: aa,
        data: centers.map((c, j) => ({x: rd.aa_series[aa][j] || 0, y: c})),
        borderColor: aaColor(aa, i), backgroundColor: 'transparent',
        borderWidth: 1.5, pointRadius: 2.5, pointHoverRadius: 4,
      }];
      if (rd.median_per_aa && rd.median_per_aa[aa] != null) sets.push({
        label: '_med_' + aa,
        data: [{x: 0, y: rd.median_per_aa[aa]}, {x: 1e9, y: rd.median_per_aa[aa]}],
        borderColor: aaColor(aa, i), backgroundColor: 'transparent',
        borderWidth: 1, borderDash: [5, 5], pointRadius: 0, pointHoverRadius: 0,
      });
      return sets;
    });
  }
  const sets = [{
    label: label(run),
    data: centers.map((c, i) => ({x: rd.counts[i], y: c})),
    borderColor: color(run), backgroundColor: 'transparent',
    borderWidth: 1.5, pointRadius: 2.5, pointHoverRadius: 4,
  }];
  if (rd.median != null) sets.push({
    label: '_med',
    data: [{x: 0, y: rd.median}, {x: 1e9, y: rd.median}],
    borderColor: color(run), backgroundColor: 'transparent',
    borderWidth: 1, borderDash: [5, 5], pointRadius: 0, pointHoverRadius: 0,
  });
  return sets;
}

// Build per-run DOM cells (canvas + export btn) for a chart grid
function _buildIonGridDOM(gridId, prefix, runsData, activeRuns) {
  const grid = document.getElementById(gridId);
  if (!grid || !runsData || !runsData.length) return;
  grid.innerHTML = '';
  activeRuns.forEach((run, idx) => {
    const rd = runsData.find(x => x.run === run);
    if (!rd) return;
    const canvasId = `${prefix}_${idx}`;
    const box = document.createElement('div');
    box.className = 'chart-box';
    box.innerHTML =
      `<h4>${label(run)}</h4>` +
      `<button class="export-btn" onclick="exportChart('${canvasId}')">PNG</button>` +
      `<div style="position:relative;height:375px"><canvas id="${canvasId}"></canvas></div>`;
    grid.appendChild(box);
  });
}

// Generic intersected precursor computation + rendering.
// cfg: { dataKey, intGridId, nrmGridId, intBadgeId, nrmBadgeId,
//         intMsgId, nrmMsgId, intPrefix, nrmPrefix, yTitle, aas }
// cfg.aas: array of AA letters to split by (enzyme-aware), or null/undefined for combined.
function _computeAndRenderIntersected(ar, cfg) {
  const precData    = DATA.ion_sampling[cfg.dataKey] || {};
  const precAaMap   = DATA.ion_sampling.prec_aa_map  || {};
  const aas         = cfg.aas || null;
  const intGrid     = document.getElementById(cfg.intGridId);
  const intMsg      = document.getElementById(cfg.intMsgId);
  const intBadge    = document.getElementById(cfg.intBadgeId);
  const nrmGrid     = document.getElementById(cfg.nrmGridId);
  const nrmMsg      = document.getElementById(cfg.nrmMsgId);
  const nrmBadge    = document.getElementById(cfg.nrmBadgeId);

  function showMsg(grid, msg, badge, msgEl, n) {
    if (grid)  grid.innerHTML = '';
    if (msgEl) { msgEl.textContent = msg; msgEl.style.display = 'block'; }
    if (badge) badge.textContent = n > 0 ? `(${n.toLocaleString()} shared)` : '';
  }
  function hideMsg(msgEl) {
    if (msgEl) { msgEl.textContent = ''; msgEl.style.display = 'none'; }
  }

  if (ar.length < 2) {
    const msg = 'Select at least 2 active runs to compute intersected precursors.';
    showMsg(intGrid, msg, intBadge, intMsg, 0);
    showMsg(nrmGrid, msg, nrmBadge, nrmMsg, 0);
    return;
  }

  const runMaps = {};
  ar.forEach(run => { runMaps[run] = new Map(Object.entries(precData[run] || {})); });

  let common = new Set(runMaps[ar[0]].keys());
  for (let i = 1; i < ar.length; i++) {
    const m = runMaps[ar[i]];
    common = new Set([...common].filter(p => m.has(p)));
  }
  const n_common = common.size;

  if (n_common < 20) {
    const msg = n_common === 0
      ? 'No shared precursors across the selected runs.'
      : `Only ${n_common.toLocaleString()} shared precursors found (minimum 20 required).`;
    showMsg(intGrid, msg, intBadge, intMsg, n_common);
    showMsg(nrmGrid, msg, nrmBadge, nrmMsg, n_common);
    return;
  }

  hideMsg(intMsg);
  hideMsg(nrmMsg);
  if (intBadge) intBadge.textContent = `(${n_common.toLocaleString()} shared)`;

  const N_BINS    = 50;
  const commonArr = [...common];

  function medianOf(sv) {
    const m = sv.length;
    if (!m) return null;
    return m % 2 ? sv[Math.floor(m / 2)] : (sv[m / 2 - 1] + sv[m / 2]) / 2;
  }
  function binVals(vals, lo, hi, bs, nBins) {
    nBins = nBins !== undefined ? nBins : N_BINS;
    const counts = new Array(nBins).fill(0);
    vals.forEach(v => {
      if (v < lo || v > hi) return;
      counts[Math.min(nBins - 1, Math.floor((v - lo) / bs))]++;
    });
    return counts;
  }

  // ── Global histogram range ────────────────────────────────────────────────
  const allLog = [];
  ar.forEach(run => { commonArr.forEach(p => allLog.push(runMaps[run].get(p))); });
  allLog.sort((a, b) => a - b);
  const n    = allLog.length;
  const gmin = allLog[Math.max(0, Math.floor(0.01 * n))];
  const gmax = allLog[Math.min(n - 1, Math.floor(0.99 * n))];
  const bsize   = (gmax - gmin) / N_BINS;
  const centers = Array.from({length: N_BINS}, (_, i) => gmin + (i + 0.5) * bsize);

  // ── Absolute histogram (enzyme-aware if cfg.aas) ─────────────────────────
  const intRunData = ar.map(run => {
    if (aas && aas.length > 0) {
      const aa_series = {}, median_per_aa = {};
      aas.forEach(aa => {
        const precs = commonArr.filter(p => precAaMap[p] === aa);
        const vals  = precs.map(p => runMaps[run].get(p)).sort((a, b) => a - b);
        aa_series[aa]     = binVals(vals, gmin, gmax, bsize);
        median_per_aa[aa] = medianOf(vals);
      });
      const max_count = Math.max(...Object.values(aa_series).flatMap(v => v), 1);
      return {run, aa_series, median_per_aa, max_count};
    }
    const vals   = commonArr.map(p => runMaps[run].get(p)).sort((a, b) => a - b);
    const counts = binVals(vals, gmin, gmax, bsize);
    return {run, counts, max_count: Math.max(...counts, 1), median: medianOf(vals)};
  });

  _buildIonGridDOM(cfg.intGridId, cfg.intPrefix, intRunData, ar);
  ar.forEach((run, idx) => {
    const rd = intRunData.find(x => x.run === run);
    if (rd) _ionChart(`${cfg.intPrefix}_${idx}`, _ionDatasets(rd, aas, centers, run),
                      gmin, gmax, rd.max_count, cfg.yTitle);
  });

  // ── Normalized (log2 FC vs first active run) — always combined ────────────
  // N_NRM = 51 (odd) ensures FC=0 lands exactly on a bin centre at [-3,+3]:
  //   centre of bin 25 = -3 + 25.5 × (6/51) = 0.00
  const LOG2_10  = 1 / Math.log10(2);
  const run1     = ar[0];
  const r1Map    = runMaps[run1];
  const N_NRM    = 51;
  const nBsize   = 6 / N_NRM;
  const nCenters = Array.from({length: N_NRM}, (_, i) => -3 + (i + 0.5) * nBsize);

  const nrmRunData = ar.map(run => {
    const fcVals = run === run1
      ? new Array(n_common).fill(0.0)
      : commonArr.map(p => (runMaps[run].get(p) - r1Map.get(p)) * LOG2_10);
    fcVals.sort((a, b) => a - b);
    const counts = binVals(fcVals, -3, 3, nBsize, N_NRM);
    return {run, counts, max_count: Math.max(...counts, 1), median: medianOf(fcVals)};
  });

  if (nrmBadge) nrmBadge.textContent = `(${n_common.toLocaleString()} shared · vs ${label(run1)})`;

  _buildIonGridDOM(cfg.nrmGridId, cfg.nrmPrefix, nrmRunData, ar);
  ar.forEach((run, idx) => {
    const rd = nrmRunData.find(x => x.run === run);
    if (rd) _ionChart(`${cfg.nrmPrefix}_${idx}`, _ionDatasets(rd, null, nCenters, run),
                      -3, 3, rd.max_count, 'log₂ FC (vs run 1)', 1);
  });
}

// Section 12: MS2 Intensity vs Precursor M/Z — per-run grid, precomputed data.
function _renderMs2VsMz(ar) {
  const mzcData = DATA.ion_sampling.ms2_vs_mz || {};
  const badge   = document.getElementById('ion-mzc-badge');
  const msg     = document.getElementById('ion-mzc-msg');

  if (!mzcData.data || !mzcData.data.length) {
    if (msg) { msg.textContent = 'Precursor.Mz or Precursor.Quantity data not available.'; msg.style.display = 'block'; }
    return;
  }
  if (msg) { msg.textContent = ''; msg.style.display = 'none'; }
  if (badge) badge.textContent = `(${ar.length} run${ar.length > 1 ? 's' : ''})`;

  const activeData = mzcData.data.filter(rd => ar.includes(rd.run));
  _buildIonGridDOM('ion-mz2-grid', 'ion_mz2', activeData, ar);

  const [xMin, xMax] = mzcData.global_range || [0, 10];
  const [yMin, yMax] = mzcData.mz_range     || [300, 1200];

  ar.forEach((run, idx) => {
    const rd = activeData.find(x => x.run === run);
    if (rd) _mzProfileChart(`ion_mz2_${idx}`, rd, mzcData.centers, [xMin, xMax], [yMin, yMax],
                             'log₁₀(Precursor Quantity)', {showDots: true, noMedian: true});
  });
}

// Sections 13 + 14: Best Fragment M/Z — per-run grids, live JS.
function _computeAndRenderFragMz(ar) {
  const fragMzMap   = DATA.ion_sampling.frag_mz_map   || {};
  const ms2PrecData = DATA.ion_sampling.ms2_prec_data  || {};
  const fmzGrid  = document.getElementById('ion-fmz-grid');
  const fnmGrid  = document.getElementById('ion-fnm-grid');
  const fmzBadge = document.getElementById('ion-fmz-badge');
  const fnmBadge = document.getElementById('ion-fnm-badge');
  const fmzMsg   = document.getElementById('ion-fmz-msg');
  const fnmMsg   = document.getElementById('ion-fnm-msg');

  function showMsg(grid, badge, msgEl, txt) {
    if (grid)  grid.innerHTML = '';
    if (badge) badge.textContent = '';
    if (msgEl) { msgEl.textContent = txt; msgEl.style.display = 'block'; }
  }
  function hideMsg(msgEl) {
    if (msgEl) { msgEl.textContent = ''; msgEl.style.display = 'none'; }
  }

  const hasMzMap = Object.keys(fragMzMap).length > 0;
  if (ar.length < 2 || !hasMzMap) {
    const txt = !hasMzMap
      ? 'Best.Fr.Mz data not available in this report.'
      : 'Select at least 2 active runs to compute intersected precursors.';
    showMsg(fmzGrid, fmzBadge, fmzMsg, txt);
    showMsg(fnmGrid, fnmBadge, fnmMsg, txt);
    return;
  }

  // Maps restricted to precursors with known fragment m/z
  const runMaps = {};
  ar.forEach(run => {
    const m = new Map();
    Object.entries(ms2PrecData[run] || {}).forEach(([pid, v]) => {
      if (pid in fragMzMap) m.set(pid, v);
    });
    runMaps[run] = m;
  });

  let common = new Set(runMaps[ar[0]].keys());
  for (let i = 1; i < ar.length; i++) {
    const m = runMaps[ar[i]];
    common = new Set([...common].filter(p => m.has(p)));
  }
  const n_common = common.size;

  if (n_common < 20) {
    const txt = n_common === 0
      ? 'No shared precursors with fragment m/z data across selected runs.'
      : `Only ${n_common.toLocaleString()} shared precursors with fragment m/z data (minimum 20 required).`;
    showMsg(fmzGrid, fmzBadge, fmzMsg, txt);
    showMsg(fnmGrid, fnmBadge, fnmMsg, txt);
    return;
  }

  hideMsg(fmzMsg);
  hideMsg(fnmMsg);
  if (fmzBadge) fmzBadge.textContent = `(${n_common.toLocaleString()} shared precursors)`;

  const commonArr = [...common];

  // M/Z range from common precursors (p1–p99)
  const allMz = commonArr.map(p => fragMzMap[p]).sort((a, b) => a - b);
  const mzMin = Math.floor(allMz[Math.floor(0.01 * allMz.length)] / 25) * 25;
  const mzMax = Math.ceil(allMz[Math.floor(0.99 * allMz.length)] / 25) * 25;
  const BIN_SIZE = 5;
  const nBins   = Math.max(1, Math.ceil((mzMax - mzMin) / BIN_SIZE));
  const centers = Array.from({length: nBins}, (_, i) => mzMin + (i + 0.5) * BIN_SIZE);

  function medianOf(sv) {
    const m = sv.length;
    if (!m) return null;
    return m % 2 ? sv[Math.floor(m / 2)] : (sv[m / 2 - 1] + sv[m / 2]) / 2;
  }

  // Per-run binned absolute data
  const runBinnedAbs = ar.map(run => {
    const bins = Array.from({length: nBins}, () => []);
    commonArr.forEach(p => {
      const mz = fragMzMap[p];
      if (mz == null || mz <= 0) return;
      const bi = Math.min(nBins - 1, Math.floor((mz - mzMin) / BIN_SIZE));
      if (bi >= 0) bins[bi].push(runMaps[run].get(p));
    });
    const values = bins.map(bv =>
      bv.length ? medianOf(bv.slice().sort((a, b) => a - b)) : null
    );
    return {run, values};
  });

  const allAbs = runBinnedAbs.flatMap(rd => rd.values).filter(v => v != null).sort((a, b) => a - b);
  const absMin = allAbs[Math.floor(0.01 * allAbs.length)] ?? 0;
  const absMax = allAbs[Math.floor(0.99 * allAbs.length)] ?? 10;

  _buildIonGridDOM('ion-fmz-grid', 'ion_fmz', ar.map(r => ({run: r})), ar);
  ar.forEach((run, idx) => {
    const rd = runBinnedAbs.find(x => x.run === run);
    if (rd) _mzProfileChart(`ion_fmz_${idx}`, rd, centers, [absMin, absMax], [mzMin, mzMax],
                             'log₁₀(Precursor Quantity)', {noMedian: true});
  });

  // Normalized: log2 FC per bin vs first active run
  const LOG2_10  = Math.log(10) / Math.log(2);
  const run1Vals = runBinnedAbs[0].values;

  const runBinnedNrm = ar.map((run, ridx) => {
    const values = runBinnedAbs[ridx].values.map((v, i) => {
      const ref = run1Vals[i];
      if (v == null || ref == null) return null;
      const fc = (v - ref) * LOG2_10;
      return isFinite(fc) ? fc : null;
    });
    return {run, values};
  });

  const allNrm = runBinnedNrm.flatMap(rd => rd.values).filter(v => v != null && isFinite(v));
  const nrmExt = allNrm.length ? Math.min(5, Math.max(2, ...allNrm.map(Math.abs))) : 3;

  if (fnmBadge) fnmBadge.textContent = `(${n_common.toLocaleString()} shared · vs ${label(ar[0])})`;
  _buildIonGridDOM('ion-fnm-grid', 'ion_fnm', ar.map(r => ({run: r})), ar);
  ar.forEach((run, idx) => {
    const rd = runBinnedNrm.find(x => x.run === run);
    if (rd) _mzProfileChart(`ion_fnm_${idx}`, rd, centers, [-nrmExt, nrmExt], [mzMin, mzMax],
                             'log₂ FC (vs run 1)', {medLabelRight: true});
  });
}

// Section D: MS2/MS1 Intensity Ratio (live JS, per-run grid histogram).
function _computeAndRenderRatio(ar) {
  const ms1PrecData = DATA.ion_sampling.prec_data      || {};
  const ms2PrecData = DATA.ion_sampling.ms2_prec_data  || {};
  const grid  = document.getElementById('ion-rto-grid');
  const badge = document.getElementById('ion-rto-badge');
  const msg   = document.getElementById('ion-rto-msg');

  function showMsg(txt) {
    if (grid)  grid.innerHTML = '';
    if (badge) badge.textContent = '';
    if (msg)   { msg.textContent = txt; msg.style.display = 'block'; }
  }
  function hideMsg() {
    if (msg) { msg.textContent = ''; msg.style.display = 'none'; }
  }

  const ms1Maps = {}, ms2Maps = {};
  ar.forEach(run => {
    ms1Maps[run] = new Map(Object.entries(ms1PrecData[run] || {}));
    ms2Maps[run] = new Map(Object.entries(ms2PrecData[run] || {}));
  });

  let common = null;
  ar.forEach(run => {
    const both = new Set([...ms1Maps[run].keys()].filter(p => ms2Maps[run].has(p)));
    common = common ? new Set([...common].filter(p => both.has(p))) : both;
  });
  const n_common = common ? common.size : 0;

  if (!common || n_common < 20) {
    showMsg(n_common === 0
      ? 'No precursors with both MS1 and MS2 quantification in all active runs.'
      : `Only ${n_common.toLocaleString()} shared precursors with both MS1 and MS2 (minimum 20 required).`);
    return;
  }

  hideMsg();
  if (badge) badge.textContent = `(${n_common.toLocaleString()} shared)`;

  const commonArr = [...common];
  const LOG2_10 = Math.log(10) / Math.log(2);
  const N_BINS  = 50;

  const runRatios = ar.map(run => ({
    run,
    ratios: commonArr
      .map(p => (ms2Maps[run].get(p) - ms1Maps[run].get(p)) * LOG2_10)
      .filter(isFinite),
  }));

  const allRatios = runRatios.flatMap(r => r.ratios).sort((a, b) => a - b);
  const n    = allRatios.length;
  const gmin = allRatios[Math.floor(0.01 * n)];
  const gmax = allRatios[Math.floor(0.99 * n)];
  const bsize   = (gmax - gmin) / N_BINS;
  const centers = Array.from({length: N_BINS}, (_, i) => gmin + (i + 0.5) * bsize);

  function binRatios(vals) {
    const counts = new Array(N_BINS).fill(0);
    vals.forEach(v => {
      if (v < gmin || v > gmax) return;
      counts[Math.min(N_BINS - 1, Math.floor((v - gmin) / bsize))]++;
    });
    return counts;
  }
  function medianOf(sv) {
    const s = sv.slice().sort((a, b) => a - b);
    const m = s.length;
    if (!m) return null;
    return Math.round((m % 2 ? s[Math.floor(m / 2)] : (s[m / 2 - 1] + s[m / 2]) / 2) * 10000) / 10000;
  }

  const runsData = runRatios.map(({run, ratios}) => {
    const counts = binRatios(ratios);
    return {run, counts, max_count: Math.max(...counts, 1), median: medianOf(ratios)};
  });

  _buildIonGridDOM('ion-rto-grid', 'ion_rto', runsData, ar);
  ar.forEach((run, idx) => {
    const rd = runsData.find(x => x.run === run);
    if (rd) _ionChart(`ion_rto_${idx}`, _ionDatasets(rd, null, centers, run),
                      gmin, gmax, rd.max_count, 'log₂(MS2/MS1)');
  });
}

// Render all ion sampling charts — called by both init and layout update
function _renderAllIonCharts(ar) {
  const ms1Data = DATA.ion_sampling.ms1_intensity;
  const ms2Data = DATA.ion_sampling.ms2_intensity;
  const fdtData = DATA.ion_sampling.ms1_fill_dist || {};
  const fgrData = DATA.ion_sampling.ms1_fill_grad || {};
  const fd2Data = DATA.ion_sampling.ms2_fill_dist || {};
  const fg2Data = DATA.ion_sampling.ms2_fill_grad || {};
  // Use ms1 aas if available (MS1 data present); fall back to ms2 aas (MS2-only, e.g. TMT)
  const aas     = ms1Data.aas || ms2Data.aas || null;

  // Section 1: MS1 intensity (enzyme-aware)
  const [ms1Min, ms1Max] = ms1Data.global_range || [0, 10];
  ar.forEach((run, idx) => {
    const rd = (ms1Data.data || []).find(x => x.run === run);
    if (rd) _ionChart(`ion_ms1_${idx}`, _ionDatasets(rd, aas, ms1Data.centers || [], run),
                      ms1Min, ms1Max, rd.max_count);
  });

  // Sections 2 & 3: MS1 live intersection (enzyme-aware absolute)
  _computeAndRenderIntersected(ar, {
    dataKey:    'prec_data',
    intGridId:  'ion-int-grid', nrmGridId: 'ion-nrm-grid',
    intBadgeId: 'ion-int-badge', nrmBadgeId: 'ion-nrm-badge',
    intMsgId:   'ion-int-msg',  nrmMsgId:  'ion-nrm-msg',
    intPrefix:  'ion_int',      nrmPrefix: 'ion_nrm',
    yTitle:     'log₁₀(Intensity)',
    aas,
  });

  // Section 4: MS1 fill time distribution
  const [fdtMin, fdtMax] = fdtData.global_range || [0, 100];
  ar.forEach((run, idx) => {
    const rd = (fdtData.data || []).find(x => x.run === run);
    if (rd) _ionChart(`ion_fdt_${idx}`, _ionDatasets(rd, null, fdtData.centers || [], run),
                      fdtMin, fdtMax, rd.max_count, 'Fill Time (ms)');
  });

  // Section 5: MS1 fill times along gradient
  ar.forEach((run, idx) => {
    const rd = (fgrData.data || []).find(x => x.run === run);
    if (rd) _fillGradChart(`ion_fgr_${idx}`, rd.bins,
                            fgrData.rt_range || [0, 120], fgrData.ft_range || [0, 100]);
  });

  // Section 6: MS2 intensity (enzyme-aware)
  const [ms2Min, ms2Max] = ms2Data.global_range || [0, 10];
  ar.forEach((run, idx) => {
    const rd = (ms2Data.data || []).find(x => x.run === run);
    if (rd) _ionChart(`ion_ms2_${idx}`, _ionDatasets(rd, aas, ms2Data.centers || [], run),
                      ms2Min, ms2Max, rd.max_count, 'log₁₀(Precursor Quantity)');
  });

  // Sections 7 & 8: MS2 live intersection (enzyme-aware absolute)
  _computeAndRenderIntersected(ar, {
    dataKey:    'ms2_prec_data',
    intGridId:  'ion-i2t-grid', nrmGridId: 'ion-n2m-grid',
    intBadgeId: 'ion-i2t-badge', nrmBadgeId: 'ion-n2m-badge',
    intMsgId:   'ion-i2t-msg',  nrmMsgId:  'ion-n2m-msg',
    intPrefix:  'ion_i2t',      nrmPrefix: 'ion_n2m',
    yTitle:     'log₁₀(Precursor Quantity)',
    aas,
  });

  // Section 9: MS2 fill time distribution
  const [fd2Min, fd2Max] = fd2Data.global_range || [0, 100];
  ar.forEach((run, idx) => {
    const rd = (fd2Data.data || []).find(x => x.run === run);
    if (rd) _ionChart(`ion_fd2_${idx}`, _ionDatasets(rd, null, fd2Data.centers || [], run),
                      fd2Min, fd2Max, rd.max_count, 'Fill Time (ms)');
  });

  // Section 10: MS2 fill times along gradient
  ar.forEach((run, idx) => {
    const rd = (fg2Data.data || []).find(x => x.run === run);
    if (rd) _fillGradChart(`ion_fg2_${idx}`, rd.bins,
                            fg2Data.rt_range || [0, 120], fg2Data.ft_range || [0, 100]);
  });

  // Section 12: MS2 Intensity vs Precursor M/Z (per-run, precomputed)
  _renderMs2VsMz(ar);

  // Sections 13 + 14: Best Fragment M/Z (live JS, per-run grids)
  _computeAndRenderFragMz(ar);

  // Section 15: MS2/MS1 Intensity Ratio (live JS, per-run grid)
  _computeAndRenderRatio(ar);
}

function initIonSampling() {
  const ar          = getActiveRuns();
  const defaultCols = Math.min(ar.length, 6) || 1;
  const colsSel     = document.getElementById('ion-cols-sel');
  if (colsSel) colsSel.value = String(defaultCols);
  const cols = `repeat(${defaultCols},1fr)`;
  ['ion-ms1-grid','ion-int-grid','ion-nrm-grid','ion-fdt-grid','ion-fgr-grid',
   'ion-ms2-grid','ion-i2t-grid','ion-n2m-grid','ion-fd2-grid','ion-fg2-grid',
   'ion-mz2-grid','ion-fmz-grid','ion-fnm-grid','ion-rto-grid'].forEach(id => {
    const g = document.getElementById(id);
    if (g) g.style.gridTemplateColumns = cols;
  });

  const ms1Data = DATA.ion_sampling.ms1_intensity;
  const ms2Data = DATA.ion_sampling.ms2_intensity;
  const fdtData = DATA.ion_sampling.ms1_fill_dist || {};
  const fgrData = DATA.ion_sampling.ms1_fill_grad || {};
  const fd2Data = DATA.ion_sampling.ms2_fill_dist || {};
  const fg2Data = DATA.ion_sampling.ms2_fill_grad || {};

  _buildIonGridDOM('ion-ms1-grid', 'ion_ms1', ms1Data.data || [], ar);
  _buildIonGridDOM('ion-ms2-grid', 'ion_ms2', ms2Data.data || [], ar);
  _buildIonGridDOM('ion-fdt-grid', 'ion_fdt', fdtData.data || [], ar);
  _buildIonGridDOM('ion-fgr-grid', 'ion_fgr', fgrData.data || [], ar);
  _buildIonGridDOM('ion-fd2-grid', 'ion_fd2', fd2Data.data || [], ar);
  _buildIonGridDOM('ion-fg2-grid', 'ion_fg2', fg2Data.data || [], ar);
  _renderAllIonCharts(ar);
}

function updateIonLayout() {
  const n    = parseInt(document.getElementById('ion-cols-sel').value);
  const cols = `repeat(${n}, 1fr)`;
  ['ion-ms1-grid','ion-int-grid','ion-nrm-grid','ion-fdt-grid','ion-fgr-grid',
   'ion-ms2-grid','ion-i2t-grid','ion-n2m-grid','ion-fd2-grid','ion-fg2-grid',
   'ion-mz2-grid','ion-fmz-grid','ion-fnm-grid','ion-rto-grid'].forEach(id => {
    const g = document.getElementById(id);
    if (g) g.style.gridTemplateColumns = cols;
  });
  Object.keys(CHARTS).forEach(id => {
    if (/^ion_(ms1|int|nrm|fdt|fgr|ms2|i2t|n2m|fd2|fg2|mz2|fmz|fnm|rto)_/.test(id)) destroyChart(id);
  });
  const ar = getActiveRuns();
  requestAnimationFrame(() => _renderAllIonCharts(ar));
}
"""
