# Created by @thibaultcolombani with the help of Claude Code (Anthropic)
"""html_tab_chromatography.py — Chromatography tab: RT distribution and FWHM histogram."""
from __future__ import annotations
import math
from typing import Dict, List
import polars as pl
from Mods.html_utils import histogram

TAB_ID    = 'chromatography'
TAB_LABEL = 'Chromatography'


def is_visible(processor) -> bool:
    return 'Retention time' in processor.report.columns


def compute(processor) -> Dict:
    df      = processor.report
    runs    = processor.runs()
    aa_keep = _enzyme_cterm_aas(getattr(processor, 'enzymes', None))
    return {
        'rt':   _rt_dist(df, runs, aa_keep),
        'fwhm': _fwhm(df, runs, aa_keep),
    }


def _enzyme_cterm_aas(enzymes):
    """Return frozenset of expected C-terminal AAs for the given enzyme list, or None if unknown."""
    if not enzymes:
        return None
    from Mods.enzyme import ENZYMES as _E
    result = set()
    for e in enzymes:
        result.update(_E.get(e, {}).get('cterm', set()))
    return frozenset(result) if result else None


def _rt_dist(df: pl.DataFrame, runs: List[str], aa_keep=None) -> Dict:
    if 'Retention time' not in df.columns or 'C_term_aa' not in df.columns:
        return {'data': [], 'centers': [], 'rt_range': [0, 120], 'aas': []}

    gmin = float(df['Retention time'].min())
    gmax = float(df['Retention time'].max())
    bsize = 0.5
    n_bins = max(1, math.ceil((gmax - gmin) / bsize))
    centers = [round(gmin + (i + 0.5) * bsize, 3) for i in range(n_bins)]

    df2 = df.with_columns(
        ((pl.col('Retention time') - gmin) / bsize).floor().cast(pl.Int32).alias('_b')
    )

    # R first, K second, then others alphabetically; filter to enzyme AAs if known
    all_aas = sorted(df['C_term_aa'].drop_nulls().unique().to_list())
    if aa_keep:
        all_aas = [a for a in all_aas if a in aa_keep]
    priority = ['R', 'K']
    ordered_aas = [a for a in priority if a in all_aas] + sorted(a for a in all_aas if a not in priority)

    def _bcounts(frame):
        if frame.height == 0:
            return [0] * n_bins
        grp = frame.group_by('_b').agg(pl.len().alias('c')).sort('_b')
        out = [0] * n_bins
        for b, c in grp.iter_rows():
            if 0 <= b < n_bins:
                out[b] = int(c)
        return out

    data = []
    for run in runs:
        rdf = df2.filter(pl.col('Raw file') == run)
        aa_series = {aa: _bcounts(rdf.filter(pl.col('C_term_aa') == aa)) for aa in ordered_aas}
        max_count = max((max(v) for v in aa_series.values() if v), default=1)
        data.append({'run': run, 'aa_series': aa_series, 'max_count': max_count})

    return {'data': data, 'centers': centers, 'rt_range': [gmin, gmax], 'aas': ordered_aas}


def _fwhm(df: pl.DataFrame, runs: List[str], aa_keep=None) -> Dict:
    if 'Retention length' not in df.columns or 'C_term_aa' not in df.columns:
        return {'data': [], 'centers': [], 'fwhm_range': [0, 30], 'aas': []}

    # Global p1–p99 range in seconds across all runs
    all_vals = (df.filter(pl.col('Retention length') > 0)['Retention length']
                  .drop_nulls().to_list())
    if not all_vals:
        return {'data': [], 'centers': [], 'fwhm_range': [0, 30], 'aas': []}
    sv = sorted(float(v) * 60 for v in all_vals)
    n  = len(sv)
    gmin = sv[max(0, int(0.01 * n))]
    gmax = sv[min(n - 1, int(0.99 * n))]
    n_bins = 50
    bsize  = (gmax - gmin) / n_bins
    centers = [round(gmin + (i + 0.5) * bsize, 4) for i in range(n_bins)]

    # Same AA ordering as _rt_dist; filter to enzyme AAs if known
    all_aas = sorted(df['C_term_aa'].drop_nulls().unique().to_list())
    if aa_keep:
        all_aas = [a for a in all_aas if a in aa_keep]
    priority = ['R', 'K']
    ordered_aas = [a for a in priority if a in all_aas] + sorted(a for a in all_aas if a not in priority)

    df2 = (df.filter(pl.col('Retention length') > 0)
             .with_columns((pl.col('Retention length') * 60).alias('_fwhm_s'))
             .with_columns(
                 ((pl.col('_fwhm_s') - gmin) / bsize).floor().cast(pl.Int32).alias('_b')
             ))

    def _bcounts(frame):
        if frame.height == 0:
            return [0] * n_bins
        grp = frame.group_by('_b').agg(pl.len().alias('c')).sort('_b')
        out = [0] * n_bins
        for b, c in grp.iter_rows():
            if 0 <= b < n_bins:
                out[b] = int(c)
        return out

    data = []
    for run in runs:
        rdf = df2.filter(pl.col('Raw file') == run)
        if rdf.is_empty():
            continue   # run has no FWHM data — omit so JS skips the chart box
        aa_series = {aa: _bcounts(rdf.filter(pl.col('C_term_aa') == aa)) for aa in ordered_aas}
        max_count = max((max(v) for v in aa_series.values() if v), default=0)
        if max_count == 0:
            continue   # all-zero — omit
        data.append({'run': run, 'aa_series': aa_series, 'max_count': max_count})

    return {'data': data, 'centers': centers, 'fwhm_range': [gmin, gmax], 'aas': ordered_aas}


def html_content(data: Dict) -> str:
    n_runs = len((data.get('rt') or {}).get('data') or [])
    default_cols = min(max(1, n_runs), 6)
    cols_opts = ''.join(
        f'<option value="{n}"{" selected" if n == default_cols else ""}>{n}</option>'
        for n in range(1, 7)
    )
    return f"""
<div class="section">
  <div class="run-selector">
    <label>Samples per row:</label>
    <select id="chrom-cols-sel" onchange="updateChromLayout()">
      {cols_opts}
    </select>
  </div>
</div>
<div class="section">
  <h3>Retention Time Distribution</h3>
  <p class="help">Precursors per 0.5-min RT bin, coloured by C-terminal amino acid. Y = RT (min), X = count.</p>
  <div id="chrom-rt-grid" class="chart-grid" style="grid-template-columns:repeat({default_cols},1fr)"></div>
</div>
<div class="section">
  <h3>Peak Width — FWHM</h3>
  <p class="help">Peak width distribution per run. Y = FWHM (s), X = precursor count.</p>
  <div id="chrom-fwhm-grid" class="chart-grid" style="grid-template-columns:repeat({default_cols},1fr)"></div>
</div>"""


def javascript() -> str:
    return """
function _niceTicks(maxVal) {
  // 4 equal intervals → 5 labels: 0, s, 2s, 3s, 4s
  if (!maxVal || maxVal <= 0) return {step: 1, axMax: 4};
  const rawStep = maxVal / 4;
  const mag = Math.pow(10, Math.floor(Math.log10(rawStep)));
  const step = Math.ceil(rawStep / mag) * mag;
  return {step, axMax: step * 4};
}

function _rtChart(canvasId, datasets, rtMin, rtMax, maxCount) {
  destroyChart(canvasId);
  const ctx = document.getElementById(canvasId);
  if (!ctx) return;
  const {step, axMax: xMax} = _niceTicks(maxCount);
  const yMin = Math.floor(rtMin / 5) * 5;
  const yMax = Math.ceil(rtMax  / 5) * 5;
  CHARTS[canvasId] = new Chart(ctx, {
    type: 'line',
    data: { datasets },
    options: {
      responsive: true,
      maintainAspectRatio: false,
      parsing: false,
      layout: { padding: { bottom: 8, right: 12 } },
      plugins: {
        legend: { display: true, position: 'top', labels: { boxWidth: 12, font: { size: 10 } } }
      },
      elements: { point: { radius: 2.5, hoverRadius: 4 }, line: { borderWidth: 1.5 } },
      scales: {
        x: { type: 'linear', min: 0, max: xMax,
             title: { display: true, text: 'Count' },
             ticks: { stepSize: step, precision: 0, maxRotation: 0, minRotation: 0,
                      font: { size: 9 } } },
        y: { type: 'linear', reverse: true, min: yMin, max: yMax,
             title: { display: true, text: 'RT (min)' },
             ticks: {
               stepSize: 0.5,
               callback: (v) => Math.abs(Math.round(v * 10) % 50) < 1 ? v : null,
             },
             grid: {
               color: (ctx) => Math.abs(Math.round(ctx.tick.value * 10) % 50) < 1
                               ? '#ddd' : 'transparent',
             } },
      }
    }
  });
}

function _fwhmChart(canvasId, datasets, fwhmMin, fwhmMax, maxCount) {
  destroyChart(canvasId);
  const ctx = document.getElementById(canvasId);
  if (!ctx) return;
  const {step: xStep, axMax: xMax} = _niceTicks(maxCount);
  const {step: yStep, axMax: yMax} = _niceTicks(fwhmMax);
  CHARTS[canvasId] = new Chart(ctx, {
    type: 'line',
    data: { datasets },
    options: {
      responsive: true,
      maintainAspectRatio: false,
      parsing: false,
      layout: { padding: { bottom: 8, right: 12 } },
      plugins: {
        legend: { display: true, position: 'top', labels: { boxWidth: 12, font: { size: 10 } } }
      },
      elements: { point: { radius: 2.5, hoverRadius: 4 }, line: { borderWidth: 1.5 } },
      scales: {
        x: { type: 'linear', min: 0, max: xMax,
             title: { display: true, text: 'Count' },
             ticks: { stepSize: xStep, precision: 0, maxRotation: 0, minRotation: 0,
                      font: { size: 9 } } },
        y: { type: 'linear', min: 0, max: yMax,
             title: { display: true, text: 'FWHM (s)' },
             ticks: { stepSize: yStep, precision: 2 } },
      }
    }
  });
}

function initChromatography() {
  const ar  = getActiveRuns();
  const defaultCols = Math.min(ar.length, 6) || 1;
  const colsSel = document.getElementById('chrom-cols-sel');
  if (colsSel) colsSel.value = String(defaultCols);
  ['chrom-rt-grid', 'chrom-fwhm-grid'].forEach(id => {
    const g = document.getElementById(id);
    if (g) g.style.gridTemplateColumns = `repeat(${defaultCols},1fr)`;
  });

  const d   = DATA.chromatography;
  const rtD = d.rt;
  const [rtMin, rtMax] = rtD.rt_range || [0, 120];

  // ── RT grid ────────────────────────────────────────────────────────────────
  const rtGrid = document.getElementById('chrom-rt-grid');
  rtGrid.innerHTML = '';
  ar.forEach((run, idx) => {
    const rd = rtD.data.find(x => x.run === run);
    if (!rd) return;
    const canvasId = `ch_rt_${idx}`;
    const box = document.createElement('div');
    box.className = 'chart-box';
    box.innerHTML =
      `<h4>${label(run)}</h4>` +
      `<button class="export-btn" onclick="exportChart('${canvasId}')">PNG</button>` +
      `<div style="position:relative;height:375px"><canvas id="${canvasId}"></canvas></div>`;
    rtGrid.appendChild(box);
  });
  ar.forEach((run, idx) => {
    const rd = rtD.data.find(x => x.run === run);
    if (!rd) return;
    const centers  = rtD.centers;
    const datasets = rtD.aas.map((aa, i) => ({
      label: aa,
      data: centers.map((c, j) => ({x: rd.aa_series[aa][j] || 0, y: c})),
      borderColor: aaColor(aa, i),
      backgroundColor: 'transparent',
      borderWidth: 1.5,
    }));
    _rtChart(`ch_rt_${idx}`, datasets, rtMin, rtMax, rd.max_count);
  });

  // ── FWHM grid ──────────────────────────────────────────────────────────────
  const fwGrid = document.getElementById('chrom-fwhm-grid');
  fwGrid.innerHTML = '';
  const fwD = d.fwhm || {};
  const fwData = fwD.data || [];
  const [fwMin, fwMax] = fwD.fwhm_range || [0, 30];
  const fwAas = fwD.aas || [];
  const fwCenters = fwD.centers || [];
  ar.forEach((run, idx) => {
    const rd = fwData.find(x => x.run === run);
    if (!rd) return;
    const canvasId = `ch_fw_${idx}`;
    const box = document.createElement('div');
    box.className = 'chart-box';
    box.innerHTML =
      `<h4>${label(run)}</h4>` +
      `<button class="export-btn" onclick="exportChart('${canvasId}')">PNG</button>` +
      `<div style="position:relative;height:375px"><canvas id="${canvasId}"></canvas></div>`;
    fwGrid.appendChild(box);
  });
  ar.forEach((run, idx) => {
    const rd = fwData.find(x => x.run === run);
    if (!rd) return;
    const datasets = fwAas.map((aa, i) => ({
      label: aa,
      data: fwCenters.map((c, j) => ({x: rd.aa_series[aa][j] || 0, y: c})),
      borderColor: aaColor(aa, i),
      backgroundColor: 'transparent',
      borderWidth: 1.5,
    }));
    _fwhmChart(`ch_fw_${idx}`, datasets, fwMin, fwMax, rd.max_count);
  });
}

function updateChromLayout() {
  const n    = parseInt(document.getElementById('chrom-cols-sel').value);
  const cols = `repeat(${n}, 1fr)`;
  ['chrom-rt-grid', 'chrom-fwhm-grid'].forEach(id => {
    const g = document.getElementById(id);
    if (g) g.style.gridTemplateColumns = cols;
  });
  Object.keys(CHARTS).forEach(id => {
    if (id.startsWith('ch_rt_') || id.startsWith('ch_fw_')) destroyChart(id);
  });
  const ar      = getActiveRuns();
  const rtD     = DATA.chromatography.rt;
  const [rtMin, rtMax] = rtD.rt_range || [0, 120];
  const fwD     = DATA.chromatography.fwhm || {};
  const fwData  = fwD.data || [];
  const [fwMin, fwMax] = fwD.fwhm_range || [0, 30];
  const fwAas   = fwD.aas || [];
  const fwCents = fwD.centers || [];
  requestAnimationFrame(() => {
    ar.forEach((run, idx) => {
      const rtRd = rtD.data.find(x => x.run === run);
      if (rtRd) {
        const datasets = rtD.aas.map((aa, i) => ({
          label: aa,
          data: rtD.centers.map((c, j) => ({x: rtRd.aa_series[aa][j] || 0, y: c})),
          borderColor: aaColor(aa, i), backgroundColor: 'transparent', borderWidth: 1.5,
        }));
        _rtChart(`ch_rt_${idx}`, datasets, rtMin, rtMax, rtRd.max_count);
      }
      const fwRd = fwData.find(x => x.run === run);
      if (fwRd) {
        const datasets = fwAas.map((aa, i) => ({
          label: aa,
          data: fwCents.map((c, j) => ({x: fwRd.aa_series[aa][j] || 0, y: c})),
          borderColor: aaColor(aa, i), backgroundColor: 'transparent', borderWidth: 1.5,
        }));
        _fwhmChart(`ch_fw_${idx}`, datasets, fwMin, fwMax, fwRd.max_count);
      }
    });
  });
}
"""
