# Created by @thibaultcolombani with the help of Claude Code (Anthropic)
"""html_tab_plex_diagnostic.py — Plex diagnostic tab. Shown only when multiplexed data is detected."""
from __future__ import annotations
from typing import Dict, List, Optional
import math
import re
import polars as pl

TAB_ID    = 'plex_diagnostic'
TAB_LABEL = 'Plex diagnostic'


def is_visible(processor) -> bool:
    plex = getattr(processor, 'plex_report', None)
    return plex is not None and not plex.is_empty()


def compute(processor) -> Dict:
    plex = getattr(processor, 'plex_report', None)
    if plex is None or plex.is_empty():
        return {}
    plex_ch = plex.filter(pl.col('Channel').is_not_null() & (pl.col('Channel') != ''))
    runs    = sorted(plex_ch['Raw file'].drop_nulls().unique().to_list())
    ref_ch  = _find_ref_channel(plex_ch)
    has_ms2  = 'Precursor.Quantity' in plex_ch.columns and plex_ch['Precursor.Quantity'].drop_nulls().len() > 0
    # MS1 is present when plex Intensity differs from Precursor.Quantity (PSMtag/DIA).
    # TMT reporter ions: Intensity == Precursor.Quantity (same value) → no separate MS1.
    _plex_has_int = ('Intensity' in plex_ch.columns and
                     (plex_ch['Intensity'].max() or 0) > 0)
    if _plex_has_int and has_ms2:
        _sample = plex_ch.filter(pl.col('Intensity') > 0).head(500)
        _n_diff = int(_sample.select(
            (pl.col('Intensity') != pl.col('Precursor.Quantity')).sum()
        )[0, 0])
        has_ms1 = _n_diff > 5
    else:
        has_ms1 = _plex_has_int and not has_ms2
    qval_data = _plex_qval(plex_ch, runs)
    has_qval  = bool(qval_data)
    return {
        'runs':              runs,
        'has_ms1':           has_ms1,
        'has_ms2':           has_ms2,
        'has_qval':          has_qval,
        'plex_ids':          _plex_ids(plex_ch, runs),
        'plex_proteins':     _plex_proteins(plex_ch, runs),
        'plex_qval':         qval_data,
        'plex_intensity':    _plex_intensity(plex_ch, runs, ref_ch),
        'plex_ms1_dist':     _plex_ms1_dist(plex_ch, runs) if has_ms1 else {},
        'plex_cv':           _plex_cv(plex_ch, runs)        if has_ms1 else {},
        'plex_ms2_dist':     _plex_ms2_dist(plex_ch, runs) if has_ms2 else {},
        'plex_cv_ms2':       _plex_cv_ms2(plex_ch, runs)   if has_ms2 else {},
        'plex_jaccard':      _plex_jaccard(plex_ch, runs, ref_ch),
        'plex_jaccard_prot': _plex_jaccard_proteins(plex_ch, runs, ref_ch),
    }


# ---------------------------------------------------------------------------
# Python helpers
# ---------------------------------------------------------------------------

def _channels(df: pl.DataFrame) -> List[str]:
    """Channel names for `df` sorted numerically; empty strings excluded."""
    chs = [c for c in df['Channel'].drop_nulls().unique().to_list() if str(c).strip() != '']
    def _key(ch: str):
        m = re.search(r'\d+', str(ch))
        return int(m.group()) if m else float('inf')
    return sorted(chs, key=_key)


def _run_channels(plex: pl.DataFrame, run: str) -> List[str]:
    """Channels present for a specific run (non-empty, non-null)."""
    return _channels(plex.filter(pl.col('Raw file') == run))


def _find_ref_channel(plex: pl.DataFrame) -> Optional[str]:
    """Highest global median intensity = carrier/reference channel."""
    if 'Intensity' not in plex.columns:
        return None
    df = plex.filter(pl.col('Intensity') > 0)
    if df.is_empty():
        return None
    ch_meds = (df.group_by('Channel')
                 .agg(pl.col('Intensity').median().alias('med'))
                 .sort('med', descending=True))
    return ch_meds['Channel'][0] if not ch_meds.is_empty() else None


def _plex_ids(plex: pl.DataFrame, runs: List[str]) -> Dict:
    all_channels = _channels(plex)
    result: Dict = {'channels': all_channels}
    for run in runs:
        rdf          = plex.filter(pl.col('Raw file') == run)
        run_chs      = _channels(rdf)
        if 'Precursor.Id' not in rdf.columns:
            result[run] = {'channels': run_chs, 'per_channel': {}, 'intersected': 0, 'union': 0}
            continue
        ch_sets: Dict[str, set] = {
            ch: set(rdf.filter(pl.col('Channel') == ch)['Precursor.Id'].drop_nulls().to_list())
            for ch in run_chs
        }
        intersected = set.intersection(*ch_sets.values()) if ch_sets else set()
        union_set   = set.union(*ch_sets.values())        if ch_sets else set()
        result[run] = {
            'channels':    run_chs,
            'per_channel': {ch: len(s) for ch, s in ch_sets.items()},
            'intersected': len(intersected),
            'union':       len(union_set),
        }
    return result


def _plex_proteins(plex: pl.DataFrame, runs: List[str]) -> Dict:
    if 'Protein.Group' not in plex.columns:
        return {}
    all_channels = _channels(plex)
    df = plex.filter((pl.col('Proteotypic') == 1) | pl.col('Proteotypic').is_null()) if 'Proteotypic' in plex.columns else plex
    result: Dict = {'channels': all_channels}
    for run in runs:
        rdf     = df.filter(pl.col('Raw file') == run)
        run_chs = _channels(rdf)
        ch_sets: Dict[str, set] = {
            ch: set(rdf.filter(pl.col('Channel') == ch)['Protein.Group'].drop_nulls().to_list())
            for ch in run_chs
        }
        intersected = set.intersection(*ch_sets.values()) if ch_sets else set()
        union_set   = set.union(*ch_sets.values())        if ch_sets else set()
        result[run] = {
            'channels':    run_chs,
            'per_channel': {ch: len(s) for ch, s in ch_sets.items()},
            'intersected': len(intersected),
            'union':       len(union_set),
        }
    return result


def _plex_qval(plex: pl.DataFrame, runs: List[str]) -> Dict:
    if 'Channel.Q.Value' not in plex.columns:
        return {}
    from Mods.html_utils import histogram
    result: Dict = {}
    for run in runs:
        rdf     = plex.filter(pl.col('Raw file') == run)
        run_chs = _channels(rdf)
        result[run] = {'channels': run_chs}
        for ch in run_chs:
            vals = (rdf.filter((pl.col('Channel') == ch)
                               & pl.col('Channel.Q.Value').is_not_null()
                               & (pl.col('Channel.Q.Value') > 0))
                       .select('Channel.Q.Value')['Channel.Q.Value'].to_list())
            result[run][ch] = histogram([float(v) for v in vals if not math.isnan(v)])
    return result


def _plex_intensity(plex: pl.DataFrame, runs: List[str],
                    ref_ch: Optional[str] = None) -> Dict:
    """Per-run log₁₀ ratio vs carrier for each sample channel: [p5,p25,p50,p75,p95]."""
    if 'Intensity' not in plex.columns or 'Precursor.Id' not in plex.columns:
        return {}
    df = plex.filter(pl.col('Intensity') > 0)
    if ref_ch is None:
        ref_ch = _find_ref_channel(plex)
    result: Dict = {'ref_channel': ref_ch}
    for run in runs:
        rdf     = df.filter(pl.col('Raw file') == run)
        run_chs = _channels(rdf)
        if len(run_chs) < 2:
            result[run] = {'channels': [], 'values': {}}
            continue
        if ref_ch not in run_chs:
            # Use highest-intensity channel in this run as reference
            local_ref = _find_ref_channel(rdf) or run_chs[-1]
        else:
            local_ref = ref_ch
        non_ref = [ch for ch in run_chs if ch != local_ref]
        ref_df  = (rdf.filter(pl.col('Channel') == local_ref)
                      .select(['Precursor.Id', pl.col('Intensity').alias('ref_int')]))
        run_data: Dict = {}
        for ch in non_ref:
            ch_df  = (rdf.filter(pl.col('Channel') == ch)
                         .select(['Precursor.Id', pl.col('Intensity').alias('ch_int')]))
            joined = ref_df.join(ch_df, on='Precursor.Id')
            if joined.is_empty():
                run_data[ch] = None
                continue
            log_df = (joined
                      .filter((pl.col('ch_int') > 0) & (pl.col('ref_int') > 0))
                      .with_columns(
                          ((pl.col('ch_int') / pl.col('ref_int')).log() / math.log(10))
                          .alias('lr')
                      )
                      .drop_nulls('lr')
                      .filter(pl.col('lr').is_finite()))
            if log_df.is_empty():
                run_data[ch] = None
                continue
            s = log_df['lr']
            run_data[ch] = [round(float(s.quantile(p)), 4)
                            for p in [0.05, 0.25, 0.50, 0.75, 0.95]]
        result[run] = {'channels': non_ref, 'ref_channel': local_ref, 'values': run_data}
    return result


def _plex_cv(plex: pl.DataFrame, runs: List[str]) -> Dict:
    """Median MS1 CV per channel per run."""
    if not all(c in plex.columns for c in ['Protein.Group', 'seqcharge', 'Intensity', 'Channel']):
        return {}
    df = plex.filter(pl.col('Intensity') > 0)
    df = df.with_columns(
        (pl.col('Intensity') / pl.col('Intensity').mean().over(['Raw file', 'seqcharge']))
        .alias('norm1')
    )
    df = df.with_columns(
        (pl.col('norm1') / pl.col('norm1').median().over(['Raw file', 'Channel']))
        .alias('norm2')
    )
    cv_df = (
        df.filter(pl.col('norm2').is_not_null() & pl.col('norm2').is_finite() & (pl.col('norm2') > 0))
          .group_by(['Raw file', 'Channel', 'Protein.Group'])
          .agg([pl.col('norm2').count().alias('n'),
                pl.col('norm2').std().alias('std'),
                pl.col('norm2').mean().alias('mean')])
          .filter((pl.col('n') >= 3) & (pl.col('mean') > 0) & pl.col('std').is_not_null())
          .with_columns((pl.col('std') / pl.col('mean')).alias('cv'))
          .group_by(['Raw file', 'Channel'])
          .agg(pl.col('cv').median().alias('median_cv'))
    )
    result: Dict = {}
    for run in runs:
        run_chs = _run_channels(plex, run)
        rdf     = cv_df.filter(pl.col('Raw file') == run)
        result[run] = {
            'channels': run_chs,
            **{ch: (round(float(rdf.filter(pl.col('Channel') == ch)['median_cv'][0]), 4)
                    if not rdf.filter(pl.col('Channel') == ch).is_empty() else None)
               for ch in run_chs}
        }
    return result


def _intensity_dist(plex: pl.DataFrame, runs: List[str], col: str) -> Dict:
    """Shared-bin histogram: global p1–p99 range; per-run channels computed independently."""
    if col not in plex.columns:
        return {}
    all_vals = [math.log10(float(v)) for v in plex.filter(pl.col(col) > 0)[col].to_list()]
    if not all_vals:
        return {}
    sv = sorted(all_vals)
    n  = len(sv)
    g_p1  = sv[max(0, int(0.01 * n))]
    g_p99 = sv[min(n - 1, int(0.99 * n))]
    n_bins = 50
    bw     = (g_p99 - g_p1) / n_bins
    centers = [round(g_p1 + (i + 0.5) * bw, 4) for i in range(n_bins)]
    result: Dict = {'centers': centers, 'range': [g_p1, g_p99]}
    for run in runs:
        rdf     = plex.filter(pl.col('Raw file') == run)
        run_chs = _channels(rdf)
        run_max = 0
        ch_data: Dict = {'channels': run_chs}
        medians: Dict = {}
        for ch in run_chs:
            vals = [math.log10(float(v))
                    for v in rdf.filter((pl.col('Channel') == ch) & (pl.col(col) > 0))[col].to_list()]
            counts = [0] * n_bins
            for v in vals:
                if g_p1 <= v <= g_p99:
                    b = min(n_bins - 1, int((v - g_p1) / bw))
                    counts[b] += 1
            run_max = max(run_max, max(counts) if counts else 0)
            ch_data[ch] = counts
            medians[ch] = round(sorted(vals)[len(vals) // 2], 4) if vals else None
        if run_max == 0:
            continue  # no data for this column in this run (e.g., Jmod for MS2)
        ch_data['_max']     = run_max
        ch_data['_medians'] = medians
        result[run] = ch_data
    return result


def _plex_ms1_dist(plex: pl.DataFrame, runs: List[str]) -> Dict:
    return _intensity_dist(plex, runs, 'Intensity')


def _plex_ms2_dist(plex: pl.DataFrame, runs: List[str]) -> Dict:
    return _intensity_dist(plex, runs, 'Precursor.Quantity')


def _plex_cv_ms2(plex: pl.DataFrame, runs: List[str]) -> Dict:
    """Median MS2 CV per channel per run using Precursor.Quantity."""
    if not all(c in plex.columns for c in ['Protein.Group', 'seqcharge', 'Precursor.Quantity', 'Channel']):
        return {}
    df = plex.filter(pl.col('Precursor.Quantity') > 0)
    df = df.with_columns(
        (pl.col('Precursor.Quantity') / pl.col('Precursor.Quantity').mean().over(['Raw file', 'seqcharge']))
        .alias('norm1')
    )
    df = df.with_columns(
        (pl.col('norm1') / pl.col('norm1').median().over(['Raw file', 'Channel']))
        .alias('norm2')
    )
    cv_df = (
        df.filter(pl.col('norm2').is_not_null() & pl.col('norm2').is_finite() & (pl.col('norm2') > 0))
          .group_by(['Raw file', 'Channel', 'Protein.Group'])
          .agg([pl.col('norm2').count().alias('n'),
                pl.col('norm2').std().alias('std'),
                pl.col('norm2').mean().alias('mean')])
          .filter((pl.col('n') >= 3) & (pl.col('mean') > 0) & pl.col('std').is_not_null())
          .with_columns((pl.col('std') / pl.col('mean')).alias('cv'))
          .group_by(['Raw file', 'Channel'])
          .agg(pl.col('cv').median().alias('median_cv'))
    )
    result: Dict = {}
    for run in runs:
        run_chs = _run_channels(plex, run)
        rdf     = cv_df.filter(pl.col('Raw file') == run)
        if rdf.is_empty():
            continue  # no MS2 data for this run (e.g., Jmod)
        result[run] = {
            'channels': run_chs,
            **{ch: (round(float(rdf.filter(pl.col('Channel') == ch)['median_cv'][0]), 4)
                    if not rdf.filter(pl.col('Channel') == ch).is_empty() else None)
               for ch in run_chs}
        }
    return result


def _plex_jaccard(plex: pl.DataFrame, runs: List[str],
                  ref_ch: Optional[str] = None) -> Dict:
    """Pairwise Jaccard index between sample channels (carrier excluded), per run."""
    if 'Precursor.Id' not in plex.columns:
        return {}
    result: Dict = {}
    for run in runs:
        rdf        = plex.filter(pl.col('Raw file') == run)
        run_chs    = _channels(rdf)
        sample_chs = [ch for ch in run_chs if ch != ref_ch]
        if len(sample_chs) < 2:
            result[run] = {'pairs': [], 'values': []}
            continue
        pairs = [(sample_chs[i], sample_chs[j])
                 for i in range(len(sample_chs))
                 for j in range(i + 1, len(sample_chs))]
        ch_prec = {
            ch: set(rdf.filter(pl.col('Channel') == ch)['Precursor.Id'].drop_nulls().to_list())
            for ch in sample_chs
        }
        jacc = []
        for a, b in pairs:
            sa, sb = ch_prec.get(a, set()), ch_prec.get(b, set())
            inter  = len(sa & sb)
            union  = len(sa | sb)
            jacc.append(round(inter / union, 4) if union else 0)
        result[run] = {'pairs': [f'{a} / {b}' for a, b in pairs], 'values': jacc}
    return result


def _plex_jaccard_proteins(plex: pl.DataFrame, runs: List[str],
                           ref_ch: Optional[str] = None) -> Dict:
    """Pairwise protein Jaccard index between sample channels (carrier excluded), per run."""
    if 'Protein.Group' not in plex.columns:
        return {}
    result: Dict = {}
    for run in runs:
        rdf        = plex.filter(pl.col('Raw file') == run)
        run_chs    = _channels(rdf)
        sample_chs = [ch for ch in run_chs if ch != ref_ch]
        if len(sample_chs) < 2:
            result[run] = {'pairs': [], 'values': []}
            continue
        pairs = [(sample_chs[i], sample_chs[j])
                 for i in range(len(sample_chs))
                 for j in range(i + 1, len(sample_chs))]
        ch_prots = {
            ch: set(rdf.filter(pl.col('Channel') == ch)['Protein.Group'].drop_nulls().to_list())
            for ch in sample_chs
        }
        jacc = []
        for a, b in pairs:
            sa, sb = ch_prots.get(a, set()), ch_prots.get(b, set())
            inter  = len(sa & sb)
            union  = len(sa | sb)
            jacc.append(round(inter / union, 4) if union else 0)
        result[run] = {'pairs': [f'{a} / {b}' for a, b in pairs], 'values': jacc}
    return result


# ---------------------------------------------------------------------------
# HTML
# ---------------------------------------------------------------------------

def html_content(data: Dict) -> str:
    has_ms1  = data.get('has_ms1', True)
    has_ms2  = data.get('has_ms2', False)
    has_qval = data.get('has_qval', False)
    ms1_sections = """
<div class="section">
  <h3>MS1 Intensity Distribution</h3>
  <p class="help">log₁₀ MS1 intensity distribution per channel per run.</p>
  <div id="plex-ms1-grid" class="chart-grid"></div>
</div>
<div class="section">
  <h3>MS1 Quantification Variability</h3>
  <p class="help">Median MS1 CV per channel per run (normalised by per-precursor mean across channels and by
  channel loading). Lower = more reproducible quantification across peptides within a protein group.</p>
  <div id="plex-cv-grid" class="chart-grid"></div>
</div>""" if has_ms1 else ""
    ms2_sections = """
<div class="section">
  <h3>MS2 Intensity Distribution</h3>
  <p class="help">log₁₀ MS2 (Precursor.Quantity) distribution per channel per run.</p>
  <div id="plex-ms2-grid" class="chart-grid"></div>
</div>
<div class="section">
  <h3>MS2 Quantification Variability</h3>
  <p class="help">Median MS2 CV per channel per run using Precursor.Quantity (same normalisation as MS1).
  Lower = more reproducible MS2-based quantification.</p>
  <div id="plex-cv2-grid" class="chart-grid"></div>
</div>""" if has_ms2 else ""

    qval_section = """
<div class="section">
  <h3>Channel Q-Value Distribution</h3>
  <p class="help">Channel Q-Value distribution per channel for each run (values ≤ 0.01 after quality filter).
  Channels concentrated near 0 indicate high-confidence identifications.</p>
  <div id="plex-qval-grid" class="chart-grid"></div>
</div>""" if has_qval else ""

    return f"""
<div class="section">
  <div class="run-selector">
    <label>Samples per row:</label>
    <select id="plex-cols-sel" onchange="updatePlexLayout()">
      <option>1</option><option>2</option><option>3</option><option>4</option><option>5</option><option>6</option>
    </select>
  </div>
</div>
<div class="section">
  <h3>Identification Counts</h3>
  <p class="help">Per-sample precursor and protein group counts by channel (filtered: PEP ≤ 0.01 and Channel Q-Value ≤ 0.01).
  <em>In all channels</em> = identified in every channel of that run (intersection).
  <em>In any channel</em> = total unique across all channels (union).
  Protein groups use proteotypic peptides only when available.</p>
  <div class="chart-grid chart-grid-2">
    <div class="chart-box">
      <h4>Precursor IDs per Channel</h4>
      <button class="export-btn" onclick="exportChart('ch_plex_ids')">PNG</button>
      <canvas id="ch_plex_ids"></canvas>
    </div>
    <div class="chart-box">
      <h4>Protein Groups per Channel</h4>
      <button class="export-btn" onclick="exportChart('ch_plex_proteins')">PNG</button>
      <canvas id="ch_plex_proteins"></canvas>
    </div>
  </div>
</div>
{qval_section}
<div class="section">
  <h3>Relative Channel Intensity</h3>
  <p class="help">log₁₀ intensity ratio of each sample channel vs the highest-intensity (carrier) channel,
  computed on intersected precursors. Box = IQR (p25–p75), centre bar = median, whiskers = p5/p95.</p>
  <div id="plex-int-grid" class="chart-grid"></div>
</div>
{ms1_sections}
{ms2_sections}
<div class="section">
  <h3>Data Completeness — Precursors</h3>
  <p class="help">Pairwise precursor overlap between sample channels (carrier excluded).
  Jaccard Index = |A∩B| / |A∪B|; higher = more shared precursors between channels.</p>
  <div id="plex-jacc-grid" class="chart-grid"></div>
</div>
<div class="section">
  <h3>Data Completeness — Protein Groups</h3>
  <p class="help">Pairwise protein group overlap between sample channels (carrier excluded).
  Jaccard Index = |A∩B| / |A∪B|; higher = more shared protein identifications between channels.</p>
  <div id="plex-jap-grid" class="chart-grid"></div>
</div>"""


# ---------------------------------------------------------------------------
# JavaScript
# ---------------------------------------------------------------------------

def javascript() -> str:
    return """
const _PLEX_CH_COLS = ['#3498db','#e74c3c','#27ae60','#f39c12','#9b59b6','#1abc9c','#e67e22'];
const _PLEX_INT_COL = '#607d8b';
const _PLEX_UNI_COL = '#90a4ae';
const _PLEX_NO_MSG  = '<p class="help" style="padding:14px 0">No multiplexed samples are currently active. Enable them in the Files tab.</p>';

// Custom plugin: draws whiskers (p5/p95) and median bar on floating IQR bars.
const _plexBoxPlugin = {
  id: '_plexBox',
  afterDatasetsDraw(chart) {
    const {ctx, scales:{x,y}} = chart;
    if (!x || !y) return;
    chart.data.datasets.forEach((ds, di) => {
      if (!ds._whiskers) return;
      const meta = chart.getDatasetMeta(di);
      if (meta.hidden) return;
      const bcs = Array.isArray(ds.borderColor) ? ds.borderColor : null;
      ctx.save();
      ds._whiskers.forEach((w, bi) => {
        if (!w) return;
        const bar = meta.data[bi];
        if (!bar) return;
        const hw  = Math.min((bar.width || 14) * 0.35, 8);
        const col = bcs ? (bcs[bi] || '#333') : (ds.borderColor || '#333');
        const p5Y  = y.getPixelForValue(w[0]);
        const p25Y = y.getPixelForValue(w[1]);
        const p50Y = y.getPixelForValue(w[2]);
        const p75Y = y.getPixelForValue(w[3]);
        const p95Y = y.getPixelForValue(w[4]);
        ctx.strokeStyle = col;
        ctx.lineWidth = 1.5;
        ctx.beginPath();
        ctx.moveTo(bar.x, p25Y); ctx.lineTo(bar.x, p5Y);
        ctx.moveTo(bar.x - hw, p5Y); ctx.lineTo(bar.x + hw, p5Y);
        ctx.moveTo(bar.x, p75Y); ctx.lineTo(bar.x, p95Y);
        ctx.moveTo(bar.x - hw, p95Y); ctx.lineTo(bar.x + hw, p95Y);
        ctx.stroke();
        ctx.lineWidth = 2.5;
        ctx.beginPath();
        ctx.moveTo(bar.x - hw * 1.2, p50Y);
        ctx.lineTo(bar.x + hw * 1.2, p50Y);
        ctx.stroke();
      });
      ctx.restore();
    });
  }
};

function _plexActiveRuns(d) {
  const plexRuns = d.runs || [];
  return getActiveRuns().filter(r => plexRuns.includes(r));
}

const _PLEX_GRIDS = ['plex-qval-grid','plex-int-grid','plex-ms1-grid','plex-cv-grid',
                     'plex-ms2-grid','plex-cv2-grid','plex-jacc-grid','plex-jap-grid'];

function _plexCols(ar) {
  const sel = document.getElementById('plex-cols-sel');
  return sel ? (parseInt(sel.value) || Math.min(ar.length, 4)) : Math.min(ar.length, 4);
}

function initPlexDiagnostic() {
  const ar = _plexActiveRuns(DATA.plex_diagnostic);
  const d  = DATA.plex_diagnostic;
  if (!d || !ar.length) {
    _PLEX_GRIDS.forEach(id => { const el=document.getElementById(id); if(el) el.innerHTML=_PLEX_NO_MSG; });
    ['ch_plex_ids','ch_plex_proteins'].forEach(id => destroyChart(id));
    return;
  }
  const defaultCols = Math.min(ar.length, 4) || 1;
  const sel = document.getElementById('plex-cols-sel');
  if (sel) sel.value = String(defaultCols);
  _drawPlexIds(ar, d);
  _drawPlexProteins(ar, d);
  _drawPlexQval(ar, d);
  _drawPlexIntensity(ar, d);
  if (d.has_ms1) { _drawPlexMs1Dist(ar, d); _drawPlexCV(ar, d); }
  if (d.has_ms2) { _drawPlexMs2Dist(ar, d); _drawPlexCV2(ar, d); }
  _drawPlexJaccard(ar, d);
  _drawPlexJaccardProt(ar, d);
}

function updatePlexLayout() {
  const n    = parseInt(document.getElementById('plex-cols-sel').value);
  const cols = `repeat(${n}, 1fr)`;
  _PLEX_GRIDS.forEach(id => { const g=document.getElementById(id); if(g) g.style.gridTemplateColumns=cols; });
}

// ── Identification counts: x = samples, color = channels ──────────────────
// Uses global channel list for cross-run bar alignment; per-run intersection/union.
function _drawPlexCountBar(canvasId, ar, pd) {
  if (!pd || !pd.channels) return;
  const allChannels = pd.channels;
  const groups      = [...allChannels, 'In all channels', 'In any channel'];
  const colOf = x => {
    if (x === 'In all channels') return _PLEX_INT_COL;
    if (x === 'In any channel')  return _PLEX_UNI_COL;
    return _PLEX_CH_COLS[allChannels.indexOf(x) % _PLEX_CH_COLS.length];
  };
  barChart(canvasId, ar.map(r => label(r)),
    groups.map(x => ({
      label: x,
      data: ar.map(r => {
        const rd = pd[r] || {};
        if (x === 'In all channels') return rd.intersected || 0;
        if (x === 'In any channel')  return rd.union       || 0;
        return (rd.per_channel || {})[x] || 0;
      }),
      backgroundColor: colOf(x) + 'cc',
      borderColor:     colOf(x),
      borderWidth: 1,
    })),
    { aspectRatio: 1.2 }
  );
}

function _drawPlexIds(ar, d)      { _drawPlexCountBar('ch_plex_ids',      ar, d.plex_ids);      }
function _drawPlexProteins(ar, d) { _drawPlexCountBar('ch_plex_proteins', ar, d.plex_proteins); }

// ── Channel Q-Value: one chart per run, using per-run channels ─────────────
function _drawPlexQval(ar, d) {
  if (!d.plex_qval) return;
  const grid = document.getElementById('plex-qval-grid');
  if (!grid) return;
  ar.forEach((_, i) => destroyChart('plx_qv_' + i));
  grid.innerHTML = '';
  grid.style.gridTemplateColumns = `repeat(${_plexCols(ar)||1}, 1fr)`;
  ar.forEach((r, ri) => {
    const box = document.createElement('div');
    box.className = 'chart-box';
    box.innerHTML = `<h4>${label(r)}</h4>` +
      `<button class="export-btn" onclick="exportChart('plx_qv_${ri}')">PNG</button>` +
      `<div style="position:relative;height:300px"><canvas id="plx_qv_${ri}"></canvas></div>`;
    grid.appendChild(box);
  });
  ar.forEach((r, ri) => {
    const runData  = d.plex_qval[r] || {};
    const channels = runData.channels || [];
    const ds = channels.map((ch, ci) => {
      const h = runData[ch] || {centers:[], counts:[]};
      return { label:ch, data:h.centers.map((c,i)=>({x:h.counts[i], y:c})),
               borderColor:_PLEX_CH_COLS[ci%_PLEX_CH_COLS.length],
               backgroundColor:'transparent', borderWidth:2, pointRadius:0 };
    }).filter(s => s.data.length > 0);
    let qvXMax = 0;
    ds.forEach(s => s.data.forEach(p => { if (p.x > qvXMax) qvXMax = p.x; }));
    destroyChart('plx_qv_' + ri);
    const ctx = document.getElementById('plx_qv_' + ri);
    if (!ctx) return;
    CHARTS['plx_qv_'+ri] = new Chart(ctx, {
      type:'line', data:{datasets:ds},
      options:{
        responsive:true, maintainAspectRatio:false, parsing:false,
        plugins:{legend:{display:true, position:'top', labels:{boxWidth:10, font:{size:9}}}},
        elements:{point:{radius:0}, line:{borderWidth:1.5}},
        scales:{
          x:{type:'linear', beginAtZero:true, max:qvXMax, title:{display:true, text:'Count', font:{size:9}}},
          y:{type:'linear', min:0, title:{display:true, text:'Channel Q-Value', font:{size:9}}},
        },
      },
    });
  });
}

// ── Relative Channel Intensity: box plot per run, per-run channels ─────────
function _drawPlexIntensity(ar, d) {
  if (!d.plex_intensity) return;
  const globalRefCh = d.plex_intensity.ref_channel || '';
  const grid = document.getElementById('plex-int-grid');
  if (!grid) return;
  ar.forEach((_, i) => destroyChart('plx_in_' + i));
  grid.innerHTML = '';
  grid.style.gridTemplateColumns = `repeat(${_plexCols(ar)||1}, 1fr)`;
  ar.forEach((r, ri) => {
    const box = document.createElement('div');
    box.className = 'chart-box';
    box.innerHTML = `<h4>${label(r)}</h4>` +
      `<button class="export-btn" onclick="exportChart('plx_in_${ri}')">PNG</button>` +
      `<div style="position:relative;height:300px"><canvas id="plx_in_${ri}"></canvas></div>`;
    grid.appendChild(box);
  });
  ar.forEach((r, ri) => {
    const runEntry  = d.plex_intensity[r] || {};
    const nonRef    = runEntry.channels || [];
    const localRef  = runEntry.ref_channel || globalRefCh;
    const runData   = runEntry.values || {};
    if (!nonRef.length) return;
    const bgCols = nonRef.map((_, i) => _PLEX_CH_COLS[i%_PLEX_CH_COLS.length] + '88');
    const bdCols = nonRef.map((_, i) => _PLEX_CH_COLS[i%_PLEX_CH_COLS.length]);
    destroyChart('plx_in_' + ri);
    const ctx = document.getElementById('plx_in_' + ri);
    if (!ctx) return;
    CHARTS['plx_in_'+ri] = new Chart(ctx, {
      type: 'bar',
      data: {
        labels: nonRef,
        datasets: [{
          label: 'IQR',
          data: nonRef.map(ch => { const w=runData[ch]; return (w&&w[1]!=null&&w[3]!=null)?[w[1],w[3]]:null; }),
          backgroundColor: bgCols,
          borderColor:     bdCols,
          borderWidth: 1.5,
          barPercentage: 0.5,
          _whiskers: nonRef.map(ch => runData[ch] || null),
        }],
      },
      options: {
        responsive:true, maintainAspectRatio:false,
        plugins:{
          legend:{display:false},
          tooltip:{callbacks:{label:(item) => {
            const w = runData[nonRef[item.dataIndex]];
            return w ? `median: ${w[2].toFixed(2)}  IQR: [${w[1].toFixed(2)}, ${w[3].toFixed(2)}]` : '';
          }}},
        },
        scales:{
          x:{ticks:{font:{size:9}, maxRotation:0, minRotation:0}},
          y:{title:{display:true, text:`log₁₀(Int / ${localRef})`, font:{size:9}}, min:-3, max:1},
        },
      },
      plugins: [_plexBoxPlugin],
    });
  });
}

// ── MS1 / MS2 Intensity distributions: per-run channels ──────────────────
function _drawPlexIntDist(gridId, prefix, distKey, ar, yTitle) {
  if (!distKey || !distKey.centers) return;
  const grid = document.getElementById(gridId);
  if (!grid) return;
  ar.forEach((_, i) => destroyChart(prefix + i));
  grid.innerHTML = '';
  grid.style.gridTemplateColumns = `repeat(${_plexCols(ar)||1}, 1fr)`;
  const centers = distKey.centers || [];
  const yRange  = distKey.range   || [0, 1];
  ar.forEach((r, ri) => {
    const box = document.createElement('div');
    box.className = 'chart-box';
    box.innerHTML = `<h4>${label(r)}</h4>` +
      `<button class="export-btn" onclick="exportChart('${prefix}${ri}')">PNG</button>` +
      `<div style="position:relative;height:300px"><canvas id="${prefix}${ri}"></canvas></div>`;
    grid.appendChild(box);
  });
  ar.forEach((r, ri) => {
    const runData  = distKey[r] || {};
    const channels = runData.channels || [];
    const runMax   = runData['_max'] || 1;
    const medians  = runData['_medians'] || {};
    const {step, axMax: xMax} = _niceTicks(runMax);
    const histDs = channels.map((ch, ci) => {
      const counts = runData[ch] || [];
      return {
        label: ch,
        data: centers.map((c, i) => ({x: counts[i] || 0, y: c})),
        borderColor: _PLEX_CH_COLS[ci % _PLEX_CH_COLS.length],
        backgroundColor: 'transparent', borderWidth: 1.5,
        pointRadius: 2.5, pointHoverRadius: 4,
      };
    }).filter(s => s.data.length > 0);
    const medDs = channels.map((ch, ci) => {
      const med = medians[ch];
      if (med == null) return null;
      return {
        label: '_med_' + ch,
        data: [{x: 0, y: med}, {x: 1e9, y: med}],
        borderColor: _PLEX_CH_COLS[ci % _PLEX_CH_COLS.length],
        backgroundColor: 'transparent', borderWidth: 1,
        borderDash: [5, 5], pointRadius: 0,
      };
    }).filter(Boolean);
    const ds = [...histDs, ...medDs];
    const _medChannels = channels.slice();
    const _medLabelPlugin = {
      id: 'medLabels_' + prefix + ri,
      afterDraw(chart) {
        const {ctx: c2, chartArea, scales} = chart;
        if (!chartArea) return;
        _medChannels.forEach((ch, ci) => {
          const med = medians[ch];
          if (med == null) return;
          const medIdx = chart.data.datasets.findIndex(d => d.label === '_med_' + ch);
          if (medIdx >= 0 && chart.getDatasetMeta(medIdx).hidden) return;
          const yPx = scales.y.getPixelForValue(med);
          c2.save();
          c2.font = 'bold 9px sans-serif';
          c2.fillStyle = _PLEX_CH_COLS[ci % _PLEX_CH_COLS.length];
          c2.textAlign = 'left';
          c2.textBaseline = 'bottom';
          c2.fillText(med.toFixed(2), chartArea.left + 3, yPx - 2);
          c2.restore();
        });
      }
    };
    destroyChart(prefix + ri);
    const ctx = document.getElementById(prefix + ri);
    if (!ctx) return;
    CHARTS[prefix+ri] = new Chart(ctx, {
      type:'line', data:{datasets:ds},
      plugins:[_medLabelPlugin],
      options:{
        responsive:true, maintainAspectRatio:false, parsing:false,
        plugins:{
          legend:{
            display:true, position:'top',
            labels:{boxWidth:10, font:{size:9}, filter: item => !item.text.startsWith('_med_')},
            onClick(evt, item, legend) {
              Chart.defaults.plugins.legend.onClick(evt, item, legend);
              const chart = legend.chart;
              const medIdx = chart.data.datasets.findIndex(d => d.label === '_med_' + item.text);
              if (medIdx >= 0) {
                chart.getDatasetMeta(medIdx).hidden = chart.getDatasetMeta(item.datasetIndex).hidden;
                chart.update();
              }
            },
          },
        },
        elements:{line:{borderWidth:1.5}},
        scales:{
          x:{type:'linear', min:0, max:xMax,
             ticks:{stepSize:step, precision:0, maxRotation:0, minRotation:0, font:{size:9}},
             title:{display:true, text:'Count', font:{size:9}}},
          y:{type:'linear', min:yRange[0], max:yRange[1],
             title:{display:true, text:yTitle, font:{size:9}}},
        },
      },
    });
  });
}

function _drawPlexMs1Dist(ar, d) {
  _drawPlexIntDist('plex-ms1-grid', 'plx_m1_', d.plex_ms1_dist, ar, 'log₁₀(MS1 Intensity)');
}
function _drawPlexMs2Dist(ar, d) {
  if (!d.plex_ms2_dist || !d.plex_ms2_dist.centers) return;
  _drawPlexIntDist('plex-ms2-grid', 'plx_m2_', d.plex_ms2_dist, ar, 'log₁₀(MS2 Intensity)');
}

// ── CV charts: one bar chart per run, per-run channels ────────────────────
function _drawPlexCVGrid(gridId, prefix, cvKey, ar, yLabel) {
  if (!cvKey) return;
  const grid = document.getElementById(gridId);
  if (!grid) return;
  ar.forEach((_, i) => destroyChart(prefix + i));
  grid.innerHTML = '';
  grid.style.gridTemplateColumns = `repeat(${_plexCols(ar)||1}, 1fr)`;
  ar.forEach((r, ri) => {
    const box = document.createElement('div');
    box.className = 'chart-box';
    box.innerHTML = `<h4>${label(r)}</h4>` +
      `<button class="export-btn" onclick="exportChart('${prefix}${ri}')">PNG</button>` +
      `<div style="position:relative;height:300px"><canvas id="${prefix}${ri}"></canvas></div>`;
    grid.appendChild(box);
  });
  ar.forEach((r, ri) => {
    const runData  = cvKey[r] || {};
    const channels = runData.channels || [];
    destroyChart(prefix + ri);
    const ctx = document.getElementById(prefix + ri);
    if (!ctx) return;
    CHARTS[prefix+ri] = new Chart(ctx, {
      type:'bar',
      data:{
        labels: channels,
        datasets:[{
          label:'Median CV',
          data: channels.map(ch => { const v=runData[ch]; return v!=null?v:null; }),
          backgroundColor: channels.map((_,i) => _PLEX_CH_COLS[i%_PLEX_CH_COLS.length]+'cc'),
          borderColor:     channels.map((_,i) => _PLEX_CH_COLS[i%_PLEX_CH_COLS.length]),
          borderWidth:1,
        }],
      },
      options:{
        responsive:true, maintainAspectRatio:false,
        plugins:{legend:{display:false}},
        scales:{
          x:{ticks:{font:{size:9}, maxRotation:0, minRotation:0}},
          y:{beginAtZero:true, max:1, title:{display:true, text:yLabel, font:{size:9}}},
        },
      },
    });
  });
}

function _drawPlexCV(ar, d)  { _drawPlexCVGrid('plex-cv-grid',  'plx_cv_', d.plex_cv,     ar, 'Median CV'); }
function _drawPlexCV2(ar, d) { _drawPlexCVGrid('plex-cv2-grid', 'plx_c2_', d.plex_cv_ms2, ar, 'Median CV'); }

// ── Data Completeness: Jaccard per run, per-run pairs ─────────────────────
function _drawPlexJaccardGrid(gridId, prefix, dataKey, ar) {
  if (!dataKey) return;
  const grid = document.getElementById(gridId);
  if (!grid) return;
  const cols = _plexCols(ar);
  ar.forEach((_, i) => destroyChart(prefix + i));
  grid.innerHTML = '';
  grid.style.gridTemplateColumns = `repeat(${cols||1}, 1fr)`;
  ar.forEach((r, ri) => {
    const box = document.createElement('div');
    box.className = 'chart-box';
    box.innerHTML = `<h4>${label(r)}</h4>` +
      `<button class="export-btn" onclick="exportChart('${prefix}${ri}')">PNG</button>` +
      `<div style="position:relative;height:300px"><canvas id="${prefix}${ri}"></canvas></div>`;
    grid.appendChild(box);
  });
  ar.forEach((r, ri) => {
    const runData = dataKey[r] || {};
    const pairs   = runData.pairs  || [];
    const vals    = runData.values || [];
    if (!pairs.length) return;
    destroyChart(prefix + ri);
    const ctx = document.getElementById(prefix + ri);
    if (!ctx) return;
    CHARTS[prefix+ri] = new Chart(ctx, {
      type:'bar',
      data:{
        labels: pairs,
        datasets:[{
          label:'Jaccard',
          data: vals,
          backgroundColor: pairs.map((_,i) => _PLEX_CH_COLS[i%_PLEX_CH_COLS.length]+'cc'),
          borderColor:     pairs.map((_,i) => _PLEX_CH_COLS[i%_PLEX_CH_COLS.length]),
          borderWidth:1,
        }],
      },
      options:{
        responsive:true, maintainAspectRatio:false,
        plugins:{legend:{display:false}},
        scales:{
          x:{ticks:{font:{size:9}, maxRotation:45, minRotation:0}},
          y:{beginAtZero:true, max:1, title:{display:true, text:'Jaccard Index', font:{size:9}}},
        },
      },
    });
  });
}

function _drawPlexJaccard(ar, d) {
  _drawPlexJaccardGrid('plex-jacc-grid', 'plx_jc_', d.plex_jaccard, ar);
}
function _drawPlexJaccardProt(ar, d) {
  _drawPlexJaccardGrid('plex-jap-grid', 'plx_jp_', d.plex_jaccard_prot, ar);
}
"""
