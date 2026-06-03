# Created by @thibaultcolombani with the help of Claude Code (Anthropic)
"""html_tab_features.py — Feature Detection tab. Shown only when features.tsv was loaded."""
from __future__ import annotations
import math
from typing import Dict
import polars as pl
from Mods.html_utils import histogram

TAB_ID    = 'features'
TAB_LABEL = 'Feature Detection'


def is_visible(processor) -> bool:
    return processor.features is not None


def compute(processor) -> Dict:
    f = processor.features
    result: Dict = {}

    if 'charge' in f.columns:
        f2 = f.with_columns(
            pl.when(pl.col('charge') > 3).then(4).otherwise(pl.col('charge')).alias('charge')
        )
        counts = f2.group_by('charge').agg(pl.len().alias('n')).sort('charge')
        result['charge_dist'] = {str(int(r[0])): int(r[1]) for r in counts.iter_rows()}

    int_col = next((c for c in ('intensitySum', 'Intensity', 'intensity') if c in f.columns), None)
    if int_col:
        vals = f.filter(pl.col(int_col) > 0)[int_col].drop_nulls().to_list()
        result['intensity_hist'] = histogram([math.log10(float(v)) for v in vals if v > 0])

    return result


def html_content(data: Dict) -> str:
    return """
<div class="section">
  <h3>Feature Detection</h3>
  <p class="help">Metrics derived from features.tsv. Shown only when the file was loaded.</p>
  <div class="chart-grid chart-grid-2">
    <div class="chart-box">
      <h4>Feature charge distribution</h4>
      <button class="export-btn" onclick="exportChart('ch_feat_charge')">PNG</button>
      <canvas id="ch_feat_charge"></canvas>
    </div>
    <div class="chart-box">
      <h4>Feature intensity log₁₀</h4>
      <button class="export-btn" onclick="exportChart('ch_feat_int')">PNG</button>
      <canvas id="ch_feat_int"></canvas>
    </div>
  </div>
</div>"""


def javascript() -> str:
    return """
function initFeatures() {
  const f = DATA.features;

  if (f.charge_dist) {
    const keys = Object.keys(f.charge_dist).sort((a,b)=>+a-+b);
    const cols  = ['#3498db','#27ae60','#f39c12','#e74c3c'];
    barChart('ch_feat_charge', keys.map(k => k==='4'?'≥4':'z='+k), [{
      label:'Features',
      data: keys.map(k=>f.charge_dist[k]),
      backgroundColor: keys.map((_,i)=>cols[i%cols.length]+'cc')
    }]);
  }

  if (f.intensity_hist) {
    const h = f.intensity_hist;
    lineChart('ch_feat_int', [{
      label:'Features',
      data: h.centers.map((c,i) => ({x:c, y:h.counts[i]})),
      borderColor:'#9b59b6', backgroundColor:'transparent', borderWidth:2
    }], { x:{title:{display:true,text:'log₁₀(Feature Intensity)'}} });
  }
}
"""
