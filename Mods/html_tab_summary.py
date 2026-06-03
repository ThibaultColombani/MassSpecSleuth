# Created by @thibaultcolombani with the help of Claude Code (Anthropic)
"""html_tab_summary.py — Summary tab: per-run seqcharges, peptides, protein groups."""
from __future__ import annotations
from typing import Dict, List
import polars as pl
from Mods.enzyme import ENZYMES

TAB_ID    = 'summary'
TAB_LABEL = 'Summary'


def is_visible(processor) -> bool:
    return True


def _nu(frame, col: str) -> int:
    """Unique count of col in frame; 0 if frame is None or col absent."""
    if frame is None or col not in frame.columns:
        return 0
    return int(frame[col].n_unique())


def _enzyme_cterm_aas(processor) -> List[str]:
    """C-terminal AAs for the active enzyme filter: R before K, then alphabetical."""
    enzymes = getattr(processor, 'enzymes', None)
    if not enzymes:
        return []
    cterm: set = set()
    for name in enzymes:
        cterm |= ENZYMES.get(name.lower(), {}).get('cterm', set())
    priority = [aa for aa in ('R', 'K') if aa in cterm]
    return priority + sorted(aa for aa in cterm if aa not in ('R', 'K'))


def compute(processor) -> Dict:
    df   = processor.report
    runs = processor.runs()

    has_cterm     = 'C_term_aa' in df.columns
    breakdown_aas = _enzyme_cterm_aas(processor) if has_cterm else []

    rows: List[Dict] = []
    for run in runs:
        rdf = df.filter(pl.col('Raw file') == run)
        row: Dict = {
            'run':        run,
            'seqcharges': _nu(rdf, 'seqcharge'),
            'peptides':   _nu(rdf, 'Modified sequence'),
            'proteins':   _nu(rdf, 'Protein.Group'),
        }
        for aa in breakdown_aas:
            adf = rdf.filter(pl.col('C_term_aa') == aa)
            row[f'seqcharges_{aa.lower()}'] = _nu(adf, 'seqcharge')
            row[f'peptides_{aa.lower()}']   = _nu(adf, 'Modified sequence')
        rows.append(row)

    return {'rows': rows, 'breakdown_aas': breakdown_aas}


def html_content(data: Dict) -> str:
    bk = data.get('breakdown_aas', [])
    if bk:
        help_text = (f"Per-run unique counts at FDR ≤ 1%. "
                     f"All = all precursors; {', '.join(bk)} = C-terminal residues from enzyme filter.")
    else:
        help_text = "Per-run unique counts at FDR ≤ 1%. Apply --enzyme to see C-terminal breakdown."

    return f"""
<div class="section">
  <h3>Experiment Overview</h3>
  <div class="summary-grid" id="kpi-grid"></div>
</div>
<div class="section">
  <h3>Per-Run Counts</h3>
  <p class="help">{help_text}</p>
  <div class="chart-grid chart-grid-3">
    <div class="chart-box">
      <h4>Unique Seqcharges (seq+z)</h4>
      <button class="export-btn" onclick="exportChart('ch_seq')">PNG</button>
      <canvas id="ch_seq"></canvas>
    </div>
    <div class="chart-box">
      <h4>Unique Peptides (modified seq)</h4>
      <button class="export-btn" onclick="exportChart('ch_sum_pep')">PNG</button>
      <canvas id="ch_sum_pep"></canvas>
    </div>
    <div class="chart-box">
      <h4>Protein Groups</h4>
      <button class="export-btn" onclick="exportChart('ch_pg')">PNG</button>
      <canvas id="ch_pg"></canvas>
    </div>
  </div>
</div>"""


def javascript() -> str:
    return """
function initSummary() {
  const ar = getActiveRuns();
  const s  = ar.map(r => DATA.summary.rows.find(d => d.run === r)).filter(Boolean);
  const bkAAs = DATA.summary.breakdown_aas;

  const _AA_C = {R:'#e74c3c',K:'#27ae60',E:'#f39c12',D:'#3498db',
                 F:'#9b59b6',Y:'#1abc9c',W:'#e67e22',L:'#2980b9',M:'#c0392b'};
  const _aaCol = aa => _AA_C[aa] || '#95a5a6';

  const maxSeq = Math.max(0, ...s.map(d => d.seqcharges));
  const maxPep = Math.max(0, ...s.map(d => d.peptides));
  const maxPG  = Math.max(0, ...s.map(d => d.proteins));

  document.getElementById('kpi-grid').innerHTML = `
    <div class="kpi"><div class="val">${ar.length}</div><div class="lbl">Active runs</div></div>
    <div class="kpi"><div class="val">${maxSeq.toLocaleString()}</div><div class="lbl">Max seqcharges</div></div>
    <div class="kpi"><div class="val">${maxPep.toLocaleString()}</div><div class="lbl">Max peptides</div></div>
    <div class="kpi"><div class="val">${maxPG.toLocaleString()}</div><div class="lbl">Max protein groups</div></div>`;

  function makeDs(keyAll, keyAA, colorAll) {
    const ds = [{ label:'All', data:s.map(d=>d[keyAll]),
                  backgroundColor:colorAll+'cc', borderColor:colorAll, borderWidth:1 }];
    bkAAs.forEach(aa => {
      const c = _aaCol(aa);
      ds.push({ label: aa, data: s.map(d => d[keyAA + '_' + aa.toLowerCase()] || 0),
                backgroundColor: c+'cc', borderColor: c, borderWidth: 1 });
    });
    return ds;
  }

  barChart('ch_seq',     s.map(d => label(d.run)), makeDs('seqcharges', 'seqcharges', '#3498db'), { aspectRatio: 0.67 });
  barChart('ch_sum_pep', s.map(d => label(d.run)), makeDs('peptides',   'peptides',   '#e67e22'), { aspectRatio: 0.67 });
  barChart('ch_pg', s.map(d => label(d.run)), [{
    label:'Protein Groups', data:s.map(d=>d.proteins),
    backgroundColor:'#9b59b6cc', borderColor:'#7d3c98', borderWidth:1
  }], { aspectRatio: 0.67 });
}
"""
