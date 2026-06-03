# Created by @thibaultcolombani with the help of Claude Code (Anthropic)
"""html_tab_identifications.py — IDs tab: PEP/Q.Value cumulative, missed cleavages, protein groups, charge.

q_col logic: PEP when available (takes priority); Q.Value when PEP absent;
'PEP or Q.Value' only for combined multi-engine reports where runs differ.
"""
from __future__ import annotations
import math
from typing import Dict, List, Optional, Set, Tuple
import polars as pl

TAB_ID    = 'identifications'
TAB_LABEL = 'IDs'

_TRYPSIN_LIKE = {'trypsin', 'trypsin/p'}  # enzymes that use the proline exception


def is_visible(processor) -> bool:
    return True


def _run_q_col(df: pl.DataFrame, run: str) -> Optional[str]:
    """Return the quality column to use for a specific run: PEP if non-null, else Q.Value."""
    rdf = df.filter(pl.col('Raw file') == run)
    if 'PEP' in rdf.columns and rdf['PEP'].drop_nulls().len() > 0:
        return 'PEP'
    if 'Q.Value' in rdf.columns and rdf['Q.Value'].drop_nulls().len() > 0:
        return 'Q.Value'
    return None


def _q_filter_expr(df: pl.DataFrame) -> Optional[pl.Expr]:
    """
    Polars filter expression: PEP < 0.01 for rows that have PEP,
    Q.Value < 0.01 for rows where PEP is null/absent.
    Handles single-engine (PEP only / Q.Value only) and mixed combined reports.
    """
    has_pep = 'PEP' in df.columns and df['PEP'].drop_nulls().len() > 0
    has_qv  = 'Q.Value' in df.columns and df['Q.Value'].drop_nulls().len() > 0
    if has_pep and has_qv:
        return (pl.col('PEP') < 0.01) | (pl.col('PEP').is_null() & (pl.col('Q.Value') < 0.01))
    if has_pep:
        return pl.col('PEP') < 0.01
    if has_qv:
        return pl.col('Q.Value') < 0.01
    return None


def compute(processor) -> Dict:
    df   = processor.report
    runs = processor.runs()

    # Determine the displayed quality column name.
    # Use PEP whenever it is available (takes priority over Q.Value).
    # 'PEP or Q.Value' only for genuinely mixed combined reports where some runs
    # use PEP and others fall back to Q.Value.
    per_run_q = [_run_q_col(df, r) for r in runs]
    unique_q  = set(c for c in per_run_q if c is not None)
    if len(unique_q) == 1:
        q_col: Optional[str] = next(iter(unique_q))
    elif unique_q:
        q_col = 'PEP or Q.Value'
    else:
        q_col = None

    has_pq = (('PG.Q.Value' in df.columns and df['PG.Q.Value'].drop_nulls().len() > 0) or
              ('Protein.Q.Value' in df.columns and df['Protein.Q.Value'].drop_nulls().len() > 0))
    return {
        'q_col':          q_col,
        'has_pq':         has_pq,
        'pep_cumul':      _pep_cumulative(df, runs),
        'miss_cleavage':  _miss_cleavage(processor, df, runs),
        'protein_groups': _protein_groups(df, runs),
        'peptide_ids':    _peptide_ids(df, runs),
    }


# ── data helpers ──────────────────────────────────────────────────────────────

def _pep_cumulative(df: pl.DataFrame, runs: List[str]) -> Dict:
    result = {}
    for run in runs:
        col = _run_q_col(df, run)
        if not col:
            result[run] = {'log_vals': [], 'counts': [], 'col': None}
            continue
        rdf = (df.filter(pl.col('Raw file') == run)
                 .select(col).drop_nulls().sort(col))
        vals = rdf[col].to_list()
        n = len(vals)
        if n == 0:
            result[run] = {'log_vals': [], 'counts': [], 'col': col}
            continue
        step = max(1, n // 300)
        log_vals, counts = [], []
        for i, v in enumerate(vals[::step]):
            if v > 0:
                log_vals.append(round(math.log10(v), 4))
                counts.append(i * step + 1)
        result[run] = {'log_vals': log_vals, 'counts': counts, 'col': col}
    return result


def _mc_rule(processor, df: pl.DataFrame) -> Tuple[Set[str], Set[str], bool, str]:
    """
    Returns (cterm_res, nterm_res, proline_exc, rule_label) for MC counting.

    With enzyme  → cleavage residues from ENZYMES catalog.
    No enzyme    → infer from C-term AA frequency (≥2% of peptides).
    proline_exc  → True only for trypsin/trypsin-p, or when K+R are both inferred.
    """
    from Mods.enzyme import ENZYMES
    enzymes = getattr(processor, 'enzymes', None)

    if enzymes:
        cterm_res: Set[str] = set()
        nterm_res: Set[str] = set()
        for name in enzymes:
            enz = ENZYMES.get(name.lower(), {})
            cterm_res |= enz.get('cterm', set())
            nterm_res |= enz.get('nterm', set())
        proline_exc = any(e.lower() in _TRYPSIN_LIKE for e in enzymes)
        label = '+'.join(enzymes)
        return cterm_res, nterm_res, proline_exc, label

    # Infer from data C-term distribution
    if 'Sequence' not in df.columns:
        return set(), set(), False, 'unknown'
    vc = df['Sequence'].str.slice(-1).drop_nulls().value_counts().sort('count', descending=True)
    total = vc['count'].sum() or 1
    cterm_res = {str(row[0]) for row in vc.iter_rows() if row[1] / total >= 0.02}
    proline_exc = {'K', 'R'}.issubset(cterm_res)
    label = 'inferred (' + ', '.join(sorted(cterm_res)) + ')'
    return cterm_res, set(), proline_exc, label


def _miss_cleavage(processor, df: pl.DataFrame, runs: List[str]) -> Dict:
    """MC count per peptide using enzyme-aware cleavage residues (quality < 0.01, Intensity > 0)."""
    if 'Sequence' not in df.columns:
        return {}

    cterm_res, nterm_res, proline_exc, rule_label = _mc_rule(processor, df)
    if not cterm_res and not nterm_res:
        return {'rule': rule_label}

    base = df
    q_expr = _q_filter_expr(df)
    if q_expr is not None:
        base = base.filter(q_expr)
    if 'Intensity' in df.columns and df['Intensity'].max() > 0:
        base = base.filter(pl.col('Intensity') > 0)
    elif 'Precursor.Quantity' in df.columns and df['Precursor.Quantity'].max() > 0:
        base = base.filter(pl.col('Precursor.Quantity') > 0)

    # Build MC expression: count cleavage residues at internal positions
    mc_expr = pl.lit(0).cast(pl.Int32)
    for aa in sorted(cterm_res):
        # aa followed by any char = aa NOT at C-terminus
        term = pl.col('Sequence').str.count_matches(f'{aa}.')
        if proline_exc:
            term = term - pl.col('Sequence').str.count_matches(f'{aa}P')
        mc_expr = mc_expr + term
    for aa in sorted(nterm_res):
        # any char followed by aa = aa NOT at N-terminus
        mc_expr = mc_expr + pl.col('Sequence').str.count_matches(f'.{aa}')

    base = base.with_columns([
        mc_expr.alias('_mc'),
        pl.col('Sequence').str.slice(-1).alias('_cterm'),
    ])

    vc = base.filter(pl.col('_cterm').is_not_null()).group_by('_cterm').agg(pl.len().alias('n'))
    cterm_counts = {r[0]: r[1] for r in vc.iter_rows()}
    total_pep = sum(cterm_counts.values()) or 1
    cterm_aas = [aa for aa, _ in sorted(cterm_counts.items(), key=lambda x: -x[1])
                 if cterm_counts[aa] / total_pep >= 0.02]

    combined, by_cterm = {}, {aa: {} for aa in cterm_aas}
    for run in runs:
        rdf = base.filter(pl.col('Raw file') == run)
        counts = rdf.group_by('_mc').agg(pl.len().alias('n')).sort('_mc')
        combined[run] = {str(int(r[0])): int(r[1]) for r in counts.iter_rows()}
        for aa in cterm_aas:
            ac = rdf.filter(pl.col('_cterm') == aa).group_by('_mc').agg(pl.len().alias('n')).sort('_mc')
            by_cterm[aa][run] = {str(int(r[0])): int(r[1]) for r in ac.iter_rows()}

    return {'combined': combined, 'by_cterm': by_cterm, 'cterm_aas': cterm_aas, 'rule': rule_label}


def _protein_groups(df: pl.DataFrame, runs: List[str]) -> List[Dict]:
    if 'Intensity' in df.columns and df['Intensity'].max() > 0:
        base = df.filter(pl.col('Intensity') > 0)
    elif 'Precursor.Quantity' in df.columns and df['Precursor.Quantity'].max() > 0:
        base = df.filter(pl.col('Precursor.Quantity') > 0)
    else:
        base = df
    rows = []
    for run in runs:
        rdf = base.filter(pl.col('Raw file') == run)
        proto = rdf.filter((pl.col('Proteotypic') == 1) | pl.col('Proteotypic').is_null()) if 'Proteotypic' in rdf.columns else rdf
        def _pids(f):
            if 'Protein.Ids' in f.columns:
                v = f['Protein.Ids'].drop_nulls()
                if v.len() > 0: return v.n_unique()
            if 'Protein.Group' in f.columns:
                return f['Protein.Group'].drop_nulls().n_unique()
            return 0
        def _pgs(f):
            return f['Protein.Group'].drop_nulls().n_unique() if 'Protein.Group' in f.columns else 0

        # proteins_q: filter by Protein.Q.Value when the run actually has values
        has_pqv = ('Protein.Q.Value' in proto.columns and
                   proto['Protein.Q.Value'].drop_nulls().len() > 0)
        proteins_q_f = proto.filter(pl.col('Protein.Q.Value') < 0.01) if has_pqv else proto

        # pg_q: prefer PG.Q.Value (DIA-NN); fall back to Protein.Q.Value for engines
        # that lack PG-level Q-values (e.g. Jmod) so pg_q == proteins_q for those runs.
        has_pgq = ('PG.Q.Value' in proto.columns and
                   proto['PG.Q.Value'].drop_nulls().len() > 0)
        if has_pgq:
            pg_q_f = proto.filter((pl.col('PG.Q.Value') < 0.01) | pl.col('PG.Q.Value').is_null())
        elif has_pqv:
            pg_q_f = proto.filter(pl.col('Protein.Q.Value') < 0.01)
        else:
            pg_q_f = proto

        rows.append({
            'run':        run,
            'proteins':   _pids(proto),
            'proteins_q': _pids(proteins_q_f),
            'pg':         _pgs(proto),
            'pg_q':       _pgs(pg_q_f),
        })
    return rows


def _peptide_ids(df: pl.DataFrame, runs: List[str]) -> Dict:
    """Seqcharges, precursors, and charge distribution at quality threshold."""
    q_expr = _q_filter_expr(df)
    base = df.filter(q_expr) if q_expr is not None else df
    counts = []
    for run in runs:
        rdf = base.filter(pl.col('Raw file') == run)
        sc = (rdf.select(['Sequence', 'Charge']).drop_nulls().n_unique()
              if all(c in rdf.columns for c in ['Sequence', 'Charge']) else 0)
        pid = rdf['Precursor.Id'].n_unique() if 'Precursor.Id' in rdf.columns else sc
        counts.append({'run': run, 'seqcharges': sc, 'precursors': pid})
    charge = {}
    if 'Charge' in base.columns:
        base2 = base.with_columns(
            pl.when(pl.col('Charge') > 3).then(4).otherwise(pl.col('Charge')).alias('Charge')
        )
        for run in runs:
            c = (base2.filter(pl.col('Raw file') == run)
                      .group_by('Charge').agg(pl.len().alias('n')).sort('Charge'))
            charge[run] = {str(int(r[0])): int(r[1]) for r in c.iter_rows()}
    return {'counts': counts, 'charge': charge}


# ── HTML / JS ────────────────────────────────────────────────────────────────

def html_content(data: Dict) -> str:
    q_col    = data.get('q_col') or 'PEP'
    mc_rule  = data.get('miss_cleavage', {}).get('rule', '')
    mc_rule_html = f' Cleavage rule: <em>{mc_rule}</em>.' if mc_rule else ''
    return f"""
<div class="section">
  <h3>Precursor IDs</h3>
  <p class="help">Precursors ranked by {q_col} (ascending). X = cumulative precursor count, Y = {q_col} (log₁₀).
  <em>{q_col} used as the quality metric (PEP when available, Q.Value otherwise).</em></p>
  <div class="chart-grid chart-grid-1">
    <div class="chart-box">
      <h4>{q_col} cumulative</h4>
      <button class="export-btn" onclick="exportChart('ch_pep')">PNG</button>
      <canvas id="ch_pep"></canvas>
    </div>
  </div>
</div>
<div class="section">
  <h3>Missed Cleavages</h3>
  <p class="help">Internal cleavage sites within peptides ({q_col} &lt; 0.01, Intensity &gt; 0).{mc_rule_html}</p>
  <div id="mc_grid" class="chart-grid"></div>
</div>
<div class="section">
  <h3>Peptide Identifications &amp; Charge Distribution</h3>
  <p class="help">Seqcharges, precursors, and charge state distribution at {q_col} &lt; 0.01.</p>
  <div class="chart-grid chart-grid-2">
    <div class="chart-box">
      <h4>Seqcharges &amp; Precursors per run</h4>
      <button class="export-btn" onclick="exportChart('ch_pep_ids')">PNG</button>
      <canvas id="ch_pep_ids"></canvas>
    </div>
    <div class="chart-box">
      <h4>Charge state distribution ({q_col} &lt; 0.01)</h4>
      <button class="export-btn" onclick="exportChart('ch_charge')">PNG</button>
      <canvas id="ch_charge"></canvas>
    </div>
  </div>
</div>
<div class="section">
  <h3>Protein Identifications</h3>
  <div class="chart-grid chart-grid-1">
    <div class="chart-box">
      <h4>Protein IDs per run (4 stringencies)</h4>
      <button class="export-btn" onclick="exportChart('ch_prot')">PNG</button>
      <canvas id="ch_prot"></canvas>
    </div>
  </div>
</div>"""


def javascript() -> str:
    return """
function initIdentifications() {
  const ar = getActiveRuns();
  const d  = DATA.identifications;

  // Quality cumulative (PEP or Q.Value per run) — x=rank, y=log10(val) linear
  const qCol = d.q_col || 'PEP';
  const isMixed = qCol === 'PEP or Q.Value';
  const pepDS = ar.map(r => {
    const c = d.pep_cumul[r] || {log_vals:[], counts:[], col:null};
    const runLabel = isMixed && c.col ? `${label(r)} [${c.col}]` : label(r);
    return {
      label: runLabel,
      data: c.counts.map((n,i) => ({x:n, y:c.log_vals[i]})),
      borderColor: color(r), backgroundColor: 'transparent', borderWidth: 1.5
    };
  });
  let pepXMax = 0, pepYMin = 0, pepYMax = 0;
  ar.forEach(r => {
    const c = d.pep_cumul[r] || {log_vals:[], counts:[]};
    if (c.counts.length) pepXMax = Math.max(pepXMax, c.counts[c.counts.length - 1]);
    if (c.log_vals.length) {
      pepYMin = Math.min(pepYMin, Math.min(...c.log_vals));
      pepYMax = Math.max(pepYMax, Math.max(...c.log_vals));
    }
  });
  destroyChart('ch_pep');
  const pepCtx = document.getElementById('ch_pep');
  if (pepCtx) {
    CHARTS['ch_pep'] = new Chart(pepCtx, {
      type: 'line',
      data: { datasets: pepDS },
      options: {
        responsive: true, maintainAspectRatio: true, aspectRatio: 2.5,
        parsing: false,
        plugins: {
          legend: { display: true, position: 'top',
            labels: { boxWidth: 12, font: { size: 10 } }
          }
        },
        elements: { point: { radius: 0 }, line: { borderWidth: 1.5 } },
        scales: {
          x: { type: 'linear', min: 0, max: pepXMax, title: { display: true, text: 'Precursors' } },
          y: { type: 'linear', min: pepYMin, max: pepYMax,
               title: { display: true, text: qCol },
               ticks: { callback: v => { if (!Number.isInteger(v)) return null; const sup={'0':'⁰','1':'¹','2':'²','3':'³','4':'⁴','5':'⁵','6':'⁶','7':'⁷','8':'⁸','9':'⁹','-':'⁻'}; return '10'+String(v).split('').map(c=>sup[c]||c).join(''); } } }
        }
      }
    });
  }

  // Missed cleavages — all + per C-terminal AA in one row
  const mc = d.miss_cleavage;
  const mcCombined = mc.combined || {};
  const ctermAAs = mc.cterm_aas || [];
  const mcCols = ['#27ae60','#f39c12','#e74c3c','#9b59b6','#3498db'];
  const mcGrid = document.getElementById('mc_grid');
  if (mcGrid) {
    mcGrid.querySelectorAll('canvas').forEach(c => destroyChart(c.id));
    const allSlots = ['_all', ...ctermAAs];
    mcGrid.style.gridTemplateColumns = `repeat(${allSlots.length}, 1fr)`;
    mcGrid.innerHTML = allSlots.map(aa => {
      const cid   = aa === '_all' ? 'ch_mc' : 'ch_mc_' + aa;
      const title = aa === '_all' ? 'All peptides' : aa + '-ending peptides';
      return `<div class="chart-box"><h4>${title}</h4>` +
             `<button class="export-btn" onclick="exportChart('${cid}')">PNG</button>` +
             `<canvas id="${cid}"></canvas></div>`;
    }).join('');

    function _mcBar(cid, data) {
      const keys = [...new Set(ar.flatMap(r => Object.keys(data[r]||{})))].sort((a,b)=>+a-+b);
      const tots = {};
      ar.forEach(r => { tots[r] = keys.reduce((s,k) => s + ((data[r]||{})[k]||0), 0); });
      barChart(cid, ar.map(r=>label(r)),
        keys.map((k,i) => ({
          label: 'MC='+k,
          data: ar.map(r => parseFloat((100*((data[r]||{})[k]||0)/(tots[r]||1)).toFixed(2))),
          backgroundColor: mcCols[i % mcCols.length]+'cc'
        })),
        { aspectRatio: 0.8,
          scales: { y: { title: { display: true, text: 'Precursors (%)' }, max: 100 } },
          plugins: { tooltip: { mode: 'index', callbacks: { label: ctx => {
            const k = keys[ctx.datasetIndex];
            const run = ar[ctx.dataIndex];
            const raw = (data[run]||{})[k] || 0;
            return `MC=${k}: ${ctx.parsed.y.toFixed(1)}% (${raw.toLocaleString()} precursors)`;
          }}}}
        }
      );
    }

    _mcBar('ch_mc', mcCombined);
    ctermAAs.forEach(aa => _mcBar('ch_mc_' + aa, (mc.by_cterm || {})[aa] || {}));
  }

  // Peptide IDs — seqcharges + precursors at PEP < 0.01
  const pi = d.peptide_ids;
  const piAr = ar.map(r => (pi.counts||[]).find(x => x.run === r)).filter(Boolean);
  barChart('ch_pep_ids', piAr.map(x=>label(x.run)), [
    {label:'Seqcharges', data:piAr.map(x=>x.seqcharges), backgroundColor:'#3498dbcc'},
    {label:'Precursors', data:piAr.map(x=>x.precursors), backgroundColor:'#27ae60cc'},
  ], { aspectRatio: 1.25 });

  // Charge distribution (PEP < 0.01) — stacked bar
  const chKeys = [...new Set(ar.flatMap(r => Object.keys((pi.charge||{})[r]||{})))].sort((a,b)=>+a-+b);
  const chCols = ['#f1c40f','#e67e22','#e74c3c','#922b21'];
  barChart('ch_charge', ar.map(r=>label(r)),
    chKeys.map((k,i) => ({
      label: k==='4' ? '≥4' : 'z='+k,
      data: ar.map(r => ((pi.charge||{})[r]||{})[k] || 0),
      backgroundColor: chCols[i % chCols.length]+'cc'
    })),
    { aspectRatio: 1.25 }
  );

  // Protein groups — grouped bar (full width)
  const pg    = ar.map(r => d.protein_groups.find(x => x.run === r)).filter(Boolean);
  const hasPq = d.has_pq;
  const pgDs  = [
    {label:'Proteins',        data:pg.map(x=>x.proteins),   backgroundColor:'#3498dbcc'},
    {label:'Protein Groups',  data:pg.map(x=>x.pg),         backgroundColor:'#9b59b6cc'},
  ];
  if (hasPq) {
    pgDs.splice(1, 0, {label:'Proteins (q<1%)',       data:pg.map(x=>x.proteins_q), backgroundColor:'#27ae60cc'});
    pgDs.push(   {label:'Protein Groups (q<1%)',       data:pg.map(x=>x.pg_q),       backgroundColor:'#e67e22cc'});
  }
  barChart('ch_prot', pg.map(x=>label(x.run)), pgDs, { aspectRatio: 2.5 });
}
"""
