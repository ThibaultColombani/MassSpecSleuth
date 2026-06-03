# Created by @thibaultcolombani with the help of Claude Code (Anthropic)
"""
html_utils.py - Shared utilities for MassSleuth HTML tab modules.

Provides: colour palette, histogram helper, shared CSS, shared Chart.js helpers.
"""

from __future__ import annotations
import math
from typing import List, Dict

_PALETTE = [
    '#3498db', '#e74c3c', '#27ae60', '#f39c12', '#9b59b6',
    '#1abc9c', '#e67e22', '#2980b9', '#c0392b', '#16a085',
]


def run_colors(runs: List[str]) -> Dict[str, str]:
    return {r: _PALETTE[i % len(_PALETTE)] for i, r in enumerate(runs)}


def histogram(values: List[float], n_bins: int = 50) -> Dict:
    """Return {centers, counts} with 1st/99th-percentile capping."""
    if not values:
        return {'centers': [], 'counts': []}
    sv = sorted(values)
    n = len(sv)
    p1  = sv[max(0, int(0.01 * n))]
    p99 = sv[min(n - 1, int(0.99 * n))]
    if p1 >= p99:
        return {'centers': [round(p1, 4)], 'counts': [n]}
    bw = (p99 - p1) / n_bins
    counts = [0] * n_bins
    for v in values:
        if v < p1 or v > p99:
            continue
        b = min(n_bins - 1, int((v - p1) / bw))
        counts[b] += 1
    return {
        'centers': [round(p1 + (i + 0.5) * bw, 4) for i in range(n_bins)],
        'counts':  counts,
    }


# ─── Shared CSS ───────────────────────────────────────────────────────────────

SHARED_CSS = """
*{box-sizing:border-box;margin:0;padding:0}
body{font-family:'Segoe UI',sans-serif;background:#f0f2f5;color:#333}
.container{max-width:1600px;margin:0 auto;background:#fff;border-radius:10px;
  box-shadow:0 4px 12px rgba(0,0,0,.12);overflow:hidden}
.header{background:linear-gradient(135deg,#1a237e 0%,#283593 60%,#3949ab 100%);
  color:#fff;padding:28px 36px}
.header h1{font-size:2em;font-weight:300;margin-bottom:6px}
.header p{opacity:.85;font-size:.95em}
.tabs{display:flex;flex-wrap:wrap;background:#1e2a3a;border-bottom:3px solid #0d1621}
.tab{background:none;border:none;color:#b0bec5;padding:14px 22px;cursor:pointer;
  font-size:14px;font-weight:500;border-bottom:3px solid transparent;transition:.2s}
.tab:hover{background:#263545;color:#fff}
.tab.active{background:#263545;color:#64b5f6;border-bottom-color:#64b5f6}
.tab-content{display:none;padding:28px}
.tab-content.active{display:block}
.section{background:#f8f9fa;border-radius:8px;padding:20px;margin-bottom:24px;
  border:1px solid #e3e6ea}
.section h3{color:#1a237e;font-size:1.1em;border-bottom:2px solid #3949ab;
  padding-bottom:8px;margin-bottom:16px}
.section p.help{color:#666;font-size:.85em;margin-top:-10px;margin-bottom:12px;
  font-style:italic}
.chart-grid{display:grid;gap:20px;grid-template-columns:repeat(auto-fit,minmax(340px,1fr))}
.chart-grid-1{grid-template-columns:1fr}
.chart-grid-2{grid-template-columns:repeat(2,1fr)}
.chart-grid-3{grid-template-columns:repeat(3,1fr)}
.chart-box{background:#fff;border-radius:8px;padding:16px;
  box-shadow:0 2px 6px rgba(0,0,0,.08);position:relative;min-width:0;overflow:hidden}
.chart-box h4{text-align:center;font-size:.475em;color:#444;margin-bottom:10px;
  display:-webkit-box;-webkit-line-clamp:2;-webkit-box-orient:vertical;line-clamp:2;overflow:hidden}
.chart-box canvas{display:block;width:100%!important}
.export-btn{position:absolute;top:8px;right:8px;background:#3949ab;color:#fff;
  border:none;padding:3px 9px;border-radius:4px;font-size:11px;cursor:pointer}
.export-btn:hover{background:#283593}
.run-selector{margin-bottom:16px;display:flex;align-items:center;gap:10px;flex-wrap:wrap}
.run-selector label{font-size:.9em;font-weight:600;color:#444}
.run-selector select{padding:6px 10px;border:1px solid #ccc;border-radius:4px;font-size:.9em}
.summary-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(160px,1fr));
  gap:12px;margin-bottom:20px}
.kpi{background:#fff;border-radius:8px;padding:14px;text-align:center;
  box-shadow:0 2px 5px rgba(0,0,0,.07);border-top:3px solid #3949ab}
.kpi .val{font-size:1.8em;font-weight:700;color:#1a237e}
.kpi .lbl{font-size:.8em;color:#666;margin-top:4px}
.file-row{display:flex;align-items:center;gap:10px;padding:8px;background:#f8f9fa;
  border:1px solid #e3e6ea;border-radius:6px;margin-bottom:6px}
.file-row:hover{border-color:#3949ab;background:#e8eaf6}
.file-orig{flex:0 0 45%;font-weight:600;color:#1a237e;font-size:.85em;word-break:break-all}
.file-rename{flex:1;padding:6px 10px;border:1px solid #ccc;border-radius:4px;font-size:.85em}
.file-toggle{padding:6px 12px;border:none;border-radius:4px;cursor:pointer;
  font-size:.82em;font-weight:500}
.active-yes{background:#c8e6c9;color:#1b5e20}
.active-no{background:#ffcdd2;color:#b71c1c}
.drag-handle{cursor:grab;color:#aaa;font-size:1.1em;padding:0 4px;flex-shrink:0;user-select:none}
.drag-handle:active{cursor:grabbing}
.file-row.dragging{opacity:.35}
.file-row.drag-over-top{border-top:2px solid #3949ab}
.file-row.drag-over-bottom{border-bottom:2px solid #3949ab}
#file-status{min-height:1.4em;text-align:center;font-size:.88em;font-weight:500;
  margin-top:8px;padding:5px 12px;border-radius:4px}
.status-success{color:#1b5e20;background:#c8e6c9}
.status-info{color:#0d47a1;background:#bbdefb}
.btn{padding:9px 20px;border:none;border-radius:5px;cursor:pointer;
  font-size:.9em;font-weight:500;margin:4px}
.btn-primary{background:#3949ab;color:#fff}
.btn-primary:hover{background:#283593}
.btn-success{background:#2e7d32;color:#fff}
.btn-success:hover{background:#1b5e20}
.btn-secondary{background:#757575;color:#fff}
.btn-secondary:hover{background:#616161}
@media(max-width:1100px){.chart-grid-3{grid-template-columns:repeat(2,1fr)}}
@media(max-width:700px){.chart-grid-3,.chart-grid-2{grid-template-columns:1fr}}
"""


# ─── Shared JavaScript utilities ──────────────────────────────────────────────

SHARED_JS = """
const CHARTS = {};

const _AA_COL = {
  R:'#e74c3c', K:'#27ae60', H:'#9b59b6', D:'#8e44ad', E:'#f39c12',
  N:'#1abc9c', Q:'#e67e22', S:'#2ecc71', T:'#00bcd4', P:'#ff5722',
  F:'#795548', Y:'#607d8b', W:'#880e4f', M:'#4caf50', C:'#ffd700',
  A:'#90a4ae', V:'#78909c', L:'#546e7a', I:'#37474f', G:'#bdbdbd'
};
const _AA_FALLBACK = ['#3498db','#c0392b','#16a085','#d35400','#2c3e50'];
function aaColor(aa, i) { return _AA_COL[aa] || _AA_FALLBACK[i % _AA_FALLBACK.length]; }

function destroyChart(id) {
  if (CHARTS[id]) { CHARTS[id].destroy(); delete CHARTS[id]; }
}

function exportChart(id) {
  const c = document.getElementById(id);
  if (!c) return;
  const a = document.createElement('a');
  const docName = document.location.pathname.split('/').pop().replace(/\.html?$/i, '');
  a.download = (docName ? docName + '_' : '') + id + '.png';
  a.href = c.toDataURL('image/png', 1.0);
  a.click();
}

function _isTall(ctx) {
  const box = ctx.closest('.chart-box');
  return box && box.classList.contains('tall');
}

function barChart(id, labels, datasets, opts={}) {
  destroyChart(id);
  const ctx = document.getElementById(id);
  if (!ctx) return;
  CHARTS[id] = new Chart(ctx, {
    type: 'bar',
    data: { labels, datasets },
    options: {
      responsive: true,
      maintainAspectRatio: true,
      aspectRatio: opts.aspectRatio != null ? opts.aspectRatio : (_isTall(ctx) ? 0.7 : 1.0),
      plugins: { legend: { display: datasets.length > 1 }, tooltip: { mode: 'index' }, ...(opts.plugins||{}) },
      scales: {
        x: { ticks: { minRotation: 90, maxRotation: 90, font: { size: 10 } }, ...(opts.scales?.x||{}) },
        y: { beginAtZero: true, ...(opts.y||{}), ...(opts.scales?.y||{}) },
      }
    }
  });
}

function lineChart(id, datasets, opts={}) {
  destroyChart(id);
  const ctx = document.getElementById(id);
  if (!ctx) return;
  CHARTS[id] = new Chart(ctx, {
    type: 'line',
    data: { datasets },
    options: {
      responsive: true,
      maintainAspectRatio: true,
      aspectRatio: opts.aspectRatio != null ? opts.aspectRatio : (_isTall(ctx) ? 0.7 : 1.0),
      parsing: false,
      plugins: { legend: { display: true, position: 'top',
        labels: { boxWidth: 12, font: { size: 10 } } } },
      elements: { point: { radius: 0 }, line: { borderWidth: 1.5 } },
      scales: {
        x: { type: 'linear', ...(opts.x||{}) },
        y: { beginAtZero: true, ...(opts.y||{}) }
      }
    }
  });
}
"""
