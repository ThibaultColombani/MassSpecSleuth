# Created by @thibaultcolombani with the help of Claude Code (Anthropic)
"""
sage_html_exporter.py - Orchestrator for the Sage interactive HTML QC report.

Same tab architecture as diann_html_exporter, jmod_html_exporter, and
maxquant_html_exporter.  Tabs that require DIA-NN/Jmod-specific data
(fill_times, features, plex_report, Precursor.Quantity) are hidden
automatically by their own is_visible() guards.

Tabs shown for a typical Sage (DDA) run:
  Files, Summary, Identifications, Chromatography
  Labeling (only when PSMtag data detected)
  Ion Sampling / Features / Plex Diagnostic — hidden (no DIA data)
"""

from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path

from Mods.sage_processor import SageProcessor
from Mods.html_utils import SHARED_CSS, SHARED_JS

# ── Register tabs in display order ───────────────────────────────────────────
import Mods.html_tab_files           as _tab_files
import Mods.html_tab_summary         as _tab_summary
import Mods.html_tab_identifications as _tab_ident
import Mods.html_tab_chromatography  as _tab_chrom
import Mods.html_tab_ion_sampling    as _tab_ion
import Mods.html_tab_labeling        as _tab_lab
import Mods.html_tab_plex_diagnostic  as _tab_plex
import Mods.html_tab_features        as _tab_feat

_ALL_TABS = [
    _tab_files,
    _tab_summary,
    _tab_ident,
    _tab_chrom,
    _tab_ion,
    _tab_lab,
    _tab_plex,
    _tab_feat,
]


# ── Assembly ──────────────────────────────────────────────────────────────────

def export_html_report(processor, output_path: str,
                       title: str = "Sage QC Report") -> None:
    """
    Build and write the Sage QC HTML report.

    Parameters
    ----------
    processor : SageProcessor | CombinedProcessor
        Must already have .normalize() called (or be a CombinedProcessor).
    output_path : str
        Destination .html file path.
    title : str
        Report title shown in the browser tab and page header.
    """
    runs      = processor.runs()
    colors    = {r: ['#3498db','#e74c3c','#27ae60','#f39c12','#9b59b6',
                     '#1abc9c','#e67e22','#2980b9','#c0392b','#16a085'][i % 10]
                 for i, r in enumerate(runs)}
    generated   = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    filter_str  = ' · '.join(processor.filters) if processor.filters else 'none'

    # ── Determine which tabs are visible ─────────────────────────────────────
    visible_tabs = [t for t in _ALL_TABS if t.is_visible(processor)]

    # ── Compute data for each visible tab ────────────────────────────────────
    print("  [html] computing tab data…")
    tab_data = {}
    for tab in visible_tabs:
        print(f"    {tab.TAB_ID}")
        tab_data[tab.TAB_ID] = tab.compute(processor)

    # ── Master JSON payload ──────────────────────────────────────────────────
    payload = {
        'runs':   runs,
        'colors': colors,
        **{t.TAB_ID: tab_data[t.TAB_ID] for t in visible_tabs},
    }

    # ── Build HTML ───────────────────────────────────────────────────────────
    print(f"  [html] assembling HTML ({len(visible_tabs)} tabs)…")

    tab_buttons = '
'.join(
        f'  <button class="tab{" active" if i == 0 else ""}" '
        f'onclick="showTab(\'{t.TAB_ID}\')">{t.TAB_LABEL}</button>'
        for i, t in enumerate(visible_tabs)
    )

    tab_contents = '
'.join(
        f'<div id="{t.TAB_ID}" class="tab-content{" active" if i == 0 else ""}">
'
        f'{t.html_content(tab_data[t.TAB_ID])}
</div>'
        for i, t in enumerate(visible_tabs)
    )

    tab_js = '
'.join(t.javascript() for t in visible_tabs)

    init_cases = '
'.join(
        f'  if (name === "{t.TAB_ID}") init{_camel(t.TAB_ID)}();'
        for t in visible_tabs
        if t.TAB_ID != 'files'
    )

    data_json = json.dumps(payload, allow_nan=False)

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width,initial-scale=1.0">
<title>MSS — {title}</title>
<script src="https://cdn.jsdelivr.net/npm/chart.js@4.4.0/dist/chart.umd.min.js"></script>
<style>
{SHARED_CSS}
</style>
</head>
<body>
<div class="container">

<div class="header">
  <h1>MSS — {title}</h1>
  <p>Generated {generated}&nbsp;|&nbsp;{len(runs)} run(s)&nbsp;|&nbsp;Filters: {filter_str}</p>
</div>

<div class="tabs">
{tab_buttons}
</div>

{tab_contents}

</div><!-- /container -->

<script src="https://cdn.jsdelivr.net/npm/chart.js@4.4.0/dist/chart.umd.min.js"></script>
<script>
const DATA = {data_json};
const RUNS = {json.dumps(runs)};
let runOrder   = [...RUNS];
let activeRuns = [...RUNS];
let renamedRuns = Object.fromEntries(RUNS.map(r => [r, r]));

function label(r)  {{ return renamedRuns[r] || r; }}
function color(r)  {{ return DATA.colors[r] || '#3498db'; }}
function getActiveRuns() {{ return runOrder.filter(r => activeRuns.includes(r)); }}

{SHARED_JS}

// ── Tab switching ─────────────────────────────────────────
function showTab(name) {{
  document.querySelectorAll('.tab-content').forEach(el => el.classList.remove('active'));
  document.querySelectorAll('.tab').forEach(el => el.classList.remove('active'));
  const content = document.getElementById(name);
  if (content) content.classList.add('active');
  if (event && event.target) event.target.classList.add('active');
  initTab(name);
}}

const _inited = {{}};
function initTab(name) {{
  if (_inited[name]) return;
  _inited[name] = true;
{init_cases}
}}

// ── Tab init functions ────────────────────────────────────
{tab_js}

// ── Boot ──────────────────────────────────────────────────
document.addEventListener('DOMContentLoaded', () => {{
  buildFileList();
  const firstContent = document.querySelector('.tab-content.active');
  if (firstContent) initTab(firstContent.id);
}});
</script>
</body>
</html>"""

    Path(output_path).write_text(html, encoding='utf-8')
    size_kb = Path(output_path).stat().st_size // 1024
    print(f"  [html] saved → {output_path}  ({size_kb} KB)")


def _camel(tab_id: str) -> str:
    return ''.join(w.capitalize() for w in tab_id.split('_'))
