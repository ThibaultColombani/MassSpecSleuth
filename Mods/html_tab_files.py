# Created by @thibaultcolombani with the help of Claude Code (Anthropic)
"""html_tab_files.py — Files tab: rename / show / hide runs."""
from __future__ import annotations
from typing import Dict

TAB_ID    = 'files'
TAB_LABEL = 'Files'


def is_visible(processor) -> bool:
    return True


def compute(processor) -> Dict:
    return {'runs': processor.runs()}


def html_content(data: Dict) -> str:
    return """
<div class="section">
  <h3>File Management</h3>
  <p class="help">Rename, hide, or reorder runs. Drag the ⠿ handle to reorder. Click "Apply" to refresh all charts.<br>
  <strong>Labeling Grouped mode:</strong> runs are grouped by the part of their display name before the first <code>-</code>.
  Example: rename runs to <code>Ctrl-1</code>, <code>Ctrl-2</code>, <code>Treat-1</code>, <code>Treat-2</code>
  → two groups <em>Ctrl</em> and <em>Treat</em>. Apply changes, then switch to Grouped in the Labeling tab.</p>
  <div id="file-list"></div>
  <div style="margin-top:14px;text-align:center">
    <button class="btn btn-success" onclick="applyRenames()">Apply Changes</button>
    <button class="btn btn-secondary" onclick="resetRenames()">Reset</button>
    <button class="btn btn-primary" onclick="saveHTML()">&#8681; Save Report</button>
  </div>
  <div id="file-status"></div>
</div>"""


def javascript() -> str:
    return """
let _dragSrc = null;

function buildFileList() {
  const el = document.getElementById('file-list');
  el.innerHTML = runOrder.map(r => `
    <div class="file-row" draggable="true" data-run="${r}"
         ondragstart="onFileDragStart(event)" ondragover="onFileDragOver(event)"
         ondragleave="onFileDragLeave(event)" ondrop="onFileDrop(event)"
         ondragend="onFileDragEnd(event)">
      <span class="drag-handle" title="Drag to reorder">⠿</span>
      <span class="file-orig">${r}</span>
      <input class="file-rename" id="inp-${r}" value="${label(r)}" placeholder="Rename…">
      <button class="file-toggle ${activeRuns.includes(r)?'active-yes':'active-no'}"
        id="tog-${r}" onclick="toggleRun('${r}')">
        ${activeRuns.includes(r)?'Active':'Hidden'}</button>
    </div>`).join('');
}

function toggleRun(r) {
  const idx = activeRuns.indexOf(r);
  if (idx >= 0) activeRuns.splice(idx, 1); else activeRuns.push(r);
  const btn = document.getElementById('tog-' + r);
  btn.textContent = activeRuns.includes(r) ? 'Active' : 'Hidden';
  btn.className = 'file-toggle ' + (activeRuns.includes(r) ? 'active-yes' : 'active-no');
}

function onFileDragStart(e) {
  _dragSrc = e.currentTarget;
  _dragSrc.classList.add('dragging');
  e.dataTransfer.effectAllowed = 'move';
  e.dataTransfer.setData('text/plain', '');
}

function onFileDragOver(e) {
  e.preventDefault();
  e.dataTransfer.dropEffect = 'move';
  const target = e.currentTarget;
  if (target === _dragSrc) return;
  const rect = target.getBoundingClientRect();
  target.classList.remove('drag-over-top', 'drag-over-bottom');
  target.classList.add(e.clientY < rect.top + rect.height / 2 ? 'drag-over-top' : 'drag-over-bottom');
}

function onFileDragLeave(e) {
  e.currentTarget.classList.remove('drag-over-top', 'drag-over-bottom');
}

function onFileDrop(e) {
  e.preventDefault();
  const target = e.currentTarget;
  if (!_dragSrc || target === _dragSrc) return;
  const rect = target.getBoundingClientRect();
  target.classList.remove('drag-over-top', 'drag-over-bottom');
  if (e.clientY < rect.top + rect.height / 2) {
    target.parentNode.insertBefore(_dragSrc, target);
  } else {
    target.parentNode.insertBefore(_dragSrc, target.nextSibling);
  }
}

function onFileDragEnd(e) {
  document.querySelectorAll('.file-row').forEach(row => {
    row.classList.remove('dragging', 'drag-over-top', 'drag-over-bottom');
  });
  _dragSrc = null;
}

function applyRenames() {
  const rows = document.querySelectorAll('#file-list .file-row');
  runOrder = Array.from(rows).map(row => row.dataset.run);
  runOrder.forEach(r => { renamedRuns[r] = document.getElementById('inp-' + r).value || r; });
  Object.keys(_inited).forEach(k => delete _inited[k]);
  const active = document.querySelector('.tab-content.active');
  if (active) initTab(active.id);
  showFileStatus('Changes applied!', 'status-success');
}

function resetRenames() {
  runOrder = [...RUNS];
  activeRuns = [...RUNS];
  RUNS.forEach(r => { renamedRuns[r] = r; });
  buildFileList();
  Object.keys(_inited).forEach(k => delete _inited[k]);
  const active = document.querySelector('.tab-content.active');
  if (active) initTab(active.id);
  showFileStatus('Report reset to original state.', 'status-info');
}

async function saveHTML() {
  let html = '<!DOCTYPE html>\\n' + document.documentElement.outerHTML;
  html = html.replace(/let runOrder += *[^;]+;/,
    'let runOrder   = ' + JSON.stringify(runOrder) + ';');
  html = html.replace(/let activeRuns += *[^;]+;/,
    'let activeRuns = ' + JSON.stringify(activeRuns) + ';');
  html = html.replace(/let renamedRuns += *[^;]+;/,
    'let renamedRuns = ' + JSON.stringify(renamedRuns) + ';');

  const blob = new Blob([html], { type: 'text/html' });
  const docFile = document.location.pathname.split('/').pop() || 'mss_report.html';

  if (window.showSaveFilePicker) {
    try {
      const handle = await window.showSaveFilePicker({
        suggestedName: docFile,
        types: [{ description: 'HTML file', accept: { 'text/html': ['.html'] } }]
      });
      const writable = await handle.createWritable();
      await writable.write(blob);
      await writable.close();
      showFileStatus('Report saved.', 'status-success');
    } catch (err) {
      if (err.name !== 'AbortError') console.error(err);
    }
  } else {
    const name = prompt('Save as:', docFile);
    if (!name) return;
    const a = document.createElement('a');
    a.href = URL.createObjectURL(blob);
    a.download = name.endsWith('.html') ? name : name + '.html';
    a.click();
    URL.revokeObjectURL(a.href);
    showFileStatus('Report saved.', 'status-success');
  }
}

function showFileStatus(msg, cls) {
  const el = document.getElementById('file-status');
  if (!el) return;
  el.textContent = msg;
  el.className = cls;
  clearTimeout(el._t);
  el._t = setTimeout(() => { el.textContent = ''; el.className = ''; }, 3000);
}
"""
