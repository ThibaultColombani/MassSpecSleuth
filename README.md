# MassSpecSleuth

MassSpecSleuth (MSS) is a Python tool that generates interactive HTML quality-control reports from mass spectrometry search engine output. It auto-detects the search engine, loads all relevant result files, and produces a self-contained HTML report you can open in any browser.

---

## Supported search engines

| Engine | Key output file(s) |
|---|---|
| **DIA-NN** | `report.parquet` |
| **Jmod** | `*_IDs.parquet` / `*_IDs.csv` |
| **MaxQuant** | `evidence.txt` |
| **Sage** | `*.sage.tsv` (+ optional `lfq.tsv`) |
| **FragPipe** | `psm.tsv` |
| **ProteomeDiscoverer** | `*_PSMs.txt` |

---

## Requirements

- **Python 3.10+**
- **[Polars](https://pola.rs/)** — the only third-party dependency

```bash
pip install -r requirements.txt
```

The generated HTML report loads **Chart.js** from the jsDelivr CDN at open time, so an internet connection is required to view it in the browser.

---

## Usage

```bash
# Single folder — engine is auto-detected
python MassSpecSleuth.py /path/to/results

# Multiple folders, same engine — merged into one report
python MassSpecSleuth.py /path/run1 /path/run2 /path/run3

# Multiple folders, mixed engines — combined multi-engine report
python MassSpecSleuth.py /path/diann_run /path/jmod_run

# Apply an enzyme specificity filter
python MassSpecSleuth.py /path/to/results --enzyme trypsin

# Custom output path
python MassSpecSleuth.py /path/to/results --output /tmp/my_report.html
```

The report is written as a single self-contained `.html` file. By default it is saved inside the first input folder and named `mss_report_<engine>_<YYYYMMDD_HHMMSS>.html` (or `mss_report_combined_…` for mixed-engine runs).

---

## Enzyme filter

Pass `--enzyme <name>` to keep only fully specific peptides for that enzyme. Multiple enzymes are combined (union of cleavage rules).

```bash
python MassSpecSleuth.py /path/to/results --enzyme trypsin lys-c
```

Run `python MassSpecSleuth.py --help` to see all 23 supported enzyme names.

---

## Project structure

```
MassSpecSleuth/
├── MassSpecSleuth.py       ← CLI entry point
├── requirements.txt
└── Mods/
    ├── data_loader.py      ← engine detection and file loading
    ├── enzyme.py           ← enzyme catalogue and specificity filter
    ├── *_processor.py      ← per-engine data normalisation
    ├── *_html_exporter.py  ← per-engine report assembly
    ├── html_utils.py       ← shared CSS/JS, chart helpers
    └── html_tab_*.py       ← one module per report tab
```

---

## Report tabs

Tabs are shown only when the data supports them (e.g. Labeling is hidden for label-free runs, Plex Diagnostic requires a multiplexed experiment).

---

### Files

Rename, hide, or reorder runs. Drag the ⠿ handle to reorder. Click "Apply" to refresh all charts. The modified report can be saved as a new self-contained HTML file.

---

### Summary

Per-run unique counts at FDR ≤ 1%. Apply `--enzyme` to see C-terminal breakdown.

<table>
<tr>
  <td align="center"><b>Unique seqcharges</b><br><img src="docs/screenshots/summary_seqcharges.png" width="190" alt="Seqcharges per run"></td>
  <td align="center"><b>Unique peptides</b><br><img src="docs/screenshots/summary_peptide_counts.png" width="190" alt="Peptide counts per run"></td>
  <td align="center"><b>Protein groups</b><br><img src="docs/screenshots/summary_protein_groups.png" width="190" alt="Protein groups per run"></td>
</tr>
</table>

---

### IDs

**Precursor IDs** — Precursors ranked by PEP (ascending). X = cumulative precursor count, Y = PEP (log₁₀). PEP used as the quality metric (PEP when available, Q.Value otherwise).

<img src="docs/screenshots/ids_cumulative_precursors.png" width="800" alt="Cumulative precursor IDs">

**Peptide Identifications & Charge Distribution** — Seqcharges, precursors, and charge state distribution at PEP < 0.01.

<table>
<tr>
  <td align="center"><b>Seqcharges &amp; Precursors per run</b><br><img src="docs/screenshots/ids_peptide_counts.png" width="290" alt="Seqcharges and precursors per run"></td>
  <td align="center"><b>Charge state distribution (PEP &lt; 0.01)</b><br><img src="docs/screenshots/ids_charge_distribution.png" width="290" alt="Charge state distribution"></td>
</tr>
</table>

**Missed Cleavages** — Internal cleavage sites within peptides (PEP < 0.01, Intensity > 0). Cleavage rule: inferred (K, R).

<table>
<tr>
  <td align="center"><b>All peptides</b><br><img src="docs/screenshots/ids_missed_cleavages.png" width="190" alt="Missed cleavages — all"></td>
  <td align="center"><b>R-ending peptides</b><br><img src="docs/screenshots/ids_missed_cleavages_R.png" width="190" alt="Missed cleavages — R-ending"></td>
  <td align="center"><b>K-ending peptides</b><br><img src="docs/screenshots/ids_missed_cleavages_K.png" width="190" alt="Missed cleavages — K-ending"></td>
</tr>
</table>

**Protein Identifications** — Protein IDs per run (4 stringencies).

<img src="docs/screenshots/ids_protein_ids.png" width="600" alt="Protein IDs">

---

### Chromatography

Per-run chromatography quality, one panel per run. A "samples per row" dropdown controls the grid layout.

**Retention Time Distribution** — Precursors per 0.5-min RT bin, coloured by C-terminal amino acid. Y = RT (min), X = count.

<img src="docs/screenshots/chrom_rt.png" width="600" alt="RT distribution">

**Peak Width — FWHM** — Peak width distribution per run. Y = FWHM (s), X = precursor count.

<img src="docs/screenshots/chrom_fwhm.png" width="600" alt="Peak width FWHM">

---

### Ion Sampling

**MS1 Intensity Distribution** — Log₁₀ MS1 precursor intensity per run (full range, 50 bins). Peaks should overlap across runs for consistent ionisation efficiency.

<img src="docs/screenshots/ion_ms1_distribution.png" width="600" alt="MS1 intensity distribution">

**MS1 Intensity – Intersected Precursors** — Log₁₀ MS1 intensity for precursors detected in every active run, coloured by C-terminal residue when an enzyme is selected. Updates when runs are hidden/shown in the Files tab. Requires ≥2 active runs with ≥20 shared precursors.

<img src="docs/screenshots/ion_ms1_intersected.png" width="600" alt="MS1 intersected">

**Normalized MS1 Intensity – Intersected Precursors** — Log₂ fold-change of MS1 intensity vs. the first active run, for intersected precursors. Axis fixed [−3, +3]. A narrow peak at 0 indicates stable ionisation. Updates live with run selection.

<img src="docs/screenshots/ion_ms1_normalized.png" width="600" alt="Normalized MS1">

**MS2 Intensity Distribution** — Log₁₀ MS2 precursor quantity per run (full range, 50 bins), coloured by C-terminal residue when an enzyme is selected.

<img src="docs/screenshots/ion_ms2_distribution.png" width="600" alt="MS2 intensity distribution">

**MS2 Intensity vs Precursor M/Z** — Median log₁₀ MS2 intensity (Precursor.Quantity) in 25-Da precursor m/z bins per run. Y-axis = precursor m/z; X-axis = log₁₀ intensity. Systematic shifts across the m/z range may reflect DIA window boundaries or instrument bias.

<img src="docs/screenshots/ion_ms2_mz_profile.png" width="600" alt="MS2 vs m/z profile">

**MS2/MS1 Intensity Ratio – Intersected Precursors** — Distribution of log₂(Precursor.Quantity / MS1 Intensity) for precursors detected in every active run with both MS1 and MS2 quantification. A narrow, stable distribution indicates consistent fragmentation efficiency. Updates live with run selection.

<img src="docs/screenshots/ion_ms2_ms1_ratio.png" width="600" alt="MS2/MS1 ratio">

---

### Labeling

Shown only for experiments using amine-reactive tags (PSMtag, TMT, mTRAQ, iTRAQ, …). Tag detection is automatic — no configuration needed.

**Tag Labeling Efficiency** — Breakdown by C-terminal residue (R, K). *Fully labeled*: all expected amine sites (N-terminus + all K side chains) carry a tag. *Partly labeled*: at least one site labeled, but not all. Bars for each residue are stacked (Full + Partly). Apply `--enzyme` to restrict to specific residues. Without `--enzyme`, residues with ≥1% of precursors are shown.

<img src="docs/screenshots/labeling_efficiency.png" width="600" alt="Labeling efficiency">

**Precursor Counts** — Total precursors and per-residue breakdown: total, fully labeled, partly labeled.

<img src="docs/screenshots/labeling_counts.png" width="600" alt="Labeling counts">

**Retention Length — Fully Labeled Peptides** — Peak width FWHM (seconds) for fully labeled precursors. Mean ± SD per run. Right panel: precursors present in all active runs only.

<table>
<tr>
  <td align="center"><b>All Fully Labeled Peptides</b><br><img src="docs/screenshots/labeling_fwhm.png" width="290" alt="FWHM — all fully labeled"></td>
  <td align="center"><b>Intersected Fully Labeled Peptides</b><br><img src="docs/screenshots/labeling_fwhm_intersected.png" width="290" alt="FWHM — intersected"></td>
</tr>
</table>

**MS1 Intensity — Fully Labeled Peptides** — log₁₀(MS1 intensity) for fully labeled precursors. Mean ± SD per run. Right panel: precursors present in all active runs only.

<table>
<tr>
  <td align="center"><b>All Fully Labeled Peptides</b><br><img src="docs/screenshots/labeling_ms1_intensity.png" width="290" alt="MS1 intensity — all fully labeled"></td>
  <td align="center"><b>Intersected Fully Labeled Peptides</b><br><img src="docs/screenshots/labeling_ms1_intensity_intersected.png" width="290" alt="MS1 intensity — intersected"></td>
</tr>
</table>

**MS1 Ratio — Intersected Fully Labeled Peptides** — log₂(sample / control) of MS1 intensity for fully labeled precursors present in all active runs. Requires ≥2 active runs.

<img src="docs/screenshots/labeling_ratio.png" width="600" alt="MS1 ratio">

---

### Plex Diagnostic

Shown only for multiplexed experiments (PSMtag, TMT, mTRAQ, …).

**Identification Counts** — Per-sample precursor and protein group counts by channel (filtered: PEP ≤ 0.01 and Channel Q-Value ≤ 0.01). *In all channels* = identified in every channel of that run (intersection). *In any channel* = total unique across all channels (union). Protein groups use proteotypic peptides only when available.

<table>
<tr>
  <td align="center"><b>Precursor IDs per Channel</b><br><img src="docs/screenshots/plex_ids.png" width="290" alt="Precursor IDs per channel"></td>
  <td align="center"><b>Protein Groups per Channel</b><br><img src="docs/screenshots/plex_proteins.png" width="290" alt="Protein groups per channel"></td>
</tr>
</table>

**Channel Q-Value Distribution** — Channel Q-Value distribution per channel for each run (values ≤ 0.01 after quality filter). Channels concentrated near 0 indicate high-confidence identifications.

<img src="docs/screenshots/plex_qval.png" width="600" alt="Channel Q-value">

**Relative Channel Intensity** — log₁₀ intensity ratio of each sample channel vs the highest-intensity (carrier) channel, computed on intersected precursors. Box = IQR (p25–p75), centre bar = median, whiskers = p5/p95.

<img src="docs/screenshots/plex_relative_intensity.png" width="600" alt="Relative channel intensity">

**MS2 Intensity Distribution** — log₁₀ MS2 (Precursor.Quantity) distribution per channel per run.

<img src="docs/screenshots/plex_ms2_distribution.png" width="600" alt="MS2 distribution per channel">

**Data Completeness — Precursors** — Pairwise precursor overlap between sample channels (carrier excluded). Jaccard Index = |A∩B| / |A∪B|; higher = more shared precursors between channels.

<img src="docs/screenshots/plex_jaccard_precursors.png" width="600" alt="Jaccard precursor overlap">
