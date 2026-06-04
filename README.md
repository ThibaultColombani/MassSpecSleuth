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

## Report tabs

Tabs are shown only when the data supports them (e.g. Labeling is hidden for label-free runs, Plex Diagnostic requires a multiplexed experiment).

---

### Files

The entry point of the report. Lists all input runs with their original filenames and lets you:

- **Rename** runs — display names propagate to all charts and chart export filenames
- **Hide / show** runs — charts update live when you toggle a run off
- **Reorder** runs by drag-and-drop — controls the color order in all charts
- **Save** the modified report as a new HTML file (preserving your renames and ordering)

---

### Summary

High-level per-run counts at FDR ≤ 1%: unique sequence-charges, peptides, and protein groups. When `--enzyme` is passed, the peptide and seqcharge bars are broken down by C-terminal amino acid (R vs K for trypsin, etc.).

<table>
<tr>
  <td align="center"><b>Unique seqcharges</b><br><img src="docs/screenshots/summary_seqcharges.png" width="190" alt="Seqcharges per run"></td>
  <td align="center"><b>Unique peptides</b><br><img src="docs/screenshots/summary_peptide_counts.png" width="190" alt="Peptide counts per run"></td>
  <td align="center"><b>Protein groups</b><br><img src="docs/screenshots/summary_protein_groups.png" width="190" alt="Protein groups per run"></td>
</tr>
</table>

---

### IDs

Deeper identification metrics, all filtered at the engine's primary quality score (PEP or Q-value ≤ 0.01).

**Precursor IDs** — precursors ranked by PEP/Q-value (ascending); cumulative count vs. score. Lets you compare score distributions across runs.

<img src="docs/screenshots/ids_cumulative_precursors.png" width="600" alt="Cumulative precursor IDs">

**Peptide identifications & charge distribution** — unique seqcharges and peptides per run side by side, plus the charge state distribution (fraction of 2+, 3+, 4+ precursors).

<table>
<tr>
  <td align="center"><b>Seqcharges &amp; Precursors per run</b><br><img src="docs/screenshots/ids_peptide_counts.png" width="290" alt="Seqcharges and precursors per run"></td>
  <td align="center"><b>Charge state distribution (PEP &lt; 0.01)</b><br><img src="docs/screenshots/ids_charge_distribution.png" width="290" alt="Charge state distribution"></td>
</tr>
</table>

**Missed cleavages** — fraction of peptides with 0, 1, 2+ internal cleavage sites (inferred from sequence, engine-agnostic). Broken down by C-terminal amino acid when an enzyme is set.

<table>
<tr>
  <td align="center"><b>All peptides</b><br><img src="docs/screenshots/ids_missed_cleavages.png" width="190" alt="Missed cleavages — all"></td>
  <td align="center"><b>R-ending peptides</b><br><img src="docs/screenshots/ids_missed_cleavages_R.png" width="190" alt="Missed cleavages — R-ending"></td>
  <td align="center"><b>K-ending peptides</b><br><img src="docs/screenshots/ids_missed_cleavages_K.png" width="190" alt="Missed cleavages — K-ending"></td>
</tr>
</table>

**Protein identifications** — unique protein groups per run at four FDR stringencies (all PSMs, Q < 0.05, Q < 0.01, Q < 0.001), making it easy to spot runs with inflated or depleted ID counts.

<img src="docs/screenshots/ids_protein_ids.png" width="600" alt="Protein IDs">

---

### Chromatography

Per-run chromatography quality, shown in individual charts (one panel per run). A "samples per row" dropdown controls the grid layout.

**Retention time distribution** — precursor count per 0.5-min RT bin, coloured by C-terminal amino acid. A flat, broad distribution indicates even sampling across the gradient; a spike at the void suggests gradient problems.

<img src="docs/screenshots/chrom_rt.png" width="600" alt="RT distribution">

**Peak width — FWHM** — distribution of chromatographic peak widths in seconds. Narrow, consistent peaks indicate good LC performance. Runs with broader or bimodal distributions flag gradient or column issues.

<img src="docs/screenshots/chrom_fwhm.png" width="600" alt="Peak width FWHM">

---

### Ion Sampling

The most information-dense tab. Covers MS1 and MS2 intensity distributions, cross-run consistency of shared precursors, m/z-dependent intensity profiles, and the MS2/MS1 ratio.

**MS1 intensity distribution** — log₁₀ MS1 area per run (full dynamic range, 50 bins). Overlapping peaks across runs indicate consistent ionisation efficiency.

<img src="docs/screenshots/ion_ms1_distribution.png" width="600" alt="MS1 intensity distribution">

**MS1 — intersected precursors** — log₁₀ MS1 intensity for precursors detected in every active run. Coloured by C-terminal residue when an enzyme filter is set. Updates live as runs are toggled in the Files tab.

<img src="docs/screenshots/ion_ms1_intersected.png" width="600" alt="MS1 intersected">

**Normalized MS1 — intersected precursors** — log₂ fold-change of MS1 intensity vs. the first active run, for shared precursors. Axis fixed at [−3, +3]. A tight peak centred at 0 means stable ionisation across runs.

<img src="docs/screenshots/ion_ms1_normalized.png" width="600" alt="Normalized MS1">

**MS2 intensity distribution** — log₁₀ MS2 precursor quantity (Precursor.Quantity) per run.

<img src="docs/screenshots/ion_ms2_distribution.png" width="600" alt="MS2 intensity distribution">

**MS2 intensity vs precursor m/z** — median log₁₀ MS2 intensity in 25-Da precursor m/z bins. Y-axis = m/z, X-axis = intensity. Systematic steps across the m/z range may reveal DIA window boundaries or instrument bias.

<img src="docs/screenshots/ion_ms2_mz_profile.png" width="600" alt="MS2 vs m/z profile">

**MS2/MS1 ratio — intersected precursors** — log₂(Precursor.Quantity / MS1 Intensity) for shared precursors. A narrow, stable distribution indicates consistent fragmentation efficiency across the run.

<img src="docs/screenshots/ion_ms2_ms1_ratio.png" width="600" alt="MS2/MS1 ratio">

---

### Labeling

Shown only for experiments using amine-reactive tags (PSMtag, TMT, mTRAQ, iTRAQ, …). Tag detection is automatic — no configuration needed.

**Tag labeling efficiency** — percentage of N-terminally and Lys-labeled peptides, broken down per C-terminal residue (R-ending vs K-ending). Stacked bars show fully labeled, partly labeled (N-term only, one or more K missed), and unlabeled fractions. Individual or grouped view available.

<img src="docs/screenshots/labeling_efficiency.png" width="600" alt="Labeling efficiency">

**Precursor counts** — total precursors per run and per residue, with fully- and partly-labeled breakdown.

<img src="docs/screenshots/labeling_counts.png" width="600" alt="Labeling counts">

**Retention length — fully labeled peptides** — peak width FWHM (seconds) for fully labeled precursors, mean ± SD per run. Right panel: precursors shared across all active runs only.

<table>
<tr>
  <td align="center"><b>All Fully Labeled Peptides</b><br><img src="docs/screenshots/labeling_fwhm.png" width="290" alt="FWHM — all fully labeled"></td>
  <td align="center"><b>Intersected Fully Labeled Peptides</b><br><img src="docs/screenshots/labeling_fwhm_intersected.png" width="290" alt="FWHM — intersected"></td>
</tr>
</table>

**MS1 intensity — fully labeled peptides** — log₁₀ MS1 intensity distribution for fully labeled precursors, mean ± SD per run. Right panel: precursors shared across all active runs only.

<table>
<tr>
  <td align="center"><b>All Fully Labeled Peptides</b><br><img src="docs/screenshots/labeling_ms1_intensity.png" width="290" alt="MS1 intensity — all fully labeled"></td>
  <td align="center"><b>Intersected Fully Labeled Peptides</b><br><img src="docs/screenshots/labeling_ms1_intensity_intersected.png" width="290" alt="MS1 intensity — intersected"></td>
</tr>
</table>

**MS1 ratio — intersected fully labeled peptides** — log₂(sample / control) MS1 intensity ratio for labeled precursors shared across runs. Useful for spotting systematic loading differences.

<img src="docs/screenshots/labeling_ratio.png" width="600" alt="MS1 ratio">

---

### Plex Diagnostic

Shown only for multiplexed experiments (PSMtag, TMT, mTRAQ, …). Covers per-channel identification counts, intensity distributions, quantification variability (CV), and data completeness (Jaccard overlap between channels).

**Identification counts** — unique precursors and protein groups per channel per run, with an "in all channels" intersection count.

<table>
<tr>
  <td align="center"><b>Precursor IDs per Channel</b><br><img src="docs/screenshots/plex_ids.png" width="290" alt="Precursor IDs per channel"></td>
  <td align="center"><b>Protein Groups per Channel</b><br><img src="docs/screenshots/plex_proteins.png" width="290" alt="Protein groups per channel"></td>
</tr>
</table>

**Channel Q-value distribution** — distribution of Channel Q-values per channel (DIA-NN / Jmod only). A tight distribution near 0 indicates confident channel-level assignments.

<img src="docs/screenshots/plex_qval.png" width="600" alt="Channel Q-value">

**Relative channel intensity** — log₁₀ intensity of each sample channel relative to the highest-intensity (carrier) channel per run. Flat bars across channels indicate balanced loading.

<img src="docs/screenshots/plex_relative_intensity.png" width="600" alt="Relative channel intensity">

**MS2 intensity distribution** — log₁₀ MS2 (Precursor.Quantity) distribution per channel per run.

<img src="docs/screenshots/plex_ms2_distribution.png" width="600" alt="MS2 distribution per channel">

**Data completeness — precursors (Jaccard)** — pairwise Jaccard overlap between sample channels (carrier excluded). High overlap means channels quantify the same set of precursors.

<img src="docs/screenshots/plex_jaccard_precursors.png" width="600" alt="Jaccard precursor overlap">

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
