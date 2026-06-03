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

| Tab | What it shows |
|---|---|
| **Summary** | Peptide, protein, and sequence-charge counts per run |
| **Identifications** | PSM-level ID counts filtered by Q-value or PEP; missed cleavages; protein groups |
| **Chromatography** | Retention time distributions and peak widths per run |
| **Features** | MS1 feature-level metrics (DIA-NN only) |
| **Ion Sampling** | MS1/MS2 intensity distributions, m/z profiles, cross-run precursor consistency |
| **Labeling** | Labeling efficiency (N-term + Lys) for TMT, PSMtag, or similar amine-reactive tags |
| **Plex Diagnostic** | Per-channel intensity distributions, CVs, Jaccard overlap, and channel Q-values |
| **Files** | Input file metadata; save / export controls |

Tabs are shown only when the data supports them (e.g. Labeling is hidden for label-free runs, Plex Diagnostic requires a multiplexed experiment).

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
