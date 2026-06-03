# Created by @thibaultcolombani with the help of Claude Code (Anthropic)
"""
MassSpecSleuth.py - Main entry point for MassSleuth

Detects the search engine from each given folder, loads all relevant files,
and routes to the appropriate processor and HTML exporter.

Usage
-----
    # Single folder (auto-detect engine)
    python MassSpecSleuth.py /path/to/results

    # Multiple folders — same engine (merged into one report)
    python MassSpecSleuth.py /path/run1 /path/run2 /path/run3

    # Multiple folders — mixed engines (combined report)
    python MassSpecSleuth.py /path/diann_run /path/jmod_run

    # With enzyme filter
    python MassSpecSleuth.py /path/to/results --enzyme trypsin

    # Custom output path
    python MassSpecSleuth.py /path/to/results --output /tmp/report.html
"""

import argparse
import sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path

import polars as pl

from Mods.data_loader import load_folder
from Mods.diann_processor import DiannProcessor
from Mods.diann_html_exporter import export_html_report as diann_export_html
from Mods.jmod_processor import JmodProcessor
from Mods.jmod_html_exporter import export_html_report as jmod_export_html
from Mods.maxquant_processor import MaxquantProcessor
from Mods.maxquant_html_exporter import export_html_report as maxquant_export_html
from Mods.sage_processor import SageProcessor
from Mods.sage_html_exporter import export_html_report as sage_export_html
from Mods.fragpipe_processor import FragpipeProcessor
from Mods.fragpipe_html_exporter import export_html_report as fragpipe_export_html
from Mods.protdisc_processor import ProtdiscProcessor
from Mods.protdisc_html_exporter import export_html_report as protdisc_export_html
from Mods.enzyme import ENZYMES


def _merge_data_dicts(datas: list) -> dict:
    """Merge multiple data dicts from load_folder() by concatenating per key."""
    merged = {}
    all_keys = {k for d in datas for k in d}
    for key in all_keys:
        frames = [d[key] for d in datas if key in d]
        merged[key] = pl.concat(frames, how='diagonal') if len(frames) > 1 else frames[0]
    return merged


def _make_processor(engine: str, data: dict):
    """Instantiate and normalize the right processor for the given engine."""
    if engine == 'diann':
        p = DiannProcessor(data)
        p.normalize()
        return p
    if engine == 'jmod':
        p = JmodProcessor(data)
        p.normalize()
        return p
    if engine == 'maxquant':
        p = MaxquantProcessor(data)
        p.normalize()
        return p
    if engine == 'sage':
        p = SageProcessor(data)
        p.normalize()
        return p
    if engine == 'fragpipe':
        p = FragpipeProcessor(data)
        p.normalize()
        return p
    if engine == 'protdisc':
        p = ProtdiscProcessor(data)
        p.normalize()
        return p
    return None


def main():
    print("MassSpecSleuth — Mass Spectrometry Analysis Tool")
    print("=" * 50)

    parser = argparse.ArgumentParser(
        description="MassSpecSleuth — mass spectrometry QC report generator.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        'folders', nargs='*',
        help="Path(s) to search-engine output folder(s). Defaults to example/DIANN.",
    )
    parser.add_argument(
        '--enzyme', nargs='+', metavar='ENZYME',
        choices=sorted(ENZYMES),
        help=(
            "Keep only fully specific peptides for this enzyme. "
            "Multiple enzymes are combined (union of cleavage rules). "
            f"Available: {', '.join(sorted(ENZYMES))}"
        ),
    )
    parser.add_argument(
        '--output', metavar='PATH',
        help="Output HTML file path. Default: first folder / massspecsleuth_report.html.",
    )

    args = parser.parse_args()

    folders = args.folders if args.folders else [
        str(Path(__file__).parent / "example" / "DIANN")
    ]

    print(f"Input: {', '.join(folders)}")
    if args.enzyme:
        print(f"Enzyme filter: {' + '.join(args.enzyme)}")

    # ── Load each folder ──────────────────────────────────────────────────────
    engine_groups: dict = defaultdict(list)   # engine → [data_dict, …]

    for folder in folders:
        engine, data = load_folder(folder)
        if engine == 'unknown' or not data:
            print(f"  Warning: no supported data found in {folder!r}, skipping.")
            continue
        engine_groups[engine].append(data)

    if not engine_groups:
        print("No supported search engine output found. Exiting.")
        sys.exit(1)

    # ── Determine output path ─────────────────────────────────────────────────
    if args.output:
        output_path = args.output
    else:
        ts = datetime.now().strftime('%Y%m%d_%H%M%S')
        suffix = 'combined' if len(engine_groups) > 1 else next(iter(engine_groups))
        output_path = str(Path(folders[0]) / f"mss_report_{suffix}_{ts}.html")

    # ── Single engine ─────────────────────────────────────────────────────────
    if len(engine_groups) == 1:
        engine, datas = next(iter(engine_groups.items()))
        data = _merge_data_dicts(datas)

        print(f"\nEngine: {engine.upper()}")
        if len(datas) > 1:
            print(f"  Merged {len(datas)} folder(s)")

        processor = _make_processor(engine, data)

        if processor is None:
            print(f"  No HTML exporter for '{engine}'. Loaded files:")
            for k, df in data.items():
                print(f"  {k}: {df.shape[0]:,} rows × {df.shape[1]} columns")
            return

        if args.enzyme:
            processor.filter_enzyme(args.enzyme)

        if engine == 'diann':
            diann_export_html(processor, output_path)
        elif engine == 'maxquant':
            maxquant_export_html(processor, output_path)
        elif engine == 'sage':
            sage_export_html(processor, output_path)
        elif engine == 'fragpipe':
            fragpipe_export_html(processor, output_path)
        elif engine == 'protdisc':
            protdisc_export_html(processor, output_path)
        else:
            jmod_export_html(processor, output_path)

    # ── Multiple engines → combined report ───────────────────────────────────
    else:
        from Mods.combined_processor import CombinedProcessor

        engines_str = ' + '.join(e.upper() for e in engine_groups)
        print(f"\nEngines: {engines_str}")

        processors = []
        for engine, datas in engine_groups.items():
            data = _merge_data_dicts(datas)
            if len(datas) > 1:
                print(f"  {engine.upper()}: merged {len(datas)} folder(s)")
            p = _make_processor(engine, data)
            if p is None:
                print(f"  No processor for '{engine}', skipping.")
                continue
            if args.enzyme:
                p.filter_enzyme(args.enzyme)
            processors.append(p)

        if not processors:
            print("No processors could be created. Exiting.")
            sys.exit(1)

        combined = CombinedProcessor(processors)
        jmod_export_html(combined, output_path, title=f"Combined Report ({engines_str})")


if __name__ == "__main__":
    main()
