"""`gwas2m-lncrna analysis`: run Analysis01-09 and save the results in three files.

    report.html          every analysis: its tables and its full printed output
    report_tables.xlsx   every table from every analysis (one sheet each)
    variants.tsv.gz      Analysis03 per-variant classification (credible-set and high-PIP variants)

Each analysis is its own file in this folder (Analysis01_John_..., ...). This module
only calls their run() functions one after another.
"""

from __future__ import annotations

import contextlib
import html
import io
import numbers
import sys
import traceback
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from gwas2m_lncrna import __version__
from gwas2m_lncrna.analysis import (
    Analysis01_John_GWAS_Selection, Analysis02_John_Finemapping_CredibleSets,
    Analysis03_John_Coding_Noncoding_Distribution, Analysis04_John_Genomic_Consequence_Distribution,
    Analysis05_John_lncRNA_Exon_Intron_Distribution, Analysis06_John_lncRNA_Locus_Representation,
    Analysis07_John_lncRNA_Enrichment, Analysis08_John_Cross_Disease_Comparison,
    Analysis09_John_Finemapping_Sensitivity,
)
from gwas2m_lncrna.analysis.helpers import ancestry_folder, diseases_for
from gwas2m_lncrna.common import Gwas2RnaError, banner
from gwas2m_lncrna.config import Config
from gwas2m_lncrna.resources import Resources

# (code, title, module) in the order they run.
ANALYSES = [
    ("A01", "GWAS selection, QC and pipeline status", Analysis01_John_GWAS_Selection),
    ("A02", "Fine-mapped loci and 95% credible sets", Analysis02_John_Finemapping_CredibleSets),
    ("A03", "Strict protein-coding CDS vs non-CDS", Analysis03_John_Coding_Noncoding_Distribution),
    ("A04", "Genomic consequence distribution", Analysis04_John_Genomic_Consequence_Distribution),
    ("A05", "lncRNA exon / intron / splice distribution", Analysis05_John_lncRNA_Exon_Intron_Distribution),
    ("A06", "lncRNA locus representation", Analysis06_John_lncRNA_Locus_Representation),
    ("A07", "Within-locus lncRNA enrichment", Analysis07_John_lncRNA_Enrichment),
    ("A08", "Cross-phenotype comparison", Analysis08_John_Cross_Disease_Comparison),
    ("A09", "Fine-mapping sensitivity and PIP weighting", Analysis09_John_Finemapping_Sensitivity),
]
VARIANT_TABLE = ("A03", "variant_classes")  # saved separately as variants.tsv.gz
HTML_MAX_ROWS = 50


class _Tee:
    def __init__(self, *streams):
        self.streams = streams

    def write(self, text):
        for stream in self.streams:
            stream.write(text)
        return len(text)

    def flush(self):
        for stream in self.streams:
            stream.flush()


def options_for(code: str, cfg: Config) -> dict:
    """Settings from gwas2m-lncrna.yaml passed to each analysis's run()."""
    a = cfg["analysis"]
    top = int(a["top"])
    legacy = bool(a["legacy_lncrna_biotypes"])
    return {
        "A01": {},
        "A02": {"top": top},
        "A03": {"high_pip": 0.5},
        "A04": {"high_pip": 0.5},
        "A05": {"high_pip": 0.5, "top": top, "include_legacy_lncrna_biotypes": legacy},
        "A06": {"top": top, "high_pip": 0.5, "include_legacy_lncrna_biotypes": legacy},
        "A07": {"permutations": int(a["permutations"]), "seed": int(a["seed"]), "bin_kb": int(a["bin_kb"]),
                "legacy_lncrna": legacy, "min_annotation_coverage": float(a["min_annotation_coverage"])},
        "A08": {"high_pip": 0.5, "top": top},
        "A09": {},
    }[code]


def run_analyses(root: Path, phenotypes: list[str], ancestry: str, cfg: Config, gtf) -> list[dict]:
    diseases = diseases_for(phenotypes)
    folder = ancestry_folder(ancestry)
    results = []
    for code, title, module in ANALYSES:
        options = options_for(code, cfg)
        if code == "A06":
            options["gtf"] = gtf
        buffer = io.StringIO()
        error = None
        with contextlib.redirect_stdout(_Tee(sys.__stdout__, buffer)):
            try:
                tables = module.run(root, diseases, folder, **options) or {}
            except Exception as exc:  # one analysis failing must not lose the others
                traceback.print_exc(file=buffer)
                error, tables = f"{type(exc).__name__}: {exc}", {}
        results.append({"code": code, "title": title, "module": module.__name__.rsplit(".", 1)[-1],
                        "tables": tables, "output": buffer.getvalue(), "error": error})
        if error:
            print(f"[{code}] FAILED: {error}", file=sys.stderr)
    return results


# =============================================================================
# FILES
# =============================================================================

CSS = """
:root{--bg:#f7f7f5;--panel:#fff;--text:#1d1f23;--muted:#5d636e;--line:#e2e2de;--accent:#2f5d8a;--head:#eef1f5;
--bad:#a3261f}
@media (prefers-color-scheme:dark){:root{--bg:#15171a;--panel:#1d2024;--text:#e6e7ea;--muted:#9aa0aa;
--line:#30343a;--accent:#7fb0e0;--head:#25292e;--bad:#f08a80}}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--text);
font:15px/1.55 system-ui,-apple-system,"Segoe UI",Roboto,sans-serif}
main{max-width:1180px;margin:0 auto;padding:32px 16px 64px}h1{font-size:26px;margin:0 0 4px}
h2{font-size:20px;margin:44px 0 6px;padding-top:10px;border-top:1px solid var(--line)}
h3{font-size:15px;margin:22px 0 6px}.meta,.desc{color:var(--muted);font-size:13px}
nav{margin:16px 0;display:flex;flex-wrap:wrap;gap:6px 14px;font-size:14px}a{color:var(--accent)}
.scroll{overflow-x:auto;background:var(--panel);border:1px solid var(--line);border-radius:8px}
table{border-collapse:collapse;font-size:12.5px;width:100%}th,td{padding:5px 8px;border-bottom:1px solid var(--line);
text-align:left;white-space:nowrap}th{background:var(--head)}td.num{text-align:right;font-variant-numeric:tabular-nums}
details{margin-top:14px}summary{cursor:pointer;color:var(--accent);font-size:14px}
pre{background:var(--panel);border:1px solid var(--line);border-radius:8px;padding:12px;overflow-x:auto;
font-size:12px;line-height:1.4;max-height:70vh}.error{color:var(--bad);font-weight:600}
"""


def _is_number(value) -> bool:
    return isinstance(value, numbers.Number) and not isinstance(value, (bool, np.bool_))


def _fmt(value):
    if value is None or (isinstance(value, numbers.Real) and pd.isna(value)):
        return "–"
    if isinstance(value, (bool, np.bool_)):
        return str(bool(value))
    if isinstance(value, numbers.Integral):
        return f"{int(value):,}"
    if isinstance(value, numbers.Real):
        value = float(value)
        return f"{value:.4g}" if abs(value) < 1e-3 and value != 0 else f"{value:,.3f}".rstrip("0").rstrip(".")
    return str(value)


def _html_table(frame: pd.DataFrame) -> str:
    head = "".join(f"<th>{html.escape(str(c))}</th>" for c in frame.columns)
    body = "".join(
        "<tr>" + "".join(f"<td{' class=num' if _is_number(v) else ''}>{html.escape(_fmt(v))}</td>" for v in row)
        + "</tr>" for row in frame.itertuples(index=False))
    return f"<div class=scroll><table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table></div>"


def sheet_name(code: str, table: str) -> str:
    return f"{code}_{table}"[:31]


def write_html(path: Path, title: str, results: list[dict], cfg: Config, files: dict) -> None:
    nav, sections = [], []
    for r in results:
        anchor = r["code"].lower()
        nav.append(f'<a href="#{anchor}">{r["code"]} {html.escape(r["title"])}</a>')
        parts = [f'<h2 id="{anchor}">{r["code"]} · {html.escape(r["title"])}</h2>',
                 f'<p class=desc>{html.escape(r["module"])}.py</p>']
        if r["error"]:
            parts.append(f'<p class=error>This analysis stopped with an error: {html.escape(r["error"])}</p>')
        for name, frame in r["tables"].items():
            if (r["code"], name) == VARIANT_TABLE:
                parts.append(f"<h3>{html.escape(name)}</h3><p class=desc>{len(frame):,} rows, saved in "
                             f"{html.escape(files['variants'].name)}.</p>")
                continue
            note = (f"Showing {HTML_MAX_ROWS} of {len(frame):,} rows; " if len(frame) > HTML_MAX_ROWS else "")
            parts.append(f"<h3>{html.escape(name)}</h3><p class=desc>{note}workbook sheet "
                         f"<b>{html.escape(sheet_name(r['code'], name))}</b></p>")
            parts.append(_html_table(frame.head(HTML_MAX_ROWS)) if not frame.empty
                         else "<p class=desc>No rows.</p>")
        parts.append(f"<details><summary>Full printed output of {html.escape(r['module'])}</summary>"
                     f"<pre>{html.escape(r['output'])}</pre></details>")
        sections.append("".join(parts))
    generated = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    path.write_text(f"""<!doctype html><html lang=en><head><meta charset=utf-8>
<meta name=viewport content="width=device-width,initial-scale=1"><title>{html.escape(title)}</title>
<style>{CSS}</style></head><body><main>
<h1>{html.escape(title)}</h1>
<p class=meta>GWAS2m-lncRNA {__version__} · generated {generated} · configuration: {html.escape(cfg.describe())}<br>
All tables: {html.escape(files['xlsx'].name)} · per-variant classes: {html.escape(files['variants'].name)}</p>
<nav>{''.join(nav)}</nav>
{''.join(sections)}
</main></body></html>""", encoding="utf-8")


def write_workbook(path: Path, results: list[dict]) -> None:
    about, sheets = [], {}
    for r in results:
        for name, frame in r["tables"].items():
            if (r["code"], name) == VARIANT_TABLE:
                continue
            sheets[sheet_name(r["code"], name)] = frame
            about.append({"SHEET": sheet_name(r["code"], name), "ANALYSIS": f"{r['code']} {r['title']}",
                          "SOURCE_FILE": r["module"] + ".py", "TABLE": name, "ROWS": len(frame)})
    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        pd.DataFrame(about).to_excel(writer, sheet_name="about", index=False)
        for name, frame in sheets.items():
            frame.to_excel(writer, sheet_name=name, index=False)
        for sheet in writer.book.worksheets:
            sheet.freeze_panes = "A2"
            for column in sheet.columns:
                width = max(len(str(c.value)) if c.value is not None else 0 for c in list(column)[:200])
                sheet.column_dimensions[column[0].column_letter].width = min(60, max(8, width + 2))


def run_analysis(cfg: Config, phenotypes: list[str], ancestry: str, out_dir: str | None = None) -> dict:
    from gwas2m_lncrna.project import Project

    projects = [Project(cfg, p, ancestry) for p in phenotypes]
    for project in projects:
        if not project.manifest_path.is_file():
            raise Gwas2RnaError(f"No results for {project.phenotype} / {project.ancestry_code}: "
                                f"{project.manifest_path} is missing. Run the pipeline first.")
    gtf = None
    try:
        resources = Resources(cfg)
        resources.ensure_gtf()
        gtf = resources.gtf
    except Gwas2RnaError as exc:
        print(f"[note] GENCODE gene annotation unavailable; Analysis06 gene-body classes are skipped: {exc}")

    banner("GWAS2M-LNCRNA ANALYSIS 01-09")
    results = run_analyses(cfg.output_dir, [p.phenotype for p in projects], ancestry, cfg, gtf)

    names = "__".join(p.slug for p in projects)
    out = Path(out_dir).resolve() if out_dir else cfg.output_dir / "reports" / f"{names}_{projects[0].ancestry_code}"
    out.mkdir(parents=True, exist_ok=True)
    files = {"html": out / "report.html", "xlsx": out / "report_tables.xlsx", "variants": out / "variants.tsv.gz"}
    variants = next((r["tables"].get(VARIANT_TABLE[1]) for r in results if r["code"] == VARIANT_TABLE[0]), None)
    (variants if variants is not None else pd.DataFrame()).to_csv(
        files["variants"], sep="\t", index=False, compression="gzip")
    write_workbook(files["xlsx"], results)
    write_html(files["html"], "GWAS2m-lncRNA report: " + ", ".join(p.phenotype for p in projects), results, cfg, files)

    banner("REPORT WRITTEN")
    for label, path in files.items():
        print(f"{label:<9}: {path}")
    failed = [r["code"] for r in results if r["error"]]
    if failed:
        print(f"Analyses that stopped with an error (see report.html): {', '.join(failed)}", file=sys.stderr)
    return files
