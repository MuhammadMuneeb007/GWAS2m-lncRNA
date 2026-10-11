"""Command-line interface: `gwas2m-lncrna <command> ...`."""

from __future__ import annotations

import argparse
import sys
import traceback

from gwas2m_lncrna import __version__
from gwas2m_lncrna.common import Gwas2RnaError, banner
from gwas2m_lncrna.config import CONFIG_FILENAME, load_config, write_template

EPILOG = """\
typical use:
  gwas2m-lncrna init                                   # optional: write gwas2m-lncrna.yaml to edit
  gwas2m-lncrna setup --ancestry EUR                   # optional: pre-download reference data
  gwas2m-lncrna run -p "coronary artery disease" -a EUR
  gwas2m-lncrna analysis -p "coronary artery disease" -a EUR

individual steps (same arguments as run):
  discover (1) download (2) qc (3) ld-reference (4) clump (5) finemap (6) annotate (7)
"""


def _common(parser: argparse.ArgumentParser, phenotype: str = "one") -> None:
    parser.add_argument("--config", help=f"configuration file (default: ./{CONFIG_FILENAME} if present)")
    parser.add_argument("--threads", type=int, help="threads per tool call (overrides config)")
    parser.add_argument("--jobs", type=int, help="studies processed in parallel (overrides config)")
    parser.add_argument("--output-dir", help="results folder (overrides config)")
    if phenotype == "one":
        parser.add_argument("-p", "--phenotype", required=True, help='e.g. "coronary artery disease"')
    elif phenotype == "many":
        parser.add_argument("-p", "--phenotype", required=True, action="append",
                            help="repeat for several phenotypes")
    parser.add_argument("-a", "--ancestry", required=phenotype != "none",
                        help="EUR, AFR, EAS, SAS or AMR")


def _step_options(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--study", action="append", help="restrict to GCST accession(s); repeatable")
    parser.add_argument("--index", type=int, action="append",
                        help="restrict to the study with this number from `status` (1..N); repeatable. "
                             "Use with $SLURM_ARRAY_TASK_ID to run one GWAS per array task")
    parser.add_argument("--force", action="store_true", help="redo work even if outputs exist")


def _discover_options(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--alias", action="append", default=[], help="extra trait search term; repeatable")
    parser.add_argument("--max-studies", type=int, help="keep at most N eligible studies (0 = all)")
    parser.add_argument("--refresh-catalog", action="store_true", help="re-download GWAS Catalog metadata")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="gwas2m-lncrna", formatter_class=argparse.RawDescriptionHelpFormatter, epilog=EPILOG,
        description="GWAS2m-lncRNA: transcript-aware lncRNA overlap in fine-mapped GWAS signals.\n"
                    "GWAS discovery -> QC -> clumping -> SuSiE-RSS fine-mapping -> VEP -> lncRNA report.")
    parser.add_argument("--version", action="version", version=f"gwas2m-lncrna {__version__}")
    sub = parser.add_subparsers(dest="command", required=True, metavar="command")

    p = sub.add_parser("init", help="write a configuration file to edit")
    p.add_argument("--path", default=CONFIG_FILENAME)
    p.add_argument("--force", action="store_true")

    p = sub.add_parser("check", help="show tools, resources and configuration")
    _common(p, phenotype="none")
    p.add_argument("--ready", action="store_true",
                   help="with --ancestry: exit 0 if the LD panel, VEP cache and genome are ready, else 3 (for scripts)")

    p = sub.add_parser("setup", help="download/build reference data (1000G LD, VEP cache, genome, GENCODE)")
    _common(p, phenotype="none")
    p.add_argument("--only", choices=["ld", "vep", "gencode"], action="append",
                   help="limit to these resources (default: all)")

    p = sub.add_parser("run", help="run steps 1-7 end to end")
    _common(p, phenotype="many")
    _step_options(p)
    _discover_options(p)
    p.add_argument("--from-step", type=int, default=1)
    p.add_argument("--to-step", type=int, default=7)
    p.add_argument("--refresh-discovery", action="store_true", help="search the GWAS Catalog again (step 1)")
    p.add_argument("--analysis", action="store_true", help="also write the report when the run finishes")

    for name, number, text in (("discover", 1, "step 1: find GWAS Catalog studies"),
                               ("download", 2, "step 2: download summary statistics"),
                               ("qc", 3, "step 3: GWASLab QC"),
                               ("ld-reference", 4, "step 4: prepare the 1000G LD panel"),
                               ("clump", 5, "step 5: LD clumping"),
                               ("finemap", 6, "step 6: SuSiE-RSS fine-mapping"),
                               ("annotate", 7, "step 7: VEP annotation")):
        p = sub.add_parser(name, help=text)
        p.set_defaults(step_number=number)
        _common(p)
        _step_options(p)
        if number == 1:
            _discover_options(p)

    p = sub.add_parser("status", help="studies of a phenotype with their index and per-step status")
    _common(p)
    p.add_argument("--count", action="store_true", help="print only the number of studies (for array sizes)")

    p = sub.add_parser("analysis", help="write report.html + report_tables.xlsx + variants.tsv.gz")
    _common(p, phenotype="many")
    p.add_argument("--out", help="report folder (default: <output_dir>/reports/<phenotypes>_<ANC>)")
    return parser


def _config(args):
    overrides = {"threads": getattr(args, "threads", None), "jobs": getattr(args, "jobs", None),
                 "output_dir": getattr(args, "output_dir", None)}
    return load_config(getattr(args, "config", None), overrides)


def _selected_studies(project, args) -> list[str] | None:
    """GCST accessions chosen with --study and/or --index (None = all studies)."""
    chosen = list(args.study or [])
    if getattr(args, "index", None):
        chosen += project.accessions_for(args.index)
    return chosen or None


def _discover_kwargs(args) -> dict:
    return {"aliases": args.alias or [], "max_studies": args.max_studies,
            "refresh_catalog": True if args.refresh_catalog else None}


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    try:
        return dispatch(args)
    except Gwas2RnaError as exc:
        print(f"\nERROR: {exc}", file=sys.stderr)
        return 2
    except KeyboardInterrupt:
        print("\nInterrupted.", file=sys.stderr)
        return 130


def dispatch(args) -> int:
    if args.command == "init":
        path = write_template(args.path, args.force)
        print(f"Wrote {path}\nEdit it, then run commands from this folder (or pass --config {path}).")
        return 0

    cfg = _config(args)
    from gwas2m_lncrna.project import Project

    if args.command == "check" and args.ready:
        from gwas2m_lncrna.common import canonical_ancestry
        from gwas2m_lncrna.resources import Resources
        if not args.ancestry:
            raise Gwas2RnaError("check --ready needs --ancestry (e.g. --ancestry EUR).")
        resources = Resources(cfg)
        code = canonical_ancestry(args.ancestry)[0]
        missing = [name for name, ok in ((f"1000G {code} LD panel", resources.ld_ready(code)),
                                         ("VEP cache", resources.vep_ready()),
                                         ("GRCh38 FASTA", resources.fasta_ready())) if not ok]
        print("READY" if not missing else "MISSING: " + ", ".join(missing))
        return 0 if not missing else 3

    if args.command == "check":
        import pandas as pd
        from gwas2m_lncrna.resources import Resources
        from gwas2m_lncrna.tools import tool_report
        banner("GWAS2M-LNCRNA CHECK")
        print(f"Configuration : {cfg.describe()}\nResults       : {cfg.output_dir}\n"
              f"Resources     : {cfg.resources_dir}\nThreads/jobs  : {cfg.threads} / {cfg.jobs}\n")
        print(pd.DataFrame(tool_report(cfg)).to_string(index=False))
        print()
        print(pd.DataFrame(Resources(cfg).status_rows()).to_string(index=False))
        return 0

    if args.command == "setup":
        from gwas2m_lncrna.common import canonical_ancestry
        from gwas2m_lncrna.resources import Resources
        only = set(args.only or ["ld", "vep", "gencode"])
        resources = Resources(cfg)
        if "ld" in only:
            if not args.ancestry:
                raise Gwas2RnaError("setup needs --ancestry for the LD panel (or use --only vep / --only gencode).")
            resources.ensure_ld(canonical_ancestry(args.ancestry)[0])
        if "vep" in only:
            resources.ensure(None, vep=True)
        if "gencode" in only:
            resources.ensure_gtf()
        banner("SETUP COMPLETE")
        for row in resources.status_rows():
            print(f"{'ready  ' if row['ready'] else 'missing'}  {row['resource']:<28} {row['location']}")
        return 0

    if args.command == "run":
        from gwas2m_lncrna.pipeline import run_pipeline
        if args.index and len(args.phenotype) > 1:
            raise Gwas2RnaError("--index refers to one phenotype's study list; give a single --phenotype.")
        failed = []
        for phenotype in args.phenotype:
            try:
                project = Project(cfg, phenotype, args.ancestry)
                selected = _selected_studies(project, args)
                run_pipeline(project, args.from_step, args.to_step, args.force,
                             args.refresh_discovery, selected, _discover_kwargs(args))
            except Gwas2RnaError as exc:
                if len(args.phenotype) == 1:
                    raise
                print(f"\nERROR ({phenotype}): {exc}\nContinuing with the next phenotype.", file=sys.stderr)
                failed.append(phenotype)
        finished = [p for p in args.phenotype if p not in failed]
        if args.analysis and finished:
            from gwas2m_lncrna.analysis.report import run_analysis
            run_analysis(cfg, finished, args.ancestry)
        if failed:
            print(f"\nPhenotypes that stopped with an error: {', '.join(failed)}", file=sys.stderr)
            return 2
        return 0

    if args.command == "status":
        project = Project(cfg, args.phenotype, args.ancestry)
        if args.count:
            print(len(project.studies()))
            return 0
        banner(f"STATUS | {project.phenotype} | {project.ancestry_code}")
        table = project.status_table()
        print(table.to_string(index=False) if not table.empty else "No studies yet.")
        return 0

    if args.command == "analysis":
        from gwas2m_lncrna.analysis.report import run_analysis
        run_analysis(cfg, args.phenotype, args.ancestry, args.out)
        return 0

    if hasattr(args, "step_number"):
        from gwas2m_lncrna.pipeline import run_step
        project = Project(cfg, args.phenotype, args.ancestry)
        options = _discover_kwargs(args) if args.step_number == 1 else None
        selected = _selected_studies(project, args) if args.step_number > 1 else None
        run_step(project, args.step_number, args.force, selected, options)
        return 0

    raise Gwas2RnaError(f"Unknown command {args.command}")


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception:
        traceback.print_exc()
        raise SystemExit(1)
