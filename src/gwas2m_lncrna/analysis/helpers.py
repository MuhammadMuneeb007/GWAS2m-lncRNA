"""Small helpers shared by the Analysis01-09 files (phenotype labels and the standalone CLI)."""

from __future__ import annotations

import argparse
from functools import lru_cache
from pathlib import Path

from gwas2m_lncrna.common import canonical_ancestry, normalize_text, slugify


def normalize_label(text) -> str:
    return normalize_text(str(text or "").replace("’", "'"))


@lru_cache(maxsize=None)
def direct_labels(disease: str) -> frozenset:
    """Normalised trait labels that count as the disease itself (name + built-in aliases)."""
    from gwas2m_lncrna.steps.Step01_Plan_GWAS import build_terms
    return frozenset(normalize_label(t) for t in build_terms(disease, []))


def diseases_for(phenotypes: list[str]) -> list[tuple[str, str]]:
    """[(display name, folder slug), ...] - the DISEASES list used inside each analysis."""
    return [(p.strip(), slugify(p)) for p in phenotypes]


def ancestry_folder(ancestry: str) -> str:
    """'EUR' -> 'european', the ancestry folder name used in the output layout."""
    return slugify(canonical_ancestry(ancestry)[1])


def analysis_main(run, description: str, argv=None, extra=()):
    """Command line for running one analysis file on its own, e.g.

    python -m gwas2m_lncrna.analysis.Analysis05_John_lncRNA_Exon_Intron_Distribution \\
        -p "coronary artery disease" -p "Parkinson's disease" -a EUR
    """
    from gwas2m_lncrna.config import load_config

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-p", "--phenotype", action="append", required=True, help="repeat for several phenotypes")
    parser.add_argument("-a", "--ancestry", default="EUR")
    parser.add_argument("--config", help="gwas2m-lncrna.yaml (its output_dir is used as --root)")
    parser.add_argument("--root", help="output folder holding 01_gwas_catalog/, 06_finemapping/, ...")
    for flag, options in extra:
        parser.add_argument(flag, **options)
    args = parser.parse_args(argv)
    root = Path(args.root) if args.root else load_config(args.config).output_dir
    if not (root / "01_gwas_catalog").is_dir():
        parser.error(f"No GWAS2m output found at {root} (expected 01_gwas_catalog/). Use --root or --config.")
    options = {flag.lstrip("-").replace("-", "_"): getattr(args, flag.lstrip("-").replace("-", "_"))
               for flag, _ in extra}
    return run(root, diseases_for(args.phenotype), ancestry_folder(args.ancestry), **options)
