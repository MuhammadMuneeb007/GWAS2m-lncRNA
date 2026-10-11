"""Locate external programs (PLINK2, Rscript, VEP, samtools)."""

from __future__ import annotations

import os
import re
import shutil
import subprocess
from pathlib import Path

from gwas2m_lncrna.common import Gwas2RnaError
from gwas2m_lncrna.config import Config

PURPOSE = {
    "plink2": "LD panels, clumping and LD matrices (steps 4-6)",
    "rscript": "SuSiE-RSS fine-mapping (step 6)",
    "vep": "variant annotation (step 7)",
    "samtools": "indexing the GRCh38 FASTA for VEP (setup)",
}
EXECUTABLE = {"plink2": "plink2", "rscript": "Rscript", "vep": "vep", "samtools": "samtools"}


def find_tool(cfg: Config, name: str, required: bool = True) -> Path | None:
    configured = cfg["tools"].get(name)
    if configured:
        path = cfg.path(configured)
        if path and path.is_file() and os.access(path, os.X_OK):
            return path
        if required:
            raise Gwas2RnaError(f"tools.{name} = {configured} is not an executable file.")
        return None
    located = shutil.which(EXECUTABLE[name])
    if located:
        return Path(located).resolve()
    if required:
        raise Gwas2RnaError(
            f"'{EXECUTABLE[name]}' was not found on PATH (needed for {PURPOSE[name]}). "
            "Activate the environment (`mamba activate gwas2m-lncrna`) or set tools."
            f"{name} in gwas2m-lncrna.yaml."
        )
    return None


def require_susier(rscript: Path) -> None:
    check = subprocess.run(
        [str(rscript), "-e", 'quit(status=ifelse(requireNamespace("susieR", quietly=TRUE),0,1))'],
        capture_output=True, text=True,
    )
    if check.returncode != 0:
        raise Gwas2RnaError("R package susieR is not installed for " + str(rscript))


def vep_release(vep: Path) -> str | None:
    try:
        out = subprocess.run([str(vep), "--help"], capture_output=True, text=True, timeout=120)
    except (OSError, subprocess.TimeoutExpired):
        return None
    match = re.search(r"ensembl-vep\s*:\s*(\d+)", out.stdout + out.stderr)
    return match.group(1) if match else None


def tool_report(cfg: Config) -> list[dict]:
    rows = []
    for name in EXECUTABLE:
        path = find_tool(cfg, name, required=False)
        rows.append({"tool": EXECUTABLE[name], "path": str(path) if path else "NOT FOUND",
                     "needed_for": PURPOSE[name]})
    rscript = find_tool(cfg, "rscript", required=False)
    if rscript:
        try:
            require_susier(rscript)
            rows.append({"tool": "R:susieR", "path": "installed", "needed_for": PURPOSE["rscript"]})
        except Gwas2RnaError:
            rows.append({"tool": "R:susieR", "path": "NOT FOUND", "needed_for": PURPOSE["rscript"]})
    try:
        import gwaslab  # noqa: F401
        rows.append({"tool": "python:gwaslab", "path": "installed", "needed_for": "summary-statistics QC (step 3)"})
    except ImportError:
        rows.append({"tool": "python:gwaslab", "path": "NOT FOUND", "needed_for": "summary-statistics QC (step 3)"})
    return rows
