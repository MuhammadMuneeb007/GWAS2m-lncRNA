"""One phenotype x ancestry run, using the original GWAS2m folder layout.

    <output_dir>/
        01_gwas_catalog/<phenotype>/<ancestry>/   GWAS_download_manifest.tsv, study_selection.tsv, ...
        02_summary_stats/<phenotype>/<ancestry>/raw/<GCST>/   downloaded files
        02_summary_stats/<phenotype>/<ancestry>/qc/<GCST>/    GWASLab QC output
        05_ld_clumping/<phenotype>/<ancestry>/<GCST>/
        06_finemapping/<phenotype>/<ancestry>/<GCST>/
        07_annotation/<phenotype>/<ancestry>/<GCST>/
        pipeline_audit/<phenotype>/<ancestry>/status/<GCST>/<Stage>.json
        pipeline_audit/<phenotype>/<ancestry>/logs/<GCST>/<Stage>.log

Because the layout is unchanged, results produced earlier by GWAS2m can be read
(and resumed) by setting output_dir to that GWAS2m folder.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from gwas2m_lncrna.common import Gwas2RnaError, canonical_ancestry, read_json, slugify, utc_now, write_json
from gwas2m_lncrna.config import Config

# short step name used in the code -> stage name used in pipeline_audit status files
STAGE_NAMES = {
    "download": "Step02_Download",
    "qc": "Step03_QC",
    "ld_reference": "Step04_LD_Reference",
    "clump": "Step05_Clumping",
    "finemap": "Step06_FineMapping",
    "annotate": "Step07_VEP",
}
STUDY_STEPS = ("download", "qc", "ld_reference", "clump", "finemap", "annotate")


class Project:
    def __init__(self, cfg: Config, phenotype: str, ancestry: str):
        phenotype = (phenotype or "").strip()
        if not phenotype:
            raise Gwas2RnaError("--phenotype cannot be empty.")
        self.cfg = cfg
        self.phenotype = phenotype
        self.slug = slugify(phenotype)
        self.ancestry_code, self.ancestry_label = canonical_ancestry(ancestry)
        self.ancestry_slug = slugify(self.ancestry_label)  # e.g. "european", as in GWAS2m
        self.root = cfg.output_dir

    def __repr__(self) -> str:
        return f"Project({self.phenotype!r}, {self.ancestry_code})"

    def step_root(self, folder: str) -> Path:
        return self.root / folder / self.slug / self.ancestry_slug

    # --------------------------------------------------------------- folders
    @property
    def discovery_dir(self) -> Path:
        return self.step_root("01_gwas_catalog")

    @property
    def sumstats_dir(self) -> Path:
        return self.step_root("02_summary_stats") / "raw"

    @property
    def qc_dir(self) -> Path:
        return self.step_root("02_summary_stats") / "qc"

    @property
    def clump_dir(self) -> Path:
        return self.step_root("05_ld_clumping")

    @property
    def finemap_dir(self) -> Path:
        return self.step_root("06_finemapping")

    @property
    def annotation_dir(self) -> Path:
        return self.step_root("07_annotation")

    @property
    def audit_dir(self) -> Path:
        return self.step_root("pipeline_audit")

    def log_file(self, step: str, accession: str) -> Path:
        return self.audit_dir / "logs" / accession / f"{STAGE_NAMES[step]}.log"

    # -------------------------------------------------------------- manifest
    @property
    def manifest_path(self) -> Path:
        return self.discovery_dir / "GWAS_download_manifest.tsv"

    def manifest(self) -> pd.DataFrame:
        if not self.manifest_path.is_file():
            raise Gwas2RnaError(
                f"No study manifest for {self.phenotype} / {self.ancestry_code} ({self.manifest_path}). "
                f"Run step 1 first: gwas2m-lncrna discover --phenotype \"{self.phenotype}\" "
                f"--ancestry {self.ancestry_code}"
            )
        return pd.read_csv(self.manifest_path, sep="\t", dtype=str).fillna("")

    def studies(self, only: list[str] | None = None) -> list[dict]:
        rows = self.manifest().to_dict("records")
        if only:
            wanted = {x.strip().upper() for x in only}
            rows = [r for r in rows if r["STUDY_ACCESSION"].upper() in wanted]
        return rows

    def index_of(self, study: dict, position: int) -> int:
        """Study number shown by `status` (ARRAY_ID in the manifest, else row position)."""
        try:
            return int(float(study.get("ARRAY_ID", "") or position))
        except ValueError:
            return position

    def accessions_for(self, indices: list[int]) -> list[str]:
        """Map study numbers (1..N, as listed by `status`) to GCST accessions."""
        by_index = {self.index_of(s, i): s["STUDY_ACCESSION"] for i, s in enumerate(self.studies(), start=1)}
        missing = [i for i in indices if i not in by_index]
        if missing:
            raise Gwas2RnaError(
                f"No study with index {', '.join(map(str, missing))} for {self.phenotype} / {self.ancestry_code}; "
                f"valid indices are 1-{len(by_index)} (see `gwas2m-lncrna status`).")
        return [by_index[i] for i in indices]

    # ---------------------------------------------------------------- status
    def status_path(self, accession: str, step: str) -> Path:
        return self.audit_dir / "status" / accession / f"{STAGE_NAMES[step]}.json"

    def status(self, accession: str) -> dict:
        """{step: record} with each record's 'extra' fields merged in."""
        out = {}
        for step in STUDY_STEPS:
            record = read_json(self.status_path(accession, step))
            if record:
                out[step] = {**record.get("extra", {}), **record}
        return out

    def step_status(self, accession: str, step: str) -> str:
        return str(read_json(self.status_path(accession, step)).get("status", "NOT_RUN")).upper()

    def set_status(self, accession: str, step: str, status: str, reason: str = "", **extra) -> None:
        log = extra.pop("log", None)
        record = {
            "phenotype": self.phenotype, "phenotype_slug": self.slug,
            "ancestry_code": self.ancestry_code, "ancestry_label": self.ancestry_label,
            "study_accession": accession, "stage": STAGE_NAMES[step], "status": status,
            "reason": reason, "finished_utc": utc_now(), "log_path": log,
            "error_message": reason if status == "FAILED" else None, "extra": extra,
        }
        write_json(self.status_path(accession, step), record)

    def status_table(self) -> pd.DataFrame:
        rows = []
        for position, study in enumerate(self.studies(), start=1):
            accession = study["STUDY_ACCESSION"]
            record = self.status(accession)
            row = {"INDEX": self.index_of(study, position), "STUDY_ACCESSION": accession,
                   "TRAIT": study.get("DISEASE_TRAIT", "")}
            for step in STUDY_STEPS:
                row[STAGE_NAMES[step]] = str(record.get(step, {}).get("status", "NOT_RUN"))
            blockers = [f"{STAGE_NAMES[s]}: {record[s].get('reason', '')}" for s in STUDY_STEPS
                        if s in record and record[s].get("status") in {"FAILED", "EXCLUDED", "PARTIAL"}]
            row["NOTE"] = blockers[0] if blockers else ""
            rows.append(row)
        return pd.DataFrame(rows)
