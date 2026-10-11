"""Step 2 - download each selected study's harmonised summary statistics."""

from __future__ import annotations

from pathlib import Path

from gwas2m_lncrna.download import download
from gwas2m_lncrna.project import Project


def target_path(project: Project, study: dict) -> Path:
    return project.sumstats_dir / study["STUDY_ACCESSION"] / study["FILE_NAME"]


def download_study(project: Project, study: dict, force: bool) -> tuple:
    target = target_path(project, study)
    if force:
        target.unlink(missing_ok=True)
    status = download(study["DOWNLOAD_URL"], target, label=study["STUDY_ACCESSION"])
    return status, "", {"FILE": str(target), "BYTES": target.stat().st_size}
