"""STEP 04 - make sure the ancestry-matched 1000 Genomes GRCh38 LD panel exists.

Modified from the original Step04_Prepare_LD_Reference.py: no SLURM array and no
separate 04_ld_reference/ copy. The panel is built once per ancestry under
resources/1000G/<ANC>/ (exactly where Steps 05 and 06 read it) and is shared by
every phenotype. Each selected study gets a Step04_LD_Reference status record.
"""

from __future__ import annotations

from gwas2m_lncrna.common import banner
from gwas2m_lncrna.project import Project
from gwas2m_lncrna.resources import Resources


def prepare_ld_reference(project: Project, studies_filter: list[str] | None = None) -> dict:
    banner(f"STEP 04 - LD REFERENCE | 1000 Genomes {project.ancestry_code}")
    resources = Resources(project.cfg)
    resources.ensure_ld(project.ancestry_code)
    n = resources.count_reference_samples(project.ancestry_code)
    ld_dir = resources.ld_dir / project.ancestry_code
    # Only the studies of this run: parallel per-GWAS jobs each record their own study.
    for study in project.studies(studies_filter):
        project.set_status(study["STUDY_ACCESSION"], "ld_reference", "COMPLETE",
                           LD_DIR=str(ld_dir), REFERENCE_SAMPLES=n)
    print(f"Panel ready: {ld_dir} ({n} samples)")
    return {"LD_DIR": str(ld_dir), "REFERENCE_SAMPLES": n}
