"""Steps 1-7 for one phenotype x ancestry, in this process."""

from __future__ import annotations

from collections import Counter

from gwas2m_lncrna.common import DONE, Gwas2RnaError, banner
from gwas2m_lncrna.project import Project
from gwas2m_lncrna.resources import Resources
from gwas2m_lncrna.steps.runner import ready_for, run_for_studies
from gwas2m_lncrna.steps.Step01_Plan_GWAS import discover
from gwas2m_lncrna.steps.Step02_Download_GWAS import download_study
from gwas2m_lncrna.steps.Step03_QC_One_GWAS import qc_study
from gwas2m_lncrna.steps.Step04_Prepare_LD_Reference import prepare_ld_reference
from gwas2m_lncrna.steps.Step05_Clump_GWAS_Loci import clump_study
from gwas2m_lncrna.steps.Step06_FineMap_GWAS_Loci import finemap_study
from gwas2m_lncrna.steps.Step07_Annotate_Finemapped_Variants import annotate_study

STEP_NAMES = {1: "discover", 2: "download", 3: "qc", 4: "ld-reference", 5: "clump", 6: "finemap", 7: "annotate"}


def run_step(project: Project, number: int, force: bool = False, studies_filter: list[str] | None = None,
             discover_options: dict | None = None) -> Counter:
    cfg = project.cfg
    if number == 1:
        result = discover(project, force=force, **(discover_options or {}))
        return Counter({"COMPLETE": result["n_eligible"]})
    if number == 4:
        prepare_ld_reference(project)
        return Counter({"COMPLETE": 1})

    studies = project.studies(studies_filter)
    if number == 2:
        workers = max(1, int(cfg["resources"]["download_workers"]))
        return run_for_studies(project, "download", "STEP 2 - DOWNLOAD SUMMARY STATISTICS", download_study,
                               studies, force, jobs=workers, threads_not_processes=True)
    if number == 3:
        return run_for_studies(project, "qc", "STEP 3 - GWASLAB QC", qc_study,
                               ready_for(project, studies, "download"), force, cfg.jobs)
    if number == 5:
        Resources(cfg).ensure_ld(project.ancestry_code)
        return run_for_studies(project, "clump", "STEP 5 - LD CLUMPING", clump_study,
                               ready_for(project, studies, "qc"), force, cfg.jobs)
    if number == 6:
        Resources(cfg).ensure_ld(project.ancestry_code)
        return run_for_studies(project, "finemap", "STEP 6 - SUSIE-RSS FINE-MAPPING", finemap_study,
                               ready_for(project, studies, "clump"), force, cfg.jobs)
    if number == 7:
        candidates = ready_for(project, studies, "finemap", DONE | {"PARTIAL"})
        if candidates:
            Resources(cfg).ensure(None, vep=True)
        return run_for_studies(project, "annotate", "STEP 7 - VEP ANNOTATION", annotate_study,
                               candidates, force, cfg.jobs)
    raise Gwas2RnaError(f"Unknown step {number}")


def run_pipeline(project: Project, from_step: int = 1, to_step: int = 7, force: bool = False,
                 refresh_discovery: bool = False, studies_filter: list[str] | None = None,
                 discover_options: dict | None = None) -> None:
    if not 1 <= from_step <= to_step <= 7:
        raise Gwas2RnaError("Steps must satisfy 1 <= --from-step <= --to-step <= 7")
    banner(f"GWAS2M-LNCRNA PIPELINE | {project.phenotype} | {project.ancestry_code} | steps {from_step}-{to_step}")
    print(f"Configuration : {project.cfg.describe()}")
    print(f"Results       : {project.root}")
    print(f"Resources     : {project.cfg.resources_dir}")
    print(f"Threads/jobs  : {project.cfg.threads} / {project.cfg.jobs}")

    # Prepare shared resources first so a missing tool/resource fails early,
    # not hours into the run.
    resources = Resources(project.cfg)
    if from_step <= 6 and to_step >= 4:
        resources.ensure_ld(project.ancestry_code)
    if to_step >= 7:
        resources.ensure(None, vep=True)

    for number in range(from_step, to_step + 1):
        step_force = refresh_discovery if number == 1 else force
        run_step(project, number, step_force, studies_filter, discover_options)

    banner("PIPELINE FINISHED - STUDY STATUS")
    table = project.status_table()
    print(table.to_string(index=False) if not table.empty else "No studies.")
    print(f"\nNext: gwas2m-lncrna analysis --phenotype \"{project.phenotype}\" --ancestry {project.ancestry_code}")
