"""Step 3 - GWASLab standardisation/QC of one study (GRCh38).

Outputs per study (03_qc/<GCST>/):
    <GCST>_GRCh38_QC.tsv.gz           full cleaned GWAS
    <GCST>_GRCh38_significant.tsv.gz  P < p_threshold
    <GCST>_QC_summary.tsv / .json     counts, available fields, readiness
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from gwas2m_lncrna.common import banner, find_column, nonempty, utc_now, write_json
from gwas2m_lncrna.project import Project

VALID_EXTENSIONS = (".tsv.gz", ".txt.gz", ".vcf.gz", ".tsv", ".txt", ".vcf")


def find_downloaded_file(project: Project, study: dict) -> Path | None:
    directory = project.sumstats_dir / study["STUDY_ACCESSION"]
    expected = directory / study.get("FILE_NAME", "")
    if study.get("FILE_NAME") and nonempty(expected):
        return expected.resolve()
    if not directory.is_dir():
        return None
    candidates = [p for p in directory.iterdir()
                  if p.is_file() and p.stat().st_size > 0 and not p.name.lower().endswith(".part")
                  and p.name.lower().endswith(VALID_EXTENSIONS)]
    return candidates[0].resolve() if len(candidates) == 1 else None


def output_paths(project: Project, accession: str) -> dict[str, Path]:
    d = project.qc_dir / accession
    return {
        "dir": d,
        "full": d / f"{accession}_GRCh38_QC.tsv.gz",
        "significant": d / f"{accession}_GRCh38_significant.tsv.gz",
        "summary_tsv": d / f"{accession}_QC_summary.tsv",
        "summary_json": d / f"{accession}_QC_summary.json",
    }


def qc_study(project: Project, study: dict, force: bool) -> tuple:
    accession = study["STUDY_ACCESSION"]
    paths = output_paths(project, accession)
    if not force and all(nonempty(paths[k]) for k in ("full", "significant", "summary_tsv")):
        summary = pd.read_csv(paths["summary_tsv"], sep="\t", dtype=str).iloc[0].to_dict()
        return "SKIPPED_ALREADY_COMPLETE", "", _extra(summary)

    input_file = find_downloaded_file(project, study)
    if input_file is None:
        return "FAILED", "Downloaded summary-statistics file not found", {}
    paths["dir"].mkdir(parents=True, exist_ok=True)
    threads = project.cfg.threads
    p_threshold = float(project.cfg["qc"]["p_threshold"])

    import gwaslab as gl  # heavy import, only when QC actually runs

    banner(f"STEP 3 - GWASLAB QC | {accession}")
    print(f"Input   : {input_file}\nThreads : {threads}")
    ss = gl.Sumstats(str(input_file), fmt="auto", build="38", verbose=True)
    n_before = len(ss.data)
    ss.basic_check(remove=True, remove_dup=True, normalize=True, threads=threads, verbose=True)
    data = ss.data.copy()
    n_after_basic = len(data)

    col = {name: find_column(data.columns, cands) for name, cands in {
        "CHR": ["CHR"], "POS": ["POS"], "P": ["P"], "EA": ["EA"], "NEA": ["NEA"],
        "BETA": ["BETA"], "SE": ["SE"], "Z": ["Z"], "EAF": ["EAF"], "N": ["N"],
        "RSID": ["rsID", "RSID"], "SNPID": ["SNPID"],
    }.items()}
    missing = [k for k in ("CHR", "POS", "P") if col[k] is None]
    if missing:
        return "FAILED", "Missing essential fields after GWASLab: " + ", ".join(missing), {}

    p_col = col["P"]
    data[p_col] = pd.to_numeric(data[p_col], errors="coerce")
    n_missing_p = int(data[p_col].isna().sum())
    data = data[data[p_col].notna()].copy()
    valid_p = (data[p_col] >= 0) & (data[p_col] <= 1)
    n_invalid_p = int((~valid_p).sum())
    data = data[valid_p].copy()

    data["GWAS_ACCESSION"] = accession
    data["ANALYSIS_PHENOTYPE"] = project.phenotype
    data["ANALYSIS_ANCESTRY"] = project.ancestry_label
    data["ANALYSIS_ANCESTRY_CODE"] = project.ancestry_code
    data["PIPELINE_BUILD"] = "GRCh38"
    data = data.sort_values([col["CHR"], col["POS"]], kind="stable")
    n_after = len(data)
    significant = data[data[p_col] < p_threshold].sort_values(p_col, kind="stable")

    has = {k: col[k] is not None for k in col}
    clumping_ready = has["CHR"] and has["POS"] and has["P"] and (
        has["RSID"] or has["SNPID"] or (has["EA"] and has["NEA"]))
    # SuSiE-RSS needs alleles plus BETA+SE or Z.
    susie_ready = (has["CHR"] and has["POS"] and has["P"] and has["EA"] and has["NEA"]
                   and ((has["BETA"] and has["SE"]) or has["Z"]))
    qc_status = "PASS" if susie_ready else "WARNING" if clumping_ready else "FAIL"

    data.to_csv(paths["full"], sep="\t", index=False, compression="gzip")
    significant.to_csv(paths["significant"], sep="\t", index=False, compression="gzip")

    summary = {
        "STUDY_ACCESSION": accession, "PHENOTYPE": project.phenotype,
        "ANCESTRY": project.ancestry_label, "ANCESTRY_CODE": project.ancestry_code, "BUILD": "GRCh38",
        "INPUT_FILE": str(input_file), "FULL_QC_FILE": str(paths["full"]),
        "SIGNIFICANT_FILE": str(paths["significant"]),
        "N_BEFORE_QC": n_before, "N_AFTER_GWASLAB_BASIC_CHECK": n_after_basic, "N_AFTER_QC": n_after,
        "N_REMOVED": n_before - n_after, "N_MISSING_P_REMOVED": n_missing_p,
        "N_INVALID_P_REMOVED": n_invalid_p, "N_GENOME_WIDE_SIGNIFICANT": len(significant),
        "P_THRESHOLD": p_threshold, "MIN_P": float(data[p_col].min()) if n_after else None,
        **{f"HAS_{k}": v for k, v in has.items()},
        "CLUMPING_READY": clumping_ready, "SUSIE_READY": susie_ready, "QC_STATUS": qc_status,
        "THREADS": threads, "QC_COMPLETED_UTC": utc_now(),
    }
    pd.DataFrame([summary]).to_csv(paths["summary_tsv"], sep="\t", index=False)
    write_json(paths["summary_json"], summary)

    banner("GWAS QC COMPLETE")
    print(f"Variants before/after QC : {n_before:,} / {n_after:,}")
    print(f"Genome-wide significant  : {len(significant):,}")
    print(f"QC status                : {qc_status}")
    return "COMPLETE", "", _extra(summary)


def _extra(summary: dict) -> dict:
    return {k: summary.get(k) for k in ("QC_STATUS", "CLUMPING_READY", "SUSIE_READY",
                                       "N_AFTER_QC", "N_GENOME_WIDE_SIGNIFICANT")}
