"""Step 1 - find eligible GWAS Catalog studies for a phenotype and ancestry.

1. Downloads/caches GWAS Catalog study and ancestry metadata.
2. Matches the phenotype (plus aliases) against reported/mapped traits.
3. Keeps studies with full summary statistics.
4. Keeps studies whose discovery ancestry is exclusively the requested one.
5. Ranks by match score, then sample size.
6. Keeps studies with a harmonised summary-statistics file on the FTP site.
7. Writes gwas_manifest.tsv plus study_selection.tsv, an audit of every
   phenotype-matched study and why it was or was not selected.
"""

from __future__ import annotations

import html
import re
from urllib.parse import urljoin

import numpy as np
import pandas as pd
import requests

from gwas2m_lncrna.common import Gwas2RnaError, banner, file_lock, normalize_text, write_json
from gwas2m_lncrna.download import download, http_session
from gwas2m_lncrna.project import Project
from gwas2m_lncrna.resources import Resources

STUDIES_URL = ("https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/"
               "gwas-catalog-download-studies-v1.0.3.1.txt")
ANCESTRY_URL = ("https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/"
                "gwas-catalog-download-ancestries-v1.0.3.1.txt")
SUMMARY_ROOT = "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics"

PHENOTYPE_ALIASES = {
    "cad": ["coronary artery disease", "coronary heart disease"],
    "chd": ["coronary heart disease", "coronary artery disease"],
    "t2d": ["type 2 diabetes", "type 2 diabetes mellitus"],
    "mdd": ["major depressive disorder", "major depression", "depression"],
    "pd": ["parkinson's disease", "parkinson disease"],
    "ad": ["alzheimer's disease", "alzheimer disease"],
}

# Catalog "BROAD ANCESTRAL CATEGORY" values for each supported code.
CATALOG_ANCESTRY = {
    "EUR": "European",
    "AFR": "African",
    "EAS": "East Asian",
    "SAS": "South Asian",
    "AMR": "Hispanic or Latin American",
}

NO_SUMSTATS = {"", "na", "n/a", "nan", "none", "no", "false", "0", "not available", "nr"}


def build_terms(phenotype: str, aliases: list[str]) -> list[str]:
    terms = [phenotype] + PHENOTYPE_ALIASES.get(normalize_text(phenotype), []) + list(aliases or [])
    final, seen = [], set()
    for term in terms:
        key = normalize_text(term)
        if key and key not in seen:
            seen.add(key)
            final.append(term)
    return final


def estimate_sample_size(text) -> int:
    """Sum of all numbers >= 100 in the catalog's free-text sample description."""
    if text is None or (not isinstance(text, str) and pd.isna(text)):
        return 0
    total = 0
    for token in re.findall(r"\d[\d,]*", str(text)):
        try:
            value = int(token.replace(",", ""))
        except ValueError:
            continue
        if value >= 100:
            total += value
    return total


def parse_number(value) -> int:
    numbers = re.findall(r"\d[\d,]*", str(value or ""))
    try:
        return int(numbers[0].replace(",", "")) if numbers else 0
    except ValueError:
        return 0


def has_full_summary_stats(value) -> bool:
    if value is None or (not isinstance(value, str) and pd.isna(value)):
        return False
    return str(value).strip().lower() not in NO_SUMSTATS


def contains_term(series: pd.Series, term: str) -> pd.Series:
    term = normalize_text(term)
    if not term:
        return pd.Series(False, index=series.index)
    if " " not in term and len(term) <= 5:  # short acronyms need word boundaries
        pattern = rf"(?<![a-z0-9]){re.escape(term)}(?![a-z0-9])"
    else:
        pattern = re.escape(term)
    return series.str.contains(pattern, regex=True, na=False)


def find_phenotype_studies(studies: pd.DataFrame, terms: list[str]) -> pd.DataFrame:
    def column(name):
        if name in studies.columns:
            return studies[name].fillna("").map(normalize_text)
        return pd.Series("", index=studies.index)

    disease, mapped, title = column("DISEASE/TRAIT"), column("MAPPED_TRAIT"), column("STUDY")
    scores = pd.Series(0, index=studies.index, dtype=int)
    for term in terms:
        t = normalize_text(term)
        for mask, score in (
            (disease == t, 100), (mapped == t, 95), (contains_term(disease, term), 80),
            (contains_term(mapped, term), 70), (contains_term(title, term), 40),
        ):
            scores.loc[mask] = np.maximum(scores.loc[mask], score)
    result = studies[scores > 0].copy()
    result["MATCH_SCORE"] = scores.loc[result.index]
    return result.drop_duplicates(subset=["STUDY ACCESSION"])


def ancestry_information(ancestry: pd.DataFrame, target: str) -> dict:
    target_norm = normalize_text(target)
    if "STAGE" in ancestry.columns:
        initial = ancestry[ancestry["STAGE"].fillna("").astype(str).str.lower().str.contains("initial")]
        if not initial.empty:
            ancestry = initial
    records = {}
    for accession, group in ancestry.groupby("STUDY ACCESSION"):
        categories = {normalize_text(x) for x in group["BROAD ANCESTRAL CATEGORY"].dropna() if normalize_text(x)}
        total_n = 0
        if "NUMBER OF INDIVIDUALS" in group.columns:
            total_n = sum(parse_number(x) for x in group["NUMBER OF INDIVIDUALS"])
        records[str(accession)] = {
            # Strict: every discovery ancestry category must be the requested one.
            "STRICT_MATCH": categories == {target_norm},
            "ANCESTRY_CATEGORIES": " | ".join(sorted(categories)),
            "ANCESTRY_N": total_n,
        }
    return records


def gwas_bucket(accession: str) -> str:
    numeric = accession.replace("GCST", "")
    number = int(numeric)
    width = max(len(numeric), 8)
    start = ((number - 1) // 1000) * 1000 + 1
    return f"GCST{start:0{width}d}-GCST{start + 999:0{width}d}"


def find_harmonised_file(session: requests.Session, accession: str) -> dict | None:
    directory = f"{SUMMARY_ROOT}/{gwas_bucket(accession)}/{accession}/harmonised/"
    try:
        response = session.get(directory, timeout=(20, 60))
    except requests.RequestException:
        return None
    if response.status_code != 200:
        return None
    files = []
    for href in re.findall(r'href=["\']([^"\']+)["\']', response.text, flags=re.I):
        name = html.unescape(href).split("/")[-1]
        lower = name.lower()
        if not name or lower.endswith((".tbi", ".log", ".yaml", ".yml", ".md5", ".json")):
            continue
        if ".h.tsv" in lower or lower.endswith((".tsv.gz", ".tsv")):
            files.append(name)
    if not files:
        return None

    def priority(name):
        lower = name.lower()
        for rank, suffix in enumerate((".h.tsv.gz", ".h.tsv", ".tsv.gz", ".tsv")):
            if lower.endswith(suffix):
                return rank
        return 10

    filename = sorted(files, key=priority)[0]
    return {"URL": urljoin(directory, filename), "FILENAME": filename}


def _selection_table(candidates, phenotype, ancestry_label, info) -> pd.DataFrame:
    rows = []
    for _, row in candidates.iterrows():
        accession = str(row["STUDY ACCESSION"])
        meta = info.get(accession, {})
        rows.append({
            "STUDY_ACCESSION": accession, "REQUESTED_PHENOTYPE": phenotype,
            "DISEASE_TRAIT": row.get("DISEASE/TRAIT", ""), "MAPPED_TRAIT": row.get("MAPPED_TRAIT", ""),
            "MATCH_SCORE": int(row["MATCH_SCORE"]), "REQUESTED_ANCESTRY": ancestry_label,
            "ANCESTRY_CATEGORIES": meta.get("ANCESTRY_CATEGORIES", ""),
            "ANCESTRY_MATCH": bool(meta.get("STRICT_MATCH", False)),
            "INITIAL_SAMPLE_SIZE": row.get("INITIAL SAMPLE SIZE", ""),
            "SUMMARY_STATS_AVAILABLE": has_full_summary_stats(row.get("FULL SUMMARY STATISTICS")),
            "HARMONISED_FILE_AVAILABLE": "NOT_CHECKED", "DOWNLOAD_URL": "",
            "ELIGIBLE": False, "EXCLUSION_REASON": "", "RANK": "",
        })
    table = pd.DataFrame(rows, columns=[
        "STUDY_ACCESSION", "REQUESTED_PHENOTYPE", "DISEASE_TRAIT", "MAPPED_TRAIT", "MATCH_SCORE",
        "REQUESTED_ANCESTRY", "ANCESTRY_CATEGORIES", "ANCESTRY_MATCH", "INITIAL_SAMPLE_SIZE",
        "SUMMARY_STATS_AVAILABLE", "HARMONISED_FILE_AVAILABLE", "DOWNLOAD_URL", "ELIGIBLE",
        "EXCLUSION_REASON", "RANK",
    ])
    table.index = table["STUDY_ACCESSION"].tolist()
    no_sumstats = ~table["SUMMARY_STATS_AVAILABLE"].astype(bool)
    table.loc[no_sumstats, "EXCLUSION_REASON"] = "NO_FULL_SUMMARY_STATISTICS"
    no_anc = (table["EXCLUSION_REASON"] == "") & ~table["ANCESTRY_MATCH"].astype(bool)
    table.loc[no_anc & (table["ANCESTRY_CATEGORIES"] == ""), "EXCLUSION_REASON"] = "ANCESTRY_NOT_REPORTED"
    table.loc[no_anc & (table["ANCESTRY_CATEGORIES"] != ""), "EXCLUSION_REASON"] = "ANCESTRY_MISMATCH"
    return table


def load_catalog(project: Project, refresh: bool) -> tuple[pd.DataFrame, pd.DataFrame]:
    cache = Resources(project.cfg).gwas_catalog_dir
    studies_file = cache / "gwas_catalog_studies.tsv"
    ancestry_file = cache / "gwas_catalog_ancestry.tsv"
    with file_lock(cache / ".catalog.lock"):
        for url, path in ((STUDIES_URL, studies_file), (ANCESTRY_URL, ancestry_file)):
            if refresh:
                path.unlink(missing_ok=True)
            if path.is_file():
                print(f"Using cached GWAS Catalog metadata: {path}")
            download(url, path)
    read = dict(sep="\t", dtype=str, low_memory=False)
    return pd.read_csv(studies_file, **read), pd.read_csv(ancestry_file, **read)


def discover(project: Project, aliases: list[str] | None = None, max_studies: int | None = None,
             refresh_catalog: bool | None = None, force: bool = False) -> dict:
    settings = project.cfg["discovery"]
    aliases = list(settings.get("aliases") or []) + list(aliases or [])
    max_studies = int(settings.get("max_studies") or 0) if max_studies is None else int(max_studies)
    refresh = bool(settings.get("refresh_catalog")) if refresh_catalog is None else refresh_catalog
    if max_studies < 0:
        raise Gwas2RnaError("max_studies must be >= 0")

    banner(f"STEP 1 - GWAS DISCOVERY | {project.phenotype} | {project.ancestry_code}")
    if project.manifest_path.is_file() and not force:
        manifest = project.manifest()
        print(f"Using existing manifest ({len(manifest)} studies): {project.manifest_path}")
        print("Use --refresh-discovery (run) or --force (discover) to search the catalog again.")
        return {"n_eligible": len(manifest), "manifest": project.manifest_path, "reused": True}

    out = project.discovery_dir
    out.mkdir(parents=True, exist_ok=True)
    catalog_label = CATALOG_ANCESTRY[project.ancestry_code]
    terms = build_terms(project.phenotype, aliases)
    state = {"REQUESTED_PHENOTYPE": project.phenotype, "REQUESTED_ANCESTRY": catalog_label,
             "ANCESTRY_CODE": project.ancestry_code, "SEARCH_TERMS": terms, "MAX_STUDIES": max_studies,
             "GWAS_CATALOG_STUDIES_URL": STUDIES_URL, "GWAS_CATALOG_ANCESTRY_URL": ANCESTRY_URL}
    print("Search terms: " + "; ".join(terms))
    selection = None
    try:
        studies, ancestry = load_catalog(project, refresh)
        state["N_CATALOG_STUDIES"] = int(len(studies))
        info = ancestry_information(ancestry, catalog_label)

        candidates = find_phenotype_studies(studies, terms)
        selection = _selection_table(candidates, project.phenotype, catalog_label, info)
        state["N_PHENOTYPE_MATCHED"] = int(len(candidates))
        print(f"Phenotype matches              : {len(candidates):,}")
        if candidates.empty:
            raise Gwas2RnaError(f"No GWAS Catalog studies match phenotype '{project.phenotype}'. "
                                "Try a fuller trait name or add --alias terms.")
        if "FULL SUMMARY STATISTICS" not in candidates.columns:
            raise Gwas2RnaError("GWAS Catalog studies file lacks 'FULL SUMMARY STATISTICS'.")

        candidates = candidates[candidates["FULL SUMMARY STATISTICS"].map(has_full_summary_stats)].copy()
        state["N_WITH_FULL_SUMMARY_STATISTICS"] = int(len(candidates))
        print(f"With full summary statistics   : {len(candidates):,}")

        strict = {acc for acc, meta in info.items() if meta["STRICT_MATCH"]}
        candidates = candidates[candidates["STUDY ACCESSION"].astype(str).isin(strict)].copy()
        state["N_FULL_SUMMARY_STATISTICS_AND_ANCESTRY_MATCHED"] = int(len(candidates))
        print(f"Exclusively {catalog_label:<19}: {len(candidates):,}")
        if candidates.empty:
            raise Gwas2RnaError(f"No {catalog_label}-only studies with full summary statistics.")

        accessions = candidates["STUDY ACCESSION"].astype(str)
        candidates["ANCESTRY_N"] = [int(info.get(a, {}).get("ANCESTRY_N", 0)) for a in accessions]
        candidates["ANCESTRY_CATEGORIES"] = [info.get(a, {}).get("ANCESTRY_CATEGORIES", "") for a in accessions]
        candidates["ESTIMATED_INITIAL_N"] = (
            candidates["INITIAL SAMPLE SIZE"].map(estimate_sample_size)
            if "INITIAL SAMPLE SIZE" in candidates.columns else 0
        )
        candidates["RANKING_N"] = candidates[["ANCESTRY_N", "ESTIMATED_INITIAL_N"]].max(axis=1)
        candidates = candidates.sort_values(["MATCH_SCORE", "RANKING_N"], ascending=[False, False])
        for rank, accession in enumerate(candidates["STUDY ACCESSION"].astype(str), start=1):
            selection.at[accession, "RANK"] = rank
        candidates.to_csv(out / "candidate_studies.tsv", sep="\t", index=False)

        banner("CHECKING HARMONISED SUMMARY-STATISTICS FILES")
        session = http_session()
        eligible = []
        for _, row in candidates.iterrows():
            accession = str(row["STUDY ACCESSION"])
            found = find_harmonised_file(session, accession)
            if found is None:
                print(f"{accession}: no harmonised file")
                selection.at[accession, "HARMONISED_FILE_AVAILABLE"] = False
                selection.at[accession, "EXCLUSION_REASON"] = "NO_HARMONISED_FILE"
                continue
            print(f"{accession}: {found['FILENAME']}")
            selection.at[accession, "HARMONISED_FILE_AVAILABLE"] = True
            selection.at[accession, "DOWNLOAD_URL"] = found["URL"]
            selection.at[accession, "ELIGIBLE"] = True
            eligible.append((row, found))
            if max_studies and len(eligible) >= max_studies:
                break
        if not eligible:
            raise Gwas2RnaError("No matching studies have a harmonised summary-statistics file.")

        manifest = pd.DataFrame([{
            "ARRAY_ID": i,
            "STUDY_ACCESSION": str(row["STUDY ACCESSION"]),
            "DOWNLOAD_URL": found["URL"], "FILE_NAME": found["FILENAME"],
            "PHENOTYPE": project.phenotype, "ANCESTRY": catalog_label,
            "ANCESTRY_CODE": project.ancestry_code, "ANCESTRY_LABEL": catalog_label,
            "RANKING_N": int(row["RANKING_N"]), "MATCH_SCORE": int(row["MATCH_SCORE"]),
            "DISEASE_TRAIT": row.get("DISEASE/TRAIT", ""), "MAPPED_TRAIT": row.get("MAPPED_TRAIT", ""),
            "FIRST_AUTHOR": row.get("FIRST AUTHOR", ""), "PUBMEDID": row.get("PUBMEDID", ""),
            "INITIAL_SAMPLE_SIZE": row.get("INITIAL SAMPLE SIZE", ""),
        } for i, (row, found) in enumerate(eligible, start=1)])
        manifest.to_csv(project.manifest_path, sep="\t", index=False)
        state["N_ELIGIBLE"] = int(len(manifest))
        state["DOWNLOAD_MANIFEST"] = str(project.manifest_path)

        banner("DISCOVERY COMPLETE")
        print(manifest[["ARRAY_ID", "STUDY_ACCESSION", "RANKING_N", "DISEASE_TRAIT"]].to_string(index=False))
        print(f"\nManifest: {project.manifest_path}")
        return {"n_eligible": len(manifest), "manifest": project.manifest_path, "reused": False}
    finally:
        if selection is not None:
            selection = selection.copy()
            unchecked = (selection["EXCLUSION_REASON"] == "") & ~selection["ELIGIBLE"].astype(bool)
            selection.loc[unchecked, "EXCLUSION_REASON"] = "NOT_CHECKED_MAX_STUDIES_LIMIT"
            selection.to_csv(out / "study_selection.tsv", sep="\t", index=False)
        write_json(out / "discovery_summary.json", state)
