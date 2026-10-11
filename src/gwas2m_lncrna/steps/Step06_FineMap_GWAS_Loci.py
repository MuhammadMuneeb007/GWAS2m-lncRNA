"""Step 6 - regional extraction + ancestry-matched LD + SuSiE-RSS for one study.

1. Resolves the GWAS sample size from GWAS Catalog metadata (effective N
   4/(1/Ncase+1/Ncontrol) for case-control studies when counts are available).
2. Expands each Step 5 lead by locus_window_kb and merges overlapping windows.
3. Streams the cleaned GWAS once, extracting all regional variants.
4. Matches to the 1000G panel and orients Z scores to the reference allele.
5. Builds a signed LD matrix with PLINK2 and runs susieR::susie_rss.
6. Writes PIPs, 95% credible sets and diagnostics; loci are resumable.

Outputs in 06_finemapping/<GCST>/ keep the original file names.
"""

from __future__ import annotations

import json
import re
import time
import traceback
from collections import defaultdict
from importlib import resources as importlib_resources
from pathlib import Path
from typing import Any
from urllib.parse import urljoin

import numpy as np
import pandas as pd
import requests

from gwas2m_lncrna.common import (
    banner, clean_allele, clean_chr, find_column, nonempty, parse_numbers, read_json,
    run_command, safe_int, utc_now, write_json,
)
from gwas2m_lncrna.download import http_session
from gwas2m_lncrna.project import Project
from gwas2m_lncrna.resources import Resources
from gwas2m_lncrna.tools import find_tool, require_susier

GWAS_API_BASE = "https://www.ebi.ac.uk/gwas/rest/api/v2"
CHUNK_SIZE = 500_000
LOCUS_COMPLETE_MARKER = "LOCUS_COMPLETE.json"
LOCUS_REQUIRED_OUTPUTS = ("susie_variants.tsv", "susie_summary.tsv", "susie_fit.rds")


def r_script_path() -> Path:
    return Path(str(importlib_resources.files("gwas2m_lncrna") / "data" / "run_susie_rss.R"))


# =============================================================================
# SAMPLE SIZE FROM OFFICIAL METADATA
# =============================================================================

def effective_n(n_cases: int | None, n_controls: int | None) -> int | None:
    if not n_cases or not n_controls or n_cases <= 0 or n_controls <= 0:
        return None
    return int(round(4.0 / ((1.0 / n_cases) + (1.0 / n_controls))))


def _fetch_json(session, url: str, destination: Path, refresh: bool):
    if destination.exists() and destination.stat().st_size > 0 and not refresh:
        try:
            return json.loads(destination.read_text(encoding="utf-8"))
        except ValueError:
            pass
    response = session.get(url, timeout=(20, 120))
    response.raise_for_status()
    payload = response.json()
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    time.sleep(0.08)  # stay below the documented API throttle
    return payload


def _fetch_text(session, urls: list[str], destination: Path, refresh: bool):
    if destination.exists() and destination.stat().st_size > 0 and not refresh:
        return destination.read_text(encoding="utf-8", errors="replace"), "cache"
    last_error = None
    for url in urls:
        try:
            response = session.get(url, timeout=(20, 120))
            if response.status_code == 200 and response.text.strip():
                destination.parent.mkdir(parents=True, exist_ok=True)
                destination.write_text(response.text, encoding="utf-8")
                return response.text, url
            last_error = f"HTTP {response.status_code}"
        except requests.RequestException as exc:
            last_error = repr(exc)
    return None, last_error


def _parse_yaml(text: str | None) -> dict:
    if not text:
        return {}
    try:
        import yaml
        loaded = yaml.safe_load(text)
        return loaded if isinstance(loaded, dict) else {}
    except Exception:
        return {}


def _recursive_objects(obj):
    if isinstance(obj, dict):
        yield obj
        for value in obj.values():
            yield from _recursive_objects(value)
    elif isinstance(obj, list):
        for value in obj:
            yield from _recursive_objects(value)


def _norm_key(key) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(key).lower())


def ancestry_sample_size_from_api(payload: Any) -> int | None:
    numbers = []
    for item in _recursive_objects(payload):
        stage_text = " ".join(str(v) for k, v in item.items() if "stage" in _norm_key(k)).lower()
        if stage_text and not ("initial" in stage_text or "discovery" in stage_text):
            continue
        for key, value in item.items():
            norm = _norm_key(key)
            if "numberofindividuals" in norm or norm in {"samplesize", "samplecount", "numberindividuals", "n"}:
                number = safe_int(value)
                if number and number >= 100:
                    numbers.append(number)
    return sum(numbers) if numbers else None


def study_sample_size_from_api(payload: Any) -> int | None:
    wanted = ["initialsamplesize", "discoverysample", "samplesize"]
    candidates = []
    for item in _recursive_objects(payload):
        for key, value in item.items():
            if any(w in _norm_key(key) for w in wanted):
                candidates.extend(parse_numbers(value))
                direct = safe_int(value)
                if direct and direct >= 100:
                    candidates.append(direct)
    return max(candidates) if candidates else None


def _yaml_from_directory(session, download_url: str, accession: str) -> str | None:
    directory = download_url.rsplit("/", 1)[0] + "/"
    try:
        response = session.get(directory, timeout=(20, 60))
        if response.status_code != 200:
            return None
        names = [h for h in re.findall(r'href=["\']([^"\']+)["\']', response.text, flags=re.I)
                 if h.split("/")[-1].lower().endswith("-meta.yaml") and accession.lower() in h.lower()]
        return urljoin(directory, sorted(names)[0]) if names else None
    except requests.RequestException:
        return None


def resolve_sample_size(project: Project, study: dict, refresh: bool = False) -> dict:
    accession = study["STUDY_ACCESSION"]
    metadata_dir = project.finemap_dir / "metadata" / accession
    metadata_dir.mkdir(parents=True, exist_ok=True)
    session = http_session(retries=5)
    study_url = f"{GWAS_API_BASE}/studies/{accession}"
    ancestries_url = f"{GWAS_API_BASE}/studies/{accession}/ancestries"
    study_payload = ancestry_payload = None
    study_error = ancestry_error = None
    try:
        study_payload = _fetch_json(session, study_url, metadata_dir / "study_v2.json", refresh)
    except (requests.RequestException, ValueError) as exc:
        study_error = repr(exc)
        try:
            study_payload = _fetch_json(session, f"{GWAS_API_BASE}/studies?accession_id={accession}",
                                        metadata_dir / "study_v2.json", refresh)
            study_error = None
        except (requests.RequestException, ValueError) as exc2:
            study_error = repr(exc2)
    try:
        ancestry_payload = _fetch_json(session, ancestries_url, metadata_dir / "ancestries_v2.json", refresh)
    except (requests.RequestException, ValueError) as exc:
        ancestry_error = repr(exc)

    download_url = study.get("DOWNLOAD_URL") or None
    yaml_text = yaml_source = None
    yaml_path = metadata_dir / "summary_statistics_metadata.yaml"
    if download_url:
        yaml_text, yaml_source = _fetch_text(session, [download_url + "-meta.yaml"], yaml_path, refresh)
        if yaml_text is None:
            extra = _yaml_from_directory(session, download_url, accession)
            if extra:
                yaml_text, yaml_source = _fetch_text(session, [extra], yaml_path, refresh)
    meta = _parse_yaml(yaml_text)

    n_cases, n_controls = safe_int(meta.get("case_count")), safe_int(meta.get("control_count"))
    yaml_n = safe_int(meta.get("sample_size"))
    neff = effective_n(n_cases, n_controls)
    api_ancestry_n = ancestry_sample_size_from_api(ancestry_payload)
    api_study_n = study_sample_size_from_api(study_payload)
    ranking_n = safe_int(study.get("RANKING_N"))
    catalog_numbers = parse_numbers(study.get("INITIAL_SAMPLE_SIZE"))
    catalog_n = sum(catalog_numbers) if catalog_numbers else None
    fallback = safe_int(project.cfg["finemapping"].get("sample_size"))

    for value, source in (
        (neff, "GWAS_CATALOG_SUMSTATS_YAML_CASE_CONTROL_EFFECTIVE_N"),
        (yaml_n, "GWAS_CATALOG_SUMSTATS_YAML_SAMPLE_SIZE"),
        (api_ancestry_n, "GWAS_CATALOG_API_V2_ANCESTRY"),
        (api_study_n, "GWAS_CATALOG_API_V2_STUDY"),
        (ranking_n, "STEP01_RANKING_N"),
        (catalog_n, "STEP01_INITIAL_SAMPLE_SIZE_PARSED"),
        (fallback, "USER_FALLBACK"),
        (None, "UNRESOLVED"),
    ):
        if value or source == "UNRESOLVED":
            n_gwas, n_source = value, source
            break

    result = {
        "STUDY_ACCESSION": accession, "GWAS_API_STUDY_URL": study_url,
        "GWAS_API_ANCESTRIES_URL": ancestries_url, "GWAS_API_STUDY_ERROR": study_error,
        "GWAS_API_ANCESTRIES_ERROR": ancestry_error, "SUMMARY_METADATA_SOURCE": yaml_source,
        "CASE_CONTROL_STUDY": meta.get("case_control_study"), "N_CASES": n_cases,
        "N_CONTROLS": n_controls, "N_EFFECTIVE_CASE_CONTROL": neff, "YAML_SAMPLE_SIZE": yaml_n,
        "API_ANCESTRY_SAMPLE_SIZE": api_ancestry_n, "API_STUDY_SAMPLE_SIZE": api_study_n,
        "STEP01_RANKING_N": ranking_n, "STEP01_INITIAL_SAMPLE_SIZE": study.get("INITIAL_SAMPLE_SIZE"),
        "RESOLVED_GWAS_N": n_gwas, "RESOLVED_GWAS_N_SOURCE": n_source, "METADATA_RESOLVED_UTC": utc_now(),
    }
    write_json(metadata_dir / "resolved_metadata.json", result)
    return result


# =============================================================================
# LOCI
# =============================================================================

def load_leads_and_build_loci(lead_file: Path, window_kb: int) -> pd.DataFrame:
    columns = ["LOCUS_ID", "CHR", "LOCUS_START", "LOCUS_END", "N_INDEPENDENT_SIGNALS",
               "LEAD_IDS", "LEAD_POSITIONS", "MIN_LEAD_P"]
    leads = pd.read_csv(lead_file, sep="\t", low_memory=False)
    if leads.empty:
        return pd.DataFrame(columns=columns)
    chr_col = find_column(leads.columns, ["CHROM", "#CHROM", "CHR", "SOURCE_CHROMOSOME"])
    pos_col = find_column(leads.columns, ["POS", "BP"])
    p_col = find_column(leads.columns, ["P", "P_VALUE"])
    id_col = find_column(leads.columns, ["ID", "SNP", "SNPID", "rsID", "RSID"])
    if chr_col is None or pos_col is None:
        raise RuntimeError(f"Lead file lacks chromosome/position columns: {lead_file}")
    work = pd.DataFrame({
        "CHR": leads[chr_col].map(clean_chr),
        "POS": pd.to_numeric(leads[pos_col], errors="coerce"),
        "P": pd.to_numeric(leads[p_col], errors="coerce") if p_col else np.nan,
        "ID": leads[id_col].astype(str) if id_col else "",
    }).dropna(subset=["CHR", "POS"])
    work["CHR"] = work["CHR"].astype(int)
    work["POS"] = work["POS"].astype(int)
    work = work.sort_values(["CHR", "POS", "P"], kind="stable")

    window = window_kb * 1000
    merged = []
    for chrom, group in work.groupby("CHR", sort=True):
        current = None
        for _, row in group.iterrows():
            start, end = max(1, int(row["POS"]) - window), int(row["POS"]) + window
            signal = {"id": str(row["ID"]), "pos": int(row["POS"]),
                      "p": float(row["P"]) if pd.notna(row["P"]) else np.nan}
            if current is None or start > current["end"]:
                if current is not None:
                    merged.append(current)
                current = {"chr": int(chrom), "start": start, "end": end, "signals": [signal]}
            else:
                current["end"] = max(current["end"], end)
                current["signals"].append(signal)
        if current is not None:
            merged.append(current)

    records = []
    for i, item in enumerate(merged, start=1):
        pvals = [x["p"] for x in item["signals"] if np.isfinite(x["p"])]
        records.append({
            "LOCUS_ID": f"L{i:03d}", "CHR": item["chr"], "LOCUS_START": item["start"],
            "LOCUS_END": item["end"], "N_INDEPENDENT_SIGNALS": len(item["signals"]),
            "LEAD_IDS": ";".join(x["id"] for x in item["signals"]),
            "LEAD_POSITIONS": ";".join(str(x["pos"]) for x in item["signals"]),
            "MIN_LEAD_P": min(pvals) if pvals else np.nan,
        })
    return pd.DataFrame(records, columns=columns)


# =============================================================================
# REGIONAL GWAS EXTRACTION
# =============================================================================

def inspect_qc_columns(path: Path) -> dict:
    cols = list(pd.read_csv(path, sep="\t", compression="infer", nrows=0).columns)
    mapping = {name: find_column(cols, cands) for name, cands in {
        "CHR": ["CHR"], "POS": ["POS"], "P": ["P"], "EA": ["EA"], "NEA": ["NEA"],
        "BETA": ["BETA"], "SE": ["SE"], "Z": ["Z"], "OR": ["OR"], "OR_95U": ["OR_95U"],
        "OR_95L": ["OR_95L"], "SNPID": ["SNPID", "rsID", "RSID"],
    }.items()}
    missing = [x for x in ("CHR", "POS", "P", "EA", "NEA") if mapping[x] is None]
    if missing:
        raise RuntimeError(f"QC GWAS missing required columns {missing}: {path}")
    effect_ready = bool(mapping["Z"] or (mapping["BETA"] and mapping["SE"]) or (
        mapping["OR"] and (mapping["SE"] or (mapping["OR_95U"] and mapping["OR_95L"]))))
    if not effect_ready:
        raise RuntimeError(f"QC GWAS needs Z, BETA+SE, OR+SE, or OR+95% CI for SuSiE-RSS: {path}")
    return mapping


def extract_regional_gwas_once(qc_file: Path, loci: pd.DataFrame, regional_dir: Path) -> dict[str, Path]:
    mapping = inspect_qc_columns(qc_file)
    usecols = list(dict.fromkeys(x for x in mapping.values() if x is not None))
    by_chr = defaultdict(list)
    for _, row in loci.iterrows():
        by_chr[int(row["CHR"])].append(row.to_dict())
    buffers: dict[str, list[pd.DataFrame]] = {str(x): [] for x in loci["LOCUS_ID"]}
    total = extracted = 0
    for chunk in pd.read_csv(qc_file, sep="\t", compression="infer", usecols=usecols,
                             chunksize=CHUNK_SIZE, low_memory=False):
        total += len(chunk)
        canon = pd.DataFrame(index=chunk.index)
        canon["CHR"] = chunk[mapping["CHR"]].map(clean_chr)
        canon["POS"] = pd.to_numeric(chunk[mapping["POS"]], errors="coerce")
        canon["P"] = pd.to_numeric(chunk[mapping["P"]], errors="coerce")
        canon["EA"] = chunk[mapping["EA"]].map(clean_allele)
        canon["NEA"] = chunk[mapping["NEA"]].map(clean_allele)
        canon["SNPID"] = chunk[mapping["SNPID"]].astype(str) if mapping["SNPID"] else ""
        # Effect statistics: Z, else BETA/SE, else log(OR)/SE, else SE from the 95% CI.
        if mapping["BETA"]:
            canon["BETA"] = pd.to_numeric(chunk[mapping["BETA"]], errors="coerce")
        elif mapping["OR"]:
            odds = pd.to_numeric(chunk[mapping["OR"]], errors="coerce")
            canon["BETA"] = np.log(odds.where(odds > 0))
        else:
            canon["BETA"] = np.nan
        if mapping["SE"]:
            canon["SE"] = pd.to_numeric(chunk[mapping["SE"]], errors="coerce")
        elif mapping["OR_95U"] and mapping["OR_95L"]:
            upper = pd.to_numeric(chunk[mapping["OR_95U"]], errors="coerce")
            lower = pd.to_numeric(chunk[mapping["OR_95L"]], errors="coerce")
            canon["SE"] = (np.log(upper.where(upper > 0)) - np.log(lower.where(lower > 0))) / (2.0 * 1.959963984540054)
        else:
            canon["SE"] = np.nan
        canon.loc[(~np.isfinite(canon["SE"])) | (canon["SE"] <= 0), "SE"] = np.nan
        canon["Z"] = (pd.to_numeric(chunk[mapping["Z"]], errors="coerce") if mapping["Z"]
                      else canon["BETA"] / canon["SE"])
        canon = canon.dropna(subset=["CHR", "POS", "P", "EA", "NEA", "Z"])
        if canon.empty:
            continue
        canon["CHR"] = canon["CHR"].astype(int)
        canon["POS"] = canon["POS"].astype(int)
        canon = canon[np.isfinite(canon["Z"])]
        for chrom, locus_list in by_chr.items():
            chrom_data = canon[canon["CHR"] == chrom]
            if chrom_data.empty:
                continue
            for locus in locus_list:
                mask = chrom_data["POS"].between(int(locus["LOCUS_START"]), int(locus["LOCUS_END"]))
                if mask.any():
                    selected = chrom_data.loc[mask].copy()
                    buffers[str(locus["LOCUS_ID"])].append(selected)
                    extracted += len(selected)
    regional_dir.mkdir(parents=True, exist_ok=True)
    paths = {}
    for locus_id in loci["LOCUS_ID"].astype(str):
        path = regional_dir / f"{locus_id}.tsv.gz"
        if buffers[locus_id]:
            frame = (pd.concat(buffers[locus_id], ignore_index=True)
                     .sort_values(["CHR", "POS", "P"], kind="stable")
                     .drop_duplicates(subset=["CHR", "POS", "EA", "NEA"], keep="first"))
        else:
            frame = pd.DataFrame(columns=["CHR", "POS", "P", "EA", "NEA", "SNPID", "BETA", "SE", "Z"])
        frame.to_csv(path, sep="\t", index=False, compression="gzip")
        paths[locus_id] = path
    print(f"GWAS rows scanned : {total:,}\nRegional rows hit : {extracted:,}")
    return paths


# =============================================================================
# REFERENCE MATCHING + LD
# =============================================================================

def create_regional_pvar(plink2: Path, prefix_ref: Path, chromosome: int, start: int, end: int,
                         locus_dir: Path) -> Path:
    prefix = locus_dir / "reference_region"
    pvar = Path(f"{prefix}.pvar")
    if nonempty(pvar):
        return pvar
    run_command([plink2, "--pfile", prefix_ref, "--chr", chromosome, "--from-bp", start, "--to-bp", end,
                 "--make-just-pvar", "--out", prefix], locus_dir / "reference_region.command.log")
    if not nonempty(pvar):
        raise RuntimeError(f"Regional reference PVAR not generated for chr{chromosome}:{start}-{end}")
    return pvar


def read_pvar(path: Path) -> pd.DataFrame:
    records, header = [], None
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                header = line.lstrip("#").rstrip("\n").split("\t")
                continue
            parts = line.rstrip("\n").split("\t")
            if header and len(parts) >= len(header):
                row = dict(zip(header, parts))
                chrom = clean_chr(row.get("CHROM", row.get("#CHROM", "")))
                pos = safe_int(row.get("POS"))
                if chrom and pos:
                    records.append({"REFERENCE_CHR": chrom, "REFERENCE_POS": pos,
                                    "REFERENCE_ID": row.get("ID", ""),
                                    "REFERENCE_REF": clean_allele(row.get("REF")),
                                    "REFERENCE_ALT": clean_allele(row.get("ALT"))})
    return pd.DataFrame(records)


def _allele_key(a, b):
    a, b = clean_allele(a), clean_allele(b)
    return None if a is None or b is None else tuple(sorted((a, b)))


def _clean_rsid(value) -> str | None:
    if value is None or (not isinstance(value, str) and pd.isna(value)):
        return None
    value = str(value).strip()
    return value if re.fullmatch(r"rs\d+", value, flags=re.I) else None


def match_gwas_to_reference(regional: pd.DataFrame, reference: pd.DataFrame) -> pd.DataFrame:
    by_id, by_position = defaultdict(list), defaultdict(list)
    for record in reference.to_dict("records"):
        by_id[str(record["REFERENCE_ID"])].append(record)
        by_position[int(record["REFERENCE_POS"])].append(record)
    results = []
    for row in regional.to_dict("records"):
        ea, nea = clean_allele(row.get("EA")), clean_allele(row.get("NEA"))
        if ea is None or nea is None:
            continue
        key = _allele_key(ea, nea)
        rsid = _clean_rsid(row.get("SNPID"))
        selected = method = None
        if rsid and rsid in by_id:
            candidates = [c for c in by_id[rsid] if _allele_key(c["REFERENCE_REF"], c["REFERENCE_ALT"]) == key]
            if len(candidates) == 1:
                selected, method = candidates[0], "rsID+alleles"
        if selected is None:
            candidates = [c for c in by_position.get(int(row["POS"]), [])
                          if _allele_key(c["REFERENCE_REF"], c["REFERENCE_ALT"]) == key]
            if len(candidates) == 1:
                selected, method = candidates[0], "POS+alleles"
        if selected is None:
            continue
        ref, alt = clean_allele(selected["REFERENCE_REF"]), clean_allele(selected["REFERENCE_ALT"])
        z = float(row["Z"])
        if not np.isfinite(z):
            continue
        if ea == ref and nea == alt:
            z_ref, orientation = z, "EA=REF"
        elif ea == alt and nea == ref:
            z_ref, orientation = -z, "EA=ALT_FLIPPED"
        else:
            continue
        row.update({"REFERENCE_ID": selected["REFERENCE_ID"], "REFERENCE_POS": selected["REFERENCE_POS"],
                    "REFERENCE_REF": ref, "REFERENCE_ALT": alt, "MATCH_METHOD": method,
                    "ORIENTATION": orientation, "Z_REF": z_ref})
        if pd.notna(row.get("BETA", np.nan)):
            row["BETA_REF"] = float(row["BETA"]) if orientation == "EA=REF" else -float(row["BETA"])
        results.append(row)
    if not results:
        return pd.DataFrame()
    matched = pd.DataFrame(results).sort_values("P", kind="stable").drop_duplicates("REFERENCE_ID", keep="first")
    return matched.reset_index(drop=True)


def _parse_lead_positions(value: Any) -> list[int]:
    positions = []
    for token in str(value or "").split(";"):
        try:
            pos = int(float(token.strip()))
        except ValueError:
            continue
        if pos > 0:
            positions.append(pos)
    return sorted(set(positions))


def _balanced_ld_subset(matched: pd.DataFrame, max_variants: int, lead_positions: list[int]) -> pd.DataFrame:
    """Deterministically cap an oversized locus while keeping lead neighbourhoods."""
    if len(matched) <= max_variants:
        return matched.copy()
    work = matched.copy().reset_index(drop=True)
    work["P"] = pd.to_numeric(work["P"], errors="coerce")
    work["POS"] = pd.to_numeric(work["POS"], errors="coerce")
    leads = sorted({int(x) for x in lead_positions if int(x) > 0})
    if not leads:
        return work.sort_values(["P", "POS", "REFERENCE_ID"], kind="stable").head(max_variants).reset_index(drop=True)
    pos = work["POS"].to_numpy(dtype=float)
    distances = np.abs(pos[:, None] - np.asarray(leads, dtype=float)[None, :])
    nearest = np.argmin(distances, axis=1)
    work["_NEAREST_LEAD_INDEX"] = nearest
    work["_LEAD_DISTANCE"] = distances[np.arange(len(work)), nearest]
    quota = max(1, max_variants // len(leads))
    selected: list[int] = []
    for lead_index in range(len(leads)):
        group = work[work["_NEAREST_LEAD_INDEX"] == lead_index]
        if not group.empty:
            group = group.sort_values(["_LEAD_DISTANCE", "P", "POS", "REFERENCE_ID"], kind="stable")
            selected.extend(group.head(quota).index.tolist())
    selected = list(dict.fromkeys(selected))[:max_variants]
    remaining = max_variants - len(selected)
    if remaining > 0:
        rest = work.drop(index=selected, errors="ignore").sort_values(
            ["_LEAD_DISTANCE", "P", "POS", "REFERENCE_ID"], kind="stable")
        selected.extend(rest.head(remaining).index.tolist())
    return work.loc[selected].sort_values(["P", "POS", "REFERENCE_ID"], kind="stable").reset_index(drop=True)


def calculate_ld(plink2: Path, prefix_ref: Path, ancestry_code: str, matched: pd.DataFrame, locus_dir: Path,
                 maf: float, geno: float, threads: int, max_variants: int,
                 lead_positions: list[int]) -> tuple[Path, Path, int]:
    all_ids = locus_dir / "matched_reference_ids_all.txt"
    matched["REFERENCE_ID"].to_csv(all_ids, index=False, header=False)
    qc_prefix = locus_dir / f"{ancestry_code}_reference_qc"
    # Reference MAF/GENO QC first, then the size cap.
    run_command([plink2, "--pfile", prefix_ref, "--extract", all_ids, "--maf", maf, "--geno", geno,
                 "--min-alleles", 2, "--max-alleles", 2, "--write-snplist", "--threads", threads,
                 "--out", qc_prefix], locus_dir / f"{ancestry_code}_reference_qc.command.log")
    snplist = Path(f"{qc_prefix}.snplist")
    if not snplist.exists():
        raise RuntimeError("PLINK2 reference-QC SNP list is missing")
    surviving = [x.strip() for x in snplist.read_text().splitlines() if x.strip()]
    order = {v: i for i, v in enumerate(surviving)}
    qc_matched = matched[matched["REFERENCE_ID"].astype(str).isin(order)].copy()
    qc_matched["_ORDER"] = qc_matched["REFERENCE_ID"].astype(str).map(order)
    qc_matched = qc_matched.sort_values(["_ORDER", "P"], kind="stable").drop(columns=["_ORDER"]).reset_index(drop=True)
    qc_matched.to_csv(locus_dir / "reference_qc_matched_variants.tsv.gz", sep="\t", index=False, compression="gzip")
    if len(qc_matched) < 2:
        raise RuntimeError(f"Only {len(qc_matched)} matched variants survive reference MAF/GENO QC")

    applied = len(qc_matched) > max_variants
    selected = _balanced_ld_subset(qc_matched, max_variants, lead_positions) if applied else qc_matched
    if applied:
        print(f"[WARNING] {len(qc_matched):,} variants after reference QC exceed the cap of "
              f"{max_variants:,}; keeping a balanced subset around the lead signals.")
    selected.to_csv(locus_dir / "ld_selected_variants.tsv.gz", sep="\t", index=False, compression="gzip")
    write_json(locus_dir / "ld_variant_selection_summary.json", {
        "N_MATCHED_BEFORE_REFERENCE_QC": int(len(matched)), "N_AFTER_REFERENCE_QC": int(len(qc_matched)),
        "MAX_LD_VARIANTS": int(max_variants), "SELECTION_APPLIED": bool(applied),
        "N_SELECTED_FOR_LD": int(len(selected)), "LEAD_POSITIONS": lead_positions,
        "SELECTION_METHOD": ("balanced_nearest_lead" if applied and lead_positions else
                             "association_rank_fallback" if applied else "all_reference_qc_variants"),
    })
    ids_file = locus_dir / "matched_reference_ids.txt"
    selected["REFERENCE_ID"].to_csv(ids_file, index=False, header=False)
    ld_prefix = locus_dir / f"{ancestry_code}_LD"
    run_command([plink2, "--pfile", prefix_ref, "--extract", ids_file, "--maf", maf, "--geno", geno,
                 "--min-alleles", 2, "--max-alleles", 2,
                 "--r-unphased", "square", "bin4", "ref-based", "yes-really",
                 "--threads", threads, "--out", ld_prefix], locus_dir / f"{ancestry_code}_LD.command.log")
    ld_file = locus_dir / f"{ancestry_code}_LD.unphased.vcor1.bin"
    vars_file = locus_dir / f"{ancestry_code}_LD.unphased.vcor1.bin.vars"
    if not ld_file.exists() or not vars_file.exists():
        raise RuntimeError("PLINK2 LD output is missing")
    p = len([x for x in vars_file.read_text().splitlines() if x.strip()])
    if p < 2 or ld_file.stat().st_size != p * p * 4:
        raise RuntimeError(f"LD matrix validation failed: p={p}, bytes={ld_file.stat().st_size}")
    if p > max_variants:
        raise RuntimeError(f"Internal error: LD matrix has {p:,} variants, above cap {max_variants:,}")
    return ld_file, vars_file, p


# =============================================================================
# PER-LOCUS RESUME
# =============================================================================

def _file_signature(path: Path) -> dict:
    stat = path.stat()
    return {"PATH": str(path), "SIZE": int(stat.st_size), "MTIME_NS": int(stat.st_mtime_ns)}


def _locus_identity(locus) -> dict:
    return {"LOCUS_ID": str(locus["LOCUS_ID"]), "CHR": int(locus["CHR"]),
            "LOCUS_START": int(locus["LOCUS_START"]), "LOCUS_END": int(locus["LOCUS_END"]),
            "N_INDEPENDENT_SIGNALS": int(locus["N_INDEPENDENT_SIGNALS"]),
            "LEAD_POSITIONS": str(locus.get("LEAD_POSITIONS", ""))}


def _valid_marker(locus_dir: Path, locus, signature: dict) -> dict | None:
    marker = locus_dir / LOCUS_COMPLETE_MARKER
    if not marker.exists() or (locus_dir / "LOCUS_FAILED.txt").exists():
        return None
    payload = read_json(marker)
    if payload.get("LOCUS") != _locus_identity(locus) or payload.get("RUN_SIGNATURE") != signature:
        return None
    if not all(nonempty(locus_dir / "susie" / name) for name in LOCUS_REQUIRED_OUTPUTS):
        return None
    return payload


def _read_locus_outputs(susie_dir: Path, accession, phenotype, ancestry_code, locus_id, chromosome,
                        start, end, n_signals, n_regional, n_matched):
    variants = pd.read_csv(susie_dir / "susie_variants.tsv", sep="\t", low_memory=False)
    for i, (name, value) in enumerate((("STUDY_ACCESSION", accession), ("PHENOTYPE", phenotype),
                                       ("ANCESTRY_CODE", ancestry_code), ("LOCUS_ID", locus_id))):
        variants.insert(i, name, value)
    cs = None
    cs_path = susie_dir / "credible_sets.tsv"
    if nonempty(cs_path):
        cs = pd.read_csv(cs_path, sep="\t", low_memory=False)
        if cs.empty:
            cs = None
        else:
            for i, (name, value) in enumerate((("STUDY_ACCESSION", accession), ("PHENOTYPE", phenotype),
                                               ("ANCESTRY_CODE", ancestry_code))):
                cs.insert(i, name, value)
    summary = pd.read_csv(susie_dir / "susie_summary.tsv", sep="\t", low_memory=False)
    summary.insert(0, "STUDY_ACCESSION", accession)
    summary["CHR"], summary["LOCUS_START"], summary["LOCUS_END"] = chromosome, start, end
    summary["N_INDEPENDENT_SIGNALS"] = n_signals
    summary["N_REGIONAL_GWAS_VARIANTS"] = n_regional
    summary["N_MATCHED_VARIANTS"] = n_matched
    diagnostics = susie_dir / "susie_diagnostics.tsv"
    if nonempty(diagnostics):  # LD-mismatch diagnostic (estimate_s_rss) for reporting
        diag = pd.read_csv(diagnostics, sep="\t", low_memory=False)
        for column in ("ESTIMATE_S_RSS", "ESTIMATE_S_STATUS", "MAX_ABS_Z"):
            if column in diag.columns and len(diag):
                summary[column] = diag[column].iloc[0]
    return variants, cs, summary


# =============================================================================
# STUDY WORKER
# =============================================================================

def finemap_study(project: Project, study: dict, force: bool) -> tuple:
    accession = study["STUDY_ACCESSION"]
    settings = project.cfg["finemapping"]
    out = project.finemap_dir / accession
    final_variants = out / f"{accession}_finemapped_variants.tsv.gz"
    final_cs = out / f"{accession}_95pct_credible_sets.tsv"
    final_summary = out / f"{accession}_finemapping_summary.tsv"
    failed_file = out / f"{accession}_failed_loci.tsv"
    summary_json = out / "finemapping_summary.json"
    if not force and final_variants.exists() and final_summary.exists() and summary_json.exists():
        previous = read_json(summary_json)
        n_loci = int(previous.get("N_LOCI", 0) or 0)
        n_failed = int(previous.get("N_LOCI_FAILED", 0) or 0)
        n_success = int(previous.get("N_LOCI_SUCCESS", 0) or 0)
        # Partial/failed studies are re-run; their finished loci are reused from markers.
        if n_failed == 0 and (n_success > 0 or n_loci == 0):
            return "SKIPPED_ALREADY_COMPLETE", previous.get("MESSAGE", ""), _extra(previous)

    lead_file = project.clump_dir / accession / "lead_variants.tsv"
    qc_file = project.qc_dir / accession / f"{accession}_GRCh38_QC.tsv.gz"
    for path in (lead_file, qc_file):
        if not nonempty(path):
            return "FAILED", f"Missing input: {path}", {}

    meta = resolve_sample_size(project, study)
    if not meta["RESOLVED_GWAS_N"]:
        return "FAILED", "GWAS sample size could not be resolved (set finemapping.sample_size)", {}
    gwas_n, gwas_n_source = int(meta["RESOLVED_GWAS_N"]), meta["RESOLVED_GWAS_N_SOURCE"]

    resources = Resources(project.cfg)
    plink2 = find_tool(project.cfg, "plink2")
    rscript = find_tool(project.cfg, "rscript")
    require_susier(rscript)
    reference_n = resources.count_reference_samples(project.ancestry_code)
    threads = project.cfg.threads
    window_kb = int(settings["locus_window_kb"])
    maf, geno = float(settings["reference_maf"]), float(settings["reference_geno"])
    min_variants, max_ld = int(settings["min_variants"]), int(settings["max_ld_variants"])
    coverage, ld_ridge = float(settings["coverage"]), float(settings["ld_ridge"])
    susie_l, max_susie_l = int(settings["susie_L"]), int(settings["max_susie_L"])
    estimate_s_max = int(settings["estimate_s_max_variants"])
    out.mkdir(parents=True, exist_ok=True)

    banner(f"STEP 6 - SUSIE-RSS FINE-MAPPING | {accession}")
    print(f"GWAS N {gwas_n} [{gwas_n_source}] | reference N {reference_n} | window {window_kb} kb")
    loci = load_leads_and_build_loci(lead_file, window_kb)
    loci.to_csv(out / "loci_definition.tsv", sep="\t", index=False)

    base = {"STUDY_ACCESSION": accession, "PHENOTYPE": project.phenotype,
            "ANCESTRY_CODE": project.ancestry_code, "ANCESTRY_LABEL": project.ancestry_label,
            "GWAS_N": gwas_n, "GWAS_N_SOURCE": gwas_n_source, "REFERENCE_N": reference_n}
    if loci.empty:
        # Zero independent genome-wide significant leads is a valid outcome, not a failure.
        pd.DataFrame(columns=["STUDY_ACCESSION", "PHENOTYPE", "ANCESTRY_CODE", "LOCUS_ID",
                              "REFERENCE_ID", "CHR", "POS", "PIP"]).to_csv(final_variants, sep="\t", index=False,
                                                                         compression="gzip")
        pd.DataFrame(columns=["STUDY_ACCESSION", "LOCUS_ID", "CREDIBLE_SET", "REFERENCE_ID", "CHR", "POS",
                              "PIP"]).to_csv(final_cs, sep="\t", index=False)
        row = {**base, "STATUS": "NO_SIGNIFICANT_LOCI", "N_LOCI": 0, "N_LOCI_SUCCESS": 0, "N_LOCI_FAILED": 0,
               "MESSAGE": "Step 5 produced zero independent lead variants; nothing to fine-map."}
        pd.DataFrame([row]).to_csv(final_summary, sep="\t", index=False)
        pd.DataFrame(columns=["STUDY_ACCESSION", "LOCUS_ID", "CHR", "LOCUS_START", "LOCUS_END",
                              "ERROR"]).to_csv(failed_file, sep="\t", index=False)
        write_json(summary_json, {**row, "COMPLETED_UTC": utc_now()})
        return "COMPLETE", row["MESSAGE"], _extra(row)

    locus_root = out / "loci"
    locus_root.mkdir(parents=True, exist_ok=True)
    signature = {
        "QC_FILE": _file_signature(qc_file), "LEAD_VARIANTS_FILE": _file_signature(lead_file),
        "GWAS_N": gwas_n, "REFERENCE_N": reference_n, "LOCUS_WINDOW_KB": window_kb,
        "REFERENCE_MAF": maf, "REFERENCE_GENO": geno, "MIN_VARIANTS": min_variants,
        "MAX_LD_VARIANTS": max_ld, "COVERAGE": coverage, "SUSIE_L": susie_l, "MAX_SUSIE_L": max_susie_l,
        "LD_RIDGE": ld_ridge, "ESTIMATE_S_MAX": estimate_s_max,
    }
    completed = {}
    if not force:
        for _, locus in loci.iterrows():
            marker = _valid_marker(locus_root / str(locus["LOCUS_ID"]), locus, signature)
            if marker is not None:
                completed[str(locus["LOCUS_ID"])] = marker
    pending = loci[~loci["LOCUS_ID"].astype(str).isin(completed)]
    print(f"Loci: {len(loci)} total, {len(completed)} reused, {len(pending)} to run")
    regional_paths = extract_regional_gwas_once(qc_file, pending, out / "regional_gwas") if not pending.empty else {}

    fine_frames, cs_frames, summary_frames, failures = [], [], [], []
    for _, locus in loci.iterrows():
        locus_id = str(locus["LOCUS_ID"])
        chromosome, start, end = int(locus["CHR"]), int(locus["LOCUS_START"]), int(locus["LOCUS_END"])
        n_signals = int(locus["N_INDEPENDENT_SIGNALS"])
        locus_dir = locus_root / locus_id
        locus_dir.mkdir(parents=True, exist_ok=True)
        prefix_ref = resources.ld_prefix(project.ancestry_code, chromosome)
        if locus_id in completed:
            marker = completed[locus_id]
            v, cs, s = _read_locus_outputs(locus_dir / "susie", accession, project.phenotype,
                                           project.ancestry_code, locus_id, chromosome, start, end, n_signals,
                                           marker["N_REGIONAL_GWAS_VARIANTS"], marker["N_MATCHED_VARIANTS"])
        else:
            (locus_dir / LOCUS_COMPLETE_MARKER).unlink(missing_ok=True)
            try:
                regional = pd.read_csv(regional_paths[locus_id], sep="\t", compression="infer", low_memory=False)
                if len(regional) < min_variants:
                    raise RuntimeError(f"Only {len(regional)} regional GWAS variants; minimum is {min_variants}")
                reference = read_pvar(create_regional_pvar(plink2, prefix_ref, chromosome, start, end, locus_dir))
                matched = match_gwas_to_reference(regional, reference)
                if len(matched) < min_variants:
                    raise RuntimeError(f"Only {len(matched)} GWAS/reference matched variants; minimum is {min_variants}")
                matched_file = locus_dir / "matched_variants.tsv"
                matched.to_csv(matched_file, sep="\t", index=False)
                ld_file, vars_file, p = calculate_ld(
                    plink2, prefix_ref, project.ancestry_code, matched, locus_dir, maf, geno, threads, max_ld,
                    _parse_lead_positions(locus.get("LEAD_POSITIONS", "")))
                if p < min_variants:
                    raise RuntimeError(f"Only {p} variants remain in LD matrix; minimum is {min_variants}")
                susie_dir = locus_dir / "susie"
                L = min(max(susie_l, n_signals), max_susie_l, p)
                run_command([rscript, r_script_path(), matched_file, ld_file, vars_file, susie_dir, gwas_n,
                             reference_n, L, coverage, ld_ridge, estimate_s_max, locus_id],
                            susie_dir / "r_command.log")
                v, cs, s = _read_locus_outputs(susie_dir, accession, project.phenotype, project.ancestry_code,
                                               locus_id, chromosome, start, end, n_signals, len(regional),
                                               len(matched))
                (locus_dir / "LOCUS_FAILED.txt").unlink(missing_ok=True)
                write_json(locus_dir / LOCUS_COMPLETE_MARKER, {
                    "LOCUS": _locus_identity(locus), "RUN_SIGNATURE": signature,
                    "N_REGIONAL_GWAS_VARIANTS": int(len(regional)), "N_MATCHED_VARIANTS": int(len(matched)),
                    "COMPLETED_UTC": utc_now()})
            except Exception as exc:
                text = f"{type(exc).__name__}: {exc}"
                failures.append({"STUDY_ACCESSION": accession, "LOCUS_ID": locus_id, "CHR": chromosome,
                                 "LOCUS_START": start, "LOCUS_END": end, "ERROR": text[:1000]})
                (locus_dir / "LOCUS_FAILED.txt").write_text(text + "\n\n" + traceback.format_exc(), encoding="utf-8")
                print(f"[FAILED] {locus_id}: {text}")
                continue
        fine_frames.append(v)
        if cs is not None:
            cs_frames.append(cs)
        summary_frames.append(s)

    (pd.concat(fine_frames, ignore_index=True) if fine_frames else pd.DataFrame()).to_csv(
        final_variants, sep="\t", index=False, compression="gzip")
    (pd.concat(cs_frames, ignore_index=True) if cs_frames else pd.DataFrame(
        columns=["STUDY_ACCESSION", "LOCUS_ID", "CREDIBLE_SET", "REFERENCE_ID", "CHR", "POS", "PIP"])).to_csv(
        final_cs, sep="\t", index=False)
    (pd.concat(summary_frames, ignore_index=True) if summary_frames else pd.DataFrame()).to_csv(
        final_summary, sep="\t", index=False)
    pd.DataFrame(failures, columns=["STUDY_ACCESSION", "LOCUS_ID", "CHR", "LOCUS_START", "LOCUS_END",
                                    "ERROR"]).to_csv(failed_file, sep="\t", index=False)
    message = (f"{len(failures)} of {len(loci)} loci failed; see {failed_file.name}" if failures else "")
    summary = {**base, "STATUS": "PARTIAL" if failures else "COMPLETE", "N_LOCI": int(len(loci)),
               "N_LOCI_SUCCESS": len(summary_frames), "N_LOCI_FAILED": len(failures), "MESSAGE": message,
               "FINEMAPPED_VARIANTS_FILE": str(final_variants), "CREDIBLE_SETS_FILE": str(final_cs),
               "SUMMARY_FILE": str(final_summary), "FAILED_LOCI_FILE": str(failed_file),
               "COMPLETED_UTC": utc_now()}
    write_json(summary_json, summary)
    banner("FINE-MAPPING FINISHED")
    print(f"Loci: {len(loci)} | success {len(summary_frames)} | failed {len(failures)}")
    if not summary_frames:
        return "FAILED", f"All {len(loci)} loci failed fine-mapping", _extra(summary)
    return ("PARTIAL" if failures else "COMPLETE"), message, _extra(summary)


def _extra(summary: dict) -> dict:
    return {k: summary.get(k) for k in ("N_LOCI", "N_LOCI_SUCCESS", "N_LOCI_FAILED", "GWAS_N",
                                       "GWAS_N_SOURCE", "STATUS")}
