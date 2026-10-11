"""Step 5 - ancestry-matched LD clumping of one study (PLINK2).

Streams the cleaned GWAS once keeping P <= p2, maps variants to the 1000G panel
by position + alleles, runs `plink2 --clump` per chromosome and writes
05_clumping/<GCST>/lead_variants.tsv.

Fixes relative to the original Step05 script:
  * a chromosome whose candidates all fail to map no longer crashes (KeyError);
  * a variant whose alleles disagree with the reference is no longer mapped to
    whatever variant sits at the same position (position-only matching is used
    only when the GWAS lacks allele columns).
"""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path

import pandas as pd

from gwas2m_lncrna.common import (
    banner, clean_allele, find_column, nonempty, run_command, utc_now, write_json,
)
from gwas2m_lncrna.project import Project
from gwas2m_lncrna.resources import Resources
from gwas2m_lncrna.tools import find_tool

MAPPED_COLUMNS = ["ID", "P", "CHR", "POS", "REF", "ALT", "EA", "NEA", "SOURCE_ID"]
CHUNK_SIZE = 500_000


def clean_chr_series(series: pd.Series) -> pd.Series:
    text = series.astype(str).str.strip().str.replace(r"^chr", "", regex=True, case=False)
    return pd.to_numeric(text, errors="coerce").astype("Int64")


def inspect_gwas_columns(path: Path) -> dict:
    columns = list(pd.read_csv(path, sep="\t", compression="infer", nrows=0).columns)
    mapping = {
        "CHR": find_column(columns, ["CHR"]), "POS": find_column(columns, ["POS"]),
        "P": find_column(columns, ["P"]), "EA": find_column(columns, ["EA"]),
        "NEA": find_column(columns, ["NEA"]), "SNPID": find_column(columns, ["SNPID", "rsID", "RSID"]),
    }
    missing = [x for x in ("CHR", "POS", "P") if mapping[x] is None]
    if missing:
        raise RuntimeError(f"GWAS file is missing required columns {missing}: {path}")
    return mapping


def stream_candidate_variants(qc_file: Path, p2: float, chunk_size: int = CHUNK_SIZE):
    mapping = inspect_gwas_columns(qc_file)
    usecols = [mapping["CHR"], mapping["POS"], mapping["P"]]
    usecols += [mapping[k] for k in ("EA", "NEA", "SNPID") if mapping[k] is not None]
    usecols = list(dict.fromkeys(usecols))
    per_chrom = defaultdict(list)
    n_total = n_valid = n_p2 = 0
    for chunk in pd.read_csv(qc_file, sep="\t", compression="infer", usecols=usecols,
                             chunksize=chunk_size, low_memory=False):
        n_total += len(chunk)
        frame = pd.DataFrame(index=chunk.index)
        frame["CHR"] = clean_chr_series(chunk[mapping["CHR"]])
        frame["POS"] = pd.to_numeric(chunk[mapping["POS"]], errors="coerce").astype("Int64")
        frame["P"] = pd.to_numeric(chunk[mapping["P"]], errors="coerce")
        frame["EA"] = chunk[mapping["EA"]].map(clean_allele) if mapping["EA"] else None
        frame["NEA"] = chunk[mapping["NEA"]].map(clean_allele) if mapping["NEA"] else None
        frame["SOURCE_ID"] = chunk[mapping["SNPID"]].astype(str) if mapping["SNPID"] else ""
        valid = (frame["CHR"].between(1, 22) & frame["POS"].notna()
                 & frame["P"].notna() & frame["P"].between(0, 1))
        n_valid += int(valid.sum())
        selected = frame.loc[valid & (frame["P"] <= p2)]
        n_p2 += len(selected)
        for chrom, part in selected.groupby("CHR", sort=False):
            per_chrom[int(chrom)].append(part)
    combined = {}
    for chrom, frames in per_chrom.items():
        df = pd.concat(frames, ignore_index=True).sort_values("P", kind="stable")
        combined[chrom] = df.drop_duplicates(subset=["POS", "P", "EA", "NEA"], keep="first")
    return combined, {"N_GWAS_ROWS": n_total, "N_VALID_ROWS": n_valid, "N_P_LE_P2": n_p2}


def load_reference_matches(pvar_path: Path, target_positions: set[int]) -> dict:
    """position -> [{ID, REF, ALT}] for the requested positions of a .pvar."""
    result = defaultdict(list)
    with open(pvar_path, "rt", encoding="utf-8", errors="replace") as handle:
        header = None
        for line in handle:
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                header = line.rstrip("\n").lstrip("#").split("\t")
                break
        if header is None:
            raise RuntimeError(f"Could not find PVAR header in {pvar_path}")
        index = {name: i for i, name in enumerate(header)}
        missing = [x for x in ("POS", "ID", "REF", "ALT") if x not in index]
        if missing:
            raise RuntimeError(f"PVAR missing columns {missing}: {pvar_path}")
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            try:
                pos = int(fields[index["POS"]])
                if pos in target_positions:
                    result[pos].append({"ID": fields[index["ID"]], "REF": fields[index["REF"]].upper(),
                                        "ALT": fields[index["ALT"]].upper()})
            except (ValueError, IndexError):
                continue
    return result


def map_gwas_to_reference(gwas: pd.DataFrame, pvar_path: Path) -> tuple[pd.DataFrame, dict]:
    stats = {"N_INPUT": int(len(gwas)), "N_MAPPED": 0, "N_UNMAPPED": 0,
             "N_ALLELE_MISMATCH": 0, "N_AMBIGUOUS": 0}
    if gwas.empty:
        return pd.DataFrame(columns=MAPPED_COLUMNS), stats
    reference = load_reference_matches(pvar_path, set(int(x) for x in gwas["POS"].dropna().astype(int).unique()))
    rows = []
    for row in gwas.itertuples(index=False):
        pos = int(row.POS)
        options = reference.get(pos, [])
        if not options:
            stats["N_UNMAPPED"] += 1
            continue
        ea, nea = clean_allele(row.EA), clean_allele(row.NEA)
        if ea and nea:
            matches = [o for o in options if {o["REF"], o["ALT"]} == {ea, nea}]
            if not matches:
                stats["N_ALLELE_MISMATCH"] += 1  # a different variant at this position
                continue
        else:
            matches = options  # no alleles in the GWAS: position must be unique
        if len(matches) != 1:
            stats["N_AMBIGUOUS"] += 1
            continue
        chosen = matches[0]
        rows.append({"ID": chosen["ID"], "P": float(row.P), "CHR": int(row.CHR), "POS": pos,
                     "REF": chosen["REF"], "ALT": chosen["ALT"], "EA": ea, "NEA": nea,
                     "SOURCE_ID": getattr(row, "SOURCE_ID", "")})
    mapped = pd.DataFrame(rows, columns=MAPPED_COLUMNS)
    if not mapped.empty:
        mapped = (mapped.sort_values("P", kind="stable")
                  .drop_duplicates(subset=["ID"], keep="first").reset_index(drop=True))
    stats["N_MAPPED"] = int(len(mapped))
    return mapped, stats


def qc_file_for(project: Project, accession: str) -> Path:
    return project.qc_dir / accession / f"{accession}_GRCh38_QC.tsv.gz"


def clump_study(project: Project, study: dict, force: bool) -> tuple:
    accession = study["STUDY_ACCESSION"]
    qc_record = project.status(accession).get("qc", {})
    if str(qc_record.get("CLUMPING_READY", "True")).strip().lower() not in {"true", "1", "yes"}:
        return "EXCLUDED", f"QC status {qc_record.get('QC_STATUS', '?')}: not clumping-ready", {}

    settings = project.cfg["clumping"]
    p1, p2, r2, kb = float(settings["p1"]), float(settings["p2"]), float(settings["r2"]), int(settings["kb"])
    out = project.clump_dir / accession
    lead_output = out / "lead_variants.tsv"
    mapped_output = out / "mapped_candidates.tsv.gz"
    summary_tsv = out / "clumping_summary.tsv"
    summary_json = out / "clumping_summary.json"
    if not force and all(p.exists() for p in (lead_output, mapped_output, summary_tsv, summary_json)):
        n = len(pd.read_csv(lead_output, sep="\t"))
        return "SKIPPED_ALREADY_COMPLETE", "", {"N_LEAD_VARIANTS": n}

    qc_file = qc_file_for(project, accession)
    if not nonempty(qc_file):
        return "FAILED", f"QC file not found: {qc_file}", {}
    resources = Resources(project.cfg)
    plink2 = find_tool(project.cfg, "plink2")
    threads = project.cfg.threads
    out.mkdir(parents=True, exist_ok=True)

    banner(f"STEP 5 - LD CLUMPING | {accession}")
    print(f"p1={p1} p2={p2} r2={r2} kb={kb} threads={threads}")
    per_chrom, stream_stats = stream_candidate_variants(qc_file, p2)
    print(f"GWAS rows {stream_stats['N_GWAS_ROWS']:,} | valid {stream_stats['N_VALID_ROWS']:,} | "
          f"P<=p2 {stream_stats['N_P_LE_P2']:,}")

    all_mapped, all_leads, chrom_rows = [], [], []
    for chrom in range(1, 23):
        chrom_dir = out / f"chr{chrom:02d}"
        chrom_dir.mkdir(parents=True, exist_ok=True)
        prefix = resources.ld_prefix(project.ancestry_code, chrom)
        gwas_chr = per_chrom.get(chrom, pd.DataFrame(columns=["CHR", "POS", "P", "EA", "NEA", "SOURCE_ID"]))
        mapped, map_stats = map_gwas_to_reference(gwas_chr, Path(f"{prefix}.pvar"))
        if not mapped.empty:
            all_mapped.append(mapped)
        report_prefix = chrom_dir / "plink"
        for suffix in (".clumps", ".clumps.missing_id", ".clumps.missing_allele"):
            Path(f"{report_prefix}{suffix}").unlink(missing_ok=True)
        clump_input = chrom_dir / "clump_input.tsv"
        mapped[["ID", "P"]].to_csv(clump_input, sep="\t", index=False)
        n_p1 = int((mapped["P"] <= p1).sum()) if not mapped.empty else 0
        if n_p1 > 0:
            run_command([plink2, "--pfile", prefix, "--clump", clump_input, "--clump-id-field", "ID",
                         "--clump-p-field", "P", "--clump-p1", p1, "--clump-p2", p2, "--clump-r2", r2,
                         "--clump-kb", kb, "--threads", threads, "--out", report_prefix])
        clumps_file = Path(f"{report_prefix}.clumps")
        n_leads = 0
        if nonempty(clumps_file):
            try:
                clumps = pd.read_csv(clumps_file, sep=r"\s+", engine="python")
            except pd.errors.EmptyDataError:
                clumps = pd.DataFrame()
            if not clumps.empty:
                for i, (name, value) in enumerate((
                    ("STUDY_ACCESSION", accession), ("PHENOTYPE", project.phenotype),
                    ("ANCESTRY_CODE", project.ancestry_code), ("ANCESTRY_LABEL", project.ancestry_label),
                    ("SOURCE_CHROMOSOME", chrom),
                )):
                    clumps.insert(i, name, value)
                all_leads.append(clumps)
                n_leads = len(clumps)
        chrom_rows.append({"CHR": chrom, "N_GWAS_P_LE_P2": map_stats["N_INPUT"],
                           "N_MAPPED": map_stats["N_MAPPED"], "N_UNMAPPED": map_stats["N_UNMAPPED"],
                           "N_ALLELE_MISMATCH": map_stats["N_ALLELE_MISMATCH"],
                           "N_AMBIGUOUS": map_stats["N_AMBIGUOUS"], "N_P_LE_P1_MAPPED": n_p1,
                           "N_LEAD_VARIANTS": n_leads})

    mapped_all = (pd.concat(all_mapped, ignore_index=True).sort_values(["CHR", "POS", "P"], kind="stable")
                  if all_mapped else pd.DataFrame(columns=MAPPED_COLUMNS))
    mapped_all.to_csv(mapped_output, sep="\t", index=False, compression="gzip")
    if all_leads:
        leads = pd.concat(all_leads, ignore_index=True)
        if "P" in leads.columns:
            leads = leads.sort_values("P", kind="stable")
    else:
        leads = pd.DataFrame(columns=["STUDY_ACCESSION", "PHENOTYPE", "ANCESTRY_CODE", "ANCESTRY_LABEL",
                                      "SOURCE_CHROMOSOME", "CHROM", "POS", "ID", "P"])
    leads.to_csv(lead_output, sep="\t", index=False)
    chrom_df = pd.DataFrame(chrom_rows)
    chrom_df.to_csv(out / "chromosome_clumping_summary.tsv", sep="\t", index=False)

    n_input = int(chrom_df["N_GWAS_P_LE_P2"].sum())
    summary = {
        "STUDY_ACCESSION": accession, "PHENOTYPE": project.phenotype,
        "ANCESTRY_CODE": project.ancestry_code, "BUILD": "GRCh38", "QC_FILE": str(qc_file),
        "REFERENCE_ROOT": str(resources.ld_dir / project.ancestry_code),
        "P1": p1, "P2": p2, "R2": r2, "KB": kb, **stream_stats,
        "N_REFERENCE_MAPPED": int(chrom_df["N_MAPPED"].sum()),
        "N_ALLELE_MISMATCH": int(chrom_df["N_ALLELE_MISMATCH"].sum()),
        "MAPPING_RATE": (int(chrom_df["N_MAPPED"].sum()) / n_input) if n_input else None,
        "N_INDEX_CANDIDATES_P_LE_P1": int(chrom_df["N_P_LE_P1_MAPPED"].sum()),
        "N_LEAD_VARIANTS": int(len(leads)), "THREADS": threads, "COMPLETED_UTC": utc_now(),
    }
    pd.DataFrame([summary]).to_csv(summary_tsv, sep="\t", index=False)
    write_json(summary_json, summary)
    banner("LD CLUMPING COMPLETE")
    print(f"Lead variants: {len(leads):,} | mapping rate: {summary['MAPPING_RATE']}")
    return "COMPLETE", "", {"N_LEAD_VARIANTS": int(len(leads))}
