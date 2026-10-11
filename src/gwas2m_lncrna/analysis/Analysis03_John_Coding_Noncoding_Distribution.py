#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GWAS2m - John Analysis 03 (complete replacement, v3.0)
STRICT protein-coding CDS versus non-CDS plus transcript-aware genomic landscape.

READS ONLY existing GWAS2m pipeline products. NEVER edits Step06/Step07 inputs.
run() returns the auditable tables (per-variant classes, study summary, primary
regions, feature flags, ...); `gwas2m-lncrna analysis` saves them in the report.
Streams Step07 all-consequence .tsv.gz once per study.

Definition: protein-coding CDS requires BOTH a protein_coding BIOTYPE and a
CDS-defining VEP Consequence on the SAME transcript row. A variant which is
protein-coding-intronic/UTR or lncRNA-intronic is not a coding-CDS variant.
A variant may overlap multiple transcript classes; the PRIMARY_REGION is an
explicitly prioritised single category, whereas other feature flags are
nonexclusive. One VEP_ID is counted once per GWAS accession (not per isoform).

Default: CS95 and PIP >= .5 only (fast; no giant regional-variant dictionaries).
Use --include-all-regional for the larger optional sensitivity analysis.

Examples:
  python -m gwas2m_lncrna.analysis.Analysis03_John_Coding_Noncoding_Distribution -p "parkinson's disease" -a EUR
  python -m gwas2m_lncrna.analysis.Analysis03_John_Coding_Noncoding_Distribution -p CAD -p "parkinson's disease"
  python -m gwas2m_lncrna.analysis.Analysis03_John_Coding_Noncoding_Distribution -p CAD --include-all-regional
  python -m gwas2m_lncrna.analysis.Analysis03_John_Coding_Noncoding_Distribution --self-test
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd

from gwas2m_lncrna.analysis.helpers import analysis_main, direct_labels, normalize_label

VERSION = "3.0.0-strict-transcript-cds"
# Set by run(): {phenotype folder: (display name, phenotype folder)} and the ancestry folder name.
DISEASES = {}
ANCESTRY = "european"
DONE = {"COMPLETE", "SKIPPED_ALREADY_COMPLETE"}
CDS = {
    "coding_sequence_variant", "missense_variant", "synonymous_variant",
    "stop_gained", "stop_lost", "start_lost", "start_retained_variant",
    "stop_retained_variant", "frameshift_variant", "inframe_insertion",
    "inframe_deletion", "protein_altering_variant",
    "incomplete_terminal_codon_variant",
}
SPLICE_CORE = {"splice_acceptor_variant", "splice_donor_variant"}
SPLICE_REGION = {
    "splice_region_variant", "splice_donor_region_variant",
    "splice_donor_5th_base_variant", "splice_polypyrimidine_tract_variant",
}
UTR5 = {"5_prime_utr_variant"}
UTR3 = {"3_prime_utr_variant"}
EXON = {"non_coding_transcript_exon_variant"}
INTRON = {"intron_variant"}
OTHER_TX = {"non_coding_transcript_variant", "transcript_ablation", "transcript_amplification"}
UPSTREAM = {"upstream_gene_variant"}
DOWNSTREAM = {"downstream_gene_variant"}
INTERGENIC = {"intergenic_variant"}
OTHER_TRANSCRIPT_TERMS = {"mature_mirna_variant", "nmd_transcript_variant", "non_coding_transcript_exon_variant"}
REGULATORY = {
    "regulatory_region_variant", "regulatory_region_ablation",
    "regulatory_region_amplification", "tf_binding_site_variant",
    "tfbs_ablation", "tfbs_amplification", "feature_elongation", "feature_truncation",
}
PSEUDO_BIOTYPES = {"pseudogene", "processed_pseudogene", "unprocessed_pseudogene",
                   "transcribed_processed_pseudogene", "transcribed_unprocessed_pseudogene",
                   "unitary_pseudogene", "polymorphic_pseudogene", "translated_processed_pseudogene",
                   "translated_unprocessed_pseudogene"}
OTHER_NCRNA_BIOTYPES = {
    "mirna", "snorna", "snrna", "rrna", "scrna", "srprna", "misc_rna", "sca_rna",
    "scarna", "trna", "ribozyme", "vault_rna", "yrna", "mt_trna", "mt_rrna",
    "tec", "non_coding", "macro_lncrna", "antisense", "processed_transcript",
    "sense_intronic", "sense_overlapping", "3prime_overlapping_ncrna",
    "bidirectional_promoter_lncrna",
}
FLAGS = (
    "PC_CDS", "PC_SPLICE_CORE", "PC_SPLICE_REGION", "PC_5UTR", "PC_3UTR",
    "PC_EXON_NONCDS", "PC_INTRON", "PC_OTHER_BODY", "PC_GENE_BODY",
    "LNCRNA_BODY", "LNCRNA_EXON", "LNCRNA_INTRON", "LNCRNA_SPLICE_CORE",
    "LNCRNA_SPLICE_REGION", "LNCRNA_OTHER_BODY", "OTHER_NCRNA_BODY",
    "PSEUDOGENE_BODY", "OTHER_GENE_BODY", "UPSTREAM", "DOWNSTREAM",
    "REGULATORY", "INTERGENIC",
)
PRIMARY_ORDER = (
    ("PROTEIN_CODING_CDS", "PC_CDS"),
    ("SPLICE_CORE", "ANY_SPLICE_CORE"),
    ("SPLICE_REGION", "ANY_SPLICE_REGION"),
    ("FIVE_PRIME_UTR", "PC_5UTR"),
    ("THREE_PRIME_UTR", "PC_3UTR"),
    ("LNCRNA_EXON", "LNCRNA_EXON"),
    ("OTHER_NCRNA_EXON", "OTHER_NCRNA_EXON"),
    ("PROTEIN_CODING_INTRON", "PC_INTRON"),
    ("LNCRNA_INTRON", "LNCRNA_INTRON"),
    ("PSEUDOGENE_BODY", "PSEUDOGENE_BODY"),
    ("OTHER_NCRNA_BODY", "OTHER_NCRNA_BODY"),
    ("PROTEIN_CODING_NONCDS_BODY", "PC_GENE_BODY"),
    ("OTHER_GENE_BODY", "OTHER_GENE_BODY"),
    ("REGULATORY", "REGULATORY"),
    ("UPSTREAM", "UPSTREAM"),
    ("DOWNSTREAM", "DOWNSTREAM"),
    ("INTERGENIC", "INTERGENIC"),
    ("OTHER_ANNOTATED", "ANY_ANNOTATED"),
)

try:
    csv.field_size_limit(sys.maxsize)
except OverflowError:
    csv.field_size_limit(2**31 - 1)


def read_json(path: Path) -> dict:
    try:
        with path.open("r", encoding="utf-8") as fh:
            x = json.load(fh)
        return x if isinstance(x, dict) else {}
    except (OSError, ValueError):
        return {}


def read_tsv(path: Path):
    with path.open("r", encoding="utf-8-sig", newline="") as fh:
        yield from csv.DictReader(fh, delimiter="\t")


def write_tsv(path: Path, records: list[dict], fields: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fields is None:
        fields = list(dict.fromkeys(k for row in records for k in row))
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(records)


def flag(x) -> bool:
    return str(x).strip().lower() in {"1", "true", "t", "yes", "y"}


def pip_float(x):
    try:
        z = float(x)
        return z if math.isfinite(z) and 0 <= z <= 1 else None
    except (ValueError, TypeError):
        return None


def terms(value) -> set[str]:
    return {t.strip().lower() for t in re.split(r"[&;,|]", str(value))
            if t.strip() and t.strip().lower() not in {"na", "nan", "-", "."}}


def transcript_number(value) -> bool:
    return bool(re.match(r"^[1-9]\d*(?:/[1-9]\d*)?$", str(value).strip()))


def normalize_allele(row) -> str:
    chrom = re.sub(r"^chr", "", str(row.get("CHR", "")).strip(), flags=re.I)
    pos = str(row.get("REFERENCE_POS", "")).strip()
    ref = str(row.get("REFERENCE_REF", "")).strip().upper()
    alt = str(row.get("REFERENCE_ALT", "")).strip().upper()
    if not chrom.isdigit() or not pos.isdigit() or not ref or not alt:
        return ""
    return f"{int(chrom)}:{int(pos)}:{ref}:{alt}"


def trait_label_type(trait: str, disease: str) -> str:
    s = " ".join(str(trait or "").lower().replace("’", "'").split())
    if any(t in s for t in ("family history", "parental", "proxy", "father", "mother", "sibling")):
        return "PROXY_FAMILY"
    if normalize_label(s) in direct_labels(DISEASES[disease][0]):
        return "DIRECT_LABEL"
    if any(t in s for t in ("mtag", "interaction", "subtype", "age at onset", "medication",
                            "drug response", "cognitive", "severity", "risk score")):
        return "RELATED_SPECIALIZED"
    return "REVIEW_TRAIT"


def study_manifest(root: Path, disease: str):
    dir_ = root / "01_gwas_catalog" / disease / ANCESTRY
    choices = [dir_ / "GWAS_download_manifest.tsv",
               root / "pipeline_audit" / disease / ANCESTRY / "study_manifest.tsv"]
    path = next((p for p in choices if p.is_file()), None)
    rows = list(read_tsv(path)) if path else []
    accession_rows = {}
    for row in rows:
        study = next((row.get(k, "").strip() for k in
                      ("STUDY_ACCESSION", "GWAS_CATALOG_ACCESSION", "ACCESSION", "GCST")
                      if row.get(k, "").strip()), "")
        if study and study not in accession_rows:
            accession_rows[study] = row
    # Discover any Step06 study that was added after an old discovery manifest.
    fine_dir = root / "06_finemapping" / disease / ANCESTRY
    if fine_dir.is_dir():
        for study_dir in fine_dir.iterdir():
            if study_dir.is_dir() and re.fullmatch(r"GCST\d+", study_dir.name):
                accession_rows.setdefault(study_dir.name, {})
    return accession_rows, path


def get_trait(row: dict) -> str:
    for key in ("TRAIT", "DISEASE_TRAIT", "DISEASE_TRAITS", "EFO_TRAIT",
                "PHENOTYPE_LABEL", "DISEASE", "TRAIT_LABEL", "GWAS_TRAIT"):
        v = str(row.get(key, "")).strip()
        if v and v.lower() not in {"nan", "none", "na"}:
            return v
    return ""


def study_readiness(root: Path, disease: str, study: str):
    audit = root / "pipeline_audit" / disease / ANCESTRY / "status" / study
    s06 = str(read_json(audit / "Step06_FineMapping.json").get("status", "NOT_RECORDED")).upper()
    s07 = str(read_json(audit / "Step07_VEP.json").get("status", "NOT_RECORDED")).upper()
    fm_dir = root / "06_finemapping" / disease / ANCESTRY / study
    ann_dir = root / "07_annotation" / disease / ANCESTRY / study
    fine = fm_dir / f"{study}_finemapped_variants.tsv.gz"
    meta_file = ann_dir / "annotation_summary.json"
    meta = read_json(meta_file)
    fm = read_json(fm_dir / "finemapping_summary.json")
    if str(fm.get("STATUS", "")).upper() == "NO_SIGNIFICANT_LOCI" or (fm and fm.get("N_LOCI") == 0):
        return None, "NO_SIGNIFICANT_LOCI"
    if s06 not in DONE:
        return None, "S06_" + s06
    if s07 not in DONE:
        return None, "S07_" + s07
    summary = ann_dir / f"{study}_VEP_variant_summary.tsv"
    all_rows = ann_dir / f"{study}_VEP_all_consequences.tsv.gz"
    if any(not p.is_file() or p.stat().st_size == 0 for p in (fine, summary, all_rows, meta_file)):
        return None, "MISSING_REQUIRED_STEP06_OR_STEP07_FILE"
    if (ann_dir / "ANNOTATION_FAILED.txt").exists():
        return None, "ANNOTATION_FAILED_MARKER"
    if (str(meta.get("STATUS", "")).upper() != "COMPLETE" or
            str(meta.get("ANNOTATION_SCOPE", "")).upper() != "ALL_FINEMAPPED_VARIANTS"):
        return None, "INVALID_STEP07_STATUS_OR_SCOPE"
    try:
        if (int(meta["FINEMAPPED_INPUT_SIZE"]) != fine.stat().st_size or
                int(meta["FINEMAPPED_INPUT_MTIME_NS"]) != fine.stat().st_mtime_ns):
            return None, "STALE_ANNOTATION_SOURCE_FINGERPRINT"
    except (KeyError, TypeError, ValueError):
        return None, "MISSING_ANNOTATION_SOURCE_FINGERPRINT"
    return (summary, all_rows, meta, fine), "READY"


def initial_variant(row: dict):
    pip = pip_float(row.get("PIP", ""))
    return {
        "STUDY_ACCESSION": "", "LOCUS_ID": str(row.get("LOCUS_ID", "")).strip(),
        "VEP_ID": str(row.get("VEP_ID", "")).strip(),
        "ALLELE": normalize_allele(row), "PIP": "" if pip is None else pip,
        "CS95": int(flag(row.get("IN_95_CREDIBLE_SET", ""))),
        "HIGH_PIP": 0, "ANY_REGIONAL": 0,
        "LEGACY_ANY_CONSEQUENCE_CDS": int(bool(terms(row.get("ALL_CONSEQUENCES", "")) & CDS)),
        "seen": False, "row_count": 0, "features": set(), "biotypes": set(),
        "pc_genes": set(), "lnc_genes": set(), "other_genes": set(),
        "coding_evidence": set(), "legacy_cds_rows": set(),
    }


def select_variants(summary_path: Path, high_pip: float, all_regional: bool):
    selected = {}
    counts = Counter()
    with summary_path.open("r", encoding="utf-8-sig", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        required = {"VEP_ID", "LOCUS_ID", "PIP", "IN_95_CREDIBLE_SET",
                    "CHR", "REFERENCE_POS", "REFERENCE_REF", "REFERENCE_ALT"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError("Step07 variant summary missing columns: " + ", ".join(sorted(missing)))
        for row in reader:
            counts["SUMMARY_ROWS"] += 1
            pip = pip_float(row.get("PIP", ""))
            cs = flag(row.get("IN_95_CREDIBLE_SET", ""))
            hp = pip is not None and pip >= high_pip
            if cs:
                counts["CS95_ROWS"] += 1
                if pip is None:
                    counts["CS95_INVALID_PIP"] += 1
            if hp:
                counts["HIGH_PIP_ROWS"] += 1
            if not (cs or hp or all_regional):
                continue
            vep_id = row.get("VEP_ID", "").strip()
            locus = row.get("LOCUS_ID", "").strip()
            allele = normalize_allele(row)
            if not vep_id or not locus or not allele:
                raise ValueError(f"Missing VEP_ID, locus or exact allele in summary: {row}")
            if vep_id in selected:
                raise ValueError(f"Duplicate VEP_ID in variant summary: {vep_id}")
            r = initial_variant(row)
            r["HIGH_PIP"] = int(hp)
            r["ANY_REGIONAL"] = 1
            selected[vep_id] = r
    if counts["CS95_INVALID_PIP"]:
        raise ValueError(f"{counts['CS95_INVALID_PIP']} credible-set variants have invalid PIP")
    return selected, counts


def annotated_biotype(bio: str):
    if bio == "protein_coding":
        return "pc"
    if bio in {"lncrna", "lincrna"}:
        return "lnc"
    if "pseudogene" in bio or bio in PSEUDO_BIOTYPES:
        return "pseudo"
    if bio in OTHER_NCRNA_BIOTYPES or bio.endswith("rna") or bio.endswith("_rna"):
        return "other_nc"
    return "other"


def classify_transcript_row(record: dict, rec: dict):
    t = terms(record.get("Consequence", ""))
    bio = str(record.get("BIOTYPE", "")).lower().strip()
    if not t and not bio:
        return
    f = rec["features"]
    if t:
        f.add("ANY_ANNOTATED")
    rec["biotypes"].add(bio or "unknown")
    category = annotated_biotype(bio)
    e = transcript_number(record.get("EXON", ""))
    i = transcript_number(record.get("INTRON", ""))
    sp_core = bool(t & SPLICE_CORE)
    sp_region = bool(t & SPLICE_REGION)
    coding = bool(t & CDS)
    tx = bool(e or i or sp_core or sp_region or (t & OTHER_TX)
              or (t & CDS) or (t & UTR5) or (t & UTR3) or (t & OTHER_TRANSCRIPT_TERMS))
    if sp_core:
        f.add("ANY_SPLICE_CORE")
    if sp_region:
        f.add("ANY_SPLICE_REGION")
    if t & UPSTREAM:
        f.add("UPSTREAM")
    if t & DOWNSTREAM:
        f.add("DOWNSTREAM")
    if t & INTERGENIC:
        f.add("INTERGENIC")
    if t & REGULATORY:
        f.add("REGULATORY")
    gene = str(record.get("Gene", "")).strip()
    symbol = str(record.get("SYMBOL", "")).strip()
    feature = str(record.get("Feature", "")).strip()
    stable = gene if gene and gene.lower() not in {"nan", "-", "none"} else (symbol if symbol and symbol != "-" else "")
    if category == "pc":
        if coding:
            f.add("PC_CDS")
            rec["coding_evidence"].add("|".join([stable or "NA", feature or "NA", "&".join(sorted(t & CDS))]))
        if sp_core:
            f.add("PC_SPLICE_CORE")
        if sp_region:
            f.add("PC_SPLICE_REGION")
        if t & UTR5:
            f.add("PC_5UTR")
        if t & UTR3:
            f.add("PC_3UTR")
        if e and not coding and not (t & UTR5) and not (t & UTR3):
            f.add("PC_EXON_NONCDS")
        if i or t & INTRON:
            f.add("PC_INTRON")
        if tx:
            f.add("PC_GENE_BODY")
            if stable:
                rec["pc_genes"].add(stable)
        if tx and not coding and not any(x in f for x in ("PC_INTRON", "PC_5UTR", "PC_3UTR")):
            f.add("PC_OTHER_BODY")
    elif category == "lnc":
        if e or t & EXON:
            f.add("LNCRNA_EXON")
        if i or t & INTRON:
            f.add("LNCRNA_INTRON")
        if sp_core:
            f.add("LNCRNA_SPLICE_CORE")
        if sp_region:
            f.add("LNCRNA_SPLICE_REGION")
        if t & OTHER_TX:
            f.add("LNCRNA_OTHER_BODY")
        if any(x in f for x in ("LNCRNA_EXON", "LNCRNA_INTRON", "LNCRNA_SPLICE_CORE", "LNCRNA_SPLICE_REGION", "LNCRNA_OTHER_BODY")):
            # IMPORTANT: gene collected only for true transcript intersection;
            # upstream/downstream-only calls must never be called lncRNA overlap.
            f.add("LNCRNA_BODY")
            if stable and (e or i or sp_core or sp_region or bool(t & EXON) or bool(t & INTRON) or bool(t & OTHER_TX)):
                rec["lnc_genes"].add(stable)
    elif category == "pseudo":
        if tx:
            f.add("PSEUDOGENE_BODY")
            if stable:
                rec["other_genes"].add(stable)
    elif category == "other_nc":
        if tx:
            f.add("OTHER_NCRNA_BODY")
            if e or t & EXON or "mature_mirna_variant" in t:
                f.add("OTHER_NCRNA_EXON")
            if stable:
                rec["other_genes"].add(stable)
    else:
        if tx:
            f.add("OTHER_GENE_BODY")
            if stable:
                rec["other_genes"].add(stable)


def finalize_variant(rec: dict) -> dict:
    f = rec["features"]
    if not rec["seen"]:
        primary = "UNRESOLVED_NO_VEP_ROWS"
        cds_category = "UNRESOLVED"
    else:
        primary = next((label for label, term in PRIMARY_ORDER if term in f), "UNRESOLVED_ANNOTATION")
        cds_category = "PROTEIN_CODING_CDS" if "PC_CDS" in f else (
            "NO_PROTEIN_CODING_CDS" if "ANY_ANNOTATED" in f else "UNRESOLVED")
    out = {k: v for k, v in rec.items() if k not in {
        "features", "biotypes", "pc_genes", "lnc_genes", "other_genes",
        "coding_evidence", "legacy_cds_rows", "seen", "row_count",
    }}
    out.update({
        "CODING_CLASS": cds_category, "PRIMARY_REGION": primary,
        "N_VEP_TRANSCRIPT_ROWS": rec["row_count"],
        "N_BIOTYPES": len(rec["biotypes"]),
        "ALL_BIOTYPES": ";".join(sorted(rec["biotypes"])),
        "CODING_TRANSCRIPT_EVIDENCE": ";".join(sorted(rec["coding_evidence"])),
        "PC_GENE_IDS": ";".join(sorted(rec["pc_genes"])),
        "LNCRNA_GENE_IDS": ";".join(sorted(rec["lnc_genes"])),
        "OTHER_GENE_IDS": ";".join(sorted(rec["other_genes"])),
        "LEGACY_VS_STRICT_CDS_DISAGREEMENT": int(bool(rec["LEGACY_ANY_CONSEQUENCE_CDS"]) != ("PC_CDS" in f)),
    })
    out.update({term: int(term in f) for term in FLAGS})
    return out


def annotate_selected(all_path: Path, selected: dict):
    vep_rows = 0
    with gzip.open(all_path, "rt", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        required = {"VEP_ID", "BIOTYPE", "Consequence", "Gene", "EXON", "INTRON"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError("VEP all-consequences missing columns: " + ", ".join(sorted(missing)))
        for row in reader:
            vep_rows += 1
            r = selected.get(row.get("VEP_ID", ""))
            if r is None:
                continue
            r["seen"] = True
            r["row_count"] += 1
            classify_transcript_row(row, r)
    no_vep = sum(not r["seen"] for r in selected.values())
    return vep_rows, no_vep


def summarize(rows: list[dict], pop: str) -> dict:
    if pop == "CREDIBLE95":
        chosen = [r for r in rows if r["CS95"]]
    elif pop == "HIGH_PIP":
        chosen = [r for r in rows if r["HIGH_PIP"]]
    else:
        chosen = rows
    n = len(chosen)
    categories = Counter(r["CODING_CLASS"] for r in chosen)
    primary = Counter(r["PRIMARY_REGION"] for r in chosen)
    flags = Counter({f: sum(r[f] for r in chosen) for f in FLAGS})
    if sum(categories.values()) != n or sum(primary.values()) != n:
        raise AssertionError("Classification totals do not reconcile to selected candidate denominator")
    return {"POPULATION": pop, "N_VARIANTS": n,
            "CODING_CDS": categories["PROTEIN_CODING_CDS"],
            "NO_CDS": categories["NO_PROTEIN_CODING_CDS"],
            "UNRESOLVED": categories["UNRESOLVED"],
            "LNCRNA_BODY": flags["LNCRNA_BODY"],
            "PC_GENE_BODY": flags["PC_GENE_BODY"],
            "BOTH_LNCRNA_PC": sum(bool(r["LNCRNA_BODY"] and r["PC_GENE_BODY"]) for r in chosen),
            "LEGACY_CDS": sum(r["LEGACY_ANY_CONSEQUENCE_CDS"] for r in chosen),
            "LEGACY_CDS_DISCORDANCE": sum(r["LEGACY_VS_STRICT_CDS_DISAGREEMENT"] for r in chosen),
            "CODING_COUNTS": categories, "REGION_COUNTS": primary, "FEATURE_COUNTS": flags}


def pstr(n: int, den: int) -> str:
    return f"{n:,} ({n/den*100:.3f}%)" if den else "0 (N/A)"


def print_summary(title: str, s: dict):
    n = s["N_VARIANTS"]
    print("\n" + "-" * 112)
    print(title + f" | n={n:,}")
    print("-" * 112)
    print(f"  Strict protein-coding CDS : {pstr(s['CODING_CDS'], n)}")
    print(f"  No protein-coding CDS     : {pstr(s['NO_CDS'], n)} (NOT equivalent to intergenic)")
    print(f"  Unresolved annotation     : {pstr(s['UNRESOLVED'], n)}")
    print(f"  Strict lncRNA transcript  : {pstr(s['LNCRNA_BODY'], n)}")
    print(f"  Protein-coding gene body  : {pstr(s['PC_GENE_BODY'], n)}")
    print(f"  Both transcript classes   : {pstr(s['BOTH_LNCRNA_PC'], n)}")
    print(f"  Legacy aggregated CDS     : {s['LEGACY_CDS']:,}; mismatch IDs: {s['LEGACY_CDS_DISCORDANCE']:,}")
    print("\n  MUTUALLY EXCLUSIVE PRIMARY GENOMIC REGION (explicit precedence):")
    for label, cnt in sorted(s["REGION_COUNTS"].items(), key=lambda kv: (-kv[1], kv[0])):
        print(f"    {label:35s} {pstr(cnt, n):>25}")
    print(f"    {'SUM':35s} {pstr(sum(s['REGION_COUNTS'].values()), n):>25}")
    print("\n  NONEXCLUSIVE SAME-ROW TRANSCRIPT / GENOMIC FLAGS:")
    for label in FLAGS:
        if s["FEATURE_COUNTS"][label]:
            print(f"    {label:35s} {pstr(s['FEATURE_COUNTS'][label], n):>25}")


def output_summary(summary: dict, **fields) -> dict:
    return {**fields, **{k: v for k, v in summary.items() if k not in {
        "CODING_COUNTS", "REGION_COUNTS", "FEATURE_COUNTS"}}}


def process_study(root: Path, disease: str, study: str, files, args):
    summary_file, all_vep, meta, fine = files
    selected, input_counts = select_variants(summary_file, args.high_pip, args.include_all_regional)
    print(f"  [STUDY] {study}: Step07 variants={input_counts['SUMMARY_ROWS']:,}; "
          f"CS95={input_counts['CS95_ROWS']:,}; PIP>={args.high_pip:g}={input_counts['HIGH_PIP_ROWS']:,}; "
          f"selected for transcript checks={len(selected):,}", flush=True)
    if not selected:
        return [], input_counts, {"VEP_ROWS_SCANNED": 0, "UNMATCHED_VEP_ID": 0}
    nrows, missing = annotate_selected(all_vep, selected)
    if missing:
        raise ValueError(f"{missing} candidate VEP_ID(s) not present in complete all-consequences file")
    records = [finalize_variant(v) for v in selected.values()]
    for r in records:
        r["STUDY_ACCESSION"] = study
        r["DISEASE"] = DISEASES[disease][0]
        r["PHENOTYPE_SLUG"] = disease
    ncs = sum(r["CS95"] for r in records)
    if ncs != input_counts["CS95_ROWS"]:
        raise AssertionError(f"Credible-set input/output mismatch {ncs} vs {input_counts['CS95_ROWS']}")
    return records, input_counts, {"VEP_ROWS_SCANNED": nrows,
                                    "UNMATCHED_VEP_ID": missing,
                                    "STEP07_VERSION": meta.get("STEP07_VERSION", ""),
                                    "VEP_FILE_BYTES": all_vep.stat().st_size,
                                    "VEP_FILE_MTIME_NS": all_vep.stat().st_mtime_ns,
                                    "STEP06_FILE_BYTES": fine.stat().st_size,
                                    "STEP06_FILE_MTIME_NS": fine.stat().st_mtime_ns}


def run(root, diseases, ancestry="european", high_pip=0.5, include_all_regional=False, **_):
    """Run Analysis 03 for [(display name, folder slug), ...]; return its tables."""
    global DISEASES, ANCESTRY
    DISEASES = {slug: (name, slug) for name, slug in diseases}
    ANCESTRY = ancestry
    args = argparse.Namespace(high_pip=float(high_pip), include_all_regional=bool(include_all_regional))
    if not 0 <= args.high_pip <= 1:
        raise ValueError("high_pip must be in [0,1]")
    disease_slugs = list(DISEASES)
    root = Path(root).expanduser().resolve()
    print("=" * 122)
    print(f"JOHN ANALYSIS 03 v{VERSION} | STRICT VEP SAME-ROW CODING CLASSIFICATION")
    print(f"GWAS2m root        : {root}")
    print(f"Phenotypes         : {', '.join(DISEASES[d][0] for d in disease_slugs)}")
    print(f"Populations        : CS95, HIGH_PIP>={args.high_pip:g}" +
          (", ALL_REGIONAL (higher resource use)" if args.include_all_regional else ""))
    print("Coding definition  : BIOTYPE=protein_coding AND CDS consequence SAME transcript row")
    print("No CDS != intergenic; region labels are 1-per-variant, other flags nonexclusive")
    print("Inputs             : READ ONLY; NO recoding of GWAS, SuSiE or VEP")
    print("=" * 122, flush=True)

    all_rows = []
    audit_rows = []
    study_summary = []
    error_count = 0
    populations = ["CREDIBLE95", "HIGH_PIP"] + (["ALL_REGIONAL"] if args.include_all_regional else [])
    for disease in disease_slugs:
        disease_name = DISEASES[disease][0]
        manifest, source_manifest = study_manifest(root, disease)
        print("\n" + "=" * 122)
        print(f"DISEASE: {disease_name}; discovered accessions: {len(manifest)}")
        print(f"Manifest: {source_manifest or 'discovered from Step06 directories'}")
        print("=" * 122, flush=True)
        disease_rows = []
        for study, meta_row in sorted(manifest.items()):
            files, state = study_readiness(root, disease, study)
            trait = get_trait(meta_row)
            phenotype_class = trait_label_type(trait, disease) if trait else "UNKNOWN_TRAIT_LABEL"
            audit = {"DISEASE": disease_name, "STUDY_ACCESSION": study,
                     "TRAIT": trait, "PHENOTYPE_CLASS": phenotype_class,
                     "READINESS": state, "N_CS95": "", "N_CDS": "", "N_STRICT_LNCRNA": ""}
            if files is None:
                audit_rows.append(audit)
                print(f"  [SKIP] {study}: {state}", flush=True)
                continue
            try:
                records, counts, info = process_study(root, disease, study, files, args)
                audit.update({"READINESS": "INCLUDED", "N_CS95": counts["CS95_ROWS"],
                              "N_CDS": sum(r["CS95"] and r["PC_CDS"] for r in records),
                              "N_STRICT_LNCRNA": sum(r["CS95"] and r["LNCRNA_BODY"] for r in records),
                              "VEP_ROWS_SCANNED": info["VEP_ROWS_SCANNED"],
                              "STEP07_VERSION": info.get("STEP07_VERSION", ""),
                              "VEP_FILE_BYTES": info.get("VEP_FILE_BYTES", ""),
                              "VEP_FILE_MTIME_NS": info.get("VEP_FILE_MTIME_NS", "")})
                for r in records:
                    r["TRAIT"] = trait
                    r["PHENOTYPE_CLASS"] = phenotype_class
                disease_rows.extend(records)
                all_rows.extend(records)
                for pop in populations:
                    s = summarize(records, pop)
                    study_summary.append(output_summary(s, DISEASE=disease_name,
                                                       STUDY_ACCESSION=study, TRAIT=trait,
                                                       PHENOTYPE_CLASS=phenotype_class))
                print(f"  [OK]   {study}: CS95={audit['N_CS95']:,}; strict CDS={audit['N_CDS']:,}; "
                      f"strict lncRNA={audit['N_STRICT_LNCRNA']:,}", flush=True)
            except (OSError, EOFError, ValueError, RuntimeError, AssertionError, csv.Error) as exc:
                error_count += 1
                audit["READINESS"] = "ERROR:" + str(exc)[:250]
                print(f"  [ERROR] {study}: {exc}", flush=True)
            audit_rows.append(audit)
        for pop in populations:
            s = summarize(disease_rows, pop)
            print_summary(f"{disease_name} | {pop}", s)
            study_summary.append(output_summary(s, DISEASE=disease_name,
                                               STUDY_ACCESSION="ALL_INCLUDED", TRAIT="",
                                               PHENOTYPE_CLASS="ALL_CLASSES"))
        print("\n  STUDY AVAILABILITY:")
        for label, n in sorted(Counter(x["READINESS"] for x in audit_rows if x["DISEASE"] == disease_name).items()):
            print(f"    {label:48s} {n:,}")

    if len(disease_slugs) > 1:
        for pop in populations:
            print_summary("ALL SELECTED DISEASES POOLED | " + pop, summarize(all_rows, pop))
    tables = {}
    if True:
        var_fields = (["DISEASE", "PHENOTYPE_SLUG", "STUDY_ACCESSION", "TRAIT", "PHENOTYPE_CLASS",
                      "LOCUS_ID", "VEP_ID", "ALLELE", "PIP", "CS95", "HIGH_PIP", "ANY_REGIONAL",
                      "CODING_CLASS", "PRIMARY_REGION", "N_VEP_TRANSCRIPT_ROWS", "ALL_BIOTYPES",
                      "PC_GENE_IDS", "LNCRNA_GENE_IDS", "OTHER_GENE_IDS",
                      "CODING_TRANSCRIPT_EVIDENCE", "LEGACY_ANY_CONSEQUENCE_CDS",
                      "LEGACY_VS_STRICT_CDS_DISAGREEMENT"] + list(FLAGS))
        tables["variant_classes"] = pd.DataFrame(all_rows, columns=var_fields)
        tables["study_summary"] = pd.DataFrame(study_summary)
        tables["study_eligibility"] = pd.DataFrame(audit_rows)
        cs_rows = [r for r in all_rows if r["CS95"]]
        discrepancy = [r for r in cs_rows if r["LEGACY_VS_STRICT_CDS_DISAGREEMENT"]]
        tables["legacy_vs_strict_cds"] = pd.DataFrame(discrepancy, columns=var_fields)
        feature_stats = []
        region_stats = []
        locus_stats = []
        for disease in disease_slugs:
            subset = [r for r in all_rows if r["PHENOTYPE_SLUG"] == disease]
            for pop in populations:
                s = summarize(subset, pop)
                for label, count in sorted(s["FEATURE_COUNTS"].items()):
                    feature_stats.append({"DISEASE": DISEASES[disease][0], "POPULATION": pop,
                                          "FEATURE": label, "N": count, "DENOMINATOR": s["N_VARIANTS"]})
                for label, count in sorted(s["REGION_COUNTS"].items()):
                    region_stats.append({"DISEASE": DISEASES[disease][0], "POPULATION": pop,
                                         "REGION": label, "N": count, "DENOMINATOR": s["N_VARIANTS"]})
            by_locus = defaultdict(list)
            for r in subset:
                if r["CS95"]:
                    by_locus[(r["STUDY_ACCESSION"], r["LOCUS_ID"])].append(r)
            for (study, locus), records in sorted(by_locus.items()):
                s = summarize(records, "CREDIBLE95")
                locus_stats.append(output_summary(s, DISEASE=DISEASES[disease][0],
                                                  STUDY_ACCESSION=study, LOCUS_ID=locus))
        tables["feature_flags"] = pd.DataFrame(feature_stats)
        tables["primary_regions"] = pd.DataFrame(region_stats)
        tables["credible_locus_summary"] = pd.DataFrame(locus_stats)
    print("\nVALIDATION NOTES:")
    print("- Strict lncRNA/CDS requires BIOTYPE+Consequence in same VEP transcript record.")
    print("- Absence of a protein-coding CDS does NOT establish absence of gene regulation.")
    print("- Do not conflate these annotation counts with causal effects or colocalization.")
    print("- Compare Analysis03 strict CDS against Analysis08; inspect the mismatch TSV.")
    if error_count:
        print(f"[WARNING] {error_count} eligible studies failed analysis; see the study_eligibility table")
    return tables


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    if "--self-test" in argv:
        return self_test()
    return analysis_main(run, "John Analysis03: strict protein-coding CDS vs non-CDS", argv, extra=[
        ("--high-pip", dict(type=float, default=0.5)),
        ("--include-all-regional", dict(action="store_true", help="Slow / high RAM: include every regional variant")),
    ])


def self_test() -> int:
    def r(bio, term, exon="", intron="", gene="ENSG1"):
        return {"BIOTYPE": bio, "Consequence": term, "EXON": exon, "INTRON": intron,
                "Gene": gene, "Feature": "ENST1", "SYMBOL": "TEST"}
    def new(legacy=0):
        x = initial_variant({"LOCUS_ID": "L001", "VEP_ID": "V1", "CHR": "1",
                             "REFERENCE_POS": "1", "REFERENCE_REF": "A", "REFERENCE_ALT": "C",
                             "PIP": "0.8", "IN_95_CREDIBLE_SET": "True",
                             "ALL_CONSEQUENCES": "missense_variant" if legacy else ""})
        x["seen"] = True
        return x
    # A coding term on a non-protein-coding transcript is NOT a strict CDS.
    v = new(legacy=1)
    classify_transcript_row(r("lncRNA", "missense_variant"), v)
    assert not finalize_variant(v)["PC_CDS"]
    assert finalize_variant(v)["LEGACY_VS_STRICT_CDS_DISAGREEMENT"]
    # A VEP row with protein_coding missense is a CDS, regardless of other transcript rows.
    v = new()
    classify_transcript_row(r("protein_coding", "missense_variant"), v)
    classify_transcript_row(r("lncRNA", "upstream_gene_variant", gene="ENSG2"), v)
    t = finalize_variant(v)
    assert t["PC_CDS"] and not t["LNCRNA_BODY"] and not t["LNCRNA_GENE_IDS"]
    # An lncRNA transcript body is detected by its own exon/intron record.
    v = new()
    classify_transcript_row(r("lncRNA", "non_coding_transcript_exon_variant", exon="2/4"), v)
    t = finalize_variant(v)
    assert t["LNCRNA_BODY"] and t["LNCRNA_EXON"] and t["CODING_CLASS"] == "NO_PROTEIN_CODING_CDS"
    # A protein-coding intron is not coding CDS.
    v = new()
    classify_transcript_row(r("protein_coding", "intron_variant", intron="2/4"), v)
    t = finalize_variant(v)
    assert t["PC_INTRON"] and not t["PC_CDS"] and t["PRIMARY_REGION"] == "PROTEIN_CODING_INTRON"
    # Multiple transcript labels may overlap, but only one primary class.
    v = new()
    classify_transcript_row(r("protein_coding", "intron_variant", intron="2/4"), v)
    classify_transcript_row(r("lncRNA", "intron_variant", intron="1/3", gene="ENSG2"), v)
    t = finalize_variant(v)
    assert t["PC_GENE_BODY"] and t["LNCRNA_BODY"] and t["PRIMARY_REGION"] == "PROTEIN_CODING_INTRON"
    assert summarize([t], "CREDIBLE95")["N_VARIANTS"] == 1
    print("[PASS] Same-row CDS, coding intron, lncRNA exon, near-only exclusion, overlap and denominator tests")
    return 0


if __name__ == "__main__":
    result = main()
    sys.exit(result if isinstance(result, int) else 0)
