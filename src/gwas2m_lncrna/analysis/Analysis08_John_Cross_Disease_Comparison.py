#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Analysis 08 - cross-disease fine-mapped lncRNA comparison.

Only reads existing GWAS2m Step06/Step07 outputs for the phenotypes passed in.
Prints to stdout; run() also returns the result tables.
"""

import re
import argparse
import csv
import gzip
import json
import math
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd

from gwas2m_lncrna.analysis.helpers import analysis_main, direct_labels, normalize_label

# Set by run(): [(display name, phenotype folder), ...] and the ancestry folder name.
DISEASES = []
ANCESTRY = "european"
TABLES = {}
DONE = {"COMPLETE", "SKIPPED_ALREADY_COMPLETE"}
CODING = {
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
IN_TRANSCRIPT = {
    "non_coding_transcript_variant", "non_coding_transcript_exon_variant",
    "intron_variant",
}
PC_IN_GENE = CODING | SPLICE_CORE | SPLICE_REGION | {
    "intron_variant", "5_prime_utr_variant", "3_prime_utr_variant",
    "non_coding_transcript_exon_variant", "non_coding_transcript_variant",
    "transcript_ablation", "transcript_amplification",
}
KEYS = (
    "n", "coding", "noncoding", "unresolved", "lnc", "lnc_exon",
    "lnc_intron", "lnc_splice_core", "lnc_splice_region", "lnc_other",
    "pc_gene", "both", "high_pip", "invalid_pip",
)

try:
    csv.field_size_limit(sys.maxsize)
except OverflowError:
    csv.field_size_limit(2 ** 31 - 1)


def load_json(path):
    try:
        with path.open("r", encoding="utf-8") as handle:
            obj = json.load(handle)
        return obj if isinstance(obj, dict) else {}
    except (OSError, ValueError, UnicodeError):
        return {}


def truth(value):
    return str(value).strip().lower() in {"1", "t", "true", "y", "yes"}


def pip_value(value):
    try:
        x = float(value)
        return x if math.isfinite(x) and 0.0 <= x <= 1.0 else None
    except (TypeError, ValueError):
        return None


def vep_terms(value):
    """Split VEP compound consequences while preserving same-row provenance."""
    return {
        x.strip().lower()
        for x in re.split(r"[,;&|]", str(value or ""))
        if x.strip() and x.strip().lower() not in {"nan", "none", "na", ".", "-"}
    }


def transcript_number(value):
    s = str(value or "").strip()
    return bool(s and s[0].isdigit() and s not in {"0", "0/0"})


def norm_chr(value):
    v = str(value).strip().lower().removeprefix("chr")
    if not v.isdigit():
        return None
    n = int(v)
    return str(n) if 1 <= n <= 22 else None


def allele_key(row):
    chrom = norm_chr(row.get("CHR", ""))
    pos = str(row.get("REFERENCE_POS", "")).strip()
    ref = str(row.get("REFERENCE_REF", "")).strip().upper()
    alt = str(row.get("REFERENCE_ALT", "")).strip().upper()
    if not chrom or not pos.isdigit() or int(pos) < 1 or not ref or not alt:
        return None
    return (chrom, int(pos), ref, alt)


def fmt_allele(key):
    return "%s:%s:%s:%s" % key


def phenotype_type(trait, name):
    v = " ".join(str(trait or "").lower().replace("\u2019", "'").split())
    if any(x in v for x in ("family history", "parental", "proxy", "father", "mother", "sibling")):
        return "PROXY_FAMILY"
    if normalize_label(v) in direct_labels(name):
        return "DIRECT_LABEL"
    if any(x in v for x in (
        "interaction", "medication", "drug response", "mtag", "severity",
        "age at onset", "risk score", "subtype", "dementia", "cognitive",
        "motor", "tremor", "neuropathology", "biomarker",
    )):
        return "RELATED_SPECIALIZED"
    return "REVIEW_TRAIT"


def status(root, slug, accession, stage):
    path = (root / "pipeline_audit" / slug / ANCESTRY /
            "status" / accession / (stage + ".json"))
    return str(load_json(path).get("status", "NOT_RECORDED")).upper()


def read_manifest(root, slug):
    path = (root / "01_gwas_catalog" / slug / ANCESTRY /
            "GWAS_download_manifest.tsv")
    if not path.is_file():
        path = root / "pipeline_audit" / slug / ANCESTRY / "study_manifest.tsv"
    if not path.is_file():
        return None
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def readiness(root, slug, accession):
    s06 = status(root, slug, accession, "Step06_FineMapping")
    if s06 not in DONE:
        return None, "STEP06_" + s06
    base = root / "06_finemapping" / slug / ANCESTRY / accession
    fine = base / (accession + "_finemapped_variants.tsv.gz")
    fm = load_json(base / "finemapping_summary.json")
    if str(fm.get("STATUS", "")).upper() == "NO_SIGNIFICANT_LOCI" or (
        fm and fm.get("N_LOCI") == 0
    ):
        return None, "NO_SIGNIFICANT_LOCI"
    s07 = status(root, slug, accession, "Step07_VEP")
    if s07 not in DONE:
        return None, "STEP07_" + s07
    directory = root / "07_annotation" / slug / ANCESTRY / accession
    meta = load_json(directory / "annotation_summary.json")
    summary = directory / (accession + "_VEP_variant_summary.tsv")
    all_vep = directory / (accession + "_VEP_all_consequences.tsv.gz")
    if not all(p.is_file() and p.stat().st_size > 0 for p in (fine, summary, all_vep)):
        return None, "MISSING_INPUT"
    if (directory / "ANNOTATION_FAILED.txt").exists():
        return None, "ANNOTATION_FAILED_MARKER"
    if (str(meta.get("STATUS", "")).upper() != "COMPLETE" or
            str(meta.get("ANNOTATION_SCOPE", "")).upper() != "ALL_FINEMAPPED_VARIANTS"):
        return None, "INCOMPLETE_ANNOTATION_SCOPE"
    try:
        if (int(meta["FINEMAPPED_INPUT_SIZE"]) != fine.stat().st_size or
                int(meta["FINEMAPPED_INPUT_MTIME_NS"]) != fine.stat().st_mtime_ns):
            return None, "STALE_ANNOTATION"
    except (KeyError, TypeError, ValueError):
        return None, "ANNOTATION_PROVENANCE_MISSING"
    return (summary, all_vep, meta), "READY"


def parse_candidates(summary, threshold):
    """Read each variant once; retain only CS95/high-PIP candidates."""
    chosen = {}
    total = 0
    invalid_pips = 0
    credible_bad_pip = 0
    credible_bad_key = 0
    with summary.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        needed = {
            "VEP_ID", "LOCUS_ID", "PIP", "IN_95_CREDIBLE_SET",
            "CHR", "REFERENCE_POS", "REFERENCE_REF", "REFERENCE_ALT",
        }
        missing = needed - set(reader.fieldnames or ())
        if missing:
            raise ValueError("Step07 summary missing columns: " + ", ".join(sorted(missing)))
        for row in reader:
            total += 1
            pip = pip_value(row.get("PIP"))
            credible = truth(row.get("IN_95_CREDIBLE_SET"))
            if pip is None:
                invalid_pips += 1
                if credible:
                    credible_bad_pip += 1
            if not (credible or (pip is not None and pip >= threshold)):
                continue
            key = allele_key(row)
            if key is None:
                if credible:
                    credible_bad_key += 1
                continue
            vid = row.get("VEP_ID", "").strip()
            locus = row.get("LOCUS_ID", "").strip()
            if not vid or not locus or vid in chosen:
                raise ValueError("Duplicate/empty VEP_ID or locus: " + repr(vid))
            chosen[vid] = {
                "allele": key, "locus": locus, "pip": pip,
                "credible": credible, "high_pip": pip is not None and pip >= threshold,
                "seen": False, "annotated": False, "coding": False,
                "pc_gene": False, "lnc_exon": False, "lnc_intron": False,
                "lnc_splice_core": False, "lnc_splice_region": False,
                "lnc_other": False, "genes": set(),
            }
    if credible_bad_pip or credible_bad_key:
        raise ValueError(
            "Invalid credible-set inputs: PIP=" + str(credible_bad_pip) +
            ", alleles=" + str(credible_bad_key)
        )
    return chosen, total, invalid_pips


def annotate_chosen(all_vep, chosen):
    """Stream all VEP rows; classify each candidate by same-row biotype/consequence."""
    required = {"VEP_ID", "Consequence", "BIOTYPE", "Gene", "EXON", "INTRON"}
    n_rows = 0
    with gzip.open(all_vep, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = required - set(reader.fieldnames or ())
        if missing:
            raise ValueError("Step07 transcript annotations missing: " + ", ".join(sorted(missing)))
        for row in reader:
            n_rows += 1
            rec = chosen.get(row.get("VEP_ID", ""))
            if rec is None:
                continue
            rec["seen"] = True
            terms = vep_terms(row.get("Consequence"))
            biotype = str(row.get("BIOTYPE", "")).strip().lower()
            exon = transcript_number(row.get("EXON"))
            intron = transcript_number(row.get("INTRON"))
            if terms and terms != {"-"} and terms != {"sequence_variant"}:
                rec["annotated"] = True
            if biotype == "protein_coding":
                if terms & CODING:
                    rec["coding"] = True
                if (exon or intron or (terms & PC_IN_GENE)):
                    rec["pc_gene"] = True
            if biotype not in {"lncrna", "lincrna"}:
                continue
            core = bool(terms & SPLICE_CORE)
            region = bool(terms & SPLICE_REGION)
            exonic = exon or "non_coding_transcript_exon_variant" in terms
            intronic = intron or "intron_variant" in terms
            other = "non_coding_transcript_variant" in terms
            if not (core or region or exonic or intronic or other):
                continue  # Nearby upstream/downstream is NOT transcript overlap.
            rec["lnc_splice_core"] |= core
            rec["lnc_splice_region"] |= region
            rec["lnc_exon"] |= exonic
            rec["lnc_intron"] |= intronic
            rec["lnc_other"] |= other
            gene = str(row.get("SYMBOL", "")).strip()
            if gene in {"", "-", "nan"}:
                gene = str(row.get("Gene", "")).strip()
            if gene not in {"", "-", "nan"}:
                rec["genes"].add(gene)
    missing = sum(1 for rec in chosen.values() if not rec["seen"])
    if missing:
        raise ValueError(str(missing) + " selected variants missing from VEP all-consequences file")
    return n_rows


def get_lnc(rec):
    return any(rec[k] for k in (
        "lnc_exon", "lnc_intron", "lnc_splice_core", "lnc_splice_region", "lnc_other"
    ))


def category(rec):
    if rec["lnc_splice_core"]:
        return "LNCRNA_SPLICE_CORE"
    if rec["lnc_splice_region"]:
        return "LNCRNA_SPLICE_REGION"
    if rec["lnc_exon"]:
        return "LNCRNA_EXON"
    if rec["lnc_intron"]:
        return "LNCRNA_INTRON"
    if rec["lnc_other"]:
        return "LNCRNA_OTHER"
    return "NO_LNCRNA_TRANSCRIPT"


def count_candidates(chosen):
    counts = Counter()
    high = Counter()
    unique = set()
    lnc_unique = set()
    locus_flags = defaultdict(lambda: [False, False])
    genes = Counter()
    for rec in chosen.values():
        is_lnc = get_lnc(rec)
        for name, enabled in (("cs", rec["credible"]), ("high", rec["high_pip"])):
            if not enabled:
                continue
            target = counts if name == "cs" else high
            target["n"] += 1
            target["coding" if rec["coding"] else "noncoding" if rec["annotated"] else "unresolved"] += 1
            if is_lnc:
                target["lnc"] += 1
                target[category(rec)] += 1
            if rec["pc_gene"]:
                target["pc_gene"] += 1
            if rec["pc_gene"] and is_lnc:
                target["both"] += 1
            for label in (
                "lnc_exon", "lnc_intron", "lnc_splice_core", "lnc_splice_region"
            ):
                if rec[label]:
                    target[label] += 1
        if rec["credible"]:
            key = rec["allele"]
            unique.add(key)
            if is_lnc:
                lnc_unique.add(key)
                for gene in rec["genes"]:
                    genes[gene] += 1
            lf = locus_flags[rec["locus"]]
            lf[0] |= is_lnc
            lf[1] |= rec["pc_gene"]
    locus_counts = Counter()
    for has_lnc, has_pc in locus_flags.values():
        label = (
            "BOTH" if has_lnc and has_pc else
            "LNCRNA_ONLY" if has_lnc else
            "PROTEIN_GENE_ONLY" if has_pc else "NEITHER"
        )
        locus_counts[label] += 1
    return counts, high, unique, lnc_unique, locus_counts, genes


def fmt(n, denominator):
    return "N/A" if denominator <= 0 else f"{n * 100.0 / denominator:6.2f}%"


def print_matrix(title, stats_by_disease, count_label="Credible-set variant-study rows"):
    print("\n" + "=" * 109)
    print(title)
    print("=" * 109)
    lines = (
        (count_label, "n"),
        ("Protein-coding CDS variants", "coding"),
        ("Noncoding variants", "noncoding"),
        ("Unresolved variant annotations", "unresolved"),
        ("lncRNA transcript overlap", "lnc"),
        ("  lncRNA exon overlap", "lnc_exon"),
        ("  lncRNA intron overlap", "lnc_intron"),
        ("  lncRNA splice core overlap", "lnc_splice_core"),
        ("  lncRNA splice region overlap", "lnc_splice_region"),
        ("Protein-coding transcript overlap", "pc_gene"),
        ("Both lncRNA and protein genes", "both"),
    )
    print(f"{'MEASURE':40} " + " ".join(f"{name[:21].upper():>22}" for name, _ in DISEASES))
    print("-" * 109)
    for label, k in lines:
        cells = []
        for name, _ in DISEASES:
            s = stats_by_disease.get(name, Counter())
            n = s[k]
            cells.append(
                f"{n:,} (100%)" if k == "n" and n else
                (f"{n:,} ({fmt(n, s['n'])})" if s["n"] else "N/A")
            )
            TABLES.setdefault("matrices", []).append({
                "SUBSET": title, "MEASURE": label, "PHENOTYPE": name, "N": n,
                "PERCENT": round(100 * n / s["n"], 3) if s["n"] else None, "DENOMINATOR": s["n"]})
        print(f"{label:40} " + " ".join(f"{c:>22}" for c in cells))


def print_loci(loci_by_disease):
    print("\n" + "=" * 109)
    print("CREDIBLE-SET TRANSCRIPT OVERLAP PER STUDY-LOCUS")
    print("=" * 109)
    print("These categories describe credible-set variants, NOT every gene within a 2-Mb locus.")
    print(f"{'CATEGORY':40} " + " ".join(f"{name[:21].upper():>22}" for name, _ in DISEASES))
    for label in ("BOTH", "LNCRNA_ONLY", "PROTEIN_GENE_ONLY", "NEITHER"):
        cells = []
        for disease, _ in DISEASES:
            data = loci_by_disease.get(disease, Counter())
            n = data[label]
            cells.append(f"{n:,} ({fmt(n, sum(data.values()))})")
            TABLES.setdefault("locus_classes", []).append({
                "PHENOTYPE": disease, "CLASS": label, "STUDY_LOCI": n,
                "PERCENT": round(100 * n / sum(data.values()), 3) if sum(data.values()) else None})
        print(f"{label:40} " + " ".join(f"{c:>22}" for c in cells))


def run(root, diseases, ancestry="european", high_pip=0.50, top=20, **_):
    """Run Analysis 08 for [(display name, folder slug), ...]; return its tables."""
    global DISEASES, ANCESTRY
    DISEASES, ANCESTRY = list(diseases), ancestry
    TABLES.clear()
    args = argparse.Namespace(high_pip=float(high_pip), top=max(1, int(top)))
    if not 0.0 <= args.high_pip <= 1.0:
        raise ValueError("high_pip must be within 0..1")
    root = Path(root).expanduser().resolve()

    print("=" * 109)
    print("JOHN MATTICK - ANALYSIS 08 - CROSS-DISEASE COMPARISON")
    print("GRCh38 | VEP transcript-specific lncRNA overlap | 95% credible sets")
    print("READ ONLY")
    print("=" * 109)

    studies_summary = defaultdict(list)
    disease_all = defaultdict(Counter)
    disease_direct = defaultdict(Counter)
    disease_high = defaultdict(Counter)
    disease_loci = defaultdict(Counter)
    unique_by_disease = defaultdict(set)
    unique_lnc_by_disease = defaultdict(set)
    allele_which = defaultdict(lambda: defaultdict(set))
    allele_lnc_genes = defaultdict(set)
    readyness = defaultdict(Counter)

    for disease, slug in DISEASES:
        print("\n" + "-" * 109)
        print("DISEASE:", disease.upper())
        print("-" * 109)
        manifest = read_manifest(root, slug)
        if manifest is None:
            print("WARNING: Manifest unavailable; no results for this disease.")
            readyness[disease]["MANIFEST_MISSING"] += 1
            continue
        print("Selected studies:", len(manifest))
        print(f"{'STUDY':15} {'TRAIT_CLASS':19} {'CS_VARIANTS':>12} {'LNCRNA':>10} {'LNCRNA%':>10}  DECISION")
        for meta in manifest:
            accession = str(meta.get("STUDY_ACCESSION", "")).strip()
            if not accession:
                continue
            trait = meta.get("DISEASE_TRAIT", "") or meta.get("MAPPED_TRAIT", "")
            ptype = phenotype_type(trait, disease)
            paths, reason = readiness(root, slug, accession)
            if paths is None:
                readyness[disease][reason] += 1
                print(f"{accession:15} {ptype:19} {'-':>12} {'-':>10} {'-':>10}  {reason}")
                continue
            try:
                chosen, n_regional, bad_pips = parse_candidates(paths[0], args.high_pip)
                expected = paths[2].get("N_SELECTED_VARIANTS")
                if expected is not None and int(expected) != n_regional:
                    raise ValueError(f"Step07 total mismatch: {n_regional} versus {expected}")
                n_raw = annotate_chosen(paths[1], chosen)
                counts, high, alleles, lnc_alleles, loci, genes = count_candidates(chosen)
                if sum(counts[k] for k in ("coding", "noncoding", "unresolved")) != counts["n"]:
                    raise ValueError("Coding/noncoding/unresolved accounting mismatch")
            except (OSError, EOFError, UnicodeError, csv.Error, TypeError, ValueError) as exc:
                readyness[disease]["DATA_ERROR"] += 1
                print(f"{accession:15} {ptype:19} {'-':>12} {'-':>10} {'-':>10}  DATA_ERROR: {exc}")
                continue

            readyness[disease]["INCLUDED"] += 1
            disease_all[disease].update(counts)
            disease_high[disease].update(high)
            disease_loci[disease].update(loci)
            if ptype == "DIRECT_LABEL":
                disease_direct[disease].update(counts)
            unique_by_disease[disease].update(alleles)
            unique_lnc_by_disease[disease].update(lnc_alleles)
            for key in alleles:
                allele_which[key][disease].add(accession)
            for vid in chosen.values():
                if vid["credible"] and get_lnc(vid):
                    allele_lnc_genes[vid["allele"]].update(vid["genes"])
            studies_summary[disease].append((accession, ptype, counts, n_regional, bad_pips, n_raw))
            TABLES.setdefault("studies", []).append({
                "PHENOTYPE": disease, "STUDY_ACCESSION": accession, "TRAIT_CLASS": ptype,
                "CS_VARIANTS": counts["n"], "CS_LNCRNA": counts["lnc"],
                "PERCENT_LNCRNA": round(100 * counts["lnc"] / counts["n"], 3) if counts["n"] else None,
                "CS_CODING": counts["coding"], "REGIONAL_VARIANTS": n_regional})
            print(f"{accession:15} {ptype:19} {counts['n']:12,} {counts['lnc']:10,} {fmt(counts['lnc'],counts['n']):>10}  INCLUDED")
            if bad_pips:
                print(f"  WARNING: {bad_pips:,} regional rows had invalid PIP; credible rows must all have valid PIP.")
            print(f"  Regional variants: {n_regional:,} | VEP transcript records streamed: {n_raw:,}", flush=True)

        print("\nInclusion audit:")
        for k, v in sorted(readyness[disease].items()):
            print(f"  {k:38} {v:,}")
        values = [(acc, 100*s['lnc']/s['n']) for acc, pt, s, *_ in studies_summary[disease] if s['n']]
        if values:
            a = [v for _, v in values]
            print("\nAcross included studies (lncRNA percent of credible-set rows):")
            print(f"  Median={statistics.median(a):.3f}% | Minimum={min(a):.3f}% | Maximum={max(a):.3f}%")
            print("  Warning: these are study-level descriptive percentages, not independent samples.")

    print_matrix("ALL INCLUDED STUDIES: 95% CREDIBLE-SET DISTRIBUTIONS", disease_all)
    print_matrix("DIRECT DISEASE LABELS ONLY: SENSITIVITY ANALYSIS", disease_direct)
    print_matrix(f"HIGH-PIP (>= {args.high_pip:g}) VARIANT DISTRIBUTIONS", disease_high, "High-PIP variant-study rows")
    print_loci(disease_loci)

    print("\n" + "=" * 109)
    print("UNIQUE CREDIBLE-SET ALLELES PER DISEASE (DEDUPLICATED ACROSS GWAS STUDIES)")
    print("=" * 109)
    for disease, _ in DISEASES:
        total = len(unique_by_disease[disease])
        n_lnc = len(unique_lnc_by_disease[disease])
        print(f"{disease:15} unique_credible={total:>8,} unique_lncRNA={n_lnc:>8,} overlap={fmt(n_lnc,total)}")
        TABLES.setdefault("unique_alleles", []).append({
            "PHENOTYPE": disease, "UNIQUE_CS_ALLELES": total, "UNIQUE_LNCRNA_ALLELES": n_lnc,
            "PERCENT": round(100 * n_lnc / total, 3) if total else None})

    print("\nPAIRWISE SHARED CREDIBLE-SET ALLELES")
    for i in range(len(DISEASES)):
        for j in range(i+1, len(DISEASES)):
            a = DISEASES[i][0]
            b = DISEASES[j][0]
            shared = unique_by_disease[a] & unique_by_disease[b]
            union = unique_by_disease[a] | unique_by_disease[b]
            lnc_shared = unique_lnc_by_disease[a] & unique_lnc_by_disease[b]
            print(f"  {a} vs {b}: shared={len(shared):,}, Jaccard={fmt(len(shared),len(union))}, lncRNA-shared={len(lnc_shared):,}")
            TABLES.setdefault("pairwise_sharing", []).append({
                "PHENOTYPE_A": a, "PHENOTYPE_B": b, "SHARED_CS_ALLELES": len(shared),
                "JACCARD_PERCENT": round(100 * len(shared) / len(union), 3) if union else None,
                "SHARED_LNCRNA_ALLELES": len(lnc_shared)})

    print(f"\nTOP {args.top} ALLELES IN CREDIBLE SETS OF TWO OR MORE DISEASES")
    cross = [(k, d) for k, d in allele_which.items() if len(d) >= 2]
    cross.sort(key=lambda item: (-len(item[1]), -sum(len(x) for x in item[1].values()), item[0]))
    if not cross:
        print("  No exact allele shared across two or more included disease sets.")
    for key, disease_map in cross[:args.top]:
        parts = [f"{name}({len(studies)} studies)" for name, studies in sorted(disease_map.items())]
        lnc_in = [name for name in disease_map if key in unique_lnc_by_disease[name]]
        names = ",".join(sorted(allele_lnc_genes[key])) or "NA"
        print(f"  {fmt_allele(key):30} | {'; '.join(parts)} | lncRNA in: {','.join(sorted(lnc_in)) or 'none'} | lncRNA genes: {names[:85]}")
        TABLES.setdefault("shared_alleles", []).append({
            "ALLELE": fmt_allele(key), "PHENOTYPES": "; ".join(parts),
            "LNCRNA_IN": ",".join(sorted(lnc_in)), "LNCRNA_GENES": names})

    print("\n" + "=" * 109)
    print("MANUSCRIPT INTERPRETATION AND LIMITATIONS")
    print("=" * 109)
    print("1. Main counts are variant-study records; cohort and locus overlap can inflate pooled counts.")
    print("2. Unique alleles collapse repeated coordinates across GWAS studies, but do not establish replication.")
    print("3. EXACT variant sharing is not LD-based sharing or evidence of shared causal effects.")
    print("4. Coding means VEP protein-coding CDS consequence on a protein_coding transcript.")
    print("5. lncRNA exon/intron/splice annotations are matched to lncRNA BIOTYPE on the same VEP row.")
    print("6. Transcript overlaps are nonexclusive; splice, exon and intron flags can overlap.")
    print("7. Locus classes refer to credible-set variant transcript overlaps, not the full GENCODE gene inventory.")
    print("8. DIRECT_LABEL is a conservative exact-label heuristic; verify source phenotypes manually.")
    print("9. Different numbers of discovery GWAS, sample sizes, QC and fine-mapping resolution limit comparison.")
    print("10. Descriptive percentages do NOT estimate variance explained, causal mediation or enrichment.")
    print("11. A formal disease comparison requires matched/independent cohorts and a justified statistical model.")
    print("12. Studies with no significant loci are selection outcomes, not annotation failures.")
    return {name: pd.DataFrame(rows) for name, rows in TABLES.items()}


def main(argv=None):
    return analysis_main(run, "John Analysis08: cross-disease lncRNA comparison", argv, extra=[
        ("--high-pip", dict(type=float, default=0.50)),
        ("--top", dict(type=int, default=20)),
    ])


if __name__ == "__main__":
    main()
