#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""John Mattick | Analysis 07: within-locus lncRNA enrichment (READ-ONLY).

Null: Conditional on each study's successfully fine-mapped locus and its
number of credible-set variants, the credible labels are exchangeable among
all *annotated* variants examined in the same locus (or locus + positional
window with --bin-kb). This is a within-locus exploratory null, NOT a
LD-aware, genome-wide or causal enrichment test.

Requires completed Step06 and Step07 VEP all-consequences outputs.
Phenotypes and ancestry are passed in; run() returns the result tables.
"""

import argparse
import csv
import gzip
import itertools
import json
import math
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

from gwas2m_lncrna.analysis.helpers import analysis_main

# Set by run(): [(display name, phenotype folder), ...] and the ancestry folder name.
DISEASES = []
ANCESTRY = "european"
DONE = {"COMPLETE", "SKIPPED_ALREADY_COMPLETE"}
FEATURES = ("LNC_TRANSCRIPT", "LNC_EXON", "LNC_INTRON", "LNC_SPLICE")
IN_TRANSCRIPT = {
    "intron_variant", "non_coding_transcript_variant",
    "non_coding_transcript_exon_variant", "splice_acceptor_variant",
    "splice_donor_variant", "splice_region_variant",
    "splice_donor_region_variant", "splice_donor_5th_base_variant",
    "splice_polypyrimidine_tract_variant", "transcript_ablation",
    "transcript_amplification", "exon_loss_variant",
}
SPLICE = {
    "splice_acceptor_variant", "splice_donor_variant", "splice_region_variant",
    "splice_donor_region_variant", "splice_donor_5th_base_variant",
    "splice_polypyrimidine_tract_variant",
}
try:
    csv.field_size_limit(sys.maxsize)
except OverflowError:
    csv.field_size_limit(2**31 - 1)


def read_json(p):
    try:
        with p.open(encoding="utf-8") as f:
            v = json.load(f)
        return v if isinstance(v, dict) else {}
    except (OSError, ValueError):
        return {}


def iter_tsv(p):
    with p.open(encoding="utf-8-sig", newline="") as f:
        yield from csv.DictReader(f, delimiter="\t")


def stage(root, slug, accession, name):
    p = root / "pipeline_audit" / slug / ANCESTRY / "status" / accession / (name + ".json")
    return str(read_json(p).get("status", "NOT_RECORDED")).upper()


def truth(x):
    return str(x).strip().lower() in {"true", "t", "1", "yes", "y"}


def split_terms(x):
    return {q.strip().lower() for q in re.split(r"[,;&|]", str(x))
            if q.strip().lower() not in {"", "nan", "none", "na", ".", "-"}}


def number(x):
    try:
        z = float(x)
        return z if math.isfinite(z) else None
    except (ValueError, TypeError):
        return None


def chrom_pos(row):
    ch = str(row.get("CHR", "")).strip().lower().removeprefix("chr")
    p = number(row.get("REFERENCE_POS"))
    if not ch or p is None or int(p) < 1:
        return None
    return ch, int(p)


def percent(n, total):
    return f"{100.0 * n / total:.3f}%" if total else "N/A"


def feature_flags(group, legacy=False):
    """Flags require lncRNA BIOTYPE and consequence on the SAME VEP row."""
    flags = set()
    any_annotation = False
    accepted = {"lncrna", "lincrna"}
    if legacy:
        accepted |= {"antisense", "processed_transcript", "sense_intronic",
                     "sense_overlapping", "3prime_overlapping_ncrna",
                     "bidirectional_promoter_lncrna", "macro_lncrna"}

    for row in group:
        cons = split_terms(row.get("Consequence", ""))
        if cons:
            any_annotation = True
        bio = str(row.get("BIOTYPE", "")).strip().lower()
        if bio not in accepted or not (cons & IN_TRANSCRIPT):
            continue  # In particular, exclude upstream/downstream-only annotations.
        flags.add("LNC_TRANSCRIPT")
        ex = str(row.get("EXON", "")).strip()
        intr = str(row.get("INTRON", "")).strip()
        if "non_coding_transcript_exon_variant" in cons or re.match(r"^\d+(?:/\d+)?$", ex):
            flags.add("LNC_EXON")
        if "intron_variant" in cons or re.match(r"^\d+(?:/\d+)?$", intr):
            flags.add("LNC_INTRON")
        if cons & SPLICE:
            flags.add("LNC_SPLICE")

    return flags, any_annotation


def locate_files(root, slug, accession):
    s6 = root / "06_finemapping" / slug / ANCESTRY / accession
    s7 = root / "07_annotation" / slug / ANCESTRY / accession
    return {
        "fine": s6 / (accession + "_finemapped_variants.tsv.gz"),
        "fm": s6 / "finemapping_summary.json",
        "anno": s7 / "annotation_summary.json",
        "all": s7 / (accession + "_VEP_all_consequences.tsv.gz"),
        "failure": s7 / "ANNOTATION_FAILED.txt",
    }


def eligibility(root, slug, accession):
    s6 = stage(root, slug, accession, "Step06_FineMapping")
    s7 = stage(root, slug, accession, "Step07_VEP")
    files = locate_files(root, slug, accession)
    fm = read_json(files["fm"])
    if s6 not in DONE:
        return None, "STEP06_" + s6
    if fm and (str(fm.get("STATUS", "")).upper() == "NO_SIGNIFICANT_LOCI"
               or number(fm.get("N_LOCI")) == 0):
        return None, "NO_SIGNIFICANT_LOCI"
    if s7 not in DONE:
        return None, "STEP07_" + s7
    if files["failure"].exists() or not files["fine"].is_file() or not files["all"].is_file():
        return None, "MISSING_OR_FAILED_SOURCE"
    ann = read_json(files["anno"])
    if (str(ann.get("STATUS", "")).upper() != "COMPLETE"
            or str(ann.get("ANNOTATION_SCOPE", "")).upper() != "ALL_FINEMAPPED_VARIANTS"):
        return None, "INVALID_ANNOTATION_SCOPE"
    try:
        current = (int(ann["FINEMAPPED_INPUT_SIZE"]) == files["fine"].stat().st_size
                   and int(ann["FINEMAPPED_INPUT_MTIME_NS"]) == files["fine"].stat().st_mtime_ns)
    except (ValueError, KeyError, TypeError):
        current = False
    if not current:
        return None, "STALE_STEP07_ANNOTATION"
    return (files, ann), "READY"


def analyze_study(files, ann, bin_kb, legacy):
    """Stream VEP once; accumulate four 2x2 tables per locus or locus+window."""
    # Feature counts per stratum: [N=all, n=credible, K=all positive, k=credible positive].
    counts = defaultdict(lambda: {f: [0, 0, 0, 0] for f in FEATURES})
    stats = Counter()
    prev_id = None
    expected = ann.get("N_SELECTED_VARIANTS")
    try:
        expected = int(expected) if expected is not None else None
    except (ValueError, TypeError):
        expected = None
    with gzip.open(files["all"], "rt", encoding="utf-8", newline="") as h:
        reader = csv.DictReader(h, delimiter="\t")
        required = {"VEP_ID", "LOCUS_ID", "CHR", "REFERENCE_POS",
                    "Consequence", "BIOTYPE", "IN_95_CREDIBLE_SET"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError("Missing Step07 columns: " + ", ".join(sorted(missing)))
        for vep, group in itertools.groupby(reader, key=lambda r: r["VEP_ID"]):
            if not vep or (prev_id is not None and vep <= prev_id):
                raise ValueError("VEP_ID ordering/uniqueness validation failed")
            prev_id = vep
            items = list(group)
            stats["INPUT_VARIANTS"] += 1
            first = items[0]
            loc = str(first.get("LOCUS_ID", "")).strip()
            chrom_and_pos = chrom_pos(first)
            if not loc or chrom_and_pos is None:
                stats["UNUSABLE_LOCATION"] += 1
                continue
            if any(str(x.get("LOCUS_ID", "")).strip() != loc for x in items):
                raise ValueError("Inconsistent LOCUS_ID within a VEP_ID")
            cred = truth(first.get("IN_95_CREDIBLE_SET"))
            pip = number(first.get("PIP"))
            if pip is None or not (0 <= pip <= 1):
                raise ValueError(f"Invalid SuSiE PIP in {vep}: {first.get('PIP')}")
            flags, annotated = feature_flags(items, legacy)
            if not annotated:
                stats["UNANNOTATED"] += 1
                if cred:
                    stats["UNANNOTATED_CS"] += 1
                continue
            stats["ANNOTATED"] += 1
            stats["CREDIBLE"] += int(cred)
            chrom, pos = chrom_and_pos
            binkey = (loc, chrom, (pos - 1) // (bin_kb * 1000)) if bin_kb else (loc,)
            for feat in FEATURES:
                x = counts[binkey][feat]
                x[0] += 1  # all annotated variants
                x[1] += int(cred)
                x[2] += int(feat in flags)  # background feature count
                x[3] += int(cred and feat in flags)  # observed credible feature count
    if expected is not None and stats["INPUT_VARIANTS"] != expected:
        raise ValueError("VEP_ID count differs from annotation_summary: "
                         f"{stats['INPUT_VARIANTS']} vs {expected}")
    cs_expected = ann.get("N_SELECTED_CREDIBLE_SET")
    if cs_expected is not None:
        try:
            cs_expected = int(cs_expected)
        except (ValueError, TypeError):
            cs_expected = None
        if cs_expected is not None and stats["CREDIBLE"] + stats["UNANNOTATED_CS"] != cs_expected:
            raise ValueError("Credible-set count disagrees with Step07 annotation_summary")
    if stats["ANNOTATED"] == 0:
        raise ValueError("No valid annotated variants")
    return counts, stats


def calculate_null(strata, feat, simulations, seed):
    """One-sided Monte Carlo conditional hypergeometric enrichment/depletion."""
    rng = np.random.default_rng(seed)
    rows = [vals[feat] for vals in strata.values() if vals[feat][1] > 0]
    if not rows:
        return None
    N = sum(x[0] for x in rows)
    n = sum(x[1] for x in rows)
    K = sum(x[2] for x in rows)
    observed = sum(x[3] for x in rows)
    expected = sum(x[1] * x[2] / x[0] for x in rows)
    variance = sum((x[1] * x[2] * (x[0] - x[2]) * (x[0] - x[1])) /
                   (x[0] ** 2 * (x[0] - 1)) for x in rows if x[0] > 1)
    valid = [x for x in rows if 0 < x[2] < x[0] and 0 < x[1] < x[0]]
    if not valid:
        return dict(N=N, n=n, K=K, observed=observed, expected=expected,
                    fold=observed / expected if expected else None,
                    p_over=None, p_under=None, null_lo=observed,
                    null_hi=observed, n_strata=len(rows), n_informative=0,
                    variance=variance)
    simulated = np.zeros(simulations, dtype=np.int64)
    constant = sum((x[1] if x[2] == x[0] else x[2] if x[1] == x[0] else 0)
                   for x in rows)
    simulated += constant
    for total, cred, positives, _ in valid:
        simulated += rng.hypergeometric(positives, total - positives, cred,
                                        size=simulations)
    p_over = (1 + np.count_nonzero(simulated >= observed)) / (simulations + 1)
    p_under = (1 + np.count_nonzero(simulated <= observed)) / (simulations + 1)
    lo, hi = np.quantile(simulated, [0.025, 0.975])
    return dict(N=N, n=n, K=K, observed=observed, expected=expected,
                fold=observed / expected if expected else None,
                p_over=p_over, p_under=p_under, null_lo=lo, null_hi=hi,
                n_strata=len(rows), n_informative=len(valid), variance=variance)


def bh(pvalues):
    """Benjamini-Hochberg correction over all valid tests in this run."""
    indexed = sorted(((p, i) for i, p in enumerate(pvalues) if p is not None))
    q = [None] * len(pvalues)
    if not indexed:
        return q
    running = 1.0
    m = len(indexed)
    for j in range(m - 1, -1, -1):
        p, idx = indexed[j]
        running = min(running, p * m / (j + 1))
        q[idx] = running
    return q


def fmt(x, digits=4):
    if x is None:
        return "NA"
    return f"{x:.{digits}g}"


def print_test(name, stat, qval):
    if stat is None:
        print(f"  {name:19} NO_CREDIBLE_VARIANTS")
        return
    obs = stat["observed"]
    cred = stat["n"]
    N = stat["N"]
    K = stat["K"]
    print(f"  {name:19} CS {obs:>7,}/{cred:<9,} ({percent(obs, cred):>9}) | "
          f"regional {K:>9,}/{N:<10,} ({percent(K, N):>9}) | "
          f"expected={stat['expected']:.2f} fold={fmt(stat['fold']):>7} | "
          f"p_enrich={fmt(stat['p_over'], 5):>9} q_BH={fmt(qval, 5):>9} | "
          f"null95=[{stat['null_lo']:.0f},{stat['null_hi']:.0f}] "
          f"strata={stat['n_informative']}/{stat['n_strata']}")


def run(root, diseases, ancestry="european", permutations=10000, seed=42, bin_kb=0, legacy_lncrna=False,
        min_annotation_coverage=0.90, **_):
    """Run Analysis 07 for [(display name, folder slug), ...]; return its tables."""
    global DISEASES, ANCESTRY
    DISEASES, ANCESTRY = list(diseases), ancestry
    args = argparse.Namespace(permutations=int(permutations), seed=int(seed), bin_kb=int(bin_kb),
                              legacy_lncrna=bool(legacy_lncrna),
                              min_annotation_coverage=float(min_annotation_coverage))
    if args.permutations < 100 or args.bin_kb < 0 or not (0 <= args.min_annotation_coverage <= 1):
        raise ValueError("permutations>=100; bin-kb>=0; coverage between 0 and 1")
    root = Path(root).expanduser().resolve()
    study_rows = []

    print("=" * 145)
    print("JOHN MATTICK - ANALYSIS 07 | LNCRNA ENRICHMENT | READ ONLY")
    print("=" * 145)
    print(f"Root: {root} | simulations={args.permutations:,} | seed={args.seed} | "
          f"strata={'LOCI' if args.bin_kb == 0 else f'LOCI + {args.bin_kb}kb bins'}")
    print("Null: credible-set assignments randomized within successfully mapped loci/strata.")
    print("NOT a genome-wide, LD-matched or causal enrichment test.")
    print("Input files are never modified.\n")

    reports = []
    details = defaultdict(list)
    for disease, slug in DISEASES:
        manifest = root / "01_gwas_catalog" / slug / ANCESTRY / "GWAS_download_manifest.tsv"
        if not manifest.is_file():
            manifest = root / "pipeline_audit" / slug / ANCESTRY / "study_manifest.tsv"
        print("\n" + "=" * 145)
        print(disease.upper())
        print("=" * 145)
        if not manifest.is_file():
            print("No manifest available")
            continue
        manifest_rows = list(iter_tsv(manifest))
        status_count = Counter()
        for entry in manifest_rows:
            study = str(entry.get("STUDY_ACCESSION", "")).strip()
            if not study:
                continue
            ready, reason = eligibility(root, slug, study)
            if ready is None:
                status_count[reason] += 1
                study_rows.append({"PHENOTYPE": disease, "STUDY_ACCESSION": study, "STATUS": reason})
                print(f"{study:<15} SKIPPED: {reason}")
                continue
            files, anno = ready
            try:
                strata, counts = analyze_study(files, anno, args.bin_kb, args.legacy_lncrna)
            except (OSError, UnicodeError, EOFError, csv.Error, ValueError) as exc:
                status_count["READ_ERROR"] += 1
                study_rows.append({"PHENOTYPE": disease, "STUDY_ACCESSION": study, "STATUS": "READ_ERROR",
                                   "DETAIL": str(exc)})
                print(f"{study:<15} INVALID: {type(exc).__name__}: {exc}")
                continue
            coverage = counts["ANNOTATED"] / max(1, counts["INPUT_VARIANTS"])
            study_rows.append({"PHENOTYPE": disease, "STUDY_ACCESSION": study,
                               "STATUS": ("INCLUDED" if coverage >= args.min_annotation_coverage
                                          else "LOW_ANNOTATION_COVERAGE"),
                               "VARIANTS": counts["INPUT_VARIANTS"], "ANNOTATED": counts["ANNOTATED"],
                               "ANNOTATION_COVERAGE": round(coverage, 4), "CREDIBLE": counts["CREDIBLE"],
                               "STRATA": len(strata)})
            if coverage < args.min_annotation_coverage:
                status_count["LOW_ANNOTATION_COVERAGE"] += 1
                print(f"{study:<15} EXCLUDED: annotation coverage {percent(counts['ANNOTATED'], counts['INPUT_VARIANTS'])}")
                continue
            status_count["INCLUDED"] += 1
            print(f"{study:<15} INCLUDED | variants={counts['INPUT_VARIANTS']:,} "
                  f"annotated={counts['ANNOTATED']:,} ({percent(counts['ANNOTATED'], counts['INPUT_VARIANTS'])}) "
                  f"credible={counts['CREDIBLE']:,} strata={len(strata):,}", flush=True)
            for i, feat in enumerate(FEATURES):
                # Different feature tests get deterministic independent streams.
                stat = calculate_null(strata, feat, args.permutations, args.seed + i + len(reports) * 7)
                reports.append({"disease": disease, "study": study, "feature": feat, "stat": stat})
                details[disease].append((study, feat, stat))
        print("\nSTUDY ACCOUNTING")
        for reason, n in sorted(status_count.items()):
            print(f"  {reason:<42} {n:>5,}")

    q_values = bh([r["stat"]["p_over"] if r["stat"] is not None else None for r in reports])
    for row, q in zip(reports, q_values):
        row["q"] = q

    for disease, _ in DISEASES:
        if disease not in details:
            continue
        print("\n" + "=" * 145)
        print(f"ENRICHMENT RESULTS: {disease.upper()}")
        print("=" * 145)
        entries = [r for r in reports if r["disease"] == disease]
        for study in dict.fromkeys(r["study"] for r in entries):
            print("\nSTUDY:", study)
            for r in entries:
                if r["study"] == study:
                    print_test(r["feature"], r["stat"], r["q"])

        print("\nDESCRIPTIVE SUM ACROSS STUDIES (NO COMBINED P-VALUE; COHORTS MAY OVERLAP)")
        for feat in FEATURES:
            selected = [r["stat"] for r in entries if r["feature"] == feat and r["stat"]]
            if not selected:
                continue
            N = sum(x["N"] for x in selected)
            n = sum(x["n"] for x in selected)
            K = sum(x["K"] for x in selected)
            o = sum(x["observed"] for x in selected)
            e = sum(x["expected"] for x in selected)
            print(f"  {feat:19} observed={o:,}/{n:,} ({percent(o, n)}), "
                  f"background={K:,}/{N:,} ({percent(K, N)}), "
                  f"expected={e:.2f}, fold={fmt(o / e if e else None)}")

    print("\n" + "=" * 145)
    print("SCIENTIFIC INTERPRETATION AND LIMITATIONS")
    print("=" * 145)
    print("1. Each test compares 95% credible-set membership against all annotated variants in the SAME fine-mapped locus.")
    print("2. Per-stratum sampling preserves the number of credible variants and lncRNA-overlapping variants.")
    print("3. A variant counts once per study and can overlap several lncRNA transcripts; flags are non-exclusive.")
    print("4. The lncRNA biotype and consequence must appear on the same VEP transcript row; nearby-only calls are excluded.")
    print("5. Monte Carlo p_enrich is one-sided; q_BH adjusts all testable feature-by-study enrichment tests in this run.")
    print("6. Null95 is a randomization interval for the COUNT under the null; it is not a confidence interval for fold enrichment.")
    print("7. This null does NOT preserve LD, allele frequency or PIP ascertainment; locus-only stratification does not match local gene density.")
    print("8. Cohort-overlapping GWAS are NOT independent replicates; combined figures are descriptive, not inferential.")
    print("9. For a publication-grade LD-aware null, add LD-/MAF-matched or genomic block resampling as a separate validation.")
    print("10. Transcript overlap, even when enriched, cannot establish causal mediation through lncRNAs.")
    tests = []
    for r in reports:
        s = r["stat"]
        row = {"PHENOTYPE": r["disease"], "STUDY_ACCESSION": r["study"], "FEATURE": r["feature"]}
        if s is None:
            row["STATUS"] = "NO_CREDIBLE_VARIANTS"
        else:
            row.update({"STATUS": "TESTED", "CS_WITH_FEATURE": s["observed"], "CS_VARIANTS": s["n"],
                        "REGIONAL_WITH_FEATURE": s["K"], "REGIONAL_VARIANTS": s["N"],
                        "EXPECTED": s["expected"], "FOLD": s["fold"], "P_ENRICH": s["p_over"],
                        "Q_BH": r["q"], "P_DEPLETE": s["p_under"], "NULL95_LOW": s["null_lo"],
                        "NULL95_HIGH": s["null_hi"], "STRATA_INFORMATIVE": s["n_informative"],
                        "STRATA": s["n_strata"]})
        tests.append(row)
    sums = []
    for disease, _ in DISEASES:
        for feat in FEATURES:
            selected = [r["stat"] for r in reports if r["disease"] == disease and r["feature"] == feat and r["stat"]]
            if selected:
                o, e = sum(x["observed"] for x in selected), sum(x["expected"] for x in selected)
                sums.append({"PHENOTYPE": disease, "FEATURE": feat, "STUDIES": len(selected), "OBSERVED": o,
                             "CS_VARIANTS": sum(x["n"] for x in selected),
                             "BACKGROUND": sum(x["K"] for x in selected),
                             "REGIONAL_VARIANTS": sum(x["N"] for x in selected), "EXPECTED": e,
                             "DESCRIPTIVE_FOLD": o / e if e else None})
    return {"studies": pd.DataFrame(study_rows), "tests": pd.DataFrame(tests),
            "descriptive_sums": pd.DataFrame(sums)}


def main(argv=None):
    return analysis_main(run, "John Analysis07: within-locus lncRNA enrichment", argv, extra=[
        ("--permutations", dict(type=int, default=10000)),
        ("--seed", dict(type=int, default=42)),
        ("--bin-kb", dict(type=int, default=0,
                          help="0=stratify by locus only; 250=stratify by locus + 250kb positional bin")),
        ("--legacy-lncrna", dict(action="store_true",
                                 help="Exploratory: include older noncoding transcript biotype names")),
        ("--min-annotation-coverage", dict(type=float, default=0.90)),
    ])


if __name__ == "__main__":
    main()
