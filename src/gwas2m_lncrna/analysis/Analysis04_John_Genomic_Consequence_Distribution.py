 
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
ANALYSIS 04 - JOHN MATTICK
Genomic consequence distribution of fine-mapped GWAS variants.

Phenotypes and ancestry are passed in (run() / --phenotype, --ancestry).

INPUT:
    Completed Step06 fine-mapping
    Completed Step07 VEP annotation

OUTPUT:
    Printed to the terminal; run() also returns the result tables.

Primary analysis:
    95% credible-set variants.

Additional analyses:
    All fine-mapped regional variants.
    High-PIP variants.
"""

import argparse
import csv
import json
import math
import re

from collections import Counter
from pathlib import Path

import pandas as pd

from gwas2m_lncrna.analysis.helpers import analysis_main


# Set by run(): [(display name, phenotype folder), ...] and the ancestry folder name.
DISEASES = []
ANCESTRY = "european"
TABLES = {}

DONE = {"COMPLETE", "SKIPPED_ALREADY_COMPLETE"}

# Multiple consequences may apply to the same variant.
# This is the priority order for one mutually exclusive
# primary category per variant.

ORDER = [
    "ESSENTIAL_SPLICE_SITE",
    "PROTEIN_CODING_CDS",
    "UTR_5_OR_3",
    "NONCODING_TRANSCRIPT_EXON",
    "SPLICE_REGION_OTHER",
    "INTRONIC",
    "REGULATORY",
    "TRANSCRIPT_ABLATION_AMPLIFICATION",
    "OTHER_NONCODING_TRANSCRIPT",
    "UPSTREAM",
    "DOWNSTREAM",
    "INTERGENIC",
    "OTHER_ANNOTATED",
    "UNANNOTATED",
]

CODING = {
    "coding_sequence_variant",
    "missense_variant",
    "synonymous_variant",
    "stop_gained",
    "stop_lost",
    "start_lost",
    "start_retained_variant",
    "stop_retained_variant",
    "frameshift_variant",
    "inframe_insertion",
    "inframe_deletion",
    "protein_altering_variant",
    "incomplete_terminal_codon_variant",
}

SPLICE_CORE = {
    "splice_acceptor_variant",
    "splice_donor_variant",
}

SPLICE_OTHER = {
    "splice_region_variant",
    "splice_donor_region_variant",
    "splice_donor_5th_base_variant",
    "splice_polypyrimidine_tract_variant",
}

UTR = {
    "5_prime_utr_variant",
    "3_prime_utr_variant",
}

REGULATORY = {
    "regulatory_region_variant",
    "tf_binding_site_variant",
    "tfbs_ablation",
    "tfbs_amplification",
    "regulatory_region_ablation",
    "regulatory_region_amplification",
}

POPS = [
    "CREDIBLE95",
    "ALL_REGIONAL",
    "HIGH_PIP",
]


def read_json(path):
    try:
        with path.open(encoding="utf-8") as f:
            x = json.load(f)
            return x if isinstance(x, dict) else {}
    except (OSError, ValueError):
        return {}


def read_status(audit, study, name):
    path = (
        audit / "status" / study / (name + ".json")
    )
    return str(
        read_json(path).get("status", "NOT_RECORDED")
    ).upper()


def truth(x):
    return str(x).strip().lower() in {
        "true", "t", "yes", "1"
    }


def pip_number(x):
    try:
        n = float(x)
        return (
            n if math.isfinite(n) and 0 <= n <= 1
            else None
        )
    except (ValueError, TypeError):
        return None


def tokens(s):
    return {
        x.strip().lower()
        for x in re.split(r"[,;&|]", str(s))
        if x.strip()
        and x.strip().lower()
        not in {"nan", "none", ".", "na"}
    }


def category(terms):
    """Assign one primary genomic category."""

    if not terms:
        return "UNANNOTATED"

    if terms & SPLICE_CORE:
        return "ESSENTIAL_SPLICE_SITE"

    if terms & CODING:
        return "PROTEIN_CODING_CDS"

    if terms & UTR:
        return "UTR_5_OR_3"

    if "non_coding_transcript_exon_variant" in terms:
        return "NONCODING_TRANSCRIPT_EXON"

    if terms & SPLICE_OTHER:
        return "SPLICE_REGION_OTHER"

    if "intron_variant" in terms:
        return "INTRONIC"

    if terms & REGULATORY:
        return "REGULATORY"

    if terms & {
        "transcript_ablation",
        "transcript_amplification",
    }:
        return "TRANSCRIPT_ABLATION_AMPLIFICATION"

    if "non_coding_transcript_variant" in terms:
        return "OTHER_NONCODING_TRANSCRIPT"

    if "upstream_gene_variant" in terms:
        return "UPSTREAM"

    if "downstream_gene_variant" in terms:
        return "DOWNSTREAM"

    if "intergenic_variant" in terms:
        return "INTERGENIC"

    return "OTHER_ANNOTATED"


def blank():
    return {
        "categories": Counter(),
        "terms": Counter(),
        "total": 0,
        "lncrna": 0,
        "protein_gene": 0,
        "both": 0,
    }


def new_groups():
    return {
        name: blank()
        for name in POPS
    }


def add(dst, src):
    dst["categories"].update(src["categories"])
    dst["terms"].update(src["terms"])

    for k in (
        "total", "lncrna", "protein_gene", "both"
    ):
        dst[k] += src[k]


def count_one(stats, label, terms, lnc, protein):
    stats["total"] += 1
    stats["categories"][label] += 1

    # Raw VEP terms are counted once per variant.
    stats["terms"].update(terms)

    stats["lncrna"] += int(lnc)
    stats["protein_gene"] += int(protein)
    stats["both"] += int(lnc and protein)


def percent(n, total):
    if not total:
        return "N/A"

    return f"{100 * n / total:.3f}%"


def print_result(title, stats, nterms):
    total = stats["total"]

    print("\n" + "=" * 90)
    print(title)
    print("=" * 90)

    print(f"Variants analyzed: {total:,}\n")

    print(
        f"{'GENOMIC CATEGORY':40}"
        f"{'COUNT':>15}"
        f"{'PERCENTAGE':>15}"
    )

    print("-" * 70)

    for key in ORDER:
        n = stats["categories"][key]

        print(
            f"{key:40}"
            f"{n:>15,}"
            f"{percent(n, total):>15}"
        )

    print("-" * 70)

    print(
        f"{'TOTAL':40}"
        f"{total:>15,}"
        f"{percent(total, total):>15}"
    )

    print("\nADDITIONAL TRANSCRIPT OVERLAPS")

    for label, key in [
        ("lncRNA transcript", "lncrna"),
        ("Protein-coding gene", "protein_gene"),
        ("Both biotypes", "both"),
    ]:
        n = stats[key]

        print(
            f"{label:35}"
            f"{n:>12,} "
            f"({percent(n, total)})"
        )

    if nterms:
        print("\nMOST COMMON VEP CONSEQUENCE TERMS")

        print(
            "These terms can overlap; percentages "
            "do not need to add to 100%."
        )

        for name, n in stats["terms"].most_common(nterms):
            print(
                f"{name:45}"
                f"{n:>12,} "
                f"{percent(n, total)}"
            )


def study_data(root, slug, study, high_pip):
    base = (
        root / "06_finemapping"
        / slug / ANCESTRY / study
    )

    anno = (
        root / "07_annotation"
        / slug / ANCESTRY / study
    )

    audit = (
        root / "pipeline_audit"
        / slug / ANCESTRY
    )

    s06 = read_status(
        audit, study, "Step06_FineMapping"
    )

    s07 = read_status(
        audit, study, "Step07_VEP"
    )

    if s06 not in DONE:
        return None, "S06_" + s06

    fm = read_json(
        base / "finemapping_summary.json"
    )

    if (
        str(fm.get("STATUS", "")).upper()
        == "NO_SIGNIFICANT_LOCI"
        or (fm and fm.get("N_LOCI") == 0)
    ):
        return None, "NO_SIGNIFICANT_LOCI"

    if s07 not in DONE:
        return None, "S07_" + s07

    meta = read_json(
        anno / "annotation_summary.json"
    )

    fine_file = (
        base / f"{study}_finemapped_variants.tsv.gz"
    )

    path = (
        anno / f"{study}_VEP_variant_summary.tsv"
    )

    if (
        not path.is_file()
        or not fine_file.is_file()
    ):
        return None, "MISSING_SOURCE_FILE"

    if (
        meta.get("STATUS") != "COMPLETE"
        or meta.get("ANNOTATION_SCOPE")
        != "ALL_FINEMAPPED_VARIANTS"
        or (anno / "ANNOTATION_FAILED.txt").exists()
    ):
        return None, "INVALID_VEP_STATUS"

    # Check Step07 was generated from current Step06.
    if (
        meta.get("FINEMAPPED_INPUT_SIZE")
        != fine_file.stat().st_size
        or meta.get("FINEMAPPED_INPUT_MTIME_NS")
        != fine_file.stat().st_mtime_ns
    ):
        return None, "STALE_VEP_ANNOTATION"

    data = new_groups()

    required = {
        "VEP_ID",
        "PIP",
        "IN_95_CREDIBLE_SET",
        "ALL_CONSEQUENCES",
        "ALL_BIOTYPES",
    }

    try:
        with path.open(
            encoding="utf-8", newline=""
        ) as f:

            reader = csv.DictReader(
                f, delimiter="\t"
            )

            missing = required - set(
                reader.fieldnames or []
            )

            if missing:
                return (
                    None,
                    "MISSING_COLUMNS:"
                    + ",".join(sorted(missing)),
                )

            for row in reader:

                consequences = tokens(
                    row.get("ALL_CONSEQUENCES", "")
                )

                biotypes = tokens(
                    row.get("ALL_BIOTYPES", "")
                )

                kind = category(consequences)

                lncrna = bool(
                    biotypes & {"lncrna", "lincrna"}
                )

                protein = (
                    "protein_coding" in biotypes
                )

                pip = pip_number(
                    row.get("PIP")
                )

                groups = ["ALL_REGIONAL"]

                if truth(
                    row.get("IN_95_CREDIBLE_SET")
                ):
                    groups.append("CREDIBLE95")

                if (
                    pip is not None
                    and pip >= high_pip
                ):
                    groups.append("HIGH_PIP")

                for p in groups:
                    count_one(
                        data[p],
                        kind,
                        consequences,
                        lncrna,
                        protein,
                    )

    except (
        OSError,
        UnicodeError,
        csv.Error,
        ValueError,
    ) as exc:
        return None, "READ_ERROR:" + str(exc)

    if data["ALL_REGIONAL"]["total"] == 0:
        return None, "EMPTY_ANNOTATION"

    return data, "OK"


def stats_rows(disease, population, stats):
    """Table rows for one population: primary category counts and transcript overlaps."""
    total = stats["total"]
    rows = [{"PHENOTYPE": disease, "POPULATION": population, "CATEGORY": key,
             "N": stats["categories"][key], "PERCENT": round(100 * stats["categories"][key] / total, 3) if total else None,
             "DENOMINATOR": total} for key in ORDER]
    for label, key in (("lncRNA transcript (any biotype row)", "lncrna"),
                       ("Protein-coding gene", "protein_gene"), ("Both biotypes", "both")):
        rows.append({"PHENOTYPE": disease, "POPULATION": population, "CATEGORY": "OVERLAP: " + label,
                     "N": stats[key], "PERCENT": round(100 * stats[key] / total, 3) if total else None,
                     "DENOMINATOR": total})
    return rows


def run(root, diseases, ancestry="european", high_pip=0.50, top_terms=12, **_):
    """Run Analysis 04 for [(display name, folder slug), ...]; return its tables."""
    global DISEASES, ANCESTRY
    DISEASES, ANCESTRY = list(diseases), ancestry
    TABLES.clear()
    args = argparse.Namespace(high_pip=float(high_pip), top_terms=int(top_terms))
    if not 0 <= args.high_pip <= 1 or args.top_terms < 0:
        raise ValueError("high_pip must be 0..1; top_terms must be >= 0")
    root = Path(root).resolve()

    print("=" * 105)
    print("JOHN MATTICK - ANALYSIS 04")
    print("GENOMIC CONSEQUENCE DISTRIBUTION")
    print("READ ONLY")
    print("=" * 105)

    combined = new_groups()
    disease_results = {}

    for disease, slug in DISEASES:

        print("\n" + "=" * 105)
        print(disease.upper())
        print("=" * 105)

        manifest = (
            root / "01_gwas_catalog"
            / slug / ANCESTRY
            / "GWAS_download_manifest.tsv"
        )

        if not manifest.is_file():
            manifest = (
                root / "pipeline_audit"
                / slug / ANCESTRY
                / "study_manifest.tsv"
            )

        if not manifest.is_file():
            print("No study manifest; skipping.")
            continue

        with manifest.open(
            encoding="utf-8-sig",
            newline=""
        ) as f:
            studies = list(
                csv.DictReader(f, delimiter="\t")
            )

        sums = new_groups()
        reasons = Counter()

        print(f"\nSelected studies: {len(studies)}")

        print(
            f"{'STUDY':15}"
            f"{'STATUS':34}"
            f"{'CS VARIANTS':>13}"
            f"{'ALL REGIONAL':>14}"
        )

        for meta in studies:

            study = str(
                meta.get("STUDY_ACCESSION", "")
            ).strip()

            if not study:
                continue

            result, decision = study_data(
                root, slug, study, args.high_pip
            )

            reasons[decision] += 1
            TABLES.setdefault("studies", []).append({
                "PHENOTYPE": disease, "STUDY_ACCESSION": study, "DECISION": decision,
                "CS_VARIANTS": result["CREDIBLE95"]["total"] if result else None,
                "ALL_REGIONAL_VARIANTS": result["ALL_REGIONAL"]["total"] if result else None,
            })

            if result is None:
                print(
                    f"{study:15}"
                    f"{decision:34}"
                    f"{'-':>13}"
                    f"{'-':>14}"
                )
                continue

            for pop in POPS:
                add(
                    sums[pop],
                    result[pop]
                )

                add(
                    combined[pop],
                    result[pop]
                )

            print(
                f"{study:15}"
                f"{'INCLUDED_STEP07_COMPLETE':34}"
                f"{result['CREDIBLE95']['total']:13,}"
                f"{result['ALL_REGIONAL']['total']:14,}",
                flush=True,
            )

        # ---------------------------------------
        # STUDY INCLUSION
        # ---------------------------------------

        print("\nSTUDY INCLUSION / EXCLUSION")

        for reason, n in sorted(reasons.items()):
            print(
                f"  {reason:43} {n:5,}"
            )

        # ---------------------------------------
        # PRIMARY RESULTS
        # ---------------------------------------

        print_result(
            "PRIMARY: 95% CREDIBLE-SET VARIANTS",
            sums["CREDIBLE95"],
            args.top_terms,
        )

        # ---------------------------------------
        # COMPARISON
        # ---------------------------------------

        print_result(
            "COMPARISON: ALL REGIONAL VARIANTS",
            sums["ALL_REGIONAL"],
            0,
        )

        # ---------------------------------------
        # HIGH-PIP SENSITIVITY
        # ---------------------------------------

        print_result(
            f"SENSITIVITY: PIP >= {args.high_pip:g}",
            sums["HIGH_PIP"],
            0,
        )

        disease_results[disease] = sums["CREDIBLE95"]
        for pop in POPS:
            TABLES.setdefault("categories", []).extend(stats_rows(disease, pop, sums[pop]))
        TABLES.setdefault("top_terms_cs95", []).extend(
            {"PHENOTYPE": disease, "VEP_TERM": name, "N": n,
             "PERCENT": round(100 * n / sums["CREDIBLE95"]["total"], 3) if sums["CREDIBLE95"]["total"] else None}
            for name, n in sums["CREDIBLE95"]["terms"].most_common(args.top_terms)
        )

    # -------------------------------------------
    # CROSS-PHENOTYPE COMPARISON
    # -------------------------------------------

    print("\n" + "=" * 105)
    print("CROSS-PHENOTYPE COMPARISON")
    print("95% CREDIBLE-SET VARIANTS")
    print("=" * 105)

    print(
        f"{'GENOMIC CATEGORY':37}"
        + "".join(f"{(name[:13] + ' %'):>17}" for name, _ in DISEASES)
    )

    for cat in ORDER:

        values = []

        for disease, _ in DISEASES:

            d = disease_results.get(
                disease, blank()
            )

            values.append(
                percent(
                    d["categories"][cat],
                    d["total"],
                )
            )

        print(
            f"{cat:37}"
            + "".join(f"{v:>17}" for v in values)
        )

    # -------------------------------------------
    # COMBINED DESCRIPTIVE COUNTS
    # -------------------------------------------

    print_result(
        "POOLED CREDIBLE-SET ROWS - ALL PHENOTYPES",
        combined["CREDIBLE95"],
        0,
    )

    print("\nSCIENTIFIC INTERPRETATION")

    print(
        "1. Every variant has one primary category."
    )

    print(
        "2. A variant can have multiple VEP "
        "consequences across transcripts."
    )

    print(
        "3. Protein-coding CDS means a direct "
        "coding-sequence consequence."
    )

    print(
        "4. Intronic, intergenic, UTR and most "
        "regulatory categories are noncoding DNA."
    )

    print(
        "5. lncRNA overlaps are additional, "
        "non-exclusive transcript associations."
    )

    print(
        "6. To identify lncRNA-specific exons "
        "and introns, Analysis05 must match "
        "consequence and biotype on the same "
        "transcript annotation row."
    )

    print(
        "7. Pooled GWAS rows are not independent "
        "causal variants."
    )

    return {name: pd.DataFrame(rows) for name, rows in TABLES.items()}


def main(argv=None):
    return analysis_main(run, "John Analysis04: genomic consequence distribution", argv, extra=[
        ("--high-pip", dict(type=float, default=0.50)),
        ("--top-terms", dict(type=int, default=12)),
    ])


if __name__ == "__main__":
    main()
 