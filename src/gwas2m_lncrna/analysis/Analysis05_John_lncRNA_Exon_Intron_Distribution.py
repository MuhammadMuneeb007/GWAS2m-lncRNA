
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
ANALYSIS 05 - JOHN MATTICK

lncRNA exon, intron, and splice-region distribution.

READ ONLY:
  Existing Step06 and Step07 results are read for the phenotypes passed in.
  Results are printed to the terminal; run() also returns the result tables.

One observation = one VEP_ID per GWAS study.
Transcript consequences are matched to the lncRNA
biotype on the SAME VEP annotation row.
"""

import argparse
import csv
import gzip
import heapq
import itertools
import json
import math
import re
import sys

from collections import Counter
from pathlib import Path

import pandas as pd

from gwas2m_lncrna.analysis.helpers import analysis_main


# Set by run(): [(display name, phenotype folder), ...] and the ancestry folder name.
DISEASES = []
ANCESTRY = "european"
TABLES = {}

DONE = {
    "COMPLETE",
    "SKIPPED_ALREADY_COMPLETE",
}

CORE = {
    "splice_acceptor_variant",
    "splice_donor_variant",
}

REGION = {
    "splice_region_variant",
    "splice_donor_region_variant",
    "splice_donor_5th_base_variant",
    "splice_polypyrimidine_tract_variant",
}

EXON = {
    "non_coding_transcript_exon_variant",
}

INTRON = {
    "intron_variant",
}

INSIDE_OTHER = {
    "non_coding_transcript_variant",
}

NEAR = {
    "upstream_gene_variant",
    "downstream_gene_variant",
}

EXTENDED = {
    "antisense",
    "processed_transcript",
    "sense_intronic",
    "sense_overlapping",
    "3prime_overlapping_ncrna",
    "bidirectional_promoter_lncrna",
    "macro_lncrna",
}

# Mutually exclusive primary category priority.
CATEGORIES = (
    "SPLICE_CORE",
    "SPLICE_REGION",
    "LNCRNA_EXON",
    "LNCRNA_INTRON",
    "OTHER_LNCRNA_TRANSCRIPT",
    "NEAR_LNCRNA_ONLY",
    "NO_LNCRNA_OVERLAP",
    "UNANNOTATED",
)

POPS = (
    "CREDIBLE95",
    "HIGH_PIP",
    "ALL_REGIONAL",
)

try:
    csv.field_size_limit(sys.maxsize)
except OverflowError:
    csv.field_size_limit(2**31 - 1)


def load_json(path):
    try:
        with path.open(encoding="utf-8") as f:
            data = json.load(f)

        return data if isinstance(data, dict) else {}

    except (OSError, ValueError):
        return {}


def manifest_rows(path):
    with path.open(
        encoding="utf-8-sig",
        newline=""
    ) as f:

        return list(
            csv.DictReader(f, delimiter="\t")
        )


def to_float(value):
    try:
        value = float(value)

        if math.isfinite(value) and 0 <= value <= 1:
            return value

    except (ValueError, TypeError):
        pass

    return None


def flag(value):
    return str(value).strip().lower() in {
        "true", "t", "yes", "1"
    }


def items(value):
    return {
        x.strip().lower()
        for x in re.split(r"[,;&|]", str(value))
        if x.strip()
        and x.strip().lower()
        not in {"nan", "none", ".", "na"}
    }


def has_number(value):
    return bool(
        re.match(
            r"^\d+(?:/\d+)?$",
            str(value).strip()
        )
    )


def blank():
    return {
        "n": 0,
        "cat": Counter(),
        "flags": Counter(),
        "lnc": 0,
        "near": 0,
        "lnc_and_pc": 0,
        "pc": 0,
    }


def new_groups():
    return {
        name: blank()
        for name in POPS
    }


def add(target, other):
    target["n"] += other["n"]

    target["cat"].update(
        other["cat"]
    )

    target["flags"].update(
        other["flags"]
    )

    for k in (
        "lnc",
        "near",
        "lnc_and_pc",
        "pc",
    ):
        target[k] += other[k]


def percent(n, den):
    if not den:
        return "N/A"

    return f"{100 * n / den:.3f}%"


def call_lncrna(records, include_legacy):
    """
    Match lncRNA BIOTYPE and consequence on
    the SAME transcript-level annotation row.

    Returns one primary category per variant.
    """

    seen = set()
    genes = set()
    transcripts = set()

    near = False
    has_annotation = False
    has_pc = False

    flags = Counter()

    for record in records:

        bio = str(
            record.get("BIOTYPE", "")
        ).strip().lower()

        terms = items(
            record.get("Consequence", "")
        )

        if terms:
            has_annotation = True

        # Protein-coding transcript overlap.
        if bio == "protein_coding" and (
            has_number(record.get("EXON", ""))
            or has_number(record.get("INTRON", ""))
            or any("splice" in x for x in terms)
            or terms & {
                "intron_variant",
                "5_prime_utr_variant",
                "3_prime_utr_variant",
                "coding_sequence_variant",
                "missense_variant",
                "synonymous_variant",
                "stop_gained",
                "stop_lost",
                "start_lost",
                "frameshift_variant",
                "inframe_deletion",
                "inframe_insertion",
                "protein_altering_variant",
            }
        ):
            has_pc = True

        is_lnc = (
            bio in {"lncrna", "lincrna"}
            or (
                include_legacy
                and bio in EXTENDED
            )
        )

        if not is_lnc:
            continue

        gene = str(
            record.get("Gene", "")
        ).strip()

        symbol = str(
            record.get("SYMBOL", "")
        ).strip()

        feature = str(
            record.get("Feature", "")
        ).strip()



        # Nearby does not automatically mean overlap.
        if terms & NEAR:
            near = True

        exon = bool(
            terms & EXON
            or has_number(
                record.get("EXON", "")
            )
        )

        intron = bool(
            terms & INTRON
            or has_number(
                record.get("INTRON", "")
            )
        )

        core = bool(
            terms & CORE
        )

        region = bool(
            terms & REGION
        )

        if core:
            seen.add("SPLICE_CORE")
            flags["SPLICE_CORE"] = 1

        if region:
            seen.add("SPLICE_REGION")
            flags["SPLICE_REGION"] = 1

        if exon:
            seen.add("LNCRNA_EXON")
            flags["LNCRNA_EXON"] = 1

        if intron:
            seen.add("LNCRNA_INTRON")
            flags["LNCRNA_INTRON"] = 1

        if terms & INSIDE_OTHER:
            seen.add(
                "OTHER_LNCRNA_TRANSCRIPT"
            )
        # Collect IDs only for actual same-row transcript intersections.
        if core or region or exon or intron or bool(terms & INSIDE_OTHER):
            if gene not in {"", "nan", "-"}:
                genes.add(
                    symbol
                    if symbol not in {"", "nan", "-"}
                    else gene
                )
            if feature not in {"", "nan", "-"}:
                transcripts.add(feature)

    # One primary category for each variant.
    primary = next(
        (
            x
            for x in CATEGORIES[:5]
            if x in seen
        ),
        None,
    )

    if primary is None:

        if near:
            primary = "NEAR_LNCRNA_ONLY"

        elif has_annotation:
            primary = "NO_LNCRNA_OVERLAP"

        else:
            primary = "UNANNOTATED"

    overlapping = bool(seen)

    return (
        primary,
        overlapping,
        has_pc,
        flags,
        genes,
        transcripts,
    )


def top_push(heap, item, limit):
    if len(heap) < limit:
        heapq.heappush(
            heap, item
        )

    elif item > heap[0]:
        heapq.heapreplace(
            heap, item
        )


def study_readiness(root, slug, study):

    base = (
        root
        / "06_finemapping"
        / slug
        / ANCESTRY
        / study
    )

    anno = (
        root
        / "07_annotation"
        / slug
        / ANCESTRY
        / study
    )

    audit = (
        root
        / "pipeline_audit"
        / slug
        / ANCESTRY
        / "status"
        / study
    )

    s06 = str(
        load_json(
            audit / "Step06_FineMapping.json"
        ).get(
            "status",
            "NOT_RECORDED"
        )
    ).upper()

    s07 = str(
        load_json(
            audit / "Step07_VEP.json"
        ).get(
            "status",
            "NOT_RECORDED"
        )
    ).upper()

    fm = load_json(
        base / "finemapping_summary.json"
    )

    if s06 not in DONE:
        return None, "S06_" + s06

    if (
        str(fm.get("STATUS", "")).upper()
        == "NO_SIGNIFICANT_LOCI"
        or (
            fm
            and fm.get("N_LOCI") == 0
        )
    ):
        return None, "NO_SIGNIFICANT_LOCI"

    if s07 not in DONE:
        return None, "S07_" + s07

    info = load_json(
        anno / "annotation_summary.json"
    )

    fine = (
        base
        / f"{study}_finemapped_variants.tsv.gz"
    )

    all_file = (
        anno
        / f"{study}_VEP_all_consequences.tsv.gz"
    )

    summary = (
        anno
        / f"{study}_VEP_variant_summary.tsv"
    )

    if not all([
        fine.is_file(),
        all_file.is_file(),
        summary.is_file(),
    ]):
        return None, "MISSING_REQUIRED_FILE"

    if (
        (anno / "ANNOTATION_FAILED.txt").exists()
        or info.get("STATUS") != "COMPLETE"
        or info.get("ANNOTATION_SCOPE")
        != "ALL_FINEMAPPED_VARIANTS"
    ):
        return None, "INVALID_ANNOTATION_STATUS"

    # Ensure Step07 used the current Step06 file.
    try:
        current = (
            int(info["FINEMAPPED_INPUT_SIZE"])
            == fine.stat().st_size
            and
            int(info["FINEMAPPED_INPUT_MTIME_NS"])
            == fine.stat().st_mtime_ns
        )

    except (KeyError, ValueError, TypeError):
        current = False

    if not current:
        return None, "STALE_ANNOTATION"

    return (
        all_file,
        info,
    ), "READY"


def process_study(
    all_file,
    info,
    study,
    high_pip,
    include_legacy,
    top_n,
):

    totals = new_groups()

    locus_all = set()
    locus_lnc = set()

    top = []

    unique_cs = set()
    unique_lnc = set()

    n_variants = 0
    n_raw = 0

    last_id = ""

    try:

        with gzip.open(
            all_file,
            "rt",
            encoding="utf-8",
            newline=""
        ) as handle:

            reader = csv.DictReader(
                handle,
                delimiter="\t"
            )

            required = {
                "VEP_ID",
                "BIOTYPE",
                "Consequence",
                "LOCUS_ID",
                "PIP",
                "IN_95_CREDIBLE_SET",
                "CHR",
                "REFERENCE_POS",
                "REFERENCE_REF",
                "REFERENCE_ALT",
            }

            missing = required - set(
                reader.fieldnames or []
            )

            if missing:
                return (
                    None,
                    "MISSING_COLUMNS:"
                    + ",".join(sorted(missing)),
                )

            # Step07 writes consecutive transcript
            # annotation records for each VEP_ID.
            for vep_id, group in itertools.groupby(
                reader,
                key=lambda row: row["VEP_ID"]
            ):

                if (
                    not vep_id
                    or (
                        last_id
                        and vep_id <= last_id
                    )
                ):
                    return (
                        None,
                        "VEP_IDS_NOT_SORTED_OR_DUPLICATED",
                    )

                last_id = vep_id

                records = list(group)

                n_raw += len(records)

                row = records[0]

                (
                    primary,
                    overlapping,
                    has_pc,
                    flags,
                    genes,
                    transcripts,
                ) = call_lncrna(
                    records,
                    include_legacy,
                )

                pip = to_float(
                    row.get("PIP")
                )

                credible = flag(
                    row.get(
                        "IN_95_CREDIBLE_SET"
                    )
                )

                groups = ["ALL_REGIONAL"]

                if credible:
                    groups.append(
                        "CREDIBLE95"
                    )

                if (
                    pip is not None
                    and pip >= high_pip
                ):
                    groups.append(
                        "HIGH_PIP"
                    )

                for name in groups:

                    t = totals[name]

                    t["n"] += 1
                    t["cat"][primary] += 1

                    t["flags"].update(flags)

                    t["lnc"] += int(
                        overlapping
                    )

                    t["near"] += int(
                        primary
                        == "NEAR_LNCRNA_ONLY"
                    )

                    t["pc"] += int(
                        has_pc
                    )

                    t["lnc_and_pc"] += int(
                        overlapping and has_pc
                    )

                locus = str(
                    row.get("LOCUS_ID", "")
                )

                if credible and locus:

                    locus_all.add(locus)

                    if overlapping:
                        locus_lnc.add(locus)

                fields = [
                    str(row.get(k, "")).strip()
                    for k in (
                        "CHR",
                        "REFERENCE_POS",
                        "REFERENCE_REF",
                        "REFERENCE_ALT",
                    )
                ]

                variant_key = ":".join(
                    fields
                )

                if (
                    credible
                    and all(fields)
                ):

                    unique_cs.add(
                        variant_key
                    )

                    if overlapping:
                        unique_lnc.add(
                            variant_key
                        )

                if (
                    credible
                    and overlapping
                    and pip is not None
                ):

                    top_push(
                        top,
                        (
                            pip,
                            study,
                            locus,
                            variant_key,
                            ",".join(
                                sorted(genes)
                            )[:65],
                            ",".join(
                                sorted(transcripts)
                            )[:55],
                            primary,
                        ),
                        top_n,
                    )

                n_variants += 1

    except (
        OSError,
        EOFError,
        UnicodeError,
        csv.Error,
        ValueError,
    ) as exc:

        return (
            None,
            "READ_ERROR:" + str(exc),
        )

    # Validate the variant count against Step07.
    try:
        expected = int(
            info.get("N_SELECTED_VARIANTS")
        )

    except (ValueError, TypeError):
        expected = None

    if (
        expected is not None
        and n_variants != expected
    ):
        return (
            None,
            f"VARIANT_COUNT_MISMATCH:"
            f"{n_variants}_vs_{expected}",
        )

    if n_variants == 0:
        return (
            None,
            "EMPTY_ALL_CONSEQUENCES",
        )

    return {
        "pop": totals,
        "loci": locus_all,
        "loci_lnc": locus_lnc,
        "top": top,
        "unique_cs": unique_cs,
        "unique_lnc": unique_lnc,
        "raw_rows": n_raw,
    }, "READY"


def show_stats(title, t):

    print("\n" + "-" * 94)
    print(title)
    print("-" * 94)

    print(
        f"Variants analyzed : {t['n']:,}"
    )

    print(
        f"lncRNA overlaps   : {t['lnc']:,} "
        f"({percent(t['lnc'], t['n'])})"
    )

    print()

    print(
        f"{'PRIMARY CATEGORY':32}"
        f"{'COUNT':>12}"
        f"{'% ALL':>13}"
        f"{'% lncRNA':>14}"
    )

    for category in CATEGORIES:

        n = t["cat"][category]

        within = (
            percent(n, t["lnc"])
            if category in CATEGORIES[:5]
            else "--"
        )

        print(
            f"{category:32}"
            f"{n:12,}"
            f"{percent(n, t['n']):>13}"
            f"{within:>14}"
        )

    print("-" * 72)

    print(
        f"{'TOTAL':32}"
        f"{t['n']:12,}"
        f"{percent(t['n'], t['n']):>13}"
    )

    print(
        "\nADDITIONAL NON-EXCLUSIVE "
        "TRANSCRIPT OVERLAPS"
    )

    for c in (
        "LNCRNA_EXON",
        "LNCRNA_INTRON",
        "SPLICE_CORE",
        "SPLICE_REGION",
    ):

        n = t["flags"][c]

        print(
            f"{c:30}: {n:>9,} "
            f"({percent(n, t['n'])} of all; "
            f"{percent(n, t['lnc'])} "
            f"of lncRNA overlaps)"
        )

    print(
        "\nlncRNA + protein-coding transcript overlaps:",
        f"{t['lnc_and_pc']:,}"
    )

    print(
        "Near lncRNA only (upstream/downstream):",
        f"{t['near']:,}"
    )


def stats_rows(disease, population, t):
    """Table rows for one population: primary lncRNA categories and overlap flags."""
    rows = []
    for category in CATEGORIES:
        n = t["cat"][category]
        rows.append({"PHENOTYPE": disease, "POPULATION": population, "MEASURE": category, "N": n,
                     "PERCENT_OF_ALL": round(100 * n / t["n"], 3) if t["n"] else None,
                     "PERCENT_OF_LNCRNA": (round(100 * n / t["lnc"], 3) if t["lnc"] and category in CATEGORIES[:5]
                                           else None),
                     "DENOMINATOR": t["n"]})
    for label, n in (("LNCRNA_TRANSCRIPT_OVERLAP", t["lnc"]), ("FLAG_LNCRNA_EXON", t["flags"]["LNCRNA_EXON"]),
                     ("FLAG_LNCRNA_INTRON", t["flags"]["LNCRNA_INTRON"]),
                     ("FLAG_SPLICE_CORE", t["flags"]["SPLICE_CORE"]),
                     ("FLAG_SPLICE_REGION", t["flags"]["SPLICE_REGION"]),
                     ("LNCRNA_AND_PROTEIN_CODING", t["lnc_and_pc"]), ("NEAR_LNCRNA_ONLY", t["near"])):
        rows.append({"PHENOTYPE": disease, "POPULATION": population, "MEASURE": label, "N": n,
                     "PERCENT_OF_ALL": round(100 * n / t["n"], 3) if t["n"] else None,
                     "PERCENT_OF_LNCRNA": round(100 * n / t["lnc"], 3) if t["lnc"] and label.startswith("FLAG") else None,
                     "DENOMINATOR": t["n"]})
    return rows


def run(root, diseases, ancestry="european", high_pip=0.50, top=15,
        include_legacy_lncrna_biotypes=False, **_):
    """Run Analysis 05 for [(display name, folder slug), ...]; return its tables."""
    global DISEASES, ANCESTRY
    DISEASES, ANCESTRY = list(diseases), ancestry
    TABLES.clear()
    args = argparse.Namespace(high_pip=float(high_pip), top=max(1, int(top)),
                              include_legacy_lncrna_biotypes=bool(include_legacy_lncrna_biotypes))
    if not 0 <= args.high_pip <= 1:
        raise ValueError("high_pip must be 0..1")
    root = Path(root).expanduser().resolve()

    print("=" * 110)
    print("JOHN MATTICK - ANALYSIS 05")
    print("lncRNA EXON / INTRON / SPLICE DISTRIBUTION")
    print("READ ONLY")
    print("=" * 110)

    print(
        "lncRNA definition:",
        (
            "EXTENDED EXPLORATORY"
            if args.include_legacy_lncrna_biotypes
            else "STRICT VEP BIOTYPE"
        ),
    )

    combined = new_groups()
    disease_results = {}

    for disease, slug in DISEASES:

        print("\n" + "=" * 110)
        print(disease.upper())
        print("=" * 110)

        manifest = (
            root
            / "01_gwas_catalog"
            / slug
            / ANCESTRY
            / "GWAS_download_manifest.tsv"
        )

        if not manifest.is_file():
            manifest = (
                root
                / "pipeline_audit"
                / slug
                / ANCESTRY
                / "study_manifest.tsv"
            )

        if not manifest.is_file():
            print("No manifest found; skipping.")
            continue

        studies = manifest_rows(
            manifest
        )

        counts = Counter()
        local = new_groups()

        study_top = []

        loci_all = set()
        loci_lnc = set()

        unique_all = set()
        unique_lnc = set()

        print(
            f"\nSelected GWAS: {len(studies)}"
        )

        print(
            f"\n{'STUDY':15}"
            f"{'STATUS':31}"
            f"{'CREDIBLE':>12}"
            f"{'lncRNA':>12}"
            f"{'ALL VARIANTS':>14}"
        )

        for meta in studies:

            study = str(
                meta.get(
                    "STUDY_ACCESSION", ""
                )
            ).strip()

            if not study:
                continue

            source, reason = study_readiness(
                root,
                slug,
                study,
            )

            if source is None:

                counts[reason] += 1

                print(
                    f"{study:15}"
                    f"{reason[:31]:31}"
                    f"{'-':>12}"
                    f"{'-':>12}"
                    f"{'-':>14}"
                )

                continue

            result, reason = process_study(
                source[0],
                source[1],
                study,
                args.high_pip,
                args.include_legacy_lncrna_biotypes,
                args.top,
            )

            if result is None:

                counts[reason] += 1

                print(
                    f"{study:15}"
                    f"{reason[:31]:31}"
                    f"{'-':>12}"
                    f"{'-':>12}"
                    f"{'-':>14}"
                )

                continue

            counts["INCLUDED"] += 1

            for key in POPS:

                add(
                    local[key],
                    result["pop"][key],
                )

                add(
                    combined[key],
                    result["pop"][key],
                )

            loci_all.update(
                (study, x)
                for x in result["loci"]
            )

            loci_lnc.update(
                (study, x)
                for x in result["loci_lnc"]
            )

            unique_all.update(
                result["unique_cs"]
            )

            unique_lnc.update(
                result["unique_lnc"]
            )

            for item in result["top"]:
                top_push(
                    study_top,
                    item,
                    args.top,
                )

            cs = result["pop"]["CREDIBLE95"]
            TABLES.setdefault("studies", []).append({
                "PHENOTYPE": disease, "STUDY_ACCESSION": study, "CS_VARIANTS": cs["n"],
                "CS_LNCRNA_OVERLAP": cs["lnc"],
                "PERCENT_LNCRNA": round(100 * cs["lnc"] / cs["n"], 3) if cs["n"] else None,
                "ALL_REGIONAL_VARIANTS": result["pop"]["ALL_REGIONAL"]["n"],
                "REGIONAL_LNCRNA_OVERLAP": result["pop"]["ALL_REGIONAL"]["lnc"],
            })

            print(
                f"{study:15}"
                f"{'INCLUDED':31}"
                f"{cs['n']:12,}"
                f"{cs['lnc']:12,}"
                f"{result['pop']['ALL_REGIONAL']['n']:14,}",
                flush=True,
            )

        # Study status summary.
        print("\nSTUDY READINESS")

        for reason, count in sorted(
            counts.items()
        ):
            print(
                f"  {reason:38} {count:,}"
            )

        # Main genomic results.
        show_stats(
            "PRIMARY: 95% CREDIBLE-SET VARIANTS",
            local["CREDIBLE95"],
        )

        show_stats(
            f"SENSITIVITY: PIP >= {args.high_pip:g}",
            local["HIGH_PIP"],
        )

        show_stats(
            "COMPARISON: ALL FINE-MAPPED REGIONAL VARIANTS",
            local["ALL_REGIONAL"],
        )

        # Independent counting units.
        print("\nLOCUS-LEVEL SUMMARY")

        print(
            "Study-locus units with credible variants:",
            f"{len(loci_all):,}"
        )

        print(
            "Study-locus units with lncRNA overlaps:",
            f"{len(loci_lnc):,}",
            f"({percent(len(loci_lnc), len(loci_all))})"
        )

        print("\nUNIQUE ALLELES ACROSS GWAS STUDIES")

        print(
            "Unique credible-set alleles:",
            f"{len(unique_all):,}"
        )

        print(
            "Unique lncRNA-overlapping alleles:",
            f"{len(unique_lnc):,}",
            f"({percent(len(unique_lnc), len(unique_all))})"
        )

        print(
            "\nTOP lncRNA CREDIBLE-SET VARIANTS BY PIP"
        )

        for (
            pip,
            study,
            locus,
            variant_key,
            genes,
            transcripts,
            primary,
        ) in sorted(
            study_top,
            reverse=True,
        ):

            print(
                f"{study:13} "
                f"{locus:7} "
                f"PIP={pip:.5f} "
                f"{variant_key:27} "
                f"{primary:23} "
                f"gene={genes or 'NA'} "
                f"transcript={transcripts or 'NA'}"
            )

        disease_results[disease] = (
            local["CREDIBLE95"]
        )

        for key in POPS:
            TABLES.setdefault("distribution", []).extend(stats_rows(disease, key, local[key]))
        TABLES.setdefault("locus_and_unique_alleles", []).append({
            "PHENOTYPE": disease, "STUDY_LOCI_WITH_CS": len(loci_all),
            "STUDY_LOCI_WITH_LNCRNA_CS": len(loci_lnc), "UNIQUE_CS_ALLELES": len(unique_all),
            "UNIQUE_LNCRNA_CS_ALLELES": len(unique_lnc),
        })
        TABLES.setdefault("top_lncrna_variants", []).extend(
            {"PHENOTYPE": disease, "STUDY_ACCESSION": s, "LOCUS_ID": loc, "PIP": p, "ALLELE": key,
             "PRIMARY": prim, "LNCRNA_GENES": genes, "TRANSCRIPTS": tx}
            for p, s, loc, key, genes, tx, prim in sorted(study_top, reverse=True)
        )

    # ------------------------------------------------
    # CROSS-DISEASE RESULTS
    # ------------------------------------------------

    print("\n" + "=" * 110)
    print("CROSS-DISEASE COMPARISON")
    print("95% CREDIBLE-SET VARIANTS")
    print("=" * 110)

    print(
        f"{'MEASURE':36}"
        + "".join(f"{name[:17].upper():>18}" for name, _ in DISEASES)
    )

    measures = [
        (
            "lncRNA transcript overlap",
            "lnc",
        ),
        (
            "lncRNA exon (primary)",
            "LNCRNA_EXON",
        ),
        (
            "lncRNA intron (primary)",
            "LNCRNA_INTRON",
        ),
        (
            "lncRNA splice core (primary)",
            "SPLICE_CORE",
        ),
        (
            "lncRNA splice region (primary)",
            "SPLICE_REGION",
        ),
    ]

    for label, key in measures:

        cells = []

        for disease, _ in DISEASES:

            s = disease_results.get(
                disease,
                blank(),
            )

            n = (
                s[key]
                if key == "lnc"
                else s["cat"][key]
            )

            cells.append(
                f"{n:,} ({percent(n, s['n'])})"
            )

        print(
            f"{label:36}"
            + "".join(f"{c:>18}" for c in cells)
        )

    show_stats(
        "POOLED 95% CREDIBLE-SET VARIANT RECORDS",
        combined["CREDIBLE95"],
    )

    print("\nSCIENTIFIC NOTES")

    print(
        "1. Consequence and biotype are matched "
        "on the same transcript annotation."
    )

    print(
        "2. One VEP_ID counts as one variant "
        "per study, not one per transcript."
    )

    print(
        "3. Primary categories are mutually "
        "exclusive; overlap flags are not."
    )

    print(
        "4. Upstream/downstream annotations "
        "are not counted as lncRNA transcript overlap."
    )

    print(
        "5. Shared GWAS variants and cohorts "
        "can affect pooled counts."
    )

    print(
        "6. lncRNA overlap does not prove "
        "lncRNA-mediated disease causality."
    )

    return {name: pd.DataFrame(rows) for name, rows in TABLES.items()}


def main(argv=None):
    return analysis_main(run, "John Analysis05: lncRNA exon/intron/splice distribution", argv, extra=[
        ("--high-pip", dict(type=float, default=0.50)),
        ("--top", dict(type=int, default=15)),
        ("--include-legacy-lncrna-biotypes", dict(action="store_true")),
    ])


if __name__ == "__main__":
    main()
