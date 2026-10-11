 
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
JOHN MATTICK - ANALYSIS 09

Fine-mapping sensitivity and PIP-weighted analysis.

Phenotypes and ancestry are passed in (run() / --phenotype, --ancestry).

Inputs:
    Step06 fine-mapping outputs
    Step06 95% credible-set membership
    Step07 transcript-level VEP annotations

Outputs:
    Printed to the terminal; run() also returns the result tables.

Scientific safeguards:
    - Only completed Step06/Step07 studies
    - Validate Step07 against current Step06 input
    - Validate credible-set membership
    - Require PIP values between 0 and 1
    - Match consequence and biotype on the same transcript
    - Count each VEP variant once per study
    - Distinguish row-weighted, PIP-weighted, study-
      balanced and locus-balanced summaries
"""

import csv
import gzip
import itertools
import json
import math
import re
import statistics
import sys

from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd

from gwas2m_lncrna.analysis.helpers import analysis_main, direct_labels, normalize_label


try:
    csv.field_size_limit(sys.maxsize)
except OverflowError:
    csv.field_size_limit(2**31 - 1)


# Set by run(): [(display name, phenotype folder), ...] and the ancestry folder name.
DISEASES = []
ANCESTRY = "european"
TABLES = {}

DONE = {
    "COMPLETE",
    "SKIPPED_ALREADY_COMPLETE",
}

POPS = (
    "REGIONAL",
    "CS95",
    "PIP_GE_0_10",
    "PIP_GE_0_50",
    "PIP_GE_0_90",
    "CS95_AND_PIP_GE_0_50",
)

FEATURES = (
    "CDS",
    "NONCODING",
    "LNCRNA_BODY",
    "LNCRNA_EXON",
    "LNCRNA_INTRON",
    "LNCRNA_SPLICE",
    "PC_GENE_BODY",
    "BOTH_LNCRNA_PC",
)

CDS = {
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

SPLICE = {
    "splice_acceptor_variant",
    "splice_donor_variant",
    "splice_region_variant",
    "splice_donor_region_variant",
    "splice_donor_5th_base_variant",
    "splice_polypyrimidine_tract_variant",
}

BODY = {
    "non_coding_transcript_exon_variant",
    "non_coding_transcript_variant",
    "intron_variant",
    "5_prime_utr_variant",
    "3_prime_utr_variant",
} | CDS | SPLICE

IGNORED = {
    "",
    "na",
    "nan",
    "none",
    ".",
    "-",
}


# ==========================================================
# FILE READING
# ==========================================================

def read_json(path):
    try:
        with path.open(encoding="utf-8") as stream:
            data = json.load(stream)

        return data if isinstance(data, dict) else {}

    except (OSError, ValueError, UnicodeError):
        return {}


def tsv_rows(path):
    with path.open(
        "r",
        encoding="utf-8-sig",
        newline=""
    ) as stream:

        yield from csv.DictReader(
            stream,
            delimiter="\t"
        )


def flag(value):
    return str(value).strip().lower() in {
        "true",
        "t",
        "1",
        "yes",
        "y",
    }


def number(value):
    try:
        n = float(str(value).strip())

        return n if math.isfinite(n) else None

    except (ValueError, TypeError):
        return None


def integer(value):
    # Important: do not convert nanosecond timestamps
    # through float, which would lose precision.
    try:
        return int(str(value).strip())

    except (ValueError, TypeError):
        return None


def terms(value):
    return {
        s.strip().lower()
        for s in re.split(r"[&,;|]", str(value))
        if s.strip().lower() not in IGNORED
    }


def location(value):
    return bool(
        re.fullmatch(
            r"\d+(?:/\d+)?",
            str(value).strip()
        )
    )


def fmt_percent(num, den):
    if not den:
        return "NA"

    return f"{100 * num / den:.3f}%"


# ==========================================================
# COUNTING STRUCTURES
# ==========================================================

def blank():
    return {
        "n": 0,
        "weight": 0.0,
        "counts": Counter(),
        "weighted": defaultdict(float),
    }


def make_pops():
    return {
        key: blank()
        for key in POPS
    }


def increment(target, features, pip):
    target["n"] += 1
    target["weight"] += pip

    for feature in FEATURES:
        if features[feature]:
            target["counts"][feature] += 1
            target["weighted"][feature] += pip


def merge(dst, src):
    dst["n"] += src["n"]
    dst["weight"] += src["weight"]

    dst["counts"].update(
        src["counts"]
    )

    for key, value in src["weighted"].items():
        dst["weighted"][key] += value


# ==========================================================
# STEP06 CREDIBLE-SET MEMBERSHIP
# ==========================================================

def read_membership(path):
    """
    Read actual Step06 credible-set membership.

    A variant can belong to more than one credible set,
    but is counted only once for the union of sets.
    """

    if not path.is_file():
        raise ValueError(
            "Missing Step06 95pct credible-sets membership TSV"
        )

    result = set()

    with path.open(
        encoding="utf-8-sig",
        newline=""
    ) as fh:

        reader = csv.DictReader(
            fh,
            delimiter="\t"
        )

        required = {
            "LOCUS_ID",
            "REFERENCE_ID",
            "CREDIBLE_SET",
        }

        if not required.issubset(
            set(reader.fieldnames or [])
        ):
            raise ValueError(
                "Step06 credible sets missing expected columns"
            )

        for row in reader:

            locus = str(
                row.get("LOCUS_ID", "")
            ).strip()

            variant = str(
                row.get("REFERENCE_ID", "")
            ).strip()

            if locus and variant:
                result.add(
                    (locus, variant)
                )

    return result


# ==========================================================
# STUDY DISCOVERY
# ==========================================================

def get_manifest(root, slug):

    options = (
        root
        / "01_gwas_catalog"
        / slug
        / ANCESTRY
        / "GWAS_download_manifest.tsv",

        root
        / "pipeline_audit"
        / slug
        / ANCESTRY
        / "study_manifest.tsv",
    )

    for path in options:
        if path.is_file():
            return list(tsv_rows(path)), path

    return [], None


def direct_trait(trait, disease):
    """
    Conservative exact-label classification.

    Related/proxy studies are still reported,
    but not treated as direct-disease replications.
    """

    s = re.sub(
        r"\s+",
        " ",
        str(trait).lower().replace("\u2019", "'")
    ).strip()

    return normalize_label(s) in direct_labels(disease)


def stage_status(audit, study, name):

    path = (
        audit
        / "status"
        / study
        / (name + ".json")
    )

    return str(
        read_json(path).get(
            "status",
            "NOT_RECORDED"
        )
    ).upper()


# ==========================================================
# TRANSCRIPT-SPECIFIC GENOMIC CLASSIFICATION
# ==========================================================

def analyze_transcripts(group):
    """
    Match transcript BIOTYPE and consequence
    on the same VEP annotation row.

    Multiple annotations can apply to the same
    variant, but each variant is counted once.
    """

    features = dict.fromkeys(
        FEATURES,
        False
    )

    coding = False
    lnc = False
    pc = False

    exon = False
    intron = False
    splice = False

    for row in group:

        bio = str(
            row.get("BIOTYPE", "")
        ).strip().lower()

        consequences = terms(
            row.get("Consequence", "")
        )

        # Protein-coding transcript annotation.
        if bio == "protein_coding":

            if consequences & CDS:
                coding = True

            if (
                consequences & BODY
                or location(row.get("EXON", ""))
                or location(row.get("INTRON", ""))
            ):
                pc = True

        # lncRNA transcript annotation.
        if bio in {"lncrna", "lincrna"}:

            row_exon = (
                "non_coding_transcript_exon_variant"
                in consequences
                or location(row.get("EXON", ""))
            )

            row_intron = (
                "intron_variant" in consequences
                or location(row.get("INTRON", ""))
            )

            row_splice = bool(
                consequences & SPLICE
            )

            if (
                row_exon
                or row_intron
                or row_splice
                or "non_coding_transcript_variant"
                in consequences
            ):
                lnc = True

            exon |= row_exon
            intron |= row_intron
            splice |= row_splice

    features["CDS"] = coding
    features["NONCODING"] = not coding

    features["LNCRNA_BODY"] = lnc
    features["LNCRNA_EXON"] = exon
    features["LNCRNA_INTRON"] = intron
    features["LNCRNA_SPLICE"] = splice

    features["PC_GENE_BODY"] = pc

    features["BOTH_LNCRNA_PC"] = (
        lnc and pc
    )

    return features


# ==========================================================
# CHECK COMPLETED STEP06 AND STEP07
# ==========================================================

def check_paths(root, slug, study):

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
    )

    s06 = stage_status(
        audit,
        study,
        "Step06_FineMapping"
    )

    s07 = stage_status(
        audit,
        study,
        "Step07_VEP"
    )

    if s06 not in DONE:
        return None, "S06_" + s06

    fm = read_json(
        base / "finemapping_summary.json"
    )

    if (
        str(fm.get("STATUS", "")).upper()
        == "NO_SIGNIFICANT_LOCI"
        or (
            fm
            and integer(fm.get("N_LOCI")) == 0
        )
    ):
        return None, "NO_SIGNIFICANT_LOCI"

    if s07 not in DONE:
        return None, "S07_" + s07

    meta = read_json(
        anno / "annotation_summary.json"
    )

    fine = (
        base
        / (study + "_finemapped_variants.tsv.gz")
    )

    all_path = (
        anno
        / (study + "_VEP_all_consequences.tsv.gz")
    )

    cs = (
        base
        / (study + "_95pct_credible_sets.tsv")
    )

    if not all(
        p.is_file() and p.stat().st_size > 0
        for p in (fine, all_path, cs)
    ):
        return None, "MISSING_FINEMAPPING_OR_VEP_INPUT"

    if (
        (anno / "ANNOTATION_FAILED.txt").exists()
        or str(meta.get("STATUS", "")).upper()
        != "COMPLETE"
        or str(meta.get("ANNOTATION_SCOPE", "")).upper()
        != "ALL_FINEMAPPED_VARIANTS"
    ):
        return None, "INVALID_STEP07_COMPLETION_OR_SCOPE"

    # Check exact input version and timestamp.
    if (
        integer(meta.get("FINEMAPPED_INPUT_SIZE"))
        != fine.stat().st_size
        or integer(meta.get("FINEMAPPED_INPUT_MTIME_NS"))
        != fine.stat().st_mtime_ns
    ):
        return None, "STALE_STEP07_INPUT"

    return (
        all_path,
        cs,
        meta,
    ), "READY"


# ==========================================================
# PROCESS INDIVIDUAL GWAS STUDY
# ==========================================================

def analyze_study(
    all_path,
    cs_path,
    meta,
    study,
):
    """
    Stream full transcript-level annotations.

    The Step06 membership file defines the primary
    credible-set variant population.

    Step07 membership flags are independently checked.
    """

    membership = read_membership(
        cs_path
    )

    seen_membership = set()

    pops = make_pops()

    # Per-locus information for a locus-balanced analysis.
    locus_pops = defaultdict(
        make_pops
    )

    mismatched_cs_flags = 0
    n_records = 0
    previous = ""

    with gzip.open(
        all_path,
        "rt",
        encoding="utf-8",
        newline=""
    ) as fh:

        reader = csv.DictReader(
            fh,
            delimiter="\t"
        )

        required = {
            "VEP_ID",
            "LOCUS_ID",
            "REFERENCE_ID",
            "BIOTYPE",
            "Consequence",
            "PIP",
            "IN_95_CREDIBLE_SET",
        }

        missing = required - set(
            reader.fieldnames or []
        )

        if missing:
            raise ValueError(
                "VEP columns missing: "
                + ", ".join(sorted(missing))
            )

        # Consequences for each VEP_ID are consecutive
        # in the Step07 all-consequences output.
        for vep_id, rows in itertools.groupby(
            reader,
            key=lambda r: r["VEP_ID"]
        ):

            if (
                not vep_id
                or (
                    previous
                    and vep_id <= previous
                )
            ):
                raise ValueError(
                    "VEP_ID order or duplicate violation"
                )

            previous = vep_id

            first = next(rows)

            locus = str(
                first.get("LOCUS_ID", "")
            ).strip()

            ref = str(
                first.get("REFERENCE_ID", "")
            ).strip()

            if not locus or not ref:
                raise ValueError(
                    "Missing LOCUS_ID or REFERENCE_ID"
                )

            pip = number(
                first.get("PIP")
            )

            if (
                pip is None
                or not 0 <= pip <= 1
            ):
                raise ValueError(
                    "Invalid PIP at "
                    + vep_id
                    + ": "
                    + str(first.get("PIP"))
                )

            key = (
                locus,
                ref,
            )

            is_cs = key in membership

            if is_cs:
                seen_membership.add(key)

            if is_cs != flag(
                first.get("IN_95_CREDIBLE_SET")
            ):
                mismatched_cs_flags += 1

            # Includes the first transcript row.
            features = analyze_transcripts(
                itertools.chain(
                    (first,),
                    rows
                )
            )

            labels = ["REGIONAL"]

            if is_cs:
                labels.append("CS95")

            if pip >= 0.10:
                labels.append(
                    "PIP_GE_0_10"
                )

            if pip >= 0.50:

                labels.append(
                    "PIP_GE_0_50"
                )

                if is_cs:
                    labels.append(
                        "CS95_AND_PIP_GE_0_50"
                    )

            if pip >= 0.90:
                labels.append(
                    "PIP_GE_0_90"
                )

            for label in labels:

                increment(
                    pops[label],
                    features,
                    pip
                )

                increment(
                    locus_pops[locus][label],
                    features,
                    pip
                )

            n_records += 1

    # ======================================================
    # INTEGRITY VALIDATION
    # ======================================================

    expected = integer(
        meta.get("N_SELECTED_VARIANTS")
    )

    if (
        n_records == 0
        or (
            expected is not None
            and n_records != expected
        )
    ):
        raise ValueError(
            "VEP variants mismatch: "
            f"observed={n_records}, expected={expected}"
        )

    missing_members = (
        membership - seen_membership
    )

    if missing_members:

        example = next(
            iter(missing_members)
        )

        raise ValueError(
            f"{len(missing_members)} Step06 credible-set "
            f"members absent from VEP; e.g. {example}"
        )

    if mismatched_cs_flags:
        raise ValueError(
            f"{mismatched_cs_flags} Step06/Step07 "
            "credible-set flag mismatches"
        )

    return (
        pops,
        locus_pops,
        len(membership),
    )


# ==========================================================
# PRINTING
# ==========================================================

def summary_line(name, data):

    n = data["n"]
    c = data["counts"]

    w = data["weight"]
    z = data["weighted"]

    print(
        f"{name:<23}"
        f"{n:>12,}"
        f"{fmt_percent(c['CDS'], n):>11}"
        f"{fmt_percent(c['LNCRNA_BODY'], n):>12}"
        f"{fmt_percent(c['LNCRNA_EXON'], n):>12}"
        f"{fmt_percent(c['LNCRNA_INTRON'], n):>12}"
        f"{fmt_percent(c['LNCRNA_SPLICE'], n):>12}"
    )

    print(
        f"{'':23}"
        f"{'PIP-weighted':>12}"
        f"{fmt_percent(z['CDS'], w):>11}"
        f"{fmt_percent(z['LNCRNA_BODY'], w):>12}"
        f"{fmt_percent(z['LNCRNA_EXON'], w):>12}"
        f"{fmt_percent(z['LNCRNA_INTRON'], w):>12}"
        f"{fmt_percent(z['LNCRNA_SPLICE'], w):>12}"
    )


def print_pops(title, pops):

    print("\n" + "-" * 107)
    print(title)
    print("-" * 107)

    print(
        f"{'POPULATION':<23}"
        f"{'VARIANTS':>12}"
        f"{'CDS %':>11}"
        f"{'lncBODY %':>12}"
        f"{'lncEXON %':>12}"
        f"{'lncINTRON %':>12}"
        f"{'lncSPLICE %':>12}"
    )

    for label in POPS:
        summary_line(
            label,
            pops[label]
        )


# ==========================================================
# LOCUS-BALANCED SENSITIVITY
# ==========================================================

def locus_balanced(
    locus_by_study,
    population,
):

    result = {
        k: []
        for k in (
            "CDS",
            "LNCRNA_BODY",
            "LNCRNA_EXON",
            "LNCRNA_INTRON",
            "LNCRNA_SPLICE",
        )
    }

    n_loci = 0

    for by_locus in locus_by_study.values():

        for stats in by_locus.values():

            d = stats[population]

            if d["weight"] <= 0:
                continue

            n_loci += 1

            for k in result:

                fraction = (
                    d["weighted"][k]
                    / d["weight"]
                )

                result[k].append(
                    fraction
                )

    return n_loci, result


def print_locus_balance(
    locus_by_study,
    label,
):

    for population in (
        "REGIONAL",
        "CS95",
    ):

        n, result = locus_balanced(
            locus_by_study,
            population
        )

        print(
            f"  {label} {population}: "
            f"contributing study-loci={n}"
        )

        for feature in (
            "CDS",
            "LNCRNA_BODY",
            "LNCRNA_EXON",
            "LNCRNA_INTRON",
            "LNCRNA_SPLICE",
        ):

            values = result[feature]

            print(
                f"    {feature:<20}"
                f"{fmt_percent(sum(values), len(values))}"
            )


# ==========================================================
# MAIN
# ==========================================================

def population_rows(disease, subset, pops):
    """Row-counted and PIP-weighted percentages for every population."""
    rows = []
    for label in POPS:
        d = pops[label]
        row = {"PHENOTYPE": disease, "SUBSET": subset, "POPULATION": label, "VARIANTS": d["n"],
               "PIP_SUM": round(d["weight"], 4)}
        for feature in FEATURES:
            row[f"{feature}_PCT"] = round(100 * d["counts"][feature] / d["n"], 3) if d["n"] else None
            row[f"{feature}_PIP_WEIGHTED_PCT"] = (round(100 * d["weighted"][feature] / d["weight"], 3)
                                                  if d["weight"] else None)
        rows.append(row)
    return rows


def run(root, diseases, ancestry="european", **_):
    """Run Analysis 09 for [(display name, folder slug), ...]; return its tables."""
    global DISEASES, ANCESTRY
    DISEASES, ANCESTRY = list(diseases), ancestry
    TABLES.clear()
    root = Path(root).expanduser().resolve()

    print("=" * 107)
    print("JOHN MATTICK - ANALYSIS 09")
    print("FINE-MAPPING SENSITIVITY AND PIP WEIGHTING")
    print("Step06 credible-set membership cross-checked")
    print("against Step07 transcript annotations")
    print("READ ONLY")
    print("=" * 107)

    pooled = make_pops()

    # ======================================================
    # PROCESS EACH DISEASE
    # ======================================================

    for disease, slug in DISEASES:

        print("\n" + "=" * 107)
        print(disease.upper())
        print("=" * 107)

        manifest, path = get_manifest(
            root,
            slug
        )

        if path is None:
            print(
                "No manifest available - "
                "disease not analyzed"
            )
            continue

        audit = (
            root
            / "pipeline_audit"
            / slug
            / ANCESTRY
        )

        disease_pops = make_pops()
        direct_pops = make_pops()

        per_study_pops = []
        per_study_loci = {}

        counters = Counter()

        print(
            f"Selected study accessions: {len(manifest)}"
        )

        print(
            f"{'STUDY':<15}"
            f"{'STATUS':<40}"
            f"{'REGIONAL':>13}"
            f"{'CS95':>11}"
            f"{'PIP>=0.5':>11}"
            f"{'CS PIP SUM':>14}"
        )

        for row in manifest:

            study = str(
                row.get("STUDY_ACCESSION", "")
            ).strip()

            if not study:
                continue

            prepared, reason = check_paths(
                root,
                slug,
                study
            )

            if prepared is None:

                counters[reason] += 1
                TABLES.setdefault("studies", []).append(
                    {"PHENOTYPE": disease, "STUDY_ACCESSION": study, "DECISION": reason})

                print(
                    f"{study:<15}"
                    f"{reason:<40}"
                    f"{'-':>13}"
                    f"{'-':>11}"
                    f"{'-':>11}"
                    f"{'-':>14}"
                )

                continue

            # ----------------------------------------------
            # PROCESS AND VALIDATE THE STUDY
            # ----------------------------------------------

            try:

                stats, loci, n_cs_members = analyze_study(
                    *prepared,
                    study
                )

            except (
                OSError,
                EOFError,
                ValueError,
                UnicodeError,
                csv.Error,
            ) as exc:

                counters["INVALID_INPUT"] += 1

                print(
                    f"{study:<15} "
                    f"INVALID_INPUT: {str(exc)[:130]}"
                )

                continue

            counters["INCLUDED"] += 1

            trait = (
                str(
                    row.get("DISEASE_TRAIT", "")
                ).strip()
                or str(
                    row.get("MAPPED_TRAIT", "")
                ).strip()
            )

            is_direct = direct_trait(
                trait,
                disease
            )

            if is_direct:

                counters[
                    "DIRECT_LABEL_INCLUDED"
                ] += 1

            else:

                counters[
                    "RELATED_OR_REVIEW_LABEL_INCLUDED"
                ] += 1

            per_study_pops.append(
                (study, stats, is_direct)
            )

            per_study_loci[study] = loci

            for population in POPS:

                merge(
                    disease_pops[population],
                    stats[population]
                )

                merge(
                    pooled[population],
                    stats[population]
                )

                if is_direct:

                    merge(
                        direct_pops[population],
                        stats[population]
                    )

            decision = (
                "INCLUDED DIRECT"
                if is_direct
                else "INCLUDED REVIEW"
            )
            TABLES.setdefault("studies", []).append({
                "PHENOTYPE": disease, "STUDY_ACCESSION": study, "DECISION": decision,
                "REGIONAL": stats["REGIONAL"]["n"], "CS95": stats["CS95"]["n"],
                "PIP_GE_0_50": stats["PIP_GE_0_50"]["n"], "CS_PIP_SUM": round(stats["CS95"]["weight"], 4)})

            print(
                f"{study:<15}"
                f"{decision:<40}"
                f"{stats['REGIONAL']['n']:>13,}"
                f"{stats['CS95']['n']:>11,}"
                f"{stats['PIP_GE_0_50']['n']:>11,}"
                f"{stats['CS95']['weight']:>14.3f}",
                flush=True,
            )

            if (
                stats["CS95"]["n"]
                != n_cs_members
            ):

                print(
                    "  NOTE: CS membership count="
                    f"{n_cs_members}; observed rows="
                    f"{stats['CS95']['n']}"
                )

        # ==================================================
        # DISEASE-LEVEL RESULTS
        # ==================================================

        print("\nSTUDY AUDIT")

        for key, value in sorted(
            counters.items()
        ):

            print(
                f"  {key:<44}: {value:,}"
            )

        print_pops(
            "ALL INCLUDED STUDIES: "
            "VARIANT-ROW-WEIGHTED AND PIP-WEIGHTED",
            disease_pops,
        )
        TABLES.setdefault("populations", []).extend(population_rows(disease, "ALL_STUDIES", disease_pops))
        if counters["DIRECT_LABEL_INCLUDED"]:
            TABLES["populations"].extend(population_rows(disease, "DIRECT_LABEL_STUDIES", direct_pops))

        # ==================================================
        # DIRECT DISEASE SENSITIVITY
        # ==================================================

        if counters["DIRECT_LABEL_INCLUDED"]:

            print_pops(
                "SENSITIVITY: EXACT DIRECT DISEASE "
                "PHENOTYPE LABELS ONLY",
                direct_pops,
            )

        else:

            print(
                "\nDIRECT-LABEL SENSITIVITY: "
                "no exact direct disease labels; "
                "review metadata"
            )

        # ==================================================
        # EQUAL-STUDY-WEIGHTED SENSITIVITY
        # ==================================================

        print(
            "\nEQUAL-STUDY-WEIGHTED SENSITIVITY"
        )

        print(
            "Median study-level fraction among "
            "95% credible-set variants"
        )

        for feature in (
            "CDS",
            "LNCRNA_BODY",
            "LNCRNA_EXON",
            "LNCRNA_INTRON",
            "LNCRNA_SPLICE",
        ):

            observations = [
                s["CS95"]["counts"][feature]
                / s["CS95"]["n"]

                for _, s, _ in per_study_pops

                if s["CS95"]["n"] > 0
            ]

            if observations:

                median = statistics.median(
                    observations
                )

                display = fmt_percent(
                    median,
                    1
                )

            else:
                display = "NA"

            print(
                f"  {feature:<23}"
                f"{display:>10}"
                f"   (n={len(observations)} studies)"
            )
            TABLES.setdefault("study_balanced", []).append({
                "PHENOTYPE": disease, "FEATURE": feature, "STUDIES": len(observations),
                "MEDIAN_STUDY_PCT": round(100 * statistics.median(observations), 3) if observations else None})

        # ==================================================
        # EQUAL-LOCUS-WEIGHTED SENSITIVITY
        # ==================================================

        print(
            "\nEQUAL-LOCUS-WEIGHTED SENSITIVITY"
        )

        print(
            "Mean of locus-level "
            "PIP-weighted fractions"
        )

        print_locus_balance(
            per_study_loci,
            disease
        )
        for population in ("REGIONAL", "CS95"):
            n_loci, result = locus_balanced(per_study_loci, population)
            TABLES.setdefault("locus_balanced", []).append({
                "PHENOTYPE": disease, "POPULATION": population, "STUDY_LOCI": n_loci,
                **{f"{k}_PCT": round(100 * sum(v) / len(v), 3) if v else None for k, v in result.items()}})

        print(
            "\nNOTE: Direct disease classification "
            "uses conservative exact phenotype labels."
        )

    # ======================================================
    # COMBINED SUMMARY
    # ======================================================

    print_pops(
        "COMBINED DESCRIPTIVE SUMMARY - ALL PHENOTYPES",
        pooled,
    )
    TABLES.setdefault("populations", []).extend(population_rows("ALL PHENOTYPES POOLED", "ALL_STUDIES", pooled))

    print("\nSCIENTIFIC INTERPRETATION")

    print(
        "1. Main denominator: Step06 95% "
        "credible-set variant membership."
    )

    print(
        "2. PIP thresholds define sensitivity "
        "subsets; high-PIP variants may also "
        "occur outside credible sets."
    )

    print(
        "3. PIP-weighted overlap equals "
        "sum(PIP of overlapping variants) / "
        "sum(PIP of all variants in population)."
    )

    print(
        "4. lncRNA exon, intron and splice "
        "annotations can overlap across transcripts."
    )

    print(
        "5. Regional variants are fine-mapping "
        "inputs, not a genome-wide null background."
    )

    print(
        "6. GWAS studies and genomic loci may "
        "overlap; pooled records are not independent."
    )

    print(
        "7. PIP weights do not represent "
        "phenotypic variance explained or "
        "proof of lncRNA causality."
    )

    return {name: pd.DataFrame(rows) for name, rows in TABLES.items()}


def main(argv=None):
    return analysis_main(run, "John Analysis09: fine-mapping sensitivity", argv)


if __name__ == "__main__":
    main()
 