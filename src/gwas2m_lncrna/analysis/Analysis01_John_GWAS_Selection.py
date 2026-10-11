 #!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
ANALYSIS 01 - JOHN MATTICK
GWAS discovery, selection, QC and analysis eligibility.

Phenotypes and ancestry are passed in: run(root, diseases, ancestry) or
--phenotype / --ancestry on the command line.

READ-ONLY:
  Reads the GWAS2m-lncRNA output folders and prints results to the terminal.
  run() also returns the result tables; `gwas2m-lncrna analysis` saves them.
"""

import csv
import json
import re
from collections import Counter
from pathlib import Path

import pandas as pd

from gwas2m_lncrna.analysis.helpers import analysis_main, direct_labels, normalize_label


# Set by run(): [(display name, phenotype folder), ...] and the ancestry folder name.
DISEASES = []
ANCESTRY = "european"
TABLES = {}

STAGES = [
    ("S02", "Step02_Download"),
    ("S03", "Step03_QC"),
    ("S04", "Step04_LD_Reference"),
    ("S05", "Step05_Clumping"),
    ("S06", "Step06_FineMapping"),
    ("S07", "Step07_VEP"),
]

DONE = {"COMPLETE", "SKIPPED_ALREADY_COMPLETE"}

SHORT = {
    "COMPLETE": "OK",
    "SKIPPED_ALREADY_COMPLETE": "OK",
    "RUNNING": "RUN",
    "PENDING": "PEND",
    "FAILED": "FAIL",
    "PARTIAL": "PART",
    "EXCLUDED": "EXCL",
    "NOT_AVAILABLE": "N/A",
    "NOT_RUN_UPSTREAM_FAILED": "--",
}


def banner(title):
    print("\n" + "=" * 120)
    print(title)
    print("=" * 120)


def clean(x):
    if x is None:
        return ""
    value = str(x).strip()
    if value.lower() in {"nan", "none", "null"}:
        return ""
    return value


def read_tsv(path):
    if not path.is_file():
        return []
    try:
        with path.open("r", encoding="utf-8-sig", newline="") as f:
            return list(csv.DictReader(f, delimiter="\t"))
    except (OSError, UnicodeError, csv.Error) as exc:
        print(f"WARNING: Cannot read {path}: {exc}")
        return []


def read_json(path):
    if not path.is_file():
        return {}
    try:
        with path.open("r", encoding="utf-8") as f:
            data = json.load(f)
        return data if isinstance(data, dict) else {}
    except (OSError, ValueError, UnicodeError):
        return {}


def get_stage(audit, study, stage):
    path = audit / "status" / study / f"{stage}.json"
    record = read_json(path)
    status = clean(record.get("status")).upper()
    return status or "NOT_RECORDED", record


def to_int(value):
    try:
        return int(float(clean(value)))
    except (ValueError, TypeError, OverflowError):
        return None


def valid_vep_table(path):
    """Check that Step07 has a real annotated variant table."""
    if not path.is_file() or path.stat().st_size < 30:
        return False

    required = {
        "VEP_ID", "LOCUS_ID", "PIP", "CHR", "REFERENCE_POS"
    }

    try:
        with path.open("r", encoding="utf-8", newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            columns = set(reader.fieldnames or [])
            return required.issubset(columns) and next(reader, None) is not None
    except (OSError, UnicodeError, csv.Error):
        return False


def classify_trait(trait, disease):
    """Conservative phenotype classification; ambiguous labels need review."""
    t = re.sub(r"\s+", " ", clean(trait).lower())
    t = t.replace("’", "'")

    if any(x in t for x in [
        "family history", "family-history", "parental",
        "mother", "father", "sibling", "proxy"
    ]):
        return "PROXY_FAMILY"

    if any(x in t for x in [
        "drug response", "medication", "treatment response",
        "levodopa"
    ]):
        return "TREATMENT"

    if any(x in t for x in [
        "interaction", "mtag", "inactivation",
        "subtype", "severity", "age at onset", "motor score"
    ]):
        return "SPECIALIZED"

    if any(x in t for x in [
        "dementia in parkinson", "parkinson's dementia",
        "cognitive impairment", "tremor", "dyskinesia"
    ]):
        return "RELATED_TRAIT"

    if normalize_label(t) in direct_labels(disease):
        return "DIRECT_DISEASE"

    return "REVIEW_TRAIT"


def classify_readiness(root, slug, study, statuses):
    s06 = statuses["S06"][0]
    s07 = statuses["S07"][0]

    fm_path = (
        root / "06_finemapping" / slug / ANCESTRY
        / study / "finemapping_summary.json"
    )

    fm = read_json(fm_path)

    # Legitimate study outcome, not a technical failure.
    if s06 in DONE:
        if clean(fm.get("STATUS")).upper() == "NO_SIGNIFICANT_LOCI":
            return "NO_SIGNIFICANT_LOCI"

        if fm and to_int(fm.get("N_LOCI")) == 0:
            return "NO_SIGNIFICANT_LOCI"

    vep_path = (
        root / "07_annotation" / slug / ANCESTRY
        / study / f"{study}_VEP_variant_summary.tsv"
    )

    if s07 in DONE:
        if valid_vep_table(vep_path):
            return "ANNOTATED_READY"
        return "ANNOTATION_INVALID"

    if s06 == "PARTIAL":
        return "PARTIAL_FINEMAPPING"

    if s07 == "FAILED":
        return "VEP_FAILED"

    if s07 in {"RUNNING", "PENDING"}:
        return "VEP_IN_PROGRESS"

    for step in ["S02", "S03", "S04", "S05", "S06"]:
        if statuses[step][0] in {
            "FAILED", "EXCLUDED", "NOT_RUN_UPSTREAM_FAILED"
        }:
            return "UPSTREAM_FAILED"

    if s06 in DONE:
        return "AWAITING_ANNOTATION"

    return "PROCESSING_OR_PENDING"


def print_table(headers, rows, widths):
    fmt = "  ".join("{:<" + str(x) + "}" for x in widths)
    print(fmt.format(*headers))
    print("-" * (sum(widths) + 2 * (len(widths) - 1)))

    for row in rows:
        print(fmt.format(*[
            str(value)[:width]
            for value, width in zip(row, widths)
        ]))


def first_blocker(statuses):
    for code, stage in STAGES:
        status, record = statuses[code]

        if status in {"FAILED", "PARTIAL", "EXCLUDED"}:
            reason = clean(record.get("reason"))
            error = clean(record.get("error_message"))

            detail = reason or error or "No detailed reason stored"

            if error and error not in detail:
                detail += "; " + error

            return (
                stage, status, detail,
                clean(record.get("log_path"))
            )

    return "", "", "", ""


def analyse_disease(root, disease, slug, brief=False):
    audit = root / "pipeline_audit" / slug / ANCESTRY
    catalog = root / "01_gwas_catalog" / slug / ANCESTRY

    manifest_file = audit / "study_manifest.tsv"
    if not manifest_file.is_file():
        manifest_file = catalog / "GWAS_download_manifest.tsv"

    selected = read_tsv(manifest_file)

    selection = read_tsv(catalog / "study_selection.tsv")
    if not selection:
        selection = read_tsv(audit / "study_selection.tsv")

    discovery = read_json(catalog / "discovery_summary.json")

    banner(f"ANALYSIS 01 | {disease.upper()}")

    print("Selected manifest :", manifest_file)
    print("Audit directory   :", audit)
    print("Selected GWAS     :", len(selected))
    print("Phenotype matches :", len(selection) if selection else "Unknown")
    print(
        "Catalog scanned   :",
        discovery.get("N_CATALOG_STUDIES", "Unknown")
    )

    if not selected:
        print("WARNING: No selected study manifest available.")
        return Counter(), []

    # -------------------------------------------------------
    # DISCOVERY EXCLUSIONS
    # -------------------------------------------------------

    print("\nSTEP01 DISCOVERY EXCLUSIONS")

    if selection:
        excluded = Counter(
            clean(row.get("EXCLUSION_REASON")) or "UNSPECIFIED"
            for row in selection
            if clean(row.get("ELIGIBLE")).lower()
            not in {"true", "1", "yes"}
        )

        print("Total excluded:", sum(excluded.values()))

        for reason, count in excluded.most_common():
            print(f"  {reason:<40} {count:>8,}")

        TABLES.setdefault("selection_exclusions", []).extend(
            {"PHENOTYPE": disease, "EXCLUSION_REASON": reason, "N_STUDIES": count}
            for reason, count in excluded.most_common()
        )

        audit_selected = sum(
            clean(row.get("ELIGIBLE")).lower() in {"true", "1", "yes"}
            for row in selection
        )

        if audit_selected != len(selected):
            print(
                "WARNING: Selection audit and manifest disagree:",
                audit_selected, "versus", len(selected)
            )
    else:
        print("Selection audit unavailable; exclusion counts not calculated.")

    # -------------------------------------------------------
    # STUDY STATUS / READINESS
    # -------------------------------------------------------

    counts = Counter()
    trait_counts = Counter()

    stage_counts = {
        code: Counter() for code, _ in STAGES
    }

    records = []
    display = []

    for metadata in selected:
        study = clean(metadata.get("STUDY_ACCESSION"))
        if not study:
            continue

        trait = (
            clean(metadata.get("DISEASE_TRAIT"))
            or clean(metadata.get("MAPPED_TRAIT"))
        )

        trait_type = classify_trait(trait, disease)

        statuses = {
            code: get_stage(audit, study, stage)
            for code, stage in STAGES
        }

        for code, _ in STAGES:
            stage_counts[code][statuses[code][0]] += 1

        readiness = classify_readiness(
            root, slug, study, statuses
        )

        counts[readiness] += 1
        trait_counts[trait_type] += 1

        record = {
            "study": study,
            "trait": trait,
            "trait_type": trait_type,
            "readiness": readiness,
            "statuses": statuses,
            "metadata": metadata,
        }

        records.append(record)

        display.append(
            [study, trait_type]
            + [
                SHORT.get(statuses[code][0], "?")
                for code, _ in STAGES
            ]
            + [readiness]
        )

    print("\nALL SELECTED GWAS: PIPELINE STATUS")

    headers = (
        ["STUDY", "TRAIT_CLASS"]
        + [code for code, _ in STAGES]
        + ["ANALYSIS_STATUS"]
    )

    widths = [15, 16] + [5] * len(STAGES) + [23]

    print_table(headers, display, widths)

    # -------------------------------------------------------
    # STUDY INCLUSION / EXCLUSION SUMMARY
    # -------------------------------------------------------

    print("\nSTUDY ANALYSIS READINESS")

    categories = [
        "ANNOTATED_READY",
        "NO_SIGNIFICANT_LOCI",
        "PARTIAL_FINEMAPPING",
        "UPSTREAM_FAILED",
        "VEP_FAILED",
        "ANNOTATION_INVALID",
        "VEP_IN_PROGRESS",
        "AWAITING_ANNOTATION",
        "PROCESSING_OR_PENDING",
    ]

    for category in categories:
        if counts[category]:
            print(f"  {category:<30} {counts[category]:>5,}")

    print("\nPHENOTYPE CLASSIFICATION")

    for category, count in trait_counts.most_common():
        print(f"  {category:<30} {count:>5,}")

    # -------------------------------------------------------
    # STAGE COMPLETION COUNTS
    # -------------------------------------------------------

    print("\nPIPELINE STAGE COMPLETION")

    for code, stage in STAGES:
        c = stage_counts[code]

        complete = sum(c[x] for x in DONE)

        print(
            f"  {code} {stage:<25}"
            f" OK={complete:>3}"
            f" RUN={c['RUNNING']:>3}"
            f" PEND={c['PENDING']:>3}"
            f" PART={c['PARTIAL']:>3}"
            f" FAIL={c['FAILED']:>3}"
            f" N/A={c['NOT_AVAILABLE']:>3}"
        )

    # -------------------------------------------------------
    # STUDIES ELIGIBLE FOR JOHN'S ANALYSIS
    # -------------------------------------------------------

    ready = [
        r for r in records
        if r["readiness"] == "ANNOTATED_READY"
    ]

    direct = [
        r for r in ready
        if r["trait_type"] == "DIRECT_DISEASE"
    ]

    other = [
        r for r in ready
        if r["trait_type"] != "DIRECT_DISEASE"
    ]

    print("\nJOHN MATTICK: GENOMIC DISTRIBUTION ANALYSIS")

    print("Annotated studies available :", len(ready))
    print("Direct disease labels       :", len(direct))
    print("Other/review trait labels   :", len(other))
    print(
        "No significant loci         :",
        counts["NO_SIGNIFICANT_LOCI"]
    )

    print("\nREADY GWAS ACCESSIONS:")

    if ready:
        for row in ready:
            print(
                f"  {row['study']:<15}"
                f" {row['trait_type']:<18}"
                f" {row['trait']}"
            )
    else:
        print("  None")

    # -------------------------------------------------------
    # FULL STUDY METADATA
    # -------------------------------------------------------

    if not brief:
        print("\nSTUDY PHENOTYPES AND SAMPLE SIZE")

        print_table(
            ["STUDY", "REPORTED PHENOTYPE", "REPORTED SAMPLE SIZE", "STATUS"],
            [
                [
                    r["study"],
                    r["trait"],
                    r["metadata"].get("INITIAL_SAMPLE_SIZE", ""),
                    r["readiness"],
                ]
                for r in records
            ],
            [15, 62, 35, 26],
        )

    # -------------------------------------------------------
    # FAILED / INCOMPLETE STUDIES
    # -------------------------------------------------------

    remaining = [
        r for r in records
        if r["readiness"] not in {
            "ANNOTATED_READY", "NO_SIGNIFICANT_LOCI"
        }
    ]

    print("\nINCOMPLETE OR FAILED STUDIES")

    if not remaining:
        print("  None")
    else:
        for row in remaining:
            stage, status, reason, log = first_blocker(
                row["statuses"]
            )

            print(
                f"\n  {row['study']}"
                f" | {row['readiness']}"
            )

            if stage:
                print("    Stage :", stage, status)
                print("    Reason:", reason)
            elif row["readiness"] == "VEP_IN_PROGRESS":
                print("    Reason: Step07 annotation not complete")
            else:
                print("    Reason: Awaiting valid upstream outputs")

            if log:
                print("    Log   :", log)

    for row in records:
        stage, status, reason, log = first_blocker(row["statuses"])
        TABLES.setdefault("studies", []).append({
            "PHENOTYPE": disease, "STUDY_ACCESSION": row["study"], "TRAIT": row["trait"],
            "TRAIT_CLASS": row["trait_type"],
            "INITIAL_SAMPLE_SIZE": row["metadata"].get("INITIAL_SAMPLE_SIZE", ""),
            **{name: row["statuses"][code][0] for code, name in STAGES},
            "ANALYSIS_STATUS": row["readiness"], "BLOCKING_STAGE": stage, "BLOCKING_REASON": reason,
        })

    print("\nINTERPRETATION")
    print(
        "  ANNOTATED_READY: Can be used for John's "
        "coding/noncoding distribution."
    )
    print(
        "  NO_SIGNIFICANT_LOCI: Valid processed study "
        "with no loci available for fine-mapping."
    )
    print(
        "  PARTIAL/FAILED: Exclude from primary "
        "annotation percentages until resolved."
    )
    print(
        "  RELATED/PROXY TRAITS: Review separately "
        "from direct disease GWAS."
    )

    return counts, records


def run(root, diseases, ancestry="european", brief=False, **_):
    """Run Analysis 01 for [(display name, folder slug), ...]; return its tables."""
    global DISEASES, ANCESTRY
    DISEASES, ANCESTRY = list(diseases), ancestry
    TABLES.clear()
    root = Path(root).expanduser().resolve()

    banner("JOHN MATTICK - ANALYSIS 01")

    print("GWAS2m root :", root)
    print("Mode        : READ ONLY")

    overall = Counter()
    all_records = []

    for disease, slug in DISEASES:
        counts, records = analyse_disease(
            root, disease, slug, brief
        )

        overall.update(counts)
        all_records.extend(records)

    # -------------------------------------------------------
    # OVERALL THREE-DISEASE SUMMARY
    # -------------------------------------------------------

    banner("OVERALL SUMMARY - ALL PHENOTYPES")

    print("Selected studies          :", len(all_records))
    print(
        "Ready for John's analysis :",
        overall["ANNOTATED_READY"]
    )
    print(
        "No significant loci       :",
        overall["NO_SIGNIFICANT_LOCI"]
    )
    print(
        "Partial fine-mapping      :",
        overall["PARTIAL_FINEMAPPING"]
    )
    print(
        "Failed upstream           :",
        overall["UPSTREAM_FAILED"]
    )

    accounted = sum(overall.values())
    print(
        "Other pending/invalid     :",
        accounted
        - overall["ANNOTATED_READY"]
        - overall["NO_SIGNIFICANT_LOCI"]
        - overall["PARTIAL_FINEMAPPING"]
        - overall["UPSTREAM_FAILED"]
    )

    print("\nSCIENTIFIC NOTES")
    print(
        "1. GWAS Catalog accessions are distinct datasets, "
        "not necessarily independent cohorts."
    )
    print(
        "2. Direct disease phenotypes and proxy phenotypes "
        "must be reviewed before pooling."
    )
    print(
        "3. No-significant-loci studies are counted in "
        "selection reporting, not the annotation denominator."
    )
    print(
        "4. This analysis does not estimate the percentage "
        "of phenotypic variation explained."
    )
    print(
        "5. Step07-complete studies provide the input "
        "for coding/noncoding variant distribution."
    )

    return {name: pd.DataFrame(rows) for name, rows in TABLES.items()}


def main(argv=None):
    return analysis_main(run, "John Analysis01: GWAS selection, QC and readiness", argv,
                         extra=[("--brief", dict(action="store_true", help="Hide detailed study metadata"))])


if __name__ == "__main__":
    main()
 