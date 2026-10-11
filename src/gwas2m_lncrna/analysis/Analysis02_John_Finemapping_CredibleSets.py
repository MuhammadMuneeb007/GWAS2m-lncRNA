 #!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
ANALYSIS 02 - JOHN MATTICK
Fine-mapped GWAS loci, SuSiE PIP and 95% credible sets.

Reads existing Step05 and Step06 files for the phenotypes passed in.
Prints all results to the terminal; run() also returns the result tables.

Default: include only studies with COMPLETE Step06 status.
Optional: --include-partial for exploratory analysis.
"""

import argparse
import csv
import gzip
import heapq
import json
import math
import statistics
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd

from gwas2m_lncrna.analysis.helpers import analysis_main


# Set by run(): [(display name, phenotype folder), ...] and the ancestry folder name.
DISEASES = []
ANCESTRY = "european"
TABLES = {}

DONE = {"COMPLETE", "SKIPPED_ALREADY_COMPLETE"}


def read_json(path):
    try:
        with path.open(encoding="utf-8") as f:
            obj = json.load(f)
        return obj if isinstance(obj, dict) else {}
    except (OSError, ValueError):
        return {}


def read_tsv(path):
    if not path.is_file():
        return None

    try:
        with path.open(
            encoding="utf-8-sig", newline=""
        ) as f:
            return list(
                csv.DictReader(f, delimiter="\t")
            )
    except (OSError, UnicodeError, csv.Error) as e:
        print(f"WARNING: Cannot read {path}: {e}")
        return None


def num(value):
    try:
        n = float(str(value).strip())
        return n if math.isfinite(n) else None
    except (ValueError, TypeError):
        return None


def integer(value, default=0):
    n = num(value)
    return int(n) if n is not None else default


def yes(value):
    return str(value).strip().lower() in {
        "true", "1", "yes", "y", "t"
    }


def stage_status(root, disease, study, name):
    path = (
        root / "pipeline_audit" / disease
        / ANCESTRY / "status" / study
        / f"{name}.json"
    )

    obj = read_json(path)
    return str(
        obj.get("status", "NOT_RECORDED")
    ).upper()


def add_top(heap, pip, row, limit):
    if pip is None:
        return

    item = (
        pip,
        row["STUDY"],
        row["LOCUS"],
        row["VARIANT"],
        row,
    )

    if len(heap) < limit:
        heapq.heappush(heap, item)
    elif item[:4] > heap[0][:4]:
        heapq.heapreplace(heap, item)


def analyse_study(root, disease, study, top_n):
    folder = (
        root / "06_finemapping" / disease
        / ANCESTRY / study
    )

    summary = read_json(
        folder / "finemapping_summary.json"
    )

    loci = read_tsv(
        folder / "loci_definition.tsv"
    )

    fits = read_tsv(
        folder / f"{study}_finemapping_summary.tsv"
    )

    cs_rows = read_tsv(
        folder / f"{study}_95pct_credible_sets.tsv"
    )

    variant_file = (
        folder / f"{study}_finemapped_variants.tsv.gz"
    )

    if (
        not summary
        or loci is None
        or fits is None
        or cs_rows is None
        or not variant_file.is_file()
    ):
        return None, "Missing Step06 output files"

    n_defined = integer(
        summary.get("N_LOCI"), len(loci)
    )

    n_success = integer(
        summary.get("N_LOCI_SUCCESS"), len(fits)
    )

    n_failed = integer(
        summary.get("N_LOCI_FAILED")
    )

    if n_success == 0:
        return None, "No successfully fine-mapped loci"

    if len(fits) != n_success:
        return None, "Inconsistent locus summary counts"

    locus_info = {
        str(r.get("LOCUS_ID", "")): r
        for r in loci
    }

    fit_info = {
        str(r.get("LOCUS_ID", "")): r
        for r in fits
    }

    if (
        len(fit_info) != n_success
        or any(k not in locus_info for k in fit_info)
    ):
        return None, "Invalid or duplicated locus IDs"

    # Each credible set has its own members.
    # A variant can belong to more than one set.
    cs_members = defaultdict(set)

    for row in cs_rows:
        locus = str(row.get("LOCUS_ID", ""))
        cs = str(row.get("CREDIBLE_SET", ""))
        variant = str(row.get("REFERENCE_ID", ""))

        if locus in fit_info and cs and variant:
            cs_members[(locus, cs)].add(variant)

    credible_lookup = defaultdict(set)

    for (locus, cs), members in cs_members.items():
        credible_lookup[locus].update(members)

    counts = Counter()
    locus_counts = defaultdict(Counter)
    top_variants = []
    observed_loci = set()
    flag_mismatches = 0

    # Stream the large compressed file without loading
    # millions of variants into RAM.
    try:
        with gzip.open(
            variant_file,
            "rt",
            encoding="utf-8",
            newline=""
        ) as f:

            reader = csv.DictReader(
                f, delimiter="\t"
            )

            required = {
                "LOCUS_ID", "REFERENCE_ID", "PIP"
            }

            if not required.issubset(
                set(reader.fieldnames or [])
            ):
                return None, "Missing required variant columns"

            for row in reader:
                locus = str(
                    row.get("LOCUS_ID", "")
                )

                if locus not in fit_info:
                    return None, (
                        "Unexpected locus in variants: "
                        + locus
                    )

                observed_loci.add(locus)

                counts["VARIANTS"] += 1
                locus_counts[locus]["VARIANTS"] += 1

                pip = num(row.get("PIP"))

                if pip is not None:
                    counts["PIP_SCORED"] += 1

                    if pip >= 0.01:
                        counts["PIP_001"] += 1

                    if pip >= 0.10:
                        counts["PIP_010"] += 1

                    if pip >= 0.50:
                        counts["PIP_050"] += 1

                    if pip >= 0.90:
                        counts["PIP_090"] += 1

                variant_id = str(
                    row.get("REFERENCE_ID", "")
                )

                in_cs = (
                    variant_id in credible_lookup[locus]
                )

                if in_cs:
                    counts["CS_VARIANTS"] += 1
                    locus_counts[locus]["CS_VARIANTS"] += 1

                if yes(
                    row.get("IN_95_CREDIBLE_SET")
                ) != in_cs:
                    flag_mismatches += 1

                record = {
                    "STUDY": study,
                    "LOCUS": locus,
                    "VARIANT": variant_id,
                    "CHR": str(row.get("CHR", "")),
                    "POS": str(
                        row.get(
                            "REFERENCE_POS",
                            row.get("POS", "")
                        )
                    ),
                    "PIP": pip,
                    "CS": in_cs,
                }

                add_top(
                    top_variants, pip, record, top_n
                )

    except (OSError, EOFError, UnicodeError, csv.Error) as e:
        return None, f"Cannot read variant file: {e}"

    if observed_loci != set(fit_info):
        return None, "Successful loci missing from variants"

    if flag_mismatches:
        print(
            f"WARNING {study}: "
            f"{flag_mismatches:,} credible-set flag "
            "disagreements; using membership TSV"
        )

    converged = Counter()
    regions = []

    for locus, fit in fit_info.items():
        cv = str(
            fit.get("CONVERGED", "")
        ).lower().strip()

        if cv in {"true", "t", "1", "yes"}:
            converged["YES"] += 1
        elif cv in {"false", "f", "0", "no"}:
            converged["NO"] += 1
        else:
            converged["UNKNOWN"] += 1

        info = locus_info[locus]

        regions.append({
            "STUDY": study,
            "LOCUS": locus,
            "CHR": str(info.get("CHR", "")),
            "START": str(info.get("LOCUS_START", "")),
            "END": str(info.get("LOCUS_END", "")),
            "LEADS": str(info.get("LEAD_IDS", "")),
            "LEAD_P": num(info.get("MIN_LEAD_P")),
            "MAX_PIP": num(fit.get("MAX_PIP")),
            "N_CS": sum(
                1 for key in cs_members
                if key[0] == locus
            ),
            "N_CS_VARIANTS":
                locus_counts[locus]["CS_VARIANTS"],
            "N_VARIANTS":
                locus_counts[locus]["VARIANTS"],
        })

    return {
        "N_DEFINED": n_defined,
        "N_SUCCESS": n_success,
        "N_FAILED": n_failed,
        "N_LEADS": sum(
            integer(r.get("N_INDEPENDENT_SIGNALS"))
            for r in loci
        ),
        "N_CS": len(cs_members),
        "CS_SIZES": [
            len(x) for x in cs_members.values()
        ],
        "COUNTS": counts,
        "CONVERGED": converged,
        "REGIONS": regions,
        "TOP": top_variants,
    }, ""


def run(root, diseases, ancestry="european", top=15, include_partial=False, **_):
    """Run Analysis 02 for [(display name, folder slug), ...]; return its tables."""
    global DISEASES, ANCESTRY
    DISEASES, ANCESTRY = list(diseases), ancestry
    TABLES.clear()
    args = argparse.Namespace(top=max(1, int(top)), include_partial=bool(include_partial))
    root = Path(root).resolve()

    print("=" * 110)
    print("JOHN MATTICK - ANALYSIS 02")
    print("FINE-MAPPED LOCI AND 95% CREDIBLE SETS")
    print("READ ONLY")
    print("=" * 110)

    grand = Counter()
    grand_cs_sizes = []

    for disease, slug in DISEASES:
        print("\n" + "=" * 110)
        print(disease.upper())
        print("=" * 110)

        manifest = read_tsv(
            root / "01_gwas_catalog" / slug
            / ANCESTRY / "GWAS_download_manifest.tsv"
        )

        if manifest is None:
            manifest = read_tsv(
                root / "pipeline_audit" / slug
                / ANCESTRY / "study_manifest.tsv"
            )

        if manifest is None:
            print("No selected GWAS manifest found")
            continue

        tally = Counter()
        totals = Counter()
        convergence = Counter()

        all_regions = []
        all_top = []
        cs_sizes = []

        print("Selected GWAS:", len(manifest))

        print(
            "\nSTUDY          S05  S06  S07  "
            "LOCI  OK  FAIL   CS   VARIANTS   "
            "CS-VARIANTS   DECISION"
        )

        for meta in manifest:
            study = str(
                meta.get("STUDY_ACCESSION", "")
            ).strip()

            if not study:
                continue

            s05 = stage_status(
                root, slug, study, "Step05_Clumping"
            )

            s06 = stage_status(
                root, slug, study, "Step06_FineMapping"
            )

            s07 = stage_status(
                root, slug, study, "Step07_VEP"
            )

            fm = read_json(
                root / "06_finemapping" / slug
                / ANCESTRY / study
                / "finemapping_summary.json"
            )

            defined = integer(fm.get("N_LOCI"))
            success = integer(fm.get("N_LOCI_SUCCESS"))
            failed = integer(fm.get("N_LOCI_FAILED"))

            no_loci = (
                s06 in DONE
                and (
                    str(fm.get("STATUS", "")).upper()
                    == "NO_SIGNIFICANT_LOCI"
                    or (fm and defined == 0)
                )
            )

            n_cs = "-"
            n_variants = "-"
            n_cs_variants = "-"

            if no_loci:
                decision = "NO_SIGNIFICANT_LOCI"
                tally["NO_LOCI"] += 1

            elif (
                s06 == "PARTIAL"
                and not args.include_partial
            ):
                decision = "PARTIAL_EXCLUDED"
                tally["PARTIAL_EXCLUDED"] += 1

            elif s06 in DONE or (
                args.include_partial and s06 == "PARTIAL"
            ):
                result, error = analyse_study(
                    root, slug, study, args.top
                )

                if result is None:
                    decision = "INVALID_OUTPUT"
                    tally["INVALID"] += 1
                    print(
                        f"WARNING {study}: {error}"
                    )

                else:
                    decision = (
                        "INCLUDED"
                        if s06 in DONE
                        else "PARTIAL_EXPLORATORY"
                    )

                    tally["INCLUDED"] += 1

                    totals.update(result["COUNTS"])

                    totals["LEADS"] += result["N_LEADS"]
                    totals["LOCI_DEFINED"] += result["N_DEFINED"]
                    totals["LOCI_SUCCESS"] += result["N_SUCCESS"]
                    totals["LOCI_FAILED"] += result["N_FAILED"]
                    totals["CREDIBLE_SETS"] += result["N_CS"]

                    convergence.update(result["CONVERGED"])
                    cs_sizes.extend(result["CS_SIZES"])

                    all_regions.extend(result["REGIONS"])

                    for item in result["TOP"]:
                        add_top(
                            all_top,
                            item[0],
                            item[-1],
                            args.top
                        )

                    n_cs = f'{result["N_CS"]:,}'
                    n_variants = (
                        f'{result["COUNTS"]["VARIANTS"]:,}'
                    )
                    n_cs_variants = (
                        f'{result["COUNTS"]["CS_VARIANTS"]:,}'
                    )

            else:
                decision = (
                    "FAILED"
                    if s06 == "FAILED"
                    else "PENDING_OR_BLOCKED"
                )

                tally[
                    "FAILED"
                    if s06 == "FAILED"
                    else "PENDING"
                ] += 1

            def short(status):
                return {
                    "COMPLETE": "OK",
                    "SKIPPED_ALREADY_COMPLETE": "OK",
                    "PARTIAL": "PART",
                    "FAILED": "FAIL",
                    "RUNNING": "RUN",
                    "PENDING": "PEND",
                    "NOT_AVAILABLE": "N/A",
                    "NOT_RUN_UPSTREAM_FAILED": "--",
                }.get(status, "?")

            print(
                f"{study:14} "
                f"{short(s05):5} "
                f"{short(s06):5} "
                f"{short(s07):5} "
                f"{defined:5} "
                f"{success:3} "
                f"{failed:5} "
                f"{n_cs:>5} "
                f"{n_variants:>10} "
                f"{n_cs_variants:>13} "
                f"  {decision}"
            )

            TABLES.setdefault("studies", []).append({
                "PHENOTYPE": disease, "STUDY_ACCESSION": study, "STEP05": s05, "STEP06": s06,
                "STEP07": s07, "LOCI_DEFINED": defined, "LOCI_SUCCESS": success, "LOCI_FAILED": failed,
                "CREDIBLE_SETS": n_cs, "VARIANTS": n_variants, "CS_VARIANTS": n_cs_variants,
                "DECISION": decision,
            })

        # ----------------------------------------------------
        # SUMMARY
        # ----------------------------------------------------

        print("\n" + "-" * 110)
        print("DISEASE-LEVEL SCIENTIFIC RESULTS")
        print("-" * 110)

        metrics = [
            ("Selected GWAS", len(manifest)),
            ("Included studies", tally["INCLUDED"]),
            ("No significant loci", tally["NO_LOCI"]),
            ("Partial excluded", tally["PARTIAL_EXCLUDED"]),
            ("Failed studies", tally["FAILED"]),
            ("Pending or blocked", tally["PENDING"]),
            ("Invalid Step06 outputs", tally["INVALID"]),
            ("Independent lead signals", totals["LEADS"]),
            ("Defined loci", totals["LOCI_DEFINED"]),
            ("Successfully fine-mapped loci", totals["LOCI_SUCCESS"]),
            ("Failed loci in included studies", totals["LOCI_FAILED"]),
            ("95% credible sets", totals["CREDIBLE_SETS"]),
            ("Fine-mapped regional variant rows", totals["VARIANTS"]),
            ("Variants in 95% credible sets", totals["CS_VARIANTS"]),
            ("Variants with PIP >= 0.01", totals["PIP_001"]),
            ("Variants with PIP >= 0.10", totals["PIP_010"]),
            ("Variants with PIP >= 0.50", totals["PIP_050"]),
            ("Variants with PIP >= 0.90", totals["PIP_090"]),
        ]

        for name, value in metrics:
            print(f"{name:<45}: {value:,}")

        summary = {"PHENOTYPE": disease, **{name: value for name, value in metrics}}
        if cs_sizes:
            summary.update({"Credible-set size median": statistics.median(cs_sizes),
                            "Credible-set size min": min(cs_sizes), "Credible-set size max": max(cs_sizes)})
        summary.update({"SuSiE converged": convergence["YES"], "SuSiE not converged": convergence["NO"],
                        "SuSiE convergence unknown": convergence["UNKNOWN"]})
        TABLES.setdefault("phenotype_summary", []).append(summary)

        if totals["VARIANTS"]:
            pct_cs = (
                100 * totals["CS_VARIANTS"]
                / totals["VARIANTS"]
            )

            pct_pip = (
                100 * totals["PIP_050"]
                / totals["VARIANTS"]
            )

            print(
                f"Percentage in credible sets             : "
                f"{pct_cs:.3f}%"
            )

            print(
                f"Percentage with PIP >= 0.50             : "
                f"{pct_pip:.3f}%"
            )

        if cs_sizes:
            print(
                "Credible-set size (median/min/max)      : "
                f"{statistics.median(cs_sizes):g} / "
                f"{min(cs_sizes)} / {max(cs_sizes)}"
            )

        print(
            "SuSiE convergence (yes/no/unknown)      : "
            f"{convergence['YES']}/"
            f"{convergence['NO']}/"
            f"{convergence['UNKNOWN']}"
        )

        # ----------------------------------------------------
        # TOP LOCI
        # ----------------------------------------------------

        print("\nTOP FINE-MAPPED GENOMIC LOCI")
        print(
            "STUDY          LOCUS  CHR:START-END          "
            "LEAD_VARIANT          GWAS_P       "
            "MAX_PIP   CS  CS_VARIANTS"
        )

        ranked = sorted(
            all_regions,
            key=lambda x: (
                x["MAX_PIP"]
                if x["MAX_PIP"] is not None
                else -1
            ),
            reverse=True
        )

        TABLES.setdefault("top_loci", []).extend(
            {"PHENOTYPE": disease, **row} for row in ranked[:args.top]
        )

        for row in ranked[:args.top]:
            region = (
                f'{row["CHR"]}:'
                f'{row["START"]}-{row["END"]}'
            )

            lead_p = (
                f'{row["LEAD_P"]:.2g}'
                if row["LEAD_P"] is not None
                else "NA"
            )

            max_pip = (
                f'{row["MAX_PIP"]:.4g}'
                if row["MAX_PIP"] is not None
                else "NA"
            )

            print(
                f'{row["STUDY"]:14} '
                f'{row["LOCUS"]:6} '
                f'{region[:22]:22} '
                f'{row["LEADS"][:20]:20} '
                f'{lead_p:>10} '
                f'{max_pip:>12} '
                f'{row["N_CS"]:3} '
                f'{row["N_CS_VARIANTS"]:11}'
            )

        # ----------------------------------------------------
        # TOP VARIANTS
        # ----------------------------------------------------

        print("\nTOP VARIANTS BY POSTERIOR INCLUSION PROBABILITY")
        print(
            "STUDY          LOCUS  CHR:POS       "
            "VARIANT_ID                  PIP      CS95"
        )

        TABLES.setdefault("top_variants", []).extend(
            {"PHENOTYPE": disease, **row} for *_, row in sorted(all_top, reverse=True)
        )

        for pip, _, _, _, row in sorted(
            all_top, reverse=True
        ):
            position = (
                f'{row["CHR"]}:{row["POS"]}'
            )

            print(
                f'{row["STUDY"]:14} '
                f'{row["LOCUS"]:6} '
                f'{position[:13]:13} '
                f'{row["VARIANT"][:26]:26} '
                f'{pip:9.5f} '
                f'{"YES" if row["CS"] else "NO"}'
            )

        grand.update(totals)
        grand["STUDIES_SELECTED"] += len(manifest)
        grand["STUDIES_INCLUDED"] += tally["INCLUDED"]
        grand["STUDIES_NO_LOCI"] += tally["NO_LOCI"]

        grand_cs_sizes.extend(cs_sizes)

    # --------------------------------------------------------
    # COMBINED SUMMARY
    # --------------------------------------------------------

    print("\n" + "=" * 110)
    print("COMBINED RESULTS - ALL PHENOTYPES")
    print("=" * 110)

    for key in [
        "STUDIES_SELECTED",
        "STUDIES_INCLUDED",
        "STUDIES_NO_LOCI",
        "LEADS",
        "LOCI_DEFINED",
        "LOCI_SUCCESS",
        "LOCI_FAILED",
        "CREDIBLE_SETS",
        "VARIANTS",
        "CS_VARIANTS",
        "PIP_010",
        "PIP_050",
        "PIP_090",
    ]:
        print(f"{key:<35}: {grand[key]:,}")

    if grand_cs_sizes:
        print(
            "Median credible-set size            : "
            f"{statistics.median(grand_cs_sizes):g}"
        )

    print("\nSCIENTIFIC INTERPRETATION")
    print(
        "1. All variant counts are study-locus-level "
        "counts, not globally unique variants."
    )
    print(
        "2. Credible-set membership indicates statistical "
        "support, not confirmed causality."
    )
    print(
        "3. Different GWAS studies may share cohorts "
        "and genomic loci."
    )
    print(
        "4. Coding/noncoding and lncRNA percentages "
        "are analysed in Analyses 03-06."
    )
    print(
        "5. Partial studies, if included, remain "
        "exploratory and incomplete."
    )

    return {name: pd.DataFrame(rows) for name, rows in TABLES.items()}


def main(argv=None):
    return analysis_main(run, "John Analysis02: fine-mapped loci and 95% credible sets", argv, extra=[
        ("--top", dict(type=int, default=15)),
        ("--include-partial", dict(action="store_true", help="Exploratory inclusion of successful partial loci")),
    ])


if __name__ == "__main__":
    main()
 