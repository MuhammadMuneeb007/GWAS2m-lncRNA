#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""John Mattick - Analysis 06: lncRNA locus representation (read only).

Outputs are printed to stdout; run() also returns the result tables.
Phenotypes and ancestry are passed in (run() / --phenotype, --ancestry).
A locus is one successful SuSiE locus within one GWAS accession.
Two distinct interpretations:
  1. GENCODE gene-body overlap with the full fine-mapping locus interval.
  2. VEP transcript overlap of 95% credible-set variants in that locus.
"""

import argparse
import bisect
import csv
import gzip
import itertools
import json
import math
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd

from gwas2m_lncrna.analysis.helpers import analysis_main

# Set by run(): [(display name, phenotype folder), ...] and the ancestry folder name.
DISEASES = []
ANCESTRY = "european"
TABLES = {}
DONE = {"COMPLETE", "SKIPPED_ALREADY_COMPLETE"}
CLASSES = ("BOTH", "LNCRNA_ONLY", "PROTEIN_CODING_ONLY", "NEITHER")
LNC_BIOTYPES = {"lncrna", "lincrna"}
LEGACY = {"antisense", "processed_transcript", "sense_intronic", "sense_overlapping", "3prime_overlapping_ncrna", "bidirectional_promoter_lncrna", "macro_lncrna"}
# A VEP transcript must be intersected, not just nearby.
IN_TRANSCRIPT = {
    "intron_variant", "non_coding_transcript_variant", "non_coding_transcript_exon_variant",
    "splice_acceptor_variant", "splice_donor_variant", "splice_region_variant",
    "splice_donor_region_variant", "splice_donor_5th_base_variant",
    "splice_polypyrimidine_tract_variant", "coding_sequence_variant",
    "missense_variant", "synonymous_variant", "stop_gained", "stop_lost",
    "start_lost", "start_retained_variant", "stop_retained_variant",
    "frameshift_variant", "inframe_insertion", "inframe_deletion",
    "protein_altering_variant", "incomplete_terminal_codon_variant",
    "5_prime_utr_variant", "3_prime_utr_variant", "transcript_ablation",
    "transcript_amplification", "exon_loss_variant", "mature_mirna_variant",
}
try:
    csv.field_size_limit(sys.maxsize)
except OverflowError:
    csv.field_size_limit(2**31 - 1)


def json_data(path):
    try:
        with path.open(encoding="utf-8") as h:
            x = json.load(h)
        return x if isinstance(x, dict) else {}
    except (OSError, ValueError):
        return {}


def tsv(path):
    with path.open(encoding="utf-8-sig", newline="") as h:
        yield from csv.DictReader(h, delimiter="\t")


def state(root, slug, study, step):
    path = root / "pipeline_audit" / slug / ANCESTRY / "status" / study / (step + ".json")
    return str(json_data(path).get("status", "NOT_RECORDED")).upper()


def number(x):
    try:
        a = float(x)
        return a if math.isfinite(a) and 0 <= a <= 1 else None
    except (ValueError, TypeError):
        return None


def truth(x):
    return str(x).strip().lower() in {"true", "t", "1", "yes", "y"}


def terms(value):
    return {v.strip().lower() for v in re.split(r"[&,;|]", str(value)) if v.strip()}


def chromosome(x):
    value = str(x).strip().removeprefix("chr").removeprefix("CHR")
    return value if value in {str(i) for i in range(1, 23)} else None


def group(lnc, pc):
    if lnc and pc:
        return "BOTH"
    if lnc:
        return "LNCRNA_ONLY"
    if pc:
        return "PROTEIN_CODING_ONLY"
    return "NEITHER"


def pct(n, total):
    return f"{n * 100 / total:.2f}%" if total else "N/A"


def show_counts(title, c, n):
    print("\n" + title)
    print("-" * 72)
    print(f"{'CATEGORY':29} {'N':>11} {'PERCENT':>12}")
    for cat in CLASSES:
        print(f"{cat:29} {c[cat]:11,} {pct(c[cat], n):>12}")
    print("-" * 72)
    print(f"{'TOTAL':29} {n:11,} {pct(n, n):>12}")


def gtf_attrs(s):
    return dict(re.findall(r'(\S+)\s+"([^"]*)"', s))


def load_gtf(gtf, include_legacy):
    """Load gene features only; far smaller than loading transcript/exon rows."""
    genes = defaultdict(list)
    n = Counter()
    opener = gzip.open if gtf.suffix == ".gz" else open
    with opener(gtf, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t", 8)
            if len(cols) != 9 or cols[2] != "gene":
                continue
            chrom = chromosome(cols[0])
            if chrom is None:
                continue
            attrs = gtf_attrs(cols[8])
            bio = attrs.get("gene_type", attrs.get("gene_biotype", "")).lower()
            kind = "LNC" if bio in LNC_BIOTYPES or (include_legacy and bio in LEGACY) else ("PC" if bio == "protein_coding" else "")
            if not kind:
                continue
            try:
                start, end = int(cols[3]), int(cols[4])
            except ValueError:
                continue
            if start > end:
                continue
            name = attrs.get("gene_name", attrs.get("gene_id", "UNKNOWN"))
            gene_id = attrs.get("gene_id", name)
            genes[chrom].append((start, end, kind, gene_id, name))
            n[kind] += 1
    result = {}
    for chrom, vals in genes.items():
        vals.sort(key=lambda x: x[0])
        result[chrom] = ([x[0] for x in vals], vals)
    return result, n


def overlap_genes(index, chrom, start, end):
    """Return sets of GENCODE gene IDs intersecting this genomic interval."""
    starts, features = index.get(chrom, ([], []))
    bound = bisect.bisect_right(starts, end)
    lnc_ids, pc_ids = set(), set()
    lnc_names, pc_names = set(), set()
    for pos in range(bound):
        gene_start, gene_end, kind, gid, name = features[pos]
        if gene_end < start:
            continue
        if kind == "LNC":
            lnc_ids.add(gid)
            lnc_names.add(name)
        elif kind == "PC":
            pc_ids.add(gid)
            pc_names.add(name)
    return lnc_ids, pc_ids, lnc_names, pc_names


def study_files(root, slug, accession):
    six = root / "06_finemapping" / slug / ANCESTRY / accession
    seven = root / "07_annotation" / slug / ANCESTRY / accession
    return six, seven, {
        "fm": six / "finemapping_summary.json",
        "loci": six / "loci_definition.tsv",
        "fits": six / (accession + "_finemapping_summary.tsv"),
        "fine": six / (accession + "_finemapped_variants.tsv.gz"),
        "all_vep": seven / (accession + "_VEP_all_consequences.tsv.gz"),
        "variant_vep": seven / (accession + "_VEP_variant_summary.tsv"),
        "anno": seven / "annotation_summary.json",
        "fail_marker": seven / "ANNOTATION_FAILED.txt",
    }


def readiness(root, slug, accession):
    s06 = state(root, slug, accession, "Step06_FineMapping")
    s07 = state(root, slug, accession, "Step07_VEP")
    six, seven, f = study_files(root, slug, accession)
    fm = json_data(f["fm"])
    if s06 not in DONE:
        return None, "STEP06_" + s06
    if str(fm.get("STATUS", "")).upper() == "NO_SIGNIFICANT_LOCI" or (fm and fm.get("N_LOCI") == 0):
        return None, "NO_SIGNIFICANT_LOCI"
    if s07 not in DONE:
        return None, "STEP07_" + s07
    anno = json_data(f["anno"])
    for key in ("loci", "fits", "fine", "all_vep", "variant_vep"):
        if not f[key].is_file() or f[key].stat().st_size == 0:
            return None, "MISSING_" + key.upper()
    if f["fail_marker"].exists() or anno.get("STATUS") != "COMPLETE":
        return None, "INVALID_ANNOTATION_STATUS"
    if str(anno.get("ANNOTATION_SCOPE", "")).upper() != "ALL_FINEMAPPED_VARIANTS":
        return None, "INCOMPLETE_ANNOTATION_SCOPE"
    try:
        current = (int(anno["FINEMAPPED_INPUT_SIZE"]) == f["fine"].stat().st_size
                   and int(anno["FINEMAPPED_INPUT_MTIME_NS"]) == f["fine"].stat().st_mtime_ns)
    except (ValueError, TypeError, KeyError):
        current = False
    if not current:
        return None, "STALE_VEP_ANNOTATION"
    return (f, fm, anno), "READY"


def summarize_study(f, fm, anno, accession, index, highpip, top):
    """Atomic study processing: caller commits results only after validation."""
    loci = list(tsv(f["loci"]))
    fits = list(tsv(f["fits"]))
    successful = {r.get("LOCUS_ID", "") for r in fits if r.get("LOCUS_ID")}
    if not successful or len(successful) != len(fits):
        raise ValueError("No unique successful locus IDs")
    if int(fm.get("N_LOCI_SUCCESS", -1)) != len(successful):
        raise ValueError("Fine-mapping success count disagrees with locus summaries")
    if int(fm.get("N_LOCI_FAILED", -1)) != 0:
        raise ValueError("Nonzero failed-locus count: main analysis requires complete studies")
    locus_meta = {r.get("LOCUS_ID"): r for r in loci}
    if len(locus_meta) != len(loci) or not successful.issubset(locus_meta):
        raise ValueError("Locus definitions do not match successful locus summaries")

    info = {}
    for lid in sorted(successful):
        d = locus_meta[lid]
        chrom = chromosome(d.get("CHR"))
        try:
            start, end = int(d["LOCUS_START"]), int(d["LOCUS_END"])
        except (ValueError, KeyError):
            raise ValueError("Invalid locus boundaries: " + lid)
        if not chrom or start > end:
            raise ValueError("Invalid chromosome or coordinates: " + lid)
        lnc_ids, pc_ids, lnc_names, pc_names = overlap_genes(index, chrom, start, end) if index is not None else (set(), set(), set(), set())
        info[lid] = {
            "chrom": chrom, "start": start, "end": end,
            "gtf_class": group(bool(lnc_ids), bool(pc_ids)) if index is not None else None,
            "gene_lnc": len(lnc_ids), "gene_pc": len(pc_ids),
            "gene_lnc_names": sorted(lnc_names), "gene_pc_names": sorted(pc_names),
            "cs_variants": 0, "cs_lnc": 0, "cs_pc": 0, "cs_both": 0,
            "high_variants": 0, "high_lnc": 0, "high_pc": 0,
            "cs_lnc_symbols": set(), "cs_lnc_genes": set(), "max_lnc_pip": None,
        }

    seen = 0
    last_id = None
    with gzip.open(f["all_vep"], "rt", encoding="utf-8", newline="") as h:
        reader = csv.DictReader(h, delimiter="\t")
        mandatory = {"VEP_ID", "BIOTYPE", "Consequence", "LOCUS_ID", "PIP", "IN_95_CREDIBLE_SET"}
        missing = mandatory - set(reader.fieldnames or [])
        if missing:
            raise ValueError("Missing VEP columns: " + ",".join(sorted(missing)))
        for vep_id, records in itertools.groupby(reader, key=lambda x: x.get("VEP_ID", "")):
            if not vep_id or (last_id is not None and vep_id <= last_id):
                raise ValueError("VEP_IDs not strictly increasing; unexpected annotation ordering")
            last_id = vep_id
            first = None
            vlnc = vpc = False
            names = set()
            genes = set()
            for rec in records:
                if first is None:
                    first = rec
                bio = str(rec.get("BIOTYPE", "")).lower().strip()
                c = terms(rec.get("Consequence", ""))
                if not (c & IN_TRANSCRIPT):
                    continue
                if bio in LNC_BIOTYPES or (top["legacy"] and bio in LEGACY):
                    vlnc = True
                    gene = str(rec.get("Gene", "")).strip()
                    symbol = str(rec.get("SYMBOL", "")).strip()
                    if gene not in {"", "-", "nan"}:
                        genes.add(gene)
                    if symbol not in {"", "-", "nan"}:
                        names.add(symbol)
                if bio == "protein_coding":
                    vpc = True
            seen += 1
            locus = str(first.get("LOCUS_ID", ""))
            if locus not in info:
                raise ValueError("VEP record references an unknown successful locus: " + locus)
            loc = info[locus]
            pip = number(first.get("PIP"))
            credible = truth(first.get("IN_95_CREDIBLE_SET"))
            if credible:
                loc["cs_variants"] += 1
                loc["cs_lnc"] += int(vlnc)
                loc["cs_pc"] += int(vpc)
                loc["cs_both"] += int(vlnc and vpc)
                loc["cs_lnc_symbols"].update(names)
                loc["cs_lnc_genes"].update(genes)
                if vlnc and pip is not None:
                    loc["max_lnc_pip"] = pip if loc["max_lnc_pip"] is None else max(pip, loc["max_lnc_pip"])
            if pip is not None and pip >= highpip:
                loc["high_variants"] += 1
                loc["high_lnc"] += int(vlnc)
                loc["high_pc"] += int(vpc)
    try:
        expected = int(anno["N_SELECTED_VARIANTS"])
    except (ValueError, TypeError, KeyError):
        expected = None
    if not seen or (expected is not None and seen != expected):
        raise ValueError(f"Incomplete transcript annotation: observed {seen}, expected {expected}")
    if any(loc["cs_variants"] > 0 and loc["cs_lnc"] > loc["cs_variants"] for loc in info.values()):
        raise ValueError("Impossible per-locus lncRNA counts")
    return info, seen


def summarize_loci(info, gtf_enabled):
    physical = Counter()
    vep = Counter()
    absent_cs = 0
    variant = Counter()
    high = Counter()
    gene_lnc = set()
    for d in info.values():
        if gtf_enabled:
            physical[d["gtf_class"]] += 1
        if d["cs_variants"] == 0:
            absent_cs += 1
        else:
            vep[group(d["cs_lnc"] > 0, d["cs_pc"] > 0)] += 1
        variant["TOTAL"] += d["cs_variants"]
        variant["LNC"] += d["cs_lnc"]
        variant["PC"] += d["cs_pc"]
        variant["BOTH"] += d["cs_both"]
        variant["LNC_ONLY"] += d["cs_lnc"] - d["cs_both"]
        variant["PC_ONLY"] += d["cs_pc"] - d["cs_both"]
        variant["NEITHER"] += d["cs_variants"] - d["cs_lnc"] - d["cs_pc"] + d["cs_both"]
        high["TOTAL"] += d["high_variants"]
        high["LNC"] += d["high_lnc"]
        high["PC"] += d["high_pc"]
        gene_lnc.update(d["cs_lnc_genes"])
    return physical, vep, absent_cs, variant, high, gene_lnc


def run(root, diseases, ancestry="european", gtf="auto", top=15, high_pip=0.50,
        include_legacy_lncrna_biotypes=False, **_):
    """Run Analysis 06 for [(display name, folder slug), ...]; return its tables.

    gtf: GENCODE GRCh38 GTF path, "auto" (resources/gencode/... under root) or None.
    """
    global DISEASES, ANCESTRY
    DISEASES, ANCESTRY = list(diseases), ancestry
    TABLES.clear()
    args = argparse.Namespace(top=max(1, int(top)), high_pip=float(high_pip),
                              include_legacy_lncrna_biotypes=bool(include_legacy_lncrna_biotypes))
    if not 0 <= args.high_pip <= 1:
        raise ValueError("high_pip must be between 0 and 1")
    root = Path(root).expanduser().resolve()

    print("=" * 105)
    print("JOHN MATTICK - ANALYSIS 06: lncRNA LOCUS REPRESENTATION")
    print("READ ONLY")
    print("=" * 105)
    default_gtf = root / "resources" / "gencode" / "release50" / "gencode.v50.primary_assembly.annotation.gtf.gz"
    if gtf in (None, "none"):
        gtf = None
    else:
        gtf = default_gtf if gtf == "auto" else Path(gtf).expanduser().resolve()
    index = None
    if gtf is not None and gtf.is_file():
        print("GENCODE GTF:", gtf)
        index, gene_counts = load_gtf(gtf, args.include_legacy_lncrna_biotypes)
        print("Reference genes: lncRNA={:,}; protein-coding={:,}".format(gene_counts["LNC"], gene_counts["PC"]))
        if not index:
            raise ValueError(f"GTF had no recognizable gene features: {gtf}")
    else:
        print("GENCODE GTF: NOT AVAILABLE")
        print("Physical locus gene-body classifications will be UNAVAILABLE, not zero.")
        print("Expected reference:", default_gtf)

    combined_physical = Counter()
    combined_vep = Counter()
    combined_variant = Counter()
    combined_high = Counter()
    combined_loci = 0
    combined_without_cs = 0
    per_disease = {}

    for disease, slug in DISEASES:
        print("\n" + "=" * 105)
        print(disease.upper())
        print("=" * 105)
        manifest = root / "01_gwas_catalog" / slug / ANCESTRY / "GWAS_download_manifest.tsv"
        if not manifest.is_file():
            manifest = root / "pipeline_audit" / slug / ANCESTRY / "study_manifest.tsv"
        if not manifest.is_file():
            print("No study manifest found: no analyses performed")
            continue
        studies = list(tsv(manifest))
        counts = Counter()
        all_info = {}
        print(f"Selected GWAS studies: {len(studies)}")
        print(f"{'STUDY':15} {'STATUS':28} {'LOCI':>6} {'CS_LOCI':>8} {'CS_LNC_LOCI':>12} {'CS_VARS':>11}")
        for rec in studies:
            study = str(rec.get("STUDY_ACCESSION", "")).strip()
            if not study:
                continue
            ready, reason = readiness(root, slug, study)
            if ready is None:
                counts[reason] += 1
                print(f"{study:15} {reason[:28]:28} {'-':>6} {'-':>8} {'-':>12} {'-':>11}")
                continue
            files, fm, anno = ready
            try:
                data, n_vep = summarize_study(files, fm, anno, study, index, args.high_pip, {"legacy": args.include_legacy_lncrna_biotypes})
            except (OSError, EOFError, UnicodeError, ValueError, csv.Error) as e:
                counts["DATA_ERROR"] += 1
                print(f"{study}: SKIPPED data-validation failure: {e}")
                continue
            counts["INCLUDED"] += 1
            all_info.update({(study, locus): value for locus, value in data.items()})
            n_cs_loci = sum(d["cs_variants"] > 0 for d in data.values())
            n_lnc_loci = sum(d["cs_lnc"] > 0 for d in data.values())
            n_cs = sum(d["cs_variants"] for d in data.values())
            print(f"{study:15} {'INCLUDED':28} {len(data):6,} {n_cs_loci:8,} {n_lnc_loci:12,} {n_cs:11,}", flush=True)

        physical, veps, missing_cs, variants, highs, genes = summarize_loci(all_info, index is not None)
        per_disease[disease] = (physical, veps, missing_cs, variants, highs, len(all_info))
        combined_physical.update(physical)
        combined_vep.update(veps)
        combined_variant.update(variants)
        combined_high.update(highs)
        combined_loci += len(all_info)
        combined_without_cs += missing_cs
        print("\nSTUDY READINESS")
        for kind, n in sorted(counts.items()):
            print(f"  {kind:35}: {n:,}")
        print("\nSuccessful study-locus units analyzed:", f"{len(all_info):,}")
        if index is not None:
            show_counts("A. ALL SUCCESSFUL LOCI: GENCODE GENE-BODY OVERLAP", physical, len(all_info))
            print("Note: A 2-Mb region may contain both gene classes even if neither is causal.")
        else:
            print("\nA. GENCODE GENE-BODY OVERLAP: NOT CALCULATED (reference missing)")
        n_loci_cs = sum(veps.values())
        show_counts("B. LOCI WITH >=1 CREDIBLE-SET VARIANT: VEP TRANSCRIPT OVERLAP", veps, n_loci_cs)
        print("Successful loci with zero credible-set variants:", f"{missing_cs:,}")
        print("Note: BOTH here means one or more credible variants hit each transcript class; not necessarily the same variant.")
        print("\nC. VARIANT-LEVEL CROSS-CHECK (95% CREDIBLE SET)")
        print("Total credible-set variant-study-locus records:", f"{variants['TOTAL']:,}")
        for name, key in (("lncRNA transcript", "LNC"), ("Protein-coding transcript", "PC"), ("Both on same variant", "BOTH"), ("lncRNA only", "LNC_ONLY"), ("Protein-coding only", "PC_ONLY"), ("Neither", "NEITHER")):
            print(f"  {name:29}: {variants[key]:11,}  ({pct(variants[key], variants['TOTAL'])})")
        print(f"High-PIP >= {args.high_pip:g} variants: {highs['TOTAL']:,}; lncRNA overlap: {highs['LNC']:,} ({pct(highs['LNC'], highs['TOTAL'])})")
        for category in CLASSES:
            TABLES.setdefault("locus_classes", []).append({
                "PHENOTYPE": disease, "CLASS": category,
                "LOCI_BY_GENCODE_GENE_BODY": physical[category] if index is not None else None,
                "PERCENT_OF_ALL_LOCI": (round(100 * physical[category] / len(all_info), 3)
                                        if index is not None and all_info else None),
                "LOCI_BY_CS_VARIANT_TRANSCRIPT": veps[category],
                "PERCENT_OF_LOCI_WITH_CS": round(100 * veps[category] / sum(veps.values()), 3) if veps else None,
            })
        TABLES.setdefault("variant_crosscheck", []).append({
            "PHENOTYPE": disease, "STUDY_LOCI": len(all_info), "LOCI_WITHOUT_CS": missing_cs,
            **{f"CS_{k}": variants[k] for k in ("TOTAL", "LNC", "PC", "BOTH", "LNC_ONLY", "PC_ONLY", "NEITHER")},
            "HIGH_PIP_TOTAL": highs["TOTAL"], "HIGH_PIP_LNC": highs["LNC"], "DISTINCT_LNCRNA_GENES": len(genes),
        })
        for (study, locus), d in sorted(all_info.items()):
            TABLES.setdefault("loci", []).append({
                "PHENOTYPE": disease, "STUDY_ACCESSION": study, "LOCUS_ID": locus,
                "CHR": d["chrom"], "START": d["start"], "END": d["end"], "GENCODE_CLASS": d["gtf_class"],
                "GENCODE_LNCRNA_GENES": d["gene_lnc"], "GENCODE_PC_GENES": d["gene_pc"],
                "CS_VARIANTS": d["cs_variants"], "CS_LNCRNA": d["cs_lnc"], "CS_PC": d["cs_pc"],
                "CS_BOTH": d["cs_both"], "MAX_LNCRNA_PIP": d["max_lnc_pip"],
                "CS_LNCRNA_SYMBOLS": ";".join(sorted(d["cs_lnc_symbols"])),
            })
        print("Distinct lncRNA gene IDs with credible variant overlap:", f"{len(genes):,}")
        print("\nTOP LOCI BY NUMBER OF lncRNA-OVERLAPPING CREDIBLE VARIANTS")
        ranked = sorted(all_info.items(), key=lambda z: (z[1]["cs_lnc"], z[1]["max_lnc_pip"] or -1), reverse=True)
        for (study, locus), d in ranked[:args.top]:
            lnc_txt = ",".join(d["gene_lnc_names"][:4]) if index is not None else "GTF_UNAVAILABLE"
            print(f"  {study} {locus:<6} chr{d['chrom']}:{d['start']}-{d['end']}  CS={d['cs_variants']:,}  CS_lnc={d['cs_lnc']:,}  CS_pc={d['cs_pc']:,}  GTF_lnc_genes={d['gene_lnc']}  examples={lnc_txt}")

    print("\n" + "=" * 105)
    print("CROSS-PHENOTYPE COMPARISON: STUDY-LOCUS DISTRIBUTIONS")
    print("=" * 105)
    for mode in ("GENCODE_LOCUS_INTERVAL", "VEP_CREDIBLE_VARIANT_TRANSCRIPT"):
        if mode == "GENCODE_LOCUS_INTERVAL" and index is None:
            print("GENCODE interval comparison unavailable (no GTF)")
            continue
        print("\n" + mode)
        print(f"{'CATEGORY':27}" + "".join(f"{name[:22].upper():>23}" for name, _ in DISEASES))
        for category in CLASSES:
            values = []
            for disease, _ in DISEASES:
                entry = per_disease.get(disease)
                counts = (entry[0] if mode == "GENCODE_LOCUS_INTERVAL" else entry[1]) if entry else Counter()
                den = sum(counts.values())
                values.append(f"{counts[category]:,} ({pct(counts[category], den)})")
            print(f"{category:27}" + "".join(f"{v:>23}" for v in values))

    print("\n" + "=" * 105)
    print("OVERALL DESCRIPTIVE SUMMARY (STUDY-LOCUS UNITS, NOT INDEPENDENT LOCI)")
    print("=" * 105)
    if index is not None:
        show_counts("A. LOCUS GENE-BODY OVERLAP", combined_physical, combined_loci)
    else:
        print("GENCODE physical locus classification: unavailable")
    show_counts("B. LOCUS CREDIBLE-VARIANT TRANSCRIPT OVERLAP", combined_vep, sum(combined_vep.values()))
    print("Loci with zero credible variants:", f"{combined_without_cs:,}")
    print("Total credible-set variant-study-locus records:", f"{combined_variant['TOTAL']:,}")
    print("Credible variants overlapping lncRNA transcripts:", f"{combined_variant['LNC']:,}", f"({pct(combined_variant['LNC'], combined_variant['TOTAL'])})")
    print("\nSCIENTIFIC CAUTIONS")
    print("1. These are study-locus counts; different GWAS may contain the same physical association.")
    print("2. GENCODE gene-body overlap does not imply expression, effect, or causality.")
    print("3. Near-gene VEP consequences are excluded from transcript overlaps.")
    print("4. Zero credible-set loci have no assessed CS variant; they are not counted as NEITHER.")
    print("5. Excluded failed/partial studies must be reported in the manuscript flow.")
    return {name: pd.DataFrame(rows) for name, rows in TABLES.items()}


def main(argv=None):
    return analysis_main(run, "John Analysis06: locus-level lncRNA representation", argv, extra=[
        ("--gtf", dict(default="auto", help="GENCODE GRCh38 GTF path, auto, or none")),
        ("--top", dict(type=int, default=15)),
        ("--high-pip", dict(type=float, default=0.50)),
        ("--include-legacy-lncrna-biotypes", dict(action="store_true")),
    ])


if __name__ == "__main__":
    main()
