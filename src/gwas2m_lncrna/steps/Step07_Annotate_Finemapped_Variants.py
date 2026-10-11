"""Step 7 - annotate every fine-mapped variant of one study with Ensembl VEP.

All Step 6 variants are annotated (not only credible-set members); credible-set
membership and PIP thresholds are kept as flags. Outputs in 07_annotation/<GCST>/
keep the original file names, including <GCST>_VEP_all_consequences.tsv.gz
(one row per transcript consequence), which the analysis reads.
"""

from __future__ import annotations

import io
from pathlib import Path

import numpy as np
import pandas as pd

from gwas2m_lncrna.common import (
    Gwas2RnaError, as_bool, banner, nonempty, read_json, run_command, utc_now, write_json,
)
from gwas2m_lncrna.project import Project
from gwas2m_lncrna.resources import Resources
from gwas2m_lncrna.tools import find_tool, vep_release

ASSEMBLY = "GRCh38"
STEP07_VERSION = "2.0.0"
ANNOTATION_SCOPE = "ALL_FINEMAPPED_VARIANTS"

VEP_FIELDS = [
    "Uploaded_variation", "Location", "Allele", "Gene", "Feature", "Feature_type",
    "Consequence", "IMPACT", "SYMBOL", "SYMBOL_SOURCE", "HGNC_ID", "BIOTYPE",
    "EXON", "INTRON", "HGVSc", "HGVSp", "cDNA_position", "CDS_position",
    "Protein_position", "Amino_acids", "Codons", "Existing_variation", "DISTANCE",
    "STRAND", "FLAGS", "VARIANT_CLASS", "CANONICAL", "MANE_SELECT",
    "MANE_PLUS_CLINICAL", "TSL", "APPRIS", "DOMAINS", "SIFT", "PolyPhen", "PICK",
]

CATEGORY_RANK = {"SPLICING": 0, "CODING": 1, "UTR": 2, "REGULATORY": 3, "NONCODING_RNA": 4,
                 "INTRONIC": 5, "UPSTREAM": 6, "DOWNSTREAM": 7, "INTERGENIC": 8, "OTHER": 9}
IMPACT_RANK = {"HIGH": 0, "MODERATE": 1, "LOW": 2, "MODIFIER": 3}


def load_finemapped_variants(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", compression="infer", low_memory=False)
    required = ["LOCUS_ID", "CHR", "REFERENCE_POS", "REFERENCE_REF", "REFERENCE_ALT", "PIP"]
    missing = [x for x in required if x not in df.columns]
    if missing:
        raise RuntimeError(f"Fine-mapped file missing columns {missing}: {path}")
    df["CHR"] = pd.to_numeric(df["CHR"].astype(str).str.replace(r"^chr", "", regex=True, case=False),
                              errors="coerce")
    df["REFERENCE_POS"] = pd.to_numeric(df["REFERENCE_POS"], errors="coerce")
    df["PIP"] = pd.to_numeric(df["PIP"], errors="coerce")
    df = df.dropna(subset=required)
    df = df[df["CHR"].between(1, 22)].copy()
    df["CHR"] = df["CHR"].astype(int)
    df["REFERENCE_POS"] = df["REFERENCE_POS"].astype(int)
    for column in ("REFERENCE_REF", "REFERENCE_ALT"):
        df[column] = df[column].astype(str).str.upper().str.strip()
    df = df[(df["REFERENCE_REF"] != "") & (df["REFERENCE_ALT"] != "")].copy()
    if df.empty:
        raise RuntimeError(f"No usable fine-mapped variants in {path}")
    return df


def select_variants(df: pd.DataFrame, min_pip: float) -> pd.DataFrame:
    if "IN_95_CREDIBLE_SET" in df.columns:
        credible = df["IN_95_CREDIBLE_SET"].map(as_bool)
    elif "CREDIBLE_SET" in df.columns:
        credible = ~df["CREDIBLE_SET"].fillna("").astype(str).str.strip().str.lower().isin({"", "na", "nan", "none", "."})
    else:
        credible = pd.Series(False, index=df.index)
    top = pd.Series(False, index=df.index)
    top.loc[df.groupby("LOCUS_ID")["PIP"].idxmax()] = True
    selected = df.copy()
    selected["SELECTED_CREDIBLE_SET"] = credible.to_numpy()
    selected["SELECTED_PIP_THRESHOLD"] = (df["PIP"] >= min_pip).to_numpy()
    selected["SELECTED_TOP_LOCUS"] = top.to_numpy()
    for threshold, name in ((0.01, "PIP_GE_0_01"), (0.10, "PIP_GE_0_10"), (0.50, "PIP_GE_0_50"), (0.90, "PIP_GE_0_90")):
        selected[name] = selected["PIP"] >= threshold
    selected["PRIORITIZED"] = (selected["SELECTED_CREDIBLE_SET"] | selected["SELECTED_PIP_THRESHOLD"]
                               | selected["SELECTED_TOP_LOCUS"])
    selected = (selected.sort_values(["CHR", "REFERENCE_POS", "PIP"], ascending=[True, True, False], kind="stable")
                .drop_duplicates(["CHR", "REFERENCE_POS", "REFERENCE_REF", "REFERENCE_ALT"], keep="first")
                .reset_index(drop=True))
    selected["VEP_ID"] = [f"GWASVEP_{i:08d}" for i in range(1, len(selected) + 1)]
    return selected


def write_vep_input(selected: pd.DataFrame, vcf_file: Path, mapping_file: Path) -> None:
    columns = [c for c in [
        "VEP_ID", "STUDY_ACCESSION", "PHENOTYPE", "ANCESTRY_CODE", "LOCUS_ID", "REFERENCE_ID", "CHR",
        "REFERENCE_POS", "REFERENCE_REF", "REFERENCE_ALT", "PIP", "IN_95_CREDIBLE_SET", "CREDIBLE_SET",
        "SELECTED_CREDIBLE_SET", "SELECTED_PIP_THRESHOLD", "SELECTED_TOP_LOCUS", "PIP_GE_0_01",
        "PIP_GE_0_10", "PIP_GE_0_50", "PIP_GE_0_90", "PRIORITIZED", "rsID", "SNPID", "P", "BETA", "SE", "Z",
        "Z_REF"] if c in selected.columns]
    selected[columns].to_csv(mapping_file, sep="\t", index=False)
    with open(vcf_file, "w", encoding="utf-8") as handle:
        handle.write("##fileformat=VCFv4.2\n##reference=GRCh38\n##source=gwas2m-lncrna_step7\n")
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for row in selected.itertuples(index=False):
            handle.write(f"{int(row.CHR)}\t{int(row.REFERENCE_POS)}\t{row.VEP_ID}\t"
                         f"{row.REFERENCE_REF}\t{row.REFERENCE_ALT}\t.\t.\t.\n")


def run_vep(vep: Path, cache: Path, fasta: Path, input_vcf: Path, output_tsv: Path, log_file: Path,
            threads: int, release: str) -> None:
    run_command([
        vep, "--input_file", input_vcf, "--output_file", output_tsv, "--format", "vcf", "--tab",
        "--fields", ",".join(VEP_FIELDS), "--species", "homo_sapiens", "--assembly", ASSEMBLY,
        "--cache", "--offline", "--dir_cache", cache, "--cache_version", release, "--fasta", fasta,
        "--fork", threads, "--force_overwrite", "--symbol", "--numbers", "--hgvs", "--canonical",
        "--mane", "--appris", "--tsl", "--biotype", "--variant_class", "--domains", "--regulatory",
        "--sift", "b", "--polyphen", "b", "--flag_pick", "--no_stats",
    ], log_file)
    if not nonempty(output_tsv):
        raise RuntimeError("VEP did not generate output")


def read_vep_output(path: Path) -> pd.DataFrame:
    header, lines = None, []
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            if line.startswith("#Uploaded_variation"):
                header = line.rstrip("\n").lstrip("#").split("\t")
            elif not line.startswith("#") and line.strip():
                lines.append(line)
    if header is None:
        raise RuntimeError("Could not find VEP #Uploaded_variation header")
    if not lines:
        raise RuntimeError("VEP output contains no annotation rows")
    df = pd.read_csv(io.StringIO("".join(lines)), sep="\t", names=header, dtype=str, low_memory=False)
    return df.rename(columns={"Uploaded_variation": "VEP_ID"})


def functional_category(row) -> str:
    consequence = str(row.get("Consequence", "")).lower()
    biotype = str(row.get("BIOTYPE", "")).lower()
    if any(x in consequence for x in ("splice_acceptor_variant", "splice_donor_variant", "splice_region_variant",
                                      "splice_donor_5th_base_variant", "splice_donor_region_variant",
                                      "splice_polypyrimidine_tract_variant")):
        return "SPLICING"
    if any(x in consequence for x in ("transcript_ablation", "stop_gained", "frameshift_variant", "stop_lost",
                                      "start_lost", "transcript_amplification", "inframe_insertion",
                                      "inframe_deletion", "missense_variant", "protein_altering_variant",
                                      "coding_sequence_variant", "synonymous_variant", "start_retained_variant",
                                      "stop_retained_variant", "incomplete_terminal_codon_variant")):
        return "CODING"
    if "5_prime_utr_variant" in consequence or "3_prime_utr_variant" in consequence:
        return "UTR"
    if "regulatory_region" in consequence or "tf_binding_site" in consequence:
        return "REGULATORY"
    if "non_coding_transcript" in consequence or any(x in biotype for x in (
            "lncrna", "lincrna", "mirna", "snrna", "snorna", "rrna", "ncrna", "antisense", "processed_transcript")):
        # NOTE: biotype-based, so it also catches upstream/downstream of an lncRNA.
        # The analysis uses the stricter same-row transcript-overlap definition instead.
        return "NONCODING_RNA"
    for term, label in (("intron_variant", "INTRONIC"), ("upstream_gene_variant", "UPSTREAM"),
                        ("downstream_gene_variant", "DOWNSTREAM"), ("intergenic_variant", "INTERGENIC")):
        if term in consequence:
            return label
    return "OTHER"


def select_best_consequence(merged: pd.DataFrame) -> pd.DataFrame:
    x = merged.copy()
    blank = pd.Series("", index=x.index)
    x["_PICK"] = np.where(x.get("PICK", blank).fillna("").astype(str).str.upper().isin({"1", "YES", "Y"}), 0, 1)
    x["_MANE"] = np.where(x.get("MANE_SELECT", blank).fillna("").astype(str).str.strip().isin({"", "-"}), 1, 0)
    x["_CANON"] = np.where(x.get("CANONICAL", blank).fillna("").astype(str).str.upper().isin({"YES", "1"}), 0, 1)
    x["_IMPACT"] = x.get("IMPACT", blank).fillna("").astype(str).str.upper().map(IMPACT_RANK).fillna(9)
    x["_CATEGORY"] = x["FUNCTIONAL_CATEGORY"].map(CATEGORY_RANK).fillna(99)
    order = ["VEP_ID", "_PICK", "_MANE", "_CANON", "_IMPACT", "_CATEGORY"]
    best = x.sort_values(order, kind="stable").drop_duplicates("VEP_ID", keep="first")
    return best.drop(columns=order[1:])


def _join_unique(series: pd.Series) -> str:
    seen = []
    for value in series.dropna().astype(str):
        for token in value.split(","):
            token = token.strip()
            if token and token not in seen:
                seen.append(token)
    return ";".join(seen)


def annotate_study(project: Project, study: dict, force: bool) -> tuple:
    accession = study["STUDY_ACCESSION"]
    settings = project.cfg["annotation"]
    fine_dir = project.finemap_dir / accession
    fine_file = fine_dir / f"{accession}_finemapped_variants.tsv.gz"
    fine_summary = read_json(fine_dir / "finemapping_summary.json")
    if str(fine_summary.get("STATUS", "")).upper() == "NO_SIGNIFICANT_LOCI" or fine_summary.get("N_LOCI") == 0:
        return "NOT_APPLICABLE", "No significant loci to annotate", {}
    if int(fine_summary.get("N_LOCI_SUCCESS", 0) or 0) == 0:
        return "EXCLUDED", "No successfully fine-mapped loci", {}
    if int(fine_summary.get("N_LOCI_FAILED", 0) or 0) > 0 and not settings.get("allow_partial"):
        return "EXCLUDED", ("Partial fine-mapping "
                            f"({fine_summary.get('N_LOCI_FAILED')} loci failed); set annotation.allow_partial"), {}
    if not nonempty(fine_file):
        return "FAILED", "Fine-mapped variant file missing/empty", {}

    out = project.annotation_dir / accession
    summary_json = out / "annotation_summary.json"
    variant_output = out / f"{accession}_VEP_variant_summary.tsv"
    previous = read_json(summary_json)
    current = (int(previous.get("FINEMAPPED_INPUT_SIZE", -1) or -1) == fine_file.stat().st_size
               and int(previous.get("FINEMAPPED_INPUT_MTIME_NS", -1) or -1) == fine_file.stat().st_mtime_ns)
    if (not force and nonempty(variant_output) and previous.get("STATUS") == "COMPLETE"
            and previous.get("ANNOTATION_SCOPE") == ANNOTATION_SCOPE and current):
        return "SKIPPED_ALREADY_COMPLETE", "", {"N_ANNOTATED_VARIANTS": previous.get("N_ANNOTATED_VARIANTS")}

    resources = Resources(project.cfg)
    vep = find_tool(project.cfg, "vep")
    installed = vep_release(vep)
    if installed and installed != resources.vep_release:
        raise Gwas2RnaError(f"VEP {installed} is installed but the cache is release {resources.vep_release}; "
                            "set annotation.vep_cache_release to match (and re-run setup) or install the matching VEP.")
    if not resources.vep_ready() or not resources.fasta_ready():
        resources.ensure(None, vep=True)
    threads = project.cfg.threads
    min_pip = float(settings["min_pip"])
    out.mkdir(parents=True, exist_ok=True)
    paths = {name: out / f"{accession}_{suffix}" for name, suffix in {
        "selected": "variants_selected_for_annotation.tsv", "mapping": "VEP_ID_mapping.tsv",
        "vcf": "variants_for_VEP.vcf", "raw": "VEP_raw.tsv", "log": "VEP.log",
        "all": "VEP_all_consequences.tsv.gz", "best": "VEP_best_consequence.tsv",
        "prioritized": "VEP_prioritized_variants.tsv", "category": "VEP_category_summary.tsv",
        "locus": "VEP_locus_summary.tsv"}.items()}

    banner(f"STEP 7 - VEP ANNOTATION | {accession}")
    fine = load_finemapped_variants(fine_file)
    fine["STUDY_ACCESSION"], fine["PHENOTYPE"], fine["ANCESTRY_CODE"] = (
        accession, project.phenotype, project.ancestry_code)
    selected = select_variants(fine, min_pip)
    selected.to_csv(paths["selected"], sep="\t", index=False)
    prioritized = selected.loc[selected["PRIORITIZED"]]
    prioritized.to_csv(paths["prioritized"], sep="\t", index=False)
    write_vep_input(selected, paths["vcf"], paths["mapping"])
    print(f"Unique variants sent to VEP: {len(selected):,} (priority subset {len(prioritized):,})")

    run_vep(vep, resources.vep_cache, resources.fasta, paths["vcf"], paths["raw"], paths["log"], threads,
            resources.vep_release)
    vep_df = read_vep_output(paths["raw"])
    merged = selected.merge(vep_df, on="VEP_ID", how="left", validate="one_to_many")
    merged["FUNCTIONAL_CATEGORY"] = merged.apply(functional_category, axis=1)
    merged.to_csv(paths["all"], sep="\t", index=False, compression="gzip")
    best = select_best_consequence(merged)
    best.to_csv(paths["best"], sep="\t", index=False)
    aggregations = {"N_TRANSCRIPT_CONSEQUENCES": ("VEP_ID", "size")}
    for name, source in (("ALL_CONSEQUENCES", "Consequence"), ("ALL_GENES", "Gene"), ("ALL_SYMBOLS", "SYMBOL"),
                         ("ALL_BIOTYPES", "BIOTYPE"), ("ALL_FUNCTIONAL_CATEGORIES", "FUNCTIONAL_CATEGORY")):
        if source in merged.columns:
            aggregations[name] = (source, _join_unique)
    grouped = merged.groupby("VEP_ID", sort=False, dropna=False).agg(**aggregations).reset_index()
    variant_summary = best.merge(grouped, on="VEP_ID", how="left", validate="one_to_one")
    variant_summary.to_csv(variant_output, sep="\t", index=False)
    (variant_summary.groupby("FUNCTIONAL_CATEGORY", dropna=False)
     .agg(N_VARIANTS=("VEP_ID", "nunique"), N_LOCI=("LOCUS_ID", "nunique"), MAX_PIP=("PIP", "max"))
     .reset_index().sort_values(["N_VARIANTS", "MAX_PIP"], ascending=[False, False])
     .to_csv(paths["category"], sep="\t", index=False))
    (variant_summary.groupby("LOCUS_ID", dropna=False)
     .agg(N_ANNOTATED_VARIANTS=("VEP_ID", "nunique"), MAX_PIP=("PIP", "max"),
          N_CREDIBLE_SELECTED=("SELECTED_CREDIBLE_SET", "sum"))
     .reset_index().to_csv(paths["locus"], sep="\t", index=False))

    n_annotated = int(vep_df["VEP_ID"].nunique())
    result = {
        "STUDY_ACCESSION": accession, "PHENOTYPE": project.phenotype, "ANCESTRY_CODE": project.ancestry_code,
        "BUILD": ASSEMBLY, "ANNOTATION_SCOPE": ANNOTATION_SCOPE, "MIN_PIP": min_pip,
        "N_FINEMAPPED_ROWS_AVAILABLE": int(len(fine)), "N_FINEMAPPED_LOCI": int(fine["LOCUS_ID"].nunique()),
        "N_SELECTED_VARIANTS": int(len(selected)),
        "N_SELECTED_CREDIBLE_SET": int(selected["SELECTED_CREDIBLE_SET"].sum()),
        "N_PRIORITIZED": int(selected["PRIORITIZED"].sum()),
        "N_VEP_CONSEQUENCE_ROWS": int(len(vep_df)), "N_ANNOTATED_VARIANTS": n_annotated,
        "N_UNANNOTATED_VARIANTS": int(len(selected) - n_annotated), "STEP07_VERSION": STEP07_VERSION,
        "FINEMAPPED_INPUT_SIZE": int(fine_file.stat().st_size),
        "FINEMAPPED_INPUT_MTIME_NS": int(fine_file.stat().st_mtime_ns),
        "VEP_EXECUTABLE": str(vep), "VEP_CACHE": str(resources.vep_cache), "VEP_CACHE_RELEASE": resources.vep_release,
        "VEP_FASTA": str(resources.fasta), "VEP_ALL_CONSEQUENCES": str(paths["all"]),
        "VEP_VARIANT_SUMMARY": str(variant_output), "STATUS": "COMPLETE", "COMPLETED_UTC": utc_now(),
    }
    write_json(summary_json, result)
    print(f"Annotated variants: {n_annotated:,} of {len(selected):,}")
    return "COMPLETE", "", {"N_ANNOTATED_VARIANTS": n_annotated}
