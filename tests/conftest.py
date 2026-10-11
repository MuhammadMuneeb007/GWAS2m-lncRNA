"""Synthetic Step 1/6/7 outputs (original GWAS2m layout) so the analyses run without PLINK/R/VEP."""

from __future__ import annotations

import gzip
import json

import pandas as pd
import pytest
import yaml

from gwas2m_lncrna.config import load_config
from gwas2m_lncrna.project import Project


@pytest.fixture
def make_config(tmp_path):
    def _make(**sections):
        data = {"output_dir": "results", "resources_dir": "resources", **sections}
        path = tmp_path / "gwas2m-lncrna.yaml"
        path.write_text(yaml.safe_dump(data), encoding="utf-8")
        return load_config(str(path))
    return _make


def vep_row(vep_id, biotype, consequence, exon="-", intron="-", gene="ENSG1", symbol="GENE1", feature="ENST1"):
    return {"VEP_ID": vep_id, "BIOTYPE": biotype, "Consequence": consequence, "EXON": exon, "INTRON": intron,
            "Gene": gene, "SYMBOL": symbol, "Feature": feature}


def build_study(project: Project, accession: str, variants: list[dict]) -> None:
    """variants: dicts with LOCUS_ID, POS, PIP, CS (bool), rows (VEP transcript rows)."""
    six = project.finemap_dir / accession
    seven = project.annotation_dir / accession
    six.mkdir(parents=True, exist_ok=True)
    seven.mkdir(parents=True, exist_ok=True)

    loci = sorted({v["LOCUS_ID"] for v in variants})
    pd.DataFrame([{"LOCUS_ID": l, "CHR": 1, "LOCUS_START": 1, "LOCUS_END": 2_001_000,
                   "N_INDEPENDENT_SIGNALS": 1, "LEAD_IDS": "1:1000:A:G", "LEAD_POSITIONS": "1000",
                   "MIN_LEAD_P": 1e-10} for l in loci]).to_csv(six / "loci_definition.tsv", sep="\t", index=False)
    pd.DataFrame([{"STUDY_ACCESSION": accession, "LOCUS_ID": l, "N_VARIANTS": sum(v["LOCUS_ID"] == l for v in variants),
                   "N_CREDIBLE_SETS": 1, "MAX_PIP": max(v["PIP"] for v in variants if v["LOCUS_ID"] == l),
                   "CONVERGED": True} for l in loci]).to_csv(
        six / f"{accession}_finemapping_summary.tsv", sep="\t", index=False)

    fine_rows, cs_rows, all_rows, summary_rows = [], [], [], []
    for i, v in enumerate(variants, start=1):
        record = {"STUDY_ACCESSION": accession, "LOCUS_ID": v["LOCUS_ID"], "REFERENCE_ID": f"1:{v['POS']}:A:G",
                  "CHR": 1, "POS": v["POS"], "REFERENCE_POS": v["POS"], "REFERENCE_REF": "A", "REFERENCE_ALT": "G",
                  "PIP": v["PIP"], "IN_95_CREDIBLE_SET": v["CS"], "CREDIBLE_SET": "L1" if v["CS"] else ""}
        fine_rows.append(record)
        if v["CS"]:
            cs_rows.append({k: record[k] for k in ("STUDY_ACCESSION", "LOCUS_ID", "CREDIBLE_SET", "REFERENCE_ID",
                                                   "CHR", "POS", "PIP")})
        vep_id = f"GWASVEP_{i:08d}"
        rows = v["rows"] or [{}]
        for r in rows:
            all_rows.append({**record, "VEP_ID": vep_id, **{k: r.get(k, "") for k in
                             ("BIOTYPE", "Consequence", "EXON", "INTRON", "Gene", "SYMBOL", "Feature")}})
        summary_rows.append({**record, "VEP_ID": vep_id,
                             "ALL_CONSEQUENCES": ";".join(r.get("Consequence", "") for r in rows if r),
                             "ALL_BIOTYPES": ";".join(r.get("BIOTYPE", "") for r in rows if r)})
    fine = six / f"{accession}_finemapped_variants.tsv.gz"
    pd.DataFrame(fine_rows).to_csv(fine, sep="\t", index=False, compression="gzip")
    pd.DataFrame(cs_rows, columns=["STUDY_ACCESSION", "LOCUS_ID", "CREDIBLE_SET", "REFERENCE_ID", "CHR", "POS",
                                   "PIP"]).to_csv(six / f"{accession}_95pct_credible_sets.tsv", sep="\t", index=False)
    (six / "finemapping_summary.json").write_text(json.dumps({
        "STUDY_ACCESSION": accession, "STATUS": "COMPLETE", "N_LOCI": len(loci), "N_LOCI_SUCCESS": len(loci),
        "N_LOCI_FAILED": 0}))
    with gzip.open(seven / f"{accession}_VEP_all_consequences.tsv.gz", "wt", encoding="utf-8") as h:
        pd.DataFrame(all_rows).to_csv(h, sep="\t", index=False)
    pd.DataFrame(summary_rows).to_csv(seven / f"{accession}_VEP_variant_summary.tsv", sep="\t", index=False)
    (seven / "annotation_summary.json").write_text(json.dumps({
        "STATUS": "COMPLETE", "ANNOTATION_SCOPE": "ALL_FINEMAPPED_VARIANTS",
        "N_SELECTED_VARIANTS": len(variants), "N_SELECTED_CREDIBLE_SET": sum(bool(v["CS"]) for v in variants),
        "FINEMAPPED_INPUT_SIZE": fine.stat().st_size, "FINEMAPPED_INPUT_MTIME_NS": fine.stat().st_mtime_ns}))
    for step in ("download", "qc", "ld_reference", "clump", "finemap", "annotate"):
        project.set_status(accession, step, "COMPLETE")


def build_manifest(project: Project, studies: list[tuple[str, str]]) -> None:
    project.discovery_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([{"ARRAY_ID": i, "STUDY_ACCESSION": acc, "DOWNLOAD_URL": "", "FILE_NAME": f"{acc}.h.tsv.gz",
                   "PHENOTYPE": project.phenotype, "ANCESTRY": project.ancestry_label,
                   "ANCESTRY_CODE": project.ancestry_code, "DISEASE_TRAIT": trait,
                   "INITIAL_SAMPLE_SIZE": "1,000 cases, 2,000 controls", "RANKING_N": 3000}
                  for i, (acc, trait) in enumerate(studies, start=1)]).to_csv(
        project.manifest_path, sep="\t", index=False)
