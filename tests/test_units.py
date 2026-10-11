from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from gwas2m_lncrna.analysis import (
    Analysis01_John_GWAS_Selection as A01, Analysis03_John_Coding_Noncoding_Distribution as A03,
    Analysis07_John_lncRNA_Enrichment as A07,
)
from gwas2m_lncrna.common import Gwas2RnaError, canonical_ancestry
from gwas2m_lncrna.config import load_config, write_template
from gwas2m_lncrna.project import Project
from gwas2m_lncrna.steps.Step01_Plan_GWAS import build_terms, estimate_sample_size, gwas_bucket
from gwas2m_lncrna.steps.Step05_Clump_GWAS_Loci import map_gwas_to_reference
from gwas2m_lncrna.steps.Step06_FineMap_GWAS_Loci import _balanced_ld_subset, effective_n, load_leads_and_build_loci


# ------------------------------------------------------------------ config
def test_defaults_and_overrides(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    cfg = load_config(None, {"threads": 8})
    assert cfg.threads == 8 and cfg["finemapping"]["coverage"] == 0.95
    assert cfg.output_dir == (tmp_path / "results").resolve()


def test_init_and_unknown_key(tmp_path):
    path = write_template(tmp_path / "gwas2m-lncrna.yaml")
    assert load_config(str(path))["clumping"]["r2"] == 0.1
    with pytest.raises(Gwas2RnaError):
        write_template(path)
    path.write_text("clumping:\n  r3: 0.2\n")
    with pytest.raises(Gwas2RnaError, match="clumping.r3"):
        load_config(str(path))


def test_ancestry():
    assert canonical_ancestry("european") == ("EUR", "European")
    with pytest.raises(Gwas2RnaError):
        canonical_ancestry("MID")


def test_original_gwas2m_layout(make_config):
    project = Project(make_config(), "Parkinson's disease", "EUR")
    root = project.cfg.output_dir
    assert project.manifest_path == root / "01_gwas_catalog/parkinson_s_disease/european/GWAS_download_manifest.tsv"
    assert project.finemap_dir == root / "06_finemapping/parkinson_s_disease/european"
    assert project.qc_dir == root / "02_summary_stats/parkinson_s_disease/european/qc"
    project.set_status("GCST1", "finemap", "PARTIAL", "1 of 3 loci failed", N_LOCI=3)
    status_file = root / "pipeline_audit/parkinson_s_disease/european/status/GCST1/Step06_FineMapping.json"
    assert status_file.is_file()
    assert project.status("GCST1")["finemap"]["N_LOCI"] == 3 and project.step_status("GCST1", "finemap") == "PARTIAL"


# --------------------------------------------------------------- discovery
def test_discovery_helpers():
    assert build_terms("CAD", ["Coronary artery disease"]) == ["CAD", "coronary artery disease", "coronary heart disease"]
    assert estimate_sample_size("10,000 European ancestry cases, 20,000 controls") == 30000
    assert gwas_bucket("GCST90012345") == "GCST90012001-GCST90013000"


# ---------------------------------------------------------------- clumping
def _pvar(tmp_path: Path, rows) -> Path:
    path = tmp_path / "ref.pvar"
    path.write_text("##x\n#CHROM\tPOS\tID\tREF\tALT\n" + "".join(f"1\t{p}\t{i}\t{r}\t{a}\n" for p, i, r, a in rows))
    return path


def test_clump_mapping_no_crash_when_nothing_maps(tmp_path):
    gwas = pd.DataFrame({"CHR": [1], "POS": [100], "P": [1e-9], "EA": ["A"], "NEA": ["G"], "SOURCE_ID": ["x"]})
    mapped, stats = map_gwas_to_reference(gwas, _pvar(tmp_path, [(500, "v1", "A", "G")]))
    assert mapped.empty and list(mapped[["ID", "P"]].columns) == ["ID", "P"]  # used to raise KeyError
    assert stats["N_UNMAPPED"] == 1


def test_clump_mapping_rejects_allele_mismatch(tmp_path):
    gwas = pd.DataFrame({"CHR": [1, 1], "POS": [500, 600], "P": [1e-9, 1e-8], "EA": ["CTTT", "G"],
                         "NEA": ["C", "A"], "SOURCE_ID": ["indel", "snv"]})
    mapped, stats = map_gwas_to_reference(gwas, _pvar(tmp_path, [(500, "v1", "A", "G"), (600, "v2", "A", "G")]))
    assert mapped["ID"].tolist() == ["v2"]  # the indel is no longer mapped to the SNV v1
    assert stats["N_ALLELE_MISMATCH"] == 1


# ------------------------------------------------------------ fine-mapping
def test_loci_merge_and_effective_n(tmp_path):
    leads = tmp_path / "lead.tsv"
    pd.DataFrame({"CHROM": [1, 1, 2], "POS": [1_000_000, 2_500_000, 5_000_000], "ID": ["a", "b", "c"],
                  "P": [1e-9, 1e-10, 1e-12]}).to_csv(leads, sep="\t", index=False)
    loci = load_leads_and_build_loci(leads, 1000)
    assert len(loci) == 2 and loci.iloc[0]["N_INDEPENDENT_SIGNALS"] == 2
    assert effective_n(10_000, 10_000) == 20_000 and effective_n(None, 5) is None


def test_balanced_ld_subset_respects_cap():
    matched = pd.DataFrame({"POS": np.arange(1000), "P": np.linspace(1e-9, 1, 1000),
                            "REFERENCE_ID": [f"v{i}" for i in range(1000)]})
    subset = _balanced_ld_subset(matched, 100, [10, 900])
    assert len(subset) == 100 and set(subset["POS"]) >= {10, 900}


# ---------------------------------------------------------------- analyses
def test_analysis03_self_test():
    assert A03.self_test() == 0


def test_direct_label_works_for_any_phenotype():
    assert A01.classify_trait("Parkinson's disease", "Parkinson's disease") == "DIRECT_DISEASE"
    assert A01.classify_trait("Coronary artery disease", "CAD") == "DIRECT_DISEASE"
    assert A01.classify_trait("Migraine", "migraine") == "DIRECT_DISEASE"
    assert A01.classify_trait("Parkinson's disease (family history)", "Parkinson's disease") == "PROXY_FAMILY"


def test_analysis07_null_and_bh():
    strata = {("L1",): {"LNC_TRANSCRIPT": [10, 2, 5, 2]}, ("L2",): {"LNC_TRANSCRIPT": [20, 4, 2, 0]}}
    a = A07.calculate_null(strata, "LNC_TRANSCRIPT", 2000, 7)
    assert a == A07.calculate_null(strata, "LNC_TRANSCRIPT", 2000, 7)
    assert a["expected"] == pytest.approx(2 * 5 / 10 + 4 * 2 / 20) and a["observed"] == 2
    q = A07.bh([0.01, None, 0.04, 0.03])
    assert q[1] is None and q[0] == pytest.approx(0.03)


def test_uq_sample_config_is_valid(monkeypatch, tmp_path):
    repo = Path(__file__).resolve().parents[1]
    monkeypatch.setenv("PWD", str(tmp_path))
    cfg = load_config(str(repo / "examples" / "uq" / "gwas2m-lncrna.uq.yaml"))
    assert cfg.output_dir == (tmp_path / "results").resolve()
    assert cfg.resources_dir == (tmp_path / "resources").resolve()
    assert cfg.threads == 8 and cfg.jobs == 1


def test_top5_example_config(monkeypatch, tmp_path):
    repo = Path(__file__).resolve().parents[1]
    monkeypatch.setenv("PWD", str(tmp_path))
    cfg = load_config(str(repo / "examples" / "config_top5.yaml"))
    assert cfg["discovery"]["max_studies"] == 5
    assert cfg.output_dir == (tmp_path / "results").resolve()
