import threading
from functools import partial
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer

import pandas as pd
import pytest

from gwas2m_lncrna.analysis.report import run_analysis
from gwas2m_lncrna.cli import main
from gwas2m_lncrna.download import download
from gwas2m_lncrna.project import Project
from gwas2m_lncrna.steps.runner import run_for_studies

from conftest import build_manifest, build_study, vep_row

LNC_INTRON = [vep_row("", "lncRNA", "intron_variant&non_coding_transcript_variant", intron="1/3", symbol="LINC01")]
PC_MISSENSE = [vep_row("", "protein_coding", "missense_variant", exon="2/9", symbol="PCG1")]
INTERGENIC = [vep_row("", "-", "intergenic_variant", gene="-", symbol="-", feature="-")]
NEAR_LNC = [vep_row("", "lncRNA", "upstream_gene_variant", symbol="LINC02")]


def _variants(locus, start, specs):
    return [{"LOCUS_ID": locus, "POS": start + i, "PIP": pip, "CS": cs, "rows": rows}
            for i, (pip, cs, rows) in enumerate(specs)]


@pytest.fixture
def two_phenotypes(make_config, tmp_path):
    gtf = tmp_path / "genes.gtf"
    gtf.write_text(
        '1\tHAVANA\tgene\t500\t3000000\t.\t+\t.\tgene_id "ENSG1"; gene_type "lncRNA"; gene_name "LINC01";\n'
        '1\tHAVANA\tgene\t900\t1200\t.\t+\t.\tgene_id "ENSG2"; gene_type "protein_coding"; gene_name "PCG1";\n')
    cfg = make_config(resources={"gencode_gtf": str(gtf), "auto_download": False},
                      analysis={"permutations": 500})
    cad = Project(cfg, "coronary artery disease", "EUR")
    build_manifest(cad, [("GCST1", "Coronary artery disease"),
                         ("GCST2", "Coronary artery disease x smoking interaction"),
                         ("GCST3", "Coronary artery disease")])
    build_study(cad, "GCST1", _variants("L001", 1000, [
        (0.6, True, LNC_INTRON), (0.3, True, PC_MISSENSE), (0.01, False, INTERGENIC), (0.02, False, LNC_INTRON)]))
    build_study(cad, "GCST2", _variants("L001", 1000, [(0.9, True, LNC_INTRON), (0.05, False, NEAR_LNC)]))
    cad.set_status("GCST3", "download", "FAILED", "HTTP 404")
    pd_project = Project(cfg, "Parkinson's disease", "EUR")
    build_manifest(pd_project, [("GCST9", "Parkinson's disease")])
    build_study(pd_project, "GCST9", _variants("L001", 1000, [(0.95, True, LNC_INTRON), (0.01, False, PC_MISSENSE)]))
    return cfg, tmp_path


def test_report_runs_all_nine_analyses(two_phenotypes):
    cfg, tmp_path = two_phenotypes
    files = run_analysis(cfg, ["coronary artery disease", "Parkinson's disease"], "EUR")
    assert sorted(p.name for p in files["html"].parent.iterdir()) == ["report.html", "report_tables.xlsx",
                                                                      "variants.tsv.gz"]
    page = files["html"].read_text(encoding="utf-8")
    assert "stopped with an error" not in page
    for n in range(1, 10):
        assert f"Analysis0{n}_John_" in page

    sheets = pd.read_excel(files["xlsx"], sheet_name=None)
    a01 = sheets["A01_studies"].set_index("STUDY_ACCESSION")
    assert a01.loc["GCST3", "ANALYSIS_STATUS"] == "UPSTREAM_FAILED"
    assert a01.loc["GCST1", "TRAIT_CLASS"] == "DIRECT_DISEASE"
    assert a01.loc["GCST2", "TRAIT_CLASS"] == "SPECIALIZED"

    a05 = sheets["A05_studies"].set_index("STUDY_ACCESSION")
    assert a05.loc["GCST1", "CS_VARIANTS"] == 2 and a05.loc["GCST1", "CS_LNCRNA_OVERLAP"] == 1
    dist = sheets["A05_distribution"]
    near = dist[(dist.PHENOTYPE == "coronary artery disease") & (dist.POPULATION == "ALL_REGIONAL")
                & (dist.MEASURE == "NEAR_LNCRNA_ONLY")].iloc[0]
    assert near.N == 1  # upstream of an lncRNA is "near", not overlap

    assert {"A02_studies", "A04_categories", "A06_locus_classes", "A07_tests", "A09_populations"} <= set(sheets)
    assert sheets["A08_pairwise_sharing"].iloc[0]["SHARED_CS_ALLELES"] == 1  # 1:1000:A:G in both
    assert set(sheets["A06_loci"]["GENCODE_CLASS"]) == {"BOTH"}

    variants = pd.read_csv(files["variants"], sep="\t")
    assert len(variants) == 4 and set(variants["CS95"]) == {1}


def test_single_analysis_file_runs_standalone(two_phenotypes, capsys):
    from gwas2m_lncrna.analysis import Analysis05_John_lncRNA_Exon_Intron_Distribution as A05
    cfg, _ = two_phenotypes
    tables = A05.main(["-p", "Parkinson's disease", "-a", "EUR", "--config", str(cfg.source)])
    assert tables["studies"].iloc[0]["CS_LNCRNA_OVERLAP"] == 1
    assert "PARKINSON'S DISEASE" in capsys.readouterr().out


def test_cli_status_and_analysis(two_phenotypes, capsys):
    cfg, tmp_path = two_phenotypes
    config = str(cfg.source)
    assert main(["status", "--config", config, "-p", "coronary artery disease", "-a", "EUR"]) == 0
    assert "GCST3" in capsys.readouterr().out
    assert main(["analysis", "--config", config, "-p", "Parkinson's disease", "-a", "EUR",
                 "--out", str(tmp_path / "rep")]) == 0
    assert (tmp_path / "rep" / "report.html").is_file()
    assert main(["analysis", "--config", config, "-p", "migraine", "-a", "EUR"]) == 2  # clean error


def _boom(project, study, force):
    if study["STUDY_ACCESSION"] == "GCST2":
        raise RuntimeError("bad file")
    return "COMPLETE", "", {"X": 1}


def test_runner_isolates_failures(make_config):
    project = Project(make_config(), "migraine", "EUR")
    studies = [{"STUDY_ACCESSION": a} for a in ("GCST1", "GCST2", "GCST3")]
    tally = run_for_studies(project, "qc", "TEST", _boom, studies)
    assert tally == {"COMPLETE": 2, "FAILED": 1}
    assert project.step_status("GCST2", "qc") == "FAILED"
    assert "bad file" in project.status("GCST2")["qc"]["reason"]
    assert "RuntimeError" in project.log_file("qc", "GCST2").read_text()


class RangeHandler(SimpleHTTPRequestHandler):
    def log_message(self, *args):
        pass

    def do_GET(self):
        data = (self.directory_path / self.path.lstrip("/")).read_bytes()
        start = 0
        if "Range" in self.headers:
            start = int(self.headers["Range"].split("=")[1].split("-")[0])
            self.send_response(206)
        else:
            self.send_response(200)
        self.send_header("Content-Length", str(len(data) - start))
        self.end_headers()
        self.wfile.write(data[start:])


def test_download_resumes(tmp_path):
    served = tmp_path / "served"
    served.mkdir()
    payload = bytes(range(256)) * 4000
    (served / "file.bin").write_bytes(payload)
    handler = type("H", (RangeHandler,), {"directory_path": served})
    server = ThreadingHTTPServer(("127.0.0.1", 0), partial(handler, directory=str(served)))
    threading.Thread(target=server.serve_forever, daemon=True).start()
    try:
        target = tmp_path / "out" / "file.bin"
        target.parent.mkdir()
        (tmp_path / "out" / "file.bin.part").write_bytes(payload[:1000])
        url = f"http://127.0.0.1:{server.server_address[1]}/file.bin"
        assert download(url, target) == "COMPLETE"
        assert target.read_bytes() == payload
        assert download(url, target) == "SKIPPED_ALREADY_COMPLETE"
    finally:
        server.shutdown()
