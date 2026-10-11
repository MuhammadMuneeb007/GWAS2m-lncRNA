"""Shared reference data: download on demand, build once, reuse across runs.

Layout under resources_dir (compatible with an existing GWAS2m `resources/`):

    1000G/raw/ALL.chr<N>....vcf.gz                 raw 1000 Genomes GRCh38 VCFs
    1000G/metadata/integrated_call_samples...panel  sample -> super-population
    1000G/<ANC>/chr<N>_<ANC>_GRCh38.{pgen,pvar,psam} ancestry LD panels
    genome/GRCh38/GRCh38.primary_assembly.genome.fa(.fai)
    gencode/release50/gencode.v50.primary_assembly.annotation.gtf.gz
    vep/cache/homo_sapiens/116_GRCh38/
    gwas_catalog/                                   GWAS Catalog metadata cache
"""

from __future__ import annotations

import gzip
import os
import shutil
import subprocess
import tarfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd

from gwas2m_lncrna.common import (
    ANCESTRY_LABELS, Gwas2RnaError, banner, file_lock, nonempty, run_command,
)
from gwas2m_lncrna.config import Config
from gwas2m_lncrna.download import download, http_session
from gwas2m_lncrna.tools import find_tool

ONEKG_BASE = (
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/"
    "1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL"
)
ONEKG_PANEL_URL = (
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"
    "integrated_call_samples_v3.20130502.ALL.panel"
)
GENCODE_RELEASE = "50"
GENCODE_BASE = f"https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{GENCODE_RELEASE}"
FASTA_URL = f"{GENCODE_BASE}/GRCh38.primary_assembly.genome.fa.gz"
GTF_URL = f"{GENCODE_BASE}/gencode.v{GENCODE_RELEASE}.primary_assembly.annotation.gtf.gz"
CHROMOSOMES = range(1, 23)


def vep_cache_url(release: int | str) -> str:
    return (f"https://ftp.ensembl.org/pub/release-{release}/variation/indexed_vep_cache/"
            f"homo_sapiens_vep_{release}_GRCh38.tar.gz")


def onekg_vcf_name(chrom: int) -> str:
    return f"ALL.chr{chrom}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"


class Resources:
    def __init__(self, cfg: Config):
        self.cfg = cfg
        self.root = cfg.resources_dir
        self.settings = cfg["resources"]

    # ------------------------------------------------------------------ paths
    @property
    def onekg_dir(self) -> Path:
        return self.root / "1000G"

    @property
    def ld_dir(self) -> Path:
        return self.cfg.path(self.settings.get("ld_reference_dir")) or self.onekg_dir

    @property
    def panel(self) -> Path:
        return self.onekg_dir / "metadata" / "integrated_call_samples_v3.20130502.ALL.panel"

    def raw_vcf(self, chrom: int) -> Path:
        return self.onekg_dir / "raw" / onekg_vcf_name(chrom)

    def ld_prefix(self, ancestry: str, chrom: int) -> Path:
        return self.ld_dir / ancestry / f"chr{chrom}_{ancestry}_GRCh38"

    @property
    def fasta(self) -> Path:
        return (self.cfg.path(self.settings.get("fasta"))
                or self.root / "genome" / "GRCh38" / "GRCh38.primary_assembly.genome.fa")

    @property
    def gtf(self) -> Path:
        return (self.cfg.path(self.settings.get("gencode_gtf"))
                or self.root / "gencode" / f"release{GENCODE_RELEASE}"
                / f"gencode.v{GENCODE_RELEASE}.primary_assembly.annotation.gtf.gz")

    @property
    def vep_cache(self) -> Path:
        return self.cfg.path(self.settings.get("vep_cache_dir")) or self.root / "vep" / "cache"

    @property
    def vep_release(self) -> str:
        return str(self.cfg["annotation"]["vep_cache_release"])

    @property
    def gwas_catalog_dir(self) -> Path:
        return self.root / "gwas_catalog"

    # ----------------------------------------------------------------- checks
    def ld_ready(self, ancestry: str) -> bool:
        return all(nonempty(Path(f"{self.ld_prefix(ancestry, c)}{ext}"))
                   for c in CHROMOSOMES for ext in (".pgen", ".pvar", ".psam"))

    def fasta_ready(self) -> bool:
        return nonempty(self.fasta) and nonempty(Path(str(self.fasta) + ".fai"))

    def vep_ready(self) -> bool:
        return (self.vep_cache / "homo_sapiens" / f"{self.vep_release}_GRCh38").is_dir()

    def gtf_ready(self) -> bool:
        return nonempty(self.gtf)

    def status_rows(self) -> list[dict]:
        rows = [{"resource": f"1000G LD panel {anc}", "ready": self.ld_ready(anc),
                 "location": str(self.ld_dir / anc)} for anc in ANCESTRY_LABELS]
        rows += [
            {"resource": "GRCh38 FASTA (+.fai)", "ready": self.fasta_ready(), "location": str(self.fasta)},
            {"resource": f"VEP cache {self.vep_release}_GRCh38", "ready": self.vep_ready(),
             "location": str(self.vep_cache)},
            {"resource": f"GENCODE v{GENCODE_RELEASE} GTF", "ready": self.gtf_ready(), "location": str(self.gtf)},
        ]
        return rows

    def _missing(self, what: str, ancestry: str | None = None) -> None:
        if not self.settings.get("auto_download", True):
            flag = f" --ancestry {ancestry}" if ancestry else ""
            raise Gwas2RnaError(
                f"{what} is missing under {self.root} and resources.auto_download is false. "
                f"Run `gwas2m-lncrna setup{flag}` (or point the paths in gwas2m-lncrna.yaml at existing copies)."
            )

    # ------------------------------------------------------------- 1000 Genomes
    def ensure_ld(self, ancestry: str) -> None:
        if self.ld_ready(ancestry):
            return
        self._missing(f"The 1000 Genomes {ancestry} LD panel", ancestry)
        if self.settings.get("ld_reference_dir"):
            raise Gwas2RnaError(f"resources.ld_reference_dir has no complete {ancestry} panel: {self.ld_dir / ancestry}")
        banner(f"BUILDING 1000 GENOMES {ancestry} LD PANEL (one-time)")
        plink2 = find_tool(self.cfg, "plink2")
        session = http_session()
        with file_lock(self.onekg_dir / ".download.lock"):
            download(ONEKG_PANEL_URL, self.panel, session)
            todo = [c for c in CHROMOSOMES if not self._chrom_ready(ancestry, c)]
            jobs = [(f"{ONEKG_BASE}/{onekg_vcf_name(c)}", self.raw_vcf(c)) for c in todo]
            with ThreadPoolExecutor(max_workers=max(1, int(self.settings["download_workers"]))) as pool:
                list(pool.map(lambda job: download(job[0], job[1], http_session()), jobs))
        samples = self._sample_list(ancestry)
        with file_lock(self.ld_dir / f".{ancestry}.build.lock"):
            for chrom in CHROMOSOMES:
                self._build_chrom(plink2, ancestry, chrom, samples)
        if not self.settings.get("keep_raw_1000g_vcf", True):
            for chrom in CHROMOSOMES:
                self.raw_vcf(chrom).unlink(missing_ok=True)
        if not self.ld_ready(ancestry):
            raise Gwas2RnaError(f"{ancestry} LD panel is still incomplete under {self.ld_dir / ancestry}")

    def _chrom_ready(self, ancestry: str, chrom: int) -> bool:
        prefix = self.ld_prefix(ancestry, chrom)
        return all(nonempty(Path(f"{prefix}{ext}")) for ext in (".pgen", ".pvar", ".psam"))

    def _sample_list(self, ancestry: str) -> Path:
        panel = pd.read_csv(self.panel, sep="\t", dtype=str)
        if "super_pop" not in panel.columns or "sample" not in panel.columns:
            raise Gwas2RnaError(f"Unexpected 1000 Genomes panel format: {self.panel}")
        selected = panel.loc[panel["super_pop"] == ancestry, "sample"].dropna()
        if selected.empty:
            raise Gwas2RnaError(f"No {ancestry} samples in {self.panel}")
        out = self.ld_dir / ancestry / "samples.txt"
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text("#IID\n" + "\n".join(selected) + "\n", encoding="utf-8")
        print(f"{ancestry} reference samples: {len(selected)}")
        return out

    def _build_chrom(self, plink2: Path, ancestry: str, chrom: int, samples: Path) -> None:
        if self._chrom_ready(ancestry, chrom):
            return
        prefix = self.ld_prefix(ancestry, chrom)
        building = Path(str(prefix) + ".building")
        run_command([
            plink2, "--vcf", self.raw_vcf(chrom), "--keep", samples,
            "--set-all-var-ids", "@:#:$r:$a", "--new-id-max-allele-len", "1000",
            "--rm-dup", "force-first", "--max-alleles", "2",
            "--make-pgen", "--threads", self.cfg.threads, "--out", building,
        ], Path(str(prefix) + ".build.log"))
        # Rename only after PLINK2 finished: existing files always mean "complete".
        for ext in (".pgen", ".pvar", ".psam"):
            os.replace(f"{building}{ext}", f"{prefix}{ext}")
        Path(f"{building}.log").unlink(missing_ok=True)

    # ------------------------------------------------------------------ genome
    def ensure_fasta(self) -> None:
        if self.fasta_ready():
            return
        self._missing("The indexed GRCh38 FASTA")
        if self.settings.get("fasta"):
            raise Gwas2RnaError(f"resources.fasta must be an existing FASTA with a .fai index: {self.fasta}")
        banner("PREPARING GRCh38 FASTA (one-time)")
        samtools = find_tool(self.cfg, "samtools")
        with file_lock(self.fasta.parent / ".fasta.lock"):
            if not nonempty(self.fasta):
                gz = Path(str(self.fasta) + ".gz")
                download(FASTA_URL, gz)
                partial = Path(str(self.fasta) + ".part")
                with gzip.open(gz, "rb") as source, open(partial, "wb") as sink:
                    shutil.copyfileobj(source, sink, length=16 * 1024 * 1024)
                os.replace(partial, self.fasta)
                gz.unlink(missing_ok=True)
            if not nonempty(Path(str(self.fasta) + ".fai")):
                run_command([samtools, "faidx", self.fasta])

    def ensure_gtf(self) -> None:
        if self.gtf_ready():
            return
        self._missing("The GENCODE gene annotation")
        if self.settings.get("gencode_gtf"):
            raise Gwas2RnaError(f"resources.gencode_gtf does not exist: {self.gtf}")
        with file_lock(self.gtf.parent / ".gtf.lock"):
            download(GTF_URL, self.gtf)

    # --------------------------------------------------------------------- VEP
    def ensure_vep_cache(self) -> None:
        if self.vep_ready():
            return
        self._missing(f"The VEP {self.vep_release} GRCh38 cache")
        if self.settings.get("vep_cache_dir"):
            raise Gwas2RnaError(
                f"resources.vep_cache_dir has no homo_sapiens/{self.vep_release}_GRCh38: {self.vep_cache}")
        banner(f"PREPARING VEP {self.vep_release} CACHE (one-time, 27.6 GB download)")
        archive = self.root / "vep" / f"homo_sapiens_vep_{self.vep_release}_GRCh38.tar.gz"
        with file_lock(self.root / "vep" / ".vep.lock"):
            if self.vep_ready():
                return
            download(vep_cache_url(self.vep_release), archive)
            self.vep_cache.mkdir(parents=True, exist_ok=True)
            print(f"Unpacking {archive.name} ...", flush=True)
            if shutil.which("tar"):
                subprocess.run(["tar", "-xzf", str(archive), "-C", str(self.vep_cache)], check=True)
            else:
                with tarfile.open(archive) as tar:
                    tar.extractall(self.vep_cache, filter="data")
            if not self.vep_ready():
                raise Gwas2RnaError(f"VEP archive unpacked but {self.vep_release}_GRCh38 not found in {self.vep_cache}")
            if self.settings.get("delete_vep_archive", True):
                archive.unlink(missing_ok=True)

    # ------------------------------------------------------------------ public
    def ensure(self, ancestry: str | None, ld: bool = False, vep: bool = False, gtf: bool = False) -> None:
        if ld:
            if ancestry is None:
                raise Gwas2RnaError("An ancestry is required to prepare an LD panel.")
            self.ensure_ld(ancestry)
        if vep:
            self.ensure_fasta()
            self.ensure_vep_cache()
        if gtf:
            self.ensure_gtf()

    def count_reference_samples(self, ancestry: str) -> int:
        psam = Path(f"{self.ld_prefix(ancestry, 1)}.psam")
        with psam.open(encoding="utf-8", errors="replace") as handle:
            count = sum(1 for line in handle if line.strip() and not line.startswith("#"))
        if count <= 0:
            raise Gwas2RnaError(f"Could not count samples in {psam}")
        return count
