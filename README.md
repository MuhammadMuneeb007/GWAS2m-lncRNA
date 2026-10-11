# GWAS2m-lncRNA

**A Transcript-Aware Framework for Characterizing Long Non-Coding RNA Overlap in Fine-Mapped GWAS Signals**

From a phenotype name to fine-mapped, transcript-annotated GWAS variants and an
lncRNA-focused report, with one command.

```
gwas2m-lncrna run      -p "coronary artery disease" -a EUR   # steps 1-7
gwas2m-lncrna analysis -p "coronary artery disease" -a EUR   # report
```

| Step | Command | What it does |
|---|---|---|
| 1 | `discover` | Finds GWAS Catalog studies for the phenotype with full, harmonised summary statistics and a discovery sample of only the requested ancestry. Ranks them by match quality and sample size. |
| 2 | `download` | Downloads each study's harmonised summary statistics (resumable). |
| 3 | `qc` | GWASLab standardisation and QC on GRCh38. |
| 4 | `ld-reference` | Builds the 1000 Genomes LD panel for the ancestry (once, shared by all phenotypes). |
| 5 | `clump` | PLINK2 LD clumping (p1 5e-8, p2 0.01, r² 0.1, 1 Mb). |
| 6 | `finemap` | SuSiE-RSS fine-mapping of ±1 Mb loci with an ancestry-matched LD matrix. |
| 7 | `annotate` | Ensembl VEP annotation of every fine-mapped variant, all transcripts. |
| – | `analysis` | Report: `report.html`, `report_tables.xlsx`, `variants.tsv.gz`. |

Each GWAS study is processed separately and never merged with other studies.
A study that fails is recorded with its reason and the run carries on with the rest.

## Results: coronary artery disease, Parkinson's disease and Alzheimer's disease

These results come from the accompanying study (manuscript in preparation). That study ran Steps 01–07 and
Analyses 01–09 on European-ancestry GWAS Catalog studies for three diseases.

**Datasets.** 58 GWAS were selected (20 CAD, 18 PD, 20 AD):

* 21 completed fine-mapping and annotation: 2 CAD, 8 PD and 11 AD.
* 24 were processed but had no genome-wide significant locus.
* 13 failed or were only partly fine-mapped, all of them CAD.

The 21 included studies gave 501 fine-mapped study-loci, 1,063 95% credible sets and 11,640 credible-set
variant–study records.

**lncRNA transcript overlap.** Overlap requires the lncRNA biotype and consequence on the same VEP transcript
row.

| | CAD | PD | AD | All |
|---|---|---|---|---|
| Credible-set variant–study records | 848 | 290 | 10,502 | 11,640 |
| Overlapping an lncRNA transcript | 88 (10.4%) | 117 (40.3%) | 2,194 (20.9%) | 2,399 (20.6%) |
| … in an lncRNA intron | 79 (9.3%) | 115 (39.7%) | 2,082 (19.8%) | 2,276 (19.6%) |
| … in an lncRNA exon | 21 (2.5%) | 8 (2.8%) | 254 (2.4%) | 283 (2.4%) |
| Protein-coding CDS (strict) | 18 (2.1%) | 2 (0.7%) | 265 (2.5%) | 285 (2.4%) |
| Within-locus observed / expected lncRNA overlap | 0.34 | 1.36 | 0.79 | – |

**Other findings:**

* **Overlap is mostly intronic.** Of the lncRNA-overlapping records, 95% are intronic and 12% exonic. A variant
  can be exonic in one isoform and intronic in another.
* **Locus windows overstate overlap.** 499 of the 501 study-loci contain both lncRNA and protein-coding gene
  bodies somewhere in their ~2 Mb window. But only 250 of the 454 loci with credible-set variants have one that
  actually lies in an lncRNA transcript.
* **Robust to PIP weighting.** PIP-weighted overlap is 23.4%, and among regional variants with PIP ≥ 0.5 it is
  23.2% (169 of 728).
* **Enrichment is study-dependent, not consistent.** The direct-labelled PD study GCST90476961 shows 58
  overlaps against 25.0 expected (fold 2.32, q = 0.0015). CAD shows fewer overlaps than expected.
* **Limited exact-allele sharing.** 239 credible-set alleles are shared between CAD and AD, 20 between PD and
  AD, and none between CAD and PD. One example is `17:45892531:A:C` in *MAPT-AS1*, found in both PD and AD.

**Caveats:**

* Counts are variant–study records, and GWAS of the same disease can share cohorts.
* Both included CAD studies are related or interaction phenotypes, not direct clinical CAD.
* The enrichment null shuffles labels within each locus. It is not LD-matched or genome-wide.
* Transcript overlap is not evidence that an lncRNA causes disease.

## Install

You need [mamba](https://github.com/conda-forge/miniforge) (or conda) on Linux. The bioinformatics tools
(PLINK2, VEP) are not built for Windows.

```bash
git clone https://github.com/MuhammadMuneeb007/GWAS2m-lncRNA.git gwas2m-lncrna
cd gwas2m-lncrna
mamba env create -f environment.yml      # Python package + PLINK2 + R/susieR + VEP 116
mamba activate gwas2m-lncrna
gwas2m-lncrna check                           # lists tools and reference data
```

## Reference data

These are downloaded automatically the first time a step needs them, so there is nothing extra to do.
To fetch them in advance (recommended on a cluster, see below):

```bash
gwas2m-lncrna setup --ancestry EUR
```

| Resource | Download | Used by |
|---|---|---|
| 1000 Genomes GRCh38 VCFs (chr1–22) and the per-ancestry PLINK2 panel | ~15 GB | steps 4–6 |
| Ensembl VEP 116 GRCh38 cache | 27.6 GB, unpacked then deleted | step 7 |
| GENCODE v50 GRCh38 genome (+ index) | 0.85 GB (3.1 GB unpacked) | step 7 |
| GENCODE v50 gene annotation | 0.12 GB | report |
| GWAS Catalog metadata | 0.15 GB | step 1 |

Plan for about 80 GB free while the VEP cache is unpacked and about 50 GB afterwards. Everything goes into
`resources_dir` (default `./resources`). If you already have a GWAS2m `resources/` folder, point
`resources_dir` at it and nothing is downloaded again.

## Configuration

Everything works with the built-in defaults. To change settings, create a configuration file and edit it:

```bash
gwas2m-lncrna init            # writes ./gwas2m-lncrna.yaml with every option documented
```

`gwas2m-lncrna` reads `./gwas2m-lncrna.yaml` automatically, or the file given with `--config` or by `$GWAS2M_LNCRNA_CONFIG`.
The file holds the output and resource folders, threads, tool paths and the scientific parameters for each
step. The phenotype and ancestry are given on the command line, so one file serves every run. Command-line
`--threads`, `--jobs` and `--output-dir` override the file.

## Running

```bash
gwas2m-lncrna run -p "Parkinson's disease" -a EUR              # all steps
gwas2m-lncrna run -p "Parkinson's disease" -a EUR --analysis   # all steps, then the report
gwas2m-lncrna run -p asthma -p migraine -a EUR                 # several phenotypes in turn
gwas2m-lncrna run -p CAD -a EUR --from-step 5                  # resume from clumping
gwas2m-lncrna finemap -p CAD -a EUR --study GCSTxxxxxxxx       # one step, one study
gwas2m-lncrna status -p CAD -a EUR                             # per-study status
```

Re-running is safe. Finished work is detected and skipped, interrupted downloads resume, and fine-mapping
reuses every locus that already finished. Use `--force` to redo work, and `--refresh-discovery` to search the
GWAS Catalog again. Otherwise the study list from the first discovery is reused, so the catalog changing does
not alter a project halfway through.

Useful discovery options are `--alias "coronary heart disease"` (extra search terms) and `--max-studies 5`.

## The report

```bash
gwas2m-lncrna analysis -p "coronary artery disease" -p "Parkinson's disease" -p "Alzheimer's disease" -a EUR
```

This runs the nine analysis files `Analysis01_John_...` to `Analysis09_John_...` in order, for the
phenotypes you give, and writes three files to `results/reports/<phenotypes>_<ANC>/`:

* **`report.html`** opens in any browser. It has one section per analysis, with its tables and its full
  printed output.
* **`report_tables.xlsx`** has every table from every analysis, one sheet each (`A01_studies`,
  `A05_distribution`, …), plus an `about` sheet saying which file produced each sheet.
* **`variants.tsv.gz`** is Analysis 03's per-variant classification of credible-set and high-PIP variants.

Each analysis file can also be run on its own, exactly as before. It prints its results and returns its tables:

```bash
python -m gwas2m_lncrna.analysis.Analysis05_John_lncRNA_Exon_Intron_Distribution     -p "coronary artery disease" -p "Parkinson's disease" -a EUR
```

Results you already produced with GWAS2m use the same folder layout, so you can analyse them directly by
pointing `output_dir` (or `--root`) at that GWAS2m folder.

## Using a cluster (SLURM, PBS, LSF)

gwas2m-lncrna never submits jobs itself. It runs in whatever process you start, so you put the normal command
inside your scheduler's job script. Ready-to-edit examples are in [`examples/`](examples/):

* `slurm_one_phenotype.sh` runs one phenotype end to end.
* `slurm_phenotype_array.sh` runs one array task per line of `phenotypes.txt`.

Run `gwas2m-lncrna setup --ancestry EUR` once before starting many jobs, so they do not all wait on the same
downloads. Within a job, `threads` sets how many cores each tool uses and `jobs` sets how many GWAS studies run
side by side. Request about `threads × jobs` cores.

## Output layout

The same layout as the original GWAS2m pipeline:

```
results/
  01_gwas_catalog/<phenotype>/<ancestry>/   GWAS_download_manifest.tsv, study_selection.tsv, discovery_summary.json
  02_summary_stats/<phenotype>/<ancestry>/  raw/<GCST>/ downloads, qc/<GCST>/ GWASLab output
  05_ld_clumping/<phenotype>/<ancestry>/<GCST>/
  06_finemapping/<phenotype>/<ancestry>/<GCST>/
  07_annotation/<phenotype>/<ancestry>/<GCST>/
  pipeline_audit/<phenotype>/<ancestry>/status/<GCST>/<Stage>.json   per-study step status
  pipeline_audit/<phenotype>/<ancestry>/logs/<GCST>/<Stage>.log      per-study step logs
  reports/                                                           analysis reports
```

## Package contents

All the code lives in `src/gwas2m_lncrna/`. Your scripts are inside the package under their original names.

```
steps/
  Step01_Plan_GWAS.py                          discover   (GWAS Catalog search)
  Step02_Download_GWAS.py                      download
  Step03_QC_One_GWAS.py                        qc         (GWASLab)
  Step04_Prepare_LD_Reference.py               ld-reference
  Step05_Clump_GWAS_Loci.py                    clump      (PLINK2)
  Step06_FineMap_GWAS_Loci.py                  finemap    (SuSiE-RSS; R code in data/run_susie_rss.R)
  Step07_Annotate_Finemapped_Variants.py       annotate   (VEP)
  runner.py                                    runs a step for every study
analysis/
  Analysis01_John_GWAS_Selection.py
  Analysis02_John_Finemapping_CredibleSets.py
  Analysis03_John_Coding_Noncoding_Distribution.py
  Analysis04_John_Genomic_Consequence_Distribution.py
  Analysis05_John_lncRNA_Exon_Intron_Distribution.py
  Analysis06_John_lncRNA_Locus_Representation.py
  Analysis07_John_lncRNA_Enrichment.py
  Analysis08_John_Cross_Disease_Comparison.py
  Analysis09_John_Finemapping_Sensitivity.py
  helpers.py                                   phenotype labels + running one analysis on its own
  report.py                                    runs Analysis01-09 and writes the three report files
cli.py                                         the gwas2m-lncrna command
pipeline.py                                    runs Steps 01-07 in order (replaces Step01_10_Run.py)
resources.py                                   downloads/builds reference data (replaces Step00 + Setup_*.sh)
config.py, project.py, download.py, tools.py, common.py
```

## Changes from the original scripts

* There is no SLURM code. Steps run in the current process, and every setting lives in one YAML file.
* Step 04 builds the 1000 Genomes panel where Steps 05 and 06 read it. Before, it wrote it to a different folder.
* Step 05 no longer crashes when no candidate on a chromosome maps to the reference panel.
* Step 05 no longer maps a GWAS variant onto a different variant at the same position when the alleles
  disagree, for example an indel onto a SNV.
* Step 06 carries the `estimate_s_rss` LD-mismatch diagnostic into its per-locus summary.
* Step 07 no longer treats VEP's "-" placeholder as a MANE transcript when choosing the best consequence.
* Each analysis file keeps its own logic, with these changes:
  * It works for whatever phenotypes and ancestry you pass, instead of being fixed to CAD, Parkinson's and
    Alzheimer's in European ancestry.
  * "Direct disease label" means the reported trait matches your phenotype or its aliases.
  * Cross-disease tables have one column per phenotype.
  * `run()` returns the result tables so they can be saved.

## Tests

```bash
pytest
```

The tests cover the configuration, the original folder layout, the step logic (including the bug fixes above),
resumable downloads, and all nine analyses run end to end on synthetic pipeline outputs.

## Interpretation

The report lists its limitations. Two matter most:

* Overlap percentages have no genome-wide baseline. The `ALL_REGIONAL` population, meaning every fine-mapped
  variant in the same loci, is the only comparison built in.
* The enrichment test is a within-locus null. It does not account for LD or allele frequency.

Transcript overlap does not show that an lncRNA mediates risk.
