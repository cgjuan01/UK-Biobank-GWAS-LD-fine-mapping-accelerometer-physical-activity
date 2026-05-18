# UK Biobank ExWAS + LD Fine-Mapping: Vigorous vs Sedentary Activity

This repository provides a two-script pipeline to:

1. Run an **Exome-Wide Association Study (ExWAS)** of vigorous physical activity in UK Biobank whole-exome sequencing (WES) data using **REGENIE**
2. **Fine-map** the resulting association loci using WES-derived LD to identify the most probable causal coding variants

The fine-mapped variants from this pipeline are used as genetic instruments in a downstream multi-omic Mendelian randomisation analysis linking vigorous activity to biological ageing.

---

## Why REGENIE

REGENIE is used exclusively (rather than PLINK2 `--glm`) for three reasons:

**1. Sample size.** The analysis involves ~91,000 UK Biobank participants with WES data and valid accelerometer-derived physical activity measurements. At this scale, simple linear regression produces inflated test statistics because UK Biobank contains cryptically related individuals whose shared genetics and shared lifestyle would be confounded without correction.

**2. Mixed model correction.** REGENIE's step 1 fits a whole-genome regression model across all common array variants simultaneously, producing a per-participant polygenic background correction. Step 2 then tests each exonic variant conditioned on this correction — accounting for population stratification, cryptic relatedness, and polygenic background in a single framework.

**3. UK Biobank standard.** REGENIE is the methodology used in flagship UK Biobank exome-wide analyses (e.g. Backman et al. 2021 *Nature*) and is the recommended approach for quantitative traits at biobank scale.

---

## What the pipeline does

### Stage 1 — ExWAS (`ExWAS`)

**Phenotype preparation:**
- Loads participant table (CSV or Parquet) containing accelerometer-derived activity columns and covariates
- Normalises vigorous activity and sedentary activity independently to [0, 1]
- Inverts vigorous activity so that higher values correspond to more vigorous / less sedentary behaviour (sign harmonisation)
- Optionally adds a wear-quarter covariate derived from accelerometer date to control for seasonal variation
- Harmonises participant IDs with the WES `.fam` file (inner join — only participants present in both datasets are retained)
- Writes a REGENIE-ready `.phe` file

**REGENIE ExWAS:**
- **Step 1:** Whole-genome regression on array genotypes. Corrects for population structure, cryptic relatedness, and polygenic background. Outputs per-participant prediction files (LOCO — Leave One Chromosome Out)
- **Step 2:** Per-chromosome WES association testing conditioned on step 1 predictions. Produces association statistics (BETA, SE, LOG10P) for every exonic variant

The outcome is vigorous activity; sedentary activity is included as a covariate so that the ExWAS explicitly contrasts the two behaviours rather than testing general physical activity.

---

### Stage 2 — Fine-mapping (`fine-mapping`)

**Locus identification:**
- Clumps REGENIE step 2 results using PLINK2 (`--clump`) with WES-derived LD
- Default thresholds: lead variant p < 5×10⁻⁸, clumping r² < 0.1, window ±500 kb
- LD computed from the same ~91,000 WES participants (not a reference panel)

**Per-locus fine-mapping:**
For each independent lead variant:
- Extracts all variants in a ±500 kb window from WES PLINK genotypes
- Computes a full LD correlation matrix from WES genotypes (PLINK2 `--r square`)
- Prepares FINEMAP input files (`.z`, `.ld`, `.snp`, `.master`)
- Optionally runs FINEMAP (`--sss` stochastic search)
- Optionally runs SuSiE via a generated R script (returns per-variant Posterior Inclusion Probabilities and 95% credible sets)

**Why WES-derived LD:**
The LD matrix is computed from study participants rather than an external reference panel because exome sequencing captures rare and low-frequency coding variants that are absent or poorly represented in reference panels such as 1000 Genomes. This also ensures the fine-mapping is internally consistent with the instrument discovery. The per-locus LD matrices are subsequently used as the instrument correlation matrix in the downstream GLS-IVW Mendelian randomisation.

---

## Inputs

| File | Description |
|---|---|
| Participant CSV/Parquet | Table with `eid`, accelerometer activity columns, and covariates |
| WES PLINK files | `ukb23158_c{CHR}_b0_v1.*` (OQFE pipeline, one set per chromosome) |
| Array PLINK prefix | Common variant array data for REGENIE step 1 |
| WES `.fam` | Participant list from any WES chromosome file |

---

## Outputs

**ExWAS:**
```
results/
  regenie_step1/
    step1_pred.list        # Prediction file index for step 2
    *.pred                 # Per-participant polygenic correction files
    samples.keep           # Participant ID list used in analysis
  regenie_step2_wes/
    chr1_<phenotype>.regenie   # Per-chromosome association results
    chr2_<phenotype>.regenie
    ...
pa_wes_vig_vs_sed.phe          # Phenotype file used in REGENIE
```

**Fine-mapping:**
```
finemap_vig_vs_sed/
  clump.lead.tsv                  # Independent lead variants
  chr{N}_{BP}_500kb/              # One directory per locus
    finemap.z                     # Variant z-scores
    finemap.ld                    # LD correlation matrix (used in MR)
    finemap.snp                   # Variant order matching LD matrix
    finemap.master                # FINEMAP index file
    finemap.susie.pip.tsv         # SuSiE posterior inclusion probabilities
    finemap.susie.fit.rds         # SuSiE model object
    finemap.susie.cs.rds          # SuSiE credible sets
  loci_manifest.tsv               # Summary of all fine-mapped loci
```

---

## Usage

### Step 1 — ExWAS

```bash
python ExWAS \
    --source accel_covariates.csv \
    --fam ukb23158_c1_b0_v1.fam \
    --eid-col eid \
    --phenotype-col vigorous_activity \
    --phenotype-invert \
    --adjust-col sedentary_activity \
    --covars-base 31,21022,22009_1,22009_2,22009_3,22009_4,22009_5 \
    --phe-out pa_wes_vig_vs_sed.phe \
    --array-prefix /data/ukb_array_merged \
    --wes-prefix /data/ukb23158_c{CHR}_b0_v1 \
    --threads 8
```

Covariates passed to `--covars-base`:
- `31` — sex
- `21022` — age at recruitment
- `22009_1` to `22009_5` — first five genetic principal components (population stratification)

Wear quarter is added automatically if `--date-col` is provided.

### Step 2 — Combine REGENIE output (before fine-mapping)

REGENIE writes one output file per chromosome. Combine them before fine-mapping:

```bash
# Combine all chromosomes into one file
head -1 results/regenie_step2_wes/chr1_vigorous_activity_norm_inv.regenie \
    > results/regenie_step2_wes/combined.tsv

for i in $(seq 1 22); do
    tail -n +2 results/regenie_step2_wes/chr${i}_vigorous_activity_norm_inv.regenie
done >> results/regenie_step2_wes/combined.tsv
```

### Step 3 — Fine-mapping

```bash
python fine-mapping \
    --sumstats results/regenie_step2_wes/combined.tsv \
    --wes-prefix /data/ukb23158_c{CHR}_b0_v1 \
    --keep results/regenie_step1/samples.keep \
    --outdir finemap_vig_vs_sed \
    --n 91184 \
    --p-thresh 5e-8 \
    --run-susie
```

The `--n 91184` is the effective sample size read from the REGENIE step 2 log.

---

## Requirements

- Python ≥ 3.10 (`pandas`, `numpy`)
- [PLINK2](https://www.cog-genomics.org/plink/2.0/) in PATH
- [REGENIE](https://rgcgithub.github.io/regenie/) in PATH
- R with `susieR`, `data.table`, `Matrix` installed (for `--run-susie`)
- FINEMAP binary (for `--finemap-exe`; optional)

---

## Downstream use

The outputs of this pipeline feed directly into the multi-omic MR-GAT analysis:

- `loci_manifest.tsv` + SuSiE PIP files → select the 87 fine-mapped VPA instruments (`fine-mapped-pa-snps_fixed.csv`)
- Per-locus `finemap.ld` matrices → combined into `GS_LD_square.csv` used as the instrument correlation matrix in the GLS-IVW MR estimator
- REGENIE step 2 betas and SEs → become `beta.exposure` and `se.exposure` in the protein/CpG/glycan/single-cell MR analyses

---

## Reference

UK Biobank whole-exome sequencing data: Karczewski et al. (2022); Backman et al. (2021) *Nature*

REGENIE: Mbatchou et al. (2021) *Nature Genetics*

SuSiE: Wang et al. (2020) *Journal of the Royal Statistical Society Series B*

FINEMAP: Benner et al. (2016) *Bioinformatics*
