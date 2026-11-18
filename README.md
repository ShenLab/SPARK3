# SPARK iWES3

Code for the analyses in:

> Qiao, et al., Rare coding variant burden contributes to phenotypic heterogeneity in autism.(2026).

This repository contains the custom code used for:

* de novo and inherited variant enrichment and transmission analysis
* population attributable fraction (PAF) of de novo and inherited variants
* the beta-binomial TADA (bbTADA) gene discovery model
* gene set overrepresentation analysis
* the rare variant composite score (RVCS) and its association with phenotypes

It also contains the scripts that draw the main figures of the paper.

Analyses were performed in R (version 4.4.0) and Python.

## Contents

```
SPARK3/
├── Analysis/                 enrichment, transmission, PAF, overrepresentation, RVCS
├── beta_binomial_TADA/       gene discovery model
│   ├── input/                variant tables, sample lists, mutation rates, constraint
│   ├── script/               model scripts (run in the order listed below)
│   └── output/               priors, gene counts, Bayes factors and FDR
├── Data/                     summary tables used by Analysis/ and PaperFigures/
└── PaperFigures/
    ├── script/               one script per main figure
    └── output/               figure PDFs
```

## Requirements

R 4.4.0 with the following packages.

CRAN:

```r
install.packages(c("dplyr", "tidyr", "data.table", "ggplot2", "patchwork", "ggrepel",
                   "ggpattern", "gggenes", "gridExtra", "RColorBrewer", "scales", "Cairo",
                   "openxlsx", "readxl", "readr", "stringr", "httr", "minpack.lm",
                   "VGAM", "epitools"))
```

Bioconductor:

```r
BiocManager::install(c("GenomicRanges", "GenomicFeatures", "rtracklayer", "Gviz", "trackViewer"))
```

All scripts use paths relative to the repository root. Run them from the top level of the cloned repository:

```bash
git clone https://github.com/ShenLab/SPARK3.git
cd SPARK3
```

## Beta-binomial TADA

bbTADA extends the TADA framework (He et al., 2013; Fu et al., 2022). It combines three sources of evidence for each gene:

* de novo variants in autistic probands, against the background mutation rate
* transmission of rare inherited variants from parents to autistic probands
* rare variants in autistic cases compared with UK Biobank controls

Evidence is computed separately for protein truncating variants (LoF) and three tiers of damaging missense variants (D-mis), defined by MisFit-S and AlphaMissense. The Bayes factors are combined per gene and converted to a posterior probability of association (PPA) and a false discovery rate (FDR).

The case-control component uses a beta-binomial model. Its overdispersion depends on the number of variants observed in a gene, and is calibrated on rare synonymous variants.

The full method is described in [beta_binomial_TADA/README.md](beta_binomial_TADA/README.md).

### Running

Run the four scripts in order:

```bash
Rscript beta_binomial_TADA/script/bbTADA_variant_count.R
Rscript beta_binomial_TADA/script/bbTADA_prior.R
Rscript beta_binomial_TADA/script/bbTADA_overdispersion_parameters.R
Rscript beta_binomial_TADA/script/bbTADA_BF_FDR.R
```

| Step | Script | What it does | Output |
|---|---|---|---|
| 1 | `bbTADA_variant_count.R` | Counts LoF and D-mis variants per gene from de novo, transmitted and case-control data, and joins the mutation rates | `bbTADA_gene_count.txt` |
| 2 | `bbTADA_prior.R` | Estimates the prior relative risk. LoF: a sigmoid fit across Shet deciles. D-mis: one value per tier | `bbTADA_prior_LoF_Dmis.txt`, `bbTADA_prior_LoF.pdf` |
| 3 | `bbTADA_overdispersion_parameters.R` | Fits the beta-binomial overdispersion ρ(m) = a / (a + b·m^c) by maximum likelihood | `ndep_rho_abc_by_class.txt`, `revised_overdispersion_phi_vs_n.pdf` |
| 4 | `bbTADA_BF_FDR.R` | Computes Bayes factors for each evidence type, then PPA, FDR and exome-wide P values | `bbTADA_BF_FDR.txt` |

`bbTADA_function.R` holds the Bayes factor functions. Step 4 loads it with `source()`.

Sample sizes are set at the top of `bbTADA_BF_FDR.R`:

| Parameter | Value | Meaning |
|---|---|---|
| `n_aff_dnv` | 34,247 | autistic probands with de novo calls |
| `n_unaff_dnv` | 10,666 | unaffected siblings with de novo calls |
| `n_aff_inh` | 18,898 | probands used for transmission |
| `n_case` | 33,319 | autistic cases in the case-control analysis |
| `n_ctrl` | 142,251 | UK Biobank controls |

### Input files

All files are tab separated. The `.gz` files are gzip compressed.

| File | Contents |
|---|---|
| `case_ctrl_VarCount_perGene.txt` | Per-gene synonymous, missense and LoF counts in cases and controls |
| `background_mutation_rate.txt` | No header. Gene, transcript, chromosome, then mutation rates for LoF and the three D-mis tiers |
| `gnomAD4_pLI_GeneBayes_Shet.txt` | Gene constraint: `gene`, `pLI`, `Shet` |

### Output

`bbTADA_BF_FDR.txt` has one row per gene. Key columns:

| Column | Meaning |
|---|---|
| `BF_dn_*`, `BF_in_*`, `BF_cc_*` | Bayes factors from de novo, inherited and case-control evidence, for `lof`, `mis1`, `mis2`, `mis3` |
| `BF` | Combined Bayes factor |
| `PPA` | Posterior probability of association |
| `FDR` | False discovery rate |
| `prior_*` | Prior relative risk used for each evidence type |

Genes at FDR ≤ 0.01 are saved in `Data/FDR0.01_bbTADA_genes.txt`.

## Other analyses

Each script runs on its own from the repository root.

| Script | Analysis | Input | Output |
|---|---|---|---|
| `Analysis/DNVs_fold_enrichment.R` | De novo variant enrichment in probands vs unaffected siblings (Poisson test), by MisFit-S tier and by sex and family history | `Data/autism_DNVs.txt`, `Data/autism_trio_samples.txt` | printed tables |
| `Analysis/Transmission_fold_enrichment.R` | Transmission disequilibrium test of inherited variants (binomial test against 50%), by MisFit-S and Shet bins | `Data/autism_transmission_variants.txt.gz` | printed tables |
| `Analysis/PAF_inherited_variants.R` | PAF of de novo and inherited LoF and D-mis burden, overall and by sex and family history | `Data/variant_rate_DNVs_CaseCtrl.xlsx` | `PaperFigures/output/PAF_*.pdf` |
| `Analysis/overrepresentation_analysis.R` | Allen Brain Atlas cell type and GO Cellular Component overrepresentation. Compares novel vs known and LoF vs missense predominant genes (Fisher's exact test), with Enrichr results for each group | `Data/FDR0.01_bbTADA_genes.txt` | `Data/novel_known_*.txt`, `Data/mis_lof_*.txt` |
| `Analysis/rare_variant_composite_score.R` | Per-individual rare variant composite score (RVCS) | see below | `<sample>_S.txt` |

`overrepresentation_analysis.R` downloads gene set libraries from the [Enrichr API](https://maayanlab.cloud/Enrichr/help#api), so it needs internet access.

`rare_variant_composite_score.R` runs on one individual at a time:

```bash
Rscript Analysis/rare_variant_composite_score.R <sample_id> <variant_dir>
```

`<variant_dir>` holds one variant file per individual, named `<sample_id>.txt`. The score is written to `<variant_dir>/<sample_id>_S.txt`.

## Figures

| Script | Figure |
|---|---|
| `PaperFigures/script/Fig1_fold_enrichment.R` | Fig. 1. De novo and inherited variant enrichment |
| `PaperFigures/script/Fig2_PAF_inherited_variants.R` | Fig. 2. Population attributable fraction |
| `PaperFigures/script/Fig3_association_evidence.R` | Fig. 3. Association evidence for known and novel genes |
| `PaperFigures/script/Fig4_overrepresentation.R` | Fig. 4. Brain region and cell type overrepresentation |
| `PaperFigures/script/Fig5_SORCS_lollipop.R` | Fig. 5. SORCS1, SORCS2 and SORCS3 variants along the protein |
| `PaperFigures/script/Fig6_RVCS.R` | Fig. 6. RVCS distribution and association with phenotypes |

Figures are written to `PaperFigures/output/`.

## Citation
The original TADA framework:
> He X, et al. Integrated model of de novo and inherited genetic variants yields greater power to identify risk genes. *PLoS Genet* 9, e1003671 (2013).
> Fu JM, et al. Rare coding variation provides insight into the genetic architecture and phenotypic context of autism. *Nat Genet* 54, 1320-1332 (2022).

