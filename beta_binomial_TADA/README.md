# Beta-binomial TADA (bbTADA)

bbTADA is a gene discovery model for rare coding variants. It extends TADA (He et al., 2013), using the version from Fu et al. (2022).

For each gene it combines three kinds of evidence:

| Evidence | Comparison | Model |
|---|---|---|
| De novo (`dn`) | probands vs siblings, against the mutation rate | Poisson-Gamma |
| Inherited (`in`) | transmitted vs untransmitted parental alleles in probands | Poisson-Gamma |
| Case-control (`cc`) | autistic cases vs UK Biobank controls | Poisson-Gamma total count, **beta-binomial** case/control split |

The new part is the case-control model. Rare synonymous variants in cases vs controls are over-dispersed relative to a binomial. bbTADA replaces the binomial case/control split with a beta-binomial. Its over-dispersion depends on how many variants the gene carries, and is calibrated on synonymous variants.

## Running

Run from the repository root, in this order:

```bash
Rscript beta_binomial_TADA/script/bbTADA_variant_count.R            # 1. per-gene counts
Rscript beta_binomial_TADA/script/bbTADA_prior.R                    # 2. priors
Rscript beta_binomial_TADA/script/bbTADA_overdispersion_parameters.R # 3. over-dispersion
Rscript beta_binomial_TADA/script/bbTADA_BF_FDR.R                   # 4. BF, PPA, FDR
```

`bbTADA_function.R` holds all model functions. Step 4 loads it with `source()`.

R packages: `ggplot2`, `minpack.lm`, `dplyr`, `scales`, `VGAM`, `GenomicRanges`, `stringr`, `openxlsx`, `readxl`.

## 1. Variant classes

Variants are grouped into one LoF class and three damaging missense (D-mis) tiers.

| Class | Definition |
|---|---|
| `lof` | frameshift, stop gained, splice donor, splice acceptor, start lost, stop lost; or synonymous with SpliceAI ≥ 0.8 |
| `mis1` | missense, MisFit-S ≥ 0.1 |
| `mis2` | missense, 0.03 ≤ MisFit-S < 0.1 |
| `mis3` | missense, MisFit-S < 0.03 and AlphaMissense ≥ 0.7 |

Each class has its own mutation rate per gene (`Mu_lgd`, `Mu_mis1`, `Mu_mis2`, `Mu_mis3`). Zero rates are set to 1e-8.

`bbTADA_variant_count.R` counts each class per gene in six groups: de novo in probands and siblings, transmitted and untransmitted alleles, cases and controls.

## 2. Priors

The prior for each gene is its expected relative risk $\bar\gamma$, if it is a risk gene.

**From enrichment to relative risk.** For a group of variants, the enrichment $E$ is the rate in probands (or cases, or transmitted alleles) divided by the rate in the comparison group. Only a fraction $\pi$ of genes are risk genes, so the per-gene relative risk is

$$
\bar\gamma = 1 + \frac{E - 1}{\pi}, \qquad \pi = 0.05 .
$$

This is computed separately for de novo, inherited and case-control data.

**LoF prior (gene specific).** Genes are split into 10 bins by gnomAD v4 GeneBayes $s_{het}$. Enrichment of high-confidence LoF variants (LOFTEE HC) is computed in each bin, floored at 1.0025, and converted to $\bar\gamma$. A sigmoid is fit on the log10 scale with `nlsLM`:

$$
\log_{10}\bar\gamma = \log_{10}(1.05) + \frac{L}{1 + e^{-k(\log_{10}s_{het} - x_0)}} .
$$

Each gene gets the value predicted from its own $s_{het}$. Genes above the top decile cutoff get the top bin's empirical value, so the curve is not extrapolated. The minimum is 1.05.

**D-mis prior (flat).** Each tier gets one value for all genes, from that tier's enrichment. The minimum is 1.02.

Output: `output/bbTADA_prior_LoF_Dmis.txt` and the fit plots in `output/bbTADA_prior_LoF.pdf`.

## 3. Bayes factors

Each evidence type gives a Bayes factor (BF) comparing $H_1$ (risk gene) with $H_0$ (not a risk gene).

### De novo

For $N$ trios and mutation rate $\mu$, the de novo count is

$$
x \sim \text{Poisson}(2N\mu\gamma), \qquad \gamma \mid H_1 \sim \text{Gamma}(\bar\gamma\beta,\ \beta), \qquad \gamma = 1 \text{ under } H_0 ,
$$

with $\beta = 0.2$. Integrating over $\gamma$ gives a negative binomial under $H_1$.

A BF is computed for probands and for unaffected siblings. Each is floored at 1. The de novo BF is

$$
BF_{dn} = \max\left(\frac{BF_{probands}}{BF_{siblings}},\ 1\right).
$$

So de novo variants in siblings reduce the evidence.

### Inherited

Transmitted alleles play the role of cases. Untransmitted alleles play the role of controls. Both use $n$ = the number of probands. For a gene with untransmitted count $x_{0}$ and transmitted count $x_{1}$:

$$
x_{0} \sim \text{Poisson}(qn), \qquad x_{1} \sim \text{Poisson}(q\gamma n), \qquad q \sim \text{Gamma}(\rho,\ \nu),
$$

with $\gamma = \bar\gamma$ under $H_1$ and $\gamma = 1$ under $H_0$. The rate $q$ is integrated out, which gives negative binomial marginal likelihoods. The hyperparameters are

$$
\nu = 5000, \qquad \rho_g = \frac{\nu \sum x_{0}}{2Gn}\cdot\frac{\mu_g}{\bar\mu},
$$

where $G$ is the number of genes. So $\rho$ scales with each gene's mutation rate.

### Case-control (beta-binomial)

The original TADA case-control model is the same Poisson-Gamma model as above, with $n_{case}$ and $n_{ctrl}$. Its likelihood splits exactly into two parts. Let $m = x_{case} + x_{ctrl}$:

$$
P(x_{case}, x_{ctrl}) = \underbrace{P(m)}_{\text{total count}} \times \underbrace{P(x_{case} \mid m)}_{\text{case/control split}} .
$$

1. **Total count.** $m \sim \text{NB}\left(\rho,\ \frac{\nu}{\nu + n_{ctrl} + \gamma n_{case}}\right)$, as in TADA.
2. **Split.** In TADA this is $\text{Binomial}(m, p)$ with

$$
p_1 = \frac{\bar\gamma\, n_{case}}{\bar\gamma\, n_{case} + n_{ctrl}} \ \ (H_1), \qquad p_0 = \frac{n_{case}}{n_{case} + n_{ctrl}} \ \ (H_0).
$$

bbTADA keeps part 1 and replaces the binomial in part 2 with a beta-binomial:

$$
x_{case} \mid m \sim \text{BetaBinomial}\big(m,\ pM,\ (1-p)M\big), \qquad M = \frac{1-\rho_g}{\rho_g}.
$$

The beta-binomial has the same mean as the binomial and a larger variance. Its intra-class correlation is

$$
\rho_g = \frac{a}{a + b\, m_{syn}^{\,c}},
$$

where $m_{syn}$ is the gene's rare synonymous count in cases plus controls. So over-dispersion is treated as a property of the gene. It is estimated where no signal is expected (synonymous variants), and the same $\rho_g$ is used for LoF and all D-mis tiers of that gene.

The log BF is

$$
\log BF_{cc} = \big[\log P(m \mid H_1) - \log P(m \mid H_0)\big] + \big[\log P(x_{case} \mid m, p_1) - \log P(x_{case} \mid m, p_0)\big].
$$

Special cases:

* Genes with no synonymous count ($m_{syn}$ = 0 or missing) use the binomial.
* As $b \to \infty$, $\rho_g \to 0$ and bbTADA returns exactly the original TADA BF.
* Genes with no variants or a zero mutation rate get BF = 1.

## 4. Over-dispersion parameters

`bbTADA_overdispersion_parameters.R` estimates $a$, $b$ and $c$.

**Model.** For each gene, $x_{case} \mid m \sim \text{BetaBinomial}(m, rM, (1-r)M)$, with $r = n_{case}/(n_{case}+n_{ctrl})$ and $\rho(m) = a/(a + b\,m^c)$.

**Identifiability.** Multiplying $a$ and $b$ by the same constant leaves $\rho(m)$ unchanged. So $a$ is fixed at 1, and only $b$ and $c$ are estimated.

**Fit.** All genes in the background list are included. Genes without any variant are added with zero counts. $b$ and $c$ are estimated by maximum likelihood (`optim`, Nelder-Mead, on the log scale).

**Calibration on synonymous variants.** The synonymous fit is then adjusted so that a one-sided beta-binomial mid-p test of each gene gives a genomic inflation factor $\lambda = 1$. The objective is

$$
10^6\,(\lambda - 1)^2 - \log L(b, c),
$$

starting from the maximum likelihood estimate.

Missense and LoF fits are also reported for comparison. Only the synonymous parameters are used in the BF.

| Class | a | b | c |
|---|---|---|---|
| **Synonymous (used)** | 1 | 1.740 | 0.818 |
| Missense (all) | 1 | 1.178 | 0.867 |
| LoF (all) | 1 | 3.269 | 0.639 |

Output: `output/ndep_rho_abc_by_class.txt`, and `output/revised_overdispersion_phi_vs_n.pdf`, which compares observed and modeled over-dispersion by variant count.

## 5. Combining evidence, PPA and FDR

**Per class.** The three evidence types are multiplied:

$$
BF_{class} = \max\left(BF_{dn} \cdot BF_{in} \cdot BF_{cc},\ 1\right), \qquad class \in \{lof, mis1, mis2, mis3\}.
$$

**Per gene.** $BF = BF_{lof} \cdot BF_{mis1} \cdot BF_{mis2} \cdot BF_{mis3}$.

**PPA.** With the prior fraction of risk genes $\pi = 0.05$:

$$
PPA = \frac{\pi\, BF}{1 - \pi + \pi\, BF}.
$$

**Bayesian FDR.** Genes are sorted by BF, largest first. The FDR of the top $k$ genes is the average posterior probability of the null among them:

$$
FDR_k = \frac{1}{k}\sum_{i=1}^{k}(1 - PPA_i).
$$

**Pseudo P value.** For comparison with frequentist results, FDR is converted back to a P value scale: $P_k = FDR_k \cdot k / N$, with $N$ = 18,698 genes.

## Parameters

| Parameter | Value | Where |
|---|---|---|
| $\pi$ (fraction of risk genes) | 0.05 | `bbTADA_prior.R` (`PI_RISK`), `bbTADA_BF_FDR.R` (`pi0 = 0.95`) |
| $\beta$ (de novo prior dispersion) | 0.2 | `bbTADA_BF_FDR.R` |
| $\nu$ (inherited and case-control prior) | 5000 | `bbTADA_function.R` |
| Probands with de novo calls | 34,247 | `bbTADA_BF_FDR.R` |
| Siblings with de novo calls | 10,666 | `bbTADA_BF_FDR.R` |
| Probands used for transmission | 18,898 | `bbTADA_BF_FDR.R` |
| Autistic cases | 33,319 | `bbTADA_BF_FDR.R` |
| UK Biobank controls | 142,251 | `bbTADA_BF_FDR.R` |
| $N$ for pseudo P | 18,698 | `bbTADA_BF_FDR.R` |

## Input files (`input/`)

| File | Used in | Contents |
|---|---|---|
| `case_ctrl_VarCount_perGene.txt` | 3, 4 | per-gene synonymous, missense and LoF counts in cases and controls |
| `background_mutation_rate.txt` | 1, 3 | no header: gene, transcript, chromosome, `Mu_lgd`, `Mu_mis1`, `Mu_mis2`, `Mu_mis3` |
| `gnomAD4_pLI_GeneBayes_Shet.txt` | 2 | `gene`, `pLI`, `Shet` |

## Output files (`output/`)

| File | Step | Contents |
|---|---|---|
| `bbTADA_gene_count.txt` | 1 | mutation rates and counts per gene and class |
| `bbTADA_prior_LoF_Dmis.txt` | 2 | prior relative risk per gene, evidence type and class |
| `bbTADA_prior_LoF.pdf` | 2 | LoF enrichment by $s_{het}$ bin and sigmoid fits |
| `ndep_rho_abc_by_class.txt` | 3 | $a$, $b$, $c$ by variant class |
| `revised_overdispersion_phi_vs_n.pdf` | 3 | observed vs modeled over-dispersion |
| `bbTADA_BF_FDR.txt` | 4 | BFs (`BF_dn_*`, `BF_in_*`, `BF_cc_*`, `BF_*`, `BF`), `PPA`, `FDR`, `pseudo_p`, priors and counts |

## References

* He X, et al. Integrated model of de novo and inherited genetic variants yields greater power to identify risk genes. *PLoS Genet* 9, e1003671 (2013).
* Fu JM, et al. Rare coding variation provides insight into the genetic architecture and phenotypic context of autism. *Nat Genet* 54, 1320-1332 (2022).
* Bayesian FDR as in the original TADA code (PMID: 19822692, Section 2.3).

Original TADA code: Xin He, Jiebiao Wang, Lambertus Klei and Jack Fu. Beta-binomial extension: Lu Qiao.
