### Over-dispersion of the case/control split, as a function of per-gene genotype count m = case_count + ctrl_count.
### Model: x_case | m ~ BetaBinom(m, r*M, (1-r)*M),  M = (1-rho)/rho,  rho(m) = 1 / (1 + b*m^c)
### a is fixed at 1: rho(m) = a/(a+b*m^c) is invariant to scaling (a,b) by a common constant, so a alone is never identifiable from data -- only b/a and c are. Fixing a=1 removes that redundant
### degree of freedom and leaves a well-posed 2-parameter MLE: a single optim() call per class, and the result does not depend on where it starts.
### b, c are fit here by MLE, separately for synonymous, missense (all) and LGD (all), from rare synonymous counts in case and control file (output/case_ctrl_VarCount_perGene.txt)

library(dplyr);library(ggplot2);library(scales);library(VGAM)  # VGAM: dbetabinom.ab
n_case <- 33319; n_ctrl <- 142251
r <- n_case / (n_case + n_ctrl)   # null case fraction, shared across all genes/classes
counts <- read.csv('beta_binomial_TADA/input/case_ctrl_VarCount_perGene.txt', header = TRUE, sep = '\t')

## previously fit constant (n-independent) beta-binomial -- used for the green reference line below
fit<-Coef(VGAM::vglm(cbind(counts$syn_count_case,counts$syn_count_ctrl)~1,betabinomial,trace=FALSE))
mu<-fit[1];rho<-fit[2]
alpha_const<-(1/rho-1)*mu;beta_const<-1/rho-1-alpha_const  ##alpha_const is shape1, beta_const is shape2
rho_const   <- 1 / (alpha_const + beta_const + 1)

## pad to the full background gene list: case_ctrl_VarCount_perGene.txt only has genes that carry at least one qualifying variant (any class) in case_var1 or ctrl_var1 -- genes with literally
## zero qualifying variants of every class, in both cases and controls, never appear there at all. Add them back in with all six counts = 0, so `counts` covers every gene in the background list.
all_genes <- read.table('beta_binomial_TADA/input/background_mutation_rate.txt', header = F);all_genes<-all_genes[,1:2];colnames(all_genes)<-c('gene','transID')
zero_genes <- setdiff(all_genes$gene, counts$gene)
zero_rows <- as.data.frame(matrix(0, nrow = length(zero_genes), ncol = ncol(counts) - 1))
colnames(zero_rows) <- setdiff(colnames(counts), 'gene')
zero_rows$gene <- zero_genes
counts <- bind_rows(counts, zero_rows[, colnames(counts)])
cat(sprintf('Added %d genes with zero qualifying variants (any class, both case & control); counts now covers %d genes\n', length(zero_genes), nrow(counts)))

## ---- MLE of b, c (a fixed at 1) for one variant class -----------------------------------------
neg_log_lik <- function(case_ct, m, r, b, c) {
  rho <- 1 / (1 + b * m^c)
  M   <- (1 - rho) / rho
  -sum(dbetabinom.ab(case_ct, m, r * M, (1 - r) * M, log = TRUE))
}
fit_ndep_bc <- function(case_ct, ctrl_ct, r) {
  m <- case_ct + ctrl_ct
  keep <- m > 0 & !is.na(m)
  case_ct <- case_ct[keep]; m <- m[keep]
  
  fit <- optim(c(0, 0), function(log_par) neg_log_lik(case_ct, m, r, exp(log_par[1]), exp(log_par[2])),
               method = 'Nelder-Mead')  # log(b)=log(c)=0, i.e. b=c=1
  c(a = 1, setNames(exp(fit$par), c('b', 'c')))
}

## ---- lambda (genomic inflation) of the one-sided mid-p test under rho(m) = 1/(1+b*m^c) --------
## Vectorized over genes (pbetabinom.ab/dbetabinom.ab accept vector args).
calc_ndep_bb_pvalue <- function(case_ct, m, r, b, c) {
  pval <- rep(1, length(m))
  pos <- m > 0
  rho   <- 1 / (1 + b * m[pos]^c)
  M     <- (1 - rho) / rho
  alpha <- r * M; beta <- (1 - r) * M
  sf  <- 1 - pbetabinom.ab(case_ct[pos] - 1, m[pos], alpha, beta)
  pmf <- dbetabinom.ab(case_ct[pos], m[pos], alpha, beta)
  pval[pos] <- pmax(pmin(sf - 0.5 * pmf, 1), .Machine$double.xmin)
  pval
}
lambda_of <- function(case_ct, m, r, b, c) {
  pv <- calc_ndep_bb_pvalue(case_ct, m, r, b, c)
  median(qnorm(1 - pv)^2, na.rm = TRUE) / qchisq(0.5, df = 1)
}

## To calibrarte the overdispersion of rare synonymous in cases vs UKBB controls per gene, we make that one-sided mid-p test's genomic inflation land at 1, not just maximize the count
## likelihood. Starting the nudge from its own unconstrained MLE optimum (not an arbitrary guess) means a single optim() call is enough -- that point is already close to lambda=1, so no
## multi-start search is needed to find this one reliably.
fit_syn_calibrated <- function(case_ct, ctrl_ct, r, lambda_weight = 1e6) {
  mle <- fit_ndep_bc(case_ct, ctrl_ct, r)
  
  m <- case_ct + ctrl_ct
  keep <- m > 0 & !is.na(m)
  case_ct <- case_ct[keep]; m <- m[keep]
  
  obj <- function(log_par) {
    b <- exp(log_par[1]); c <- exp(log_par[2])
    lambda_weight * (lambda_of(case_ct, m, r, b, c) - 1)^2 + neg_log_lik(case_ct, m, r, b, c)
  }
  fit <- optim(log(c(mle['b'], mle['c'])), obj, method = 'Nelder-Mead')
  b <- unname(exp(fit$par[1])); c <- unname(exp(fit$par[2]))
  c(a = 1, b = b, c = c, lambda = lambda_of(case_ct, m, r, b, c))
}

syn_cal <- fit_syn_calibrated(counts$syn_count_case, counts$syn_count_ctrl, r)
abc_syn <- syn_cal[c('a', 'b', 'c')]
cat(sprintf('Synonymous calibrated: a=%.4f b=%.4f c=%.4f (lambda = %.4f)\n', abc_syn['a'], abc_syn['b'], abc_syn['c'], syn_cal['lambda']))
abc_mis <- fit_ndep_bc(counts$mis_count_case, counts$mis_count_ctrl, r)
abc_lgd <- fit_ndep_bc(counts$LGD_count_case, counts$LGD_count_ctrl, r)

abc_table <- rbind(Synonymous = abc_syn, `Missense (all)` = abc_mis, `LGD (all)` = abc_lgd)
cat('n-dependent rho(m) = a / (a + b * m^c), fit separately by class (a fixed at 1):\n')
print(round(abc_table, 4))
write.table(data.frame(class = rownames(abc_table), abc_table), 'beta_binomial_TADA/output/ndep_rho_abc_by_class.txt', col.names = TRUE, row.names = FALSE, sep = '\t', quote = FALSE)

## ---- empirical phi(n): observed variance / binomial variance, by n-bin ----------------------
## standardize each gene's case count assuming the null binomial (mean m*r, var m*r*(1-r)); the variance of these standardized residuals within an n-bin IS phi by definition (phi=1 if the
## data really were binomial; phi>1 reflects the over-dispersion rho(m) is modeling).
## Bins are equal-GENE-COUNT windows (via ntile on rank of m). ntile() guarantees ~n_genes/n_bins genes per window regardless of how skewed m is, which is what makes each point's phi estimate
## comparably stable across the x-axis.
empirical_phi <- function(case_ct, ctrl_ct, r, n_bins = 25) {
  m <- case_ct + ctrl_ct
  keep <- m > 0 & !is.na(m)
  case_ct <- case_ct[keep]; m <- m[keep]
  z <- (case_ct - m * r) / sqrt(m * r * (1 - r))
  bin <- ntile(m, n_bins)
  
  data.frame(m = m, z = z, bin = bin) %>%
    group_by(bin) %>%
    summarise(n_genes = n(), m_lo = min(m), m_hi = max(m), m_mid = exp(mean(log(m))),
              phi = var(z), .groups = 'drop')
}

## one observed series per class, matching the three fitted curves
obs_syn <- empirical_phi(counts$syn_count_case, counts$syn_count_ctrl, r)
obs_mis <- empirical_phi(counts$mis_count_case, counts$mis_count_ctrl, r)
obs_lgd <- empirical_phi(counts$LGD_count_case, counts$LGD_count_ctrl, r)
obs_syn$group <- 'Observed synonymous'
obs_mis$group <- 'Observed missense'
obs_lgd$group <- 'Observed LGD'
obs_all <- bind_rows(obs_syn, obs_mis, obs_lgd)

## ---- fitted curves, evaluated on a dense grid for plotting -----------------------------------
## one n-dependent curve per class, each using that class's own fitted a,b,c (rather than the single synonymous-calibrated curve used previously)
n_grid <- exp(seq(log(1), log(max(c(obs_syn$m_mid, obs_mis$m_mid)) * 1.3), length.out = 400))
rho_of <- function(abc, n) abc['a'] / (abc['a'] + abc['b'] * n ^ abc['c'])

curve_long <- bind_rows(
  data.frame(n = n_grid, y = 1,
             series = 'Binomial: 1'),
  data.frame(n = n_grid, y = 1 + (n_grid - 1) * rho_const,
             series = sprintf('Constant beta-binomial: 1+(n-1)·%.4f', rho_const)),
  data.frame(n = n_grid, y = 1 + (n_grid - 1) * rho_of(abc_syn, n_grid),
             series = 'n-dependent BB, synonymous'),
  data.frame(n = n_grid, y = 1 + (n_grid - 1) * rho_of(abc_mis, n_grid),
             series = 'n-dependent BB, missense')
)
series_levels <- unique(curve_long$series)
curve_long$series <- factor(curve_long$series, levels = series_levels)
line_colors <- setNames(c('#E0772C', '#2FA86B', '#3B7DD8', '#7A52A6'), series_levels)

## ---- plot: observed vs. modeled over-dispersion, by variant count ---------------------------
p <- ggplot() +
  geom_line(data = curve_long, aes(x = n, y = y, color = series), linewidth = 1) +
  geom_point(data = obs_all, aes(x = m_mid, y = phi, shape = group),
             color = 'black', size = 2.3, stroke = 1) +
  scale_color_manual(name = NULL, values = line_colors) +
  scale_shape_manual(name = NULL, values = c('Observed synonymous' = 16, 'Observed missense' = 0)) +
  scale_x_log10(breaks = trans_breaks('log10', function(x) 10^x),
                labels = trans_format('log10', math_format(10^.x))) +
  scale_y_log10(breaks = c(1, 2, 3, 5, 10, 20, 40)) +
  labs(x = 'Variants in the gene, n (cases + controls)',
       y = 'Variance / binomial variance',
       title = NULL) +
  theme_classic(base_size = 18) +
  theme(axis.title = element_text(size = 20), axis.text = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.position = c(0.02, 0.98), legend.justification = c(0, 1),
        legend.direction = 'vertical', legend.box = 'vertical',
        legend.background = element_rect(fill = alpha('white', 0.7), color = NA))
print(p)

pdf('beta_binomial_TADA/output/revised_overdispersion_phi_vs_n.pdf', width = 8, height = 7)
print(p)
dev.off()

