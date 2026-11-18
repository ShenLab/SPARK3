# this script is used: Beta-binomial TADA priors for LoF and damaging-missense (D-mis) variants, combining de novo, transmission (TDT), and case-control enrichment.

suppressPackageStartupMessages({ library(ggplot2); library(minpack.lm) })

IN_DIR  <- "beta_binomial_TADA/input"
OUT_DIR <- "beta_binomial_TADA/output"
# Fraction of genes assumed to be true autism risk genes; converts population-level enrichment into a per-gene relative risk: RR = 1 + (enrichment - 1) / PI_RISK.
PI_RISK <- 0.05
unfactor <- function(s) as.numeric(as.character(s))

# LGD regex match takes precedence; only synonymous-annotated variants are promoted to LGD by a high SpliceAI score (missense is left as is).
classify_cq <- function(consequence, splice_ai) {
  splice_ai <- unfactor(splice_ai)
  cq <- rep(".", length(consequence))
  cq[grepl("synonymous", consequence)] <- "syn"
  cq[grepl("missense", consequence)]   <- "mis"
  cq[grepl("frameshift|stop_gained|splice_donor|splice_acceptor", consequence)] <- "LGD"
  cq[grepl("synonymous", consequence) & !is.na(splice_ai) & splice_ai >= 0.8] <- "LGD"
  cq
}

# input Data
dnv   <- read.delim(file.path(IN_DIR, "autism_DNVs.txt"))
trios <- read.delim(file.path(IN_DIR, "autism_trio_samples.txt"))
aff_dnv_n   <- sum(trios$Pheno == "Affected")
unaff_dnv_n <- sum(trios$Pheno == "Unaffected")
dnv_aff   <- dnv[dnv$Pheno == "Affected", ]
dnv_unaff <- dnv[dnv$Pheno == "Unaffected", ]

inh <- read.delim(gzfile(file.path(IN_DIR, "inherited_variants_nonEUR_relatedEUR.txt.gz")))
inh$cq <- classify_cq(inh$consequence, inh$SpliceAI_max)
inh <- inh[inh$i_asd == 2, ]
trio_tdt <- read.delim(file.path(IN_DIR, "nonEUR_relatedEUR_trio_samples.txt"))
duo_tdt  <- read.delim(file.path(IN_DIR, "nonEUR_relatedEUR_duo_samples.txt"))
aff_inh_n <- sum(trio_tdt$asd == 2) + sum(duo_tdt$asd == 2)

pli <- read.delim(file.path(IN_DIR, "gnomAD4_pLI_GeneBayes_Shet.txt"))
pli <- pli[order(pli$Shet, decreasing = TRUE), ]

case_var <- read.delim(gzfile(file.path(IN_DIR, "cases_rmDNVs_rare_variants.txt.gz")))
ctrl_var <- read.delim(gzfile(file.path(IN_DIR, "UKBB_controls_rare_variants.txt.gz")))
case_n <- 33319; ctrl_n <- 142251
case_var$cq <- classify_cq(case_var$canon_Consequence, case_var$SpliceAI_max)
ctrl_var$cq <- classify_cq(ctrl_var$canon_Consequence, ctrl_var$SpliceAI_max)
case_var <- merge(case_var, pli, by.x = "canon_gene", by.y = "gene", all.x = TRUE)
ctrl_var <- merge(ctrl_var, pli, by.x = "canon_gene", by.y = "gene", all.x = TRUE)

# ---------------------------------------------------------------------------
# LoF prior: enrichment across Shet deciles, fit as a sigmoid in log-log space
# ---------------------------------------------------------------------------
shet_cuts <- sapply(1:9, function(i) pli$Shet[floor(nrow(pli) * i / 10)]) # Decile cutpoints, top (90th percentile) to bottom (10th percentile).
SHET_BIN_LABELS <- c("<10%", "10%-20%", "20%-30%", "30%-40%", "40%-50%", "50%-60%", "60%-70%", "70%-80%", "80%-90%", ">90%")

# Logical masks for 10 bins: >=cuts[1], [cuts[2],cuts[1]), ..., <cuts[9].
bin_masks <- function(score, cuts) {
  score <- unfactor(score)
  lo <- c(cuts, -Inf); hi <- c(Inf, cuts)
  lapply(seq_along(lo), function(i) !is.na(score) & score >= lo[i] & score < hi[i])
}
# Enrichment ratio (rate in group A / rate in group B) for each of the 10 bins.
enrichment_by_bin <- function(masks_a, n_a, masks_b, n_b) {
  mapply(function(ma, mb) (sum(ma) / n_a) / (sum(mb) / n_b), masks_a, masks_b)
}
lof_base <- function(df) df$cq == "LGD" & !is.na(df$LoF) & df$LoF == "HC"

dnv_lof <- enrichment_by_bin(lapply(bin_masks(dnv_aff$Shet, shet_cuts), `&`, lof_base(dnv_aff)), aff_dnv_n,
                              lapply(bin_masks(dnv_unaff$Shet, shet_cuts), `&`, lof_base(dnv_unaff)), unaff_dnv_n)
inh_lof <- enrichment_by_bin(lapply(bin_masks(inh$Shet, shet_cuts), `&`, lof_base(inh) & inh$iGT == 1), aff_inh_n,
                              lapply(bin_masks(inh$Shet, shet_cuts), `&`, lof_base(inh) & inh$iGT == 0), aff_inh_n)
cc_lof  <- enrichment_by_bin(lapply(bin_masks(case_var$Shet, shet_cuts), `&`, lof_base(case_var)), case_n,
                              lapply(bin_masks(ctrl_var$Shet, shet_cuts), `&`, lof_base(ctrl_var)), ctrl_n)

dnv_lof <- pmax(dnv_lof, 1.0025); inh_lof <- pmax(inh_lof, 1.0025); cc_lof <- pmax(cc_lof, 1.0025)

lof_long <- data.frame(
  Enrichment  = c(dnv_lof, inh_lof, cc_lof),
  Shet        = rep(c(shet_cuts, min(pli$Shet, na.rm = TRUE)), 3),
  Inheritance = factor(rep(c("DNV", "Inherited", "Case-control"), each = 10), levels = c("DNV", "Inherited", "Case-control")),
  Percent     = factor(rep(SHET_BIN_LABELS, 3), levels = SHET_BIN_LABELS)
)
lof_long$RR <- 1 + (lof_long$Enrichment - 1) / PI_RISK

sigmoid <- function(x, L, k, x0) log10(1.05) + L / (1 + exp(-k * (x - x0)))

fit_rr_curve <- function(shet, rr) {
  d <- data.frame(log_Shet = log10(shet), log_RR = log10(rr))
  nlsLM(log_RR ~ sigmoid(log_Shet, L, k, x0), data = d, start = list(L = 10, k = 0.1, x0 = -2))
}

# Sigmoid prediction per gene, clamped to the top bin's empirical RR above the top decile cutoff (no extrapolation past the fitted range) and floored at 1.05.
predict_gene_rr <- function(fit, gene_shet, top_cut, top_rr, floor_rr = 1.05) {
  gene_shet <- unfactor(gene_shet)
  rr <- 10 ^ predict(fit, newdata = data.frame(log_Shet = log10(gene_shet)))
  rr[!is.na(gene_shet) & gene_shet >= top_cut] <- top_rr
  rr[!is.na(rr) & rr <= floor_rr] <- floor_rr
  rr
}

plot_rr_fit <- function(shet, rr, fit, ylab) {
  d <- data.frame(log_Shet = log10(shet), log_RR = log10(rr))
  grid <- seq(min(d$log_Shet), max(d$log_Shet), length.out = 100)
  pred <- data.frame(log_Shet = grid,
                      predicted_log_RR = predict(fit, newdata = data.frame(log_Shet = grid)))
  ggplot(d, aes(log_Shet, log_RR)) +
    geom_point(size = 10, color = "blue") +
    geom_line(data = pred, aes(log_Shet, predicted_log_RR), color = "black", linewidth = 8) +
    labs(x = "Log10 (Shet)", y = ylab) + theme_minimal() +
    theme(text = element_text(size = 15), axis.text = element_text(size = 15))
}

dnv_rows <- lof_long[lof_long$Inheritance == "DNV", ]
inh_rows <- lof_long[lof_long$Inheritance == "Inherited", ]
cc_rows  <- lof_long[lof_long$Inheritance == "Case-control", ]

dnv_fit <- fit_rr_curve(dnv_rows$Shet, dnv_rows$RR)
inh_fit <- fit_rr_curve(inh_rows$Shet, inh_rows$RR)
cc_fit  <- fit_rr_curve(cc_rows$Shet,  cc_rows$RR)

pli$predicted_dnv_RR <- predict_gene_rr(dnv_fit, pli$Shet, shet_cuts[1], dnv_rows$RR[1])
pli$predicted_inh_RR <- predict_gene_rr(inh_fit, pli$Shet, shet_cuts[1], inh_rows$RR[1])
pli$predicted_cc_RR  <- predict_gene_rr(cc_fit,  pli$Shet, shet_cuts[1], cc_rows$RR[1])

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
pdf(file.path(OUT_DIR, "bbTADA_prior_LoF.pdf"), width = 10, height = 5)
lof_long$Percent <- factor(lof_long$Percent, levels = rev(levels(lof_long$Percent)))
print(
  ggplot(lof_long, aes(Percent, Enrichment, group = Inheritance, color = Inheritance)) +
    geom_line(linewidth = 1.5) + geom_point(size = 3) +
    scale_color_manual(values = c(DNV = "red", Inherited = "blue", `Case-control` = "green")) +
    labs(x = "window of % of genes ranking on Shet", y = "Enrichment") + theme_minimal() +
    theme(text = element_text(size = 20), axis.text = element_text(size = 10),
          axis.title = element_text(size = 22), legend.text = element_text(size = 18),
          legend.title = element_text(size = 20))
)
print(plot_rr_fit(dnv_rows$Shet, dnv_rows$RR, dnv_fit, "RR prior of LoF (de novo) - log10 scale"))
print(plot_rr_fit(inh_rows$Shet, inh_rows$RR, inh_fit, "RR prior of LoF (trans vs nontrans) - log10 scale"))
print(plot_rr_fit(cc_rows$Shet,  cc_rows$RR,  cc_fit,  "RR prior of LoF (case vs control) - log10 scale"))
dev.off()

# ---------------------------------------------------------------------------
# D-mis prior: a flat relative risk per MisFit-S / AlphaMissense tier (unlike LoF, not gene-specific -- the same value applies to every gene)
# ---------------------------------------------------------------------------

MISSENSE_BINS <- list(
  ">=0.1"              = function(df) df$cq == "mis" & !is.na(unfactor(df$MisFit1.5_S)) & unfactor(df$MisFit1.5_S) >= 0.1,
  "[0.03,0.1)"         = function(df) df$cq == "mis" & !is.na(unfactor(df$MisFit1.5_S)) &
                                       unfactor(df$MisFit1.5_S) < 0.1 & unfactor(df$MisFit1.5_S) >= 0.03,
  "MF1.5<0.03&AM>=0.7" = function(df) df$cq == "mis" & !is.na(unfactor(df$MisFit1.5_S)) & unfactor(df$MisFit1.5_S) < 0.03 &
                                       !is.na(unfactor(df$AlphaMissense)) & unfactor(df$AlphaMissense) >= 0.7
)

mis_enrichment <- function(df_a, n_a, df_b, n_b) {
  sapply(MISSENSE_BINS, function(f) (sum(f(df_a)) / n_a) / (sum(f(df_b)) / n_b))
}
mis_rr <- function(enrichment) pmax(1 + (enrichment - 1) / PI_RISK, 1.02)

dnv_mis_rr <- mis_rr(mis_enrichment(dnv_aff, aff_dnv_n, dnv_unaff, unaff_dnv_n))
inh_mis_rr <- mis_rr(mis_enrichment(inh[inh$iGT == 1, ], aff_inh_n, inh[inh$iGT == 0, ], aff_inh_n))
cc_mis_rr  <- mis_rr(mis_enrichment(case_var, case_n, ctrl_var, ctrl_n))

# ---------------------------------------------------------------------------
# Combined output
# ---------------------------------------------------------------------------

prior <- data.frame(
  gene = pli$gene, pLI = pli$pLI, Shet = pli$Shet,
  prior_dn_lof = pli$predicted_dnv_RR,
  prior_dn_mis1 = dnv_mis_rr[1], prior_dn_mis2 = dnv_mis_rr[2], prior_dn_mis3 = dnv_mis_rr[3],
  prior_inh_lof = pli$predicted_inh_RR,
  prior_inh_mis1 = inh_mis_rr[1], prior_inh_mis2 = inh_mis_rr[2], prior_inh_mis3 = inh_mis_rr[3],
  prior_cc_lof = pli$predicted_cc_RR,
  prior_cc_mis1 = cc_mis_rr[1], prior_cc_mis2 = cc_mis_rr[2], prior_cc_mis3 = cc_mis_rr[3]
)
write.table(prior, file.path(OUT_DIR, "bbTADA_prior_LoF_Dmis.txt"),sep = "\t", quote = FALSE, row.names = FALSE)
