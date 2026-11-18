# Per-gene variant counts (LoF / D-mis tiers) from de novo, transmission (TDT), and case-control data, combined with the background mutation rate table for the beta-binomial TADA model.

IN_DIR  <- "beta_binomial_TADA/input"
OUT_DIR <- "beta_binomial_TADA/output"

unfactor <- function(s) as.numeric(as.character(s))
classify_cq <- function(consequence, splice_ai) {
  splice_ai <- unfactor(splice_ai)
  cq <- rep(".", length(consequence))
  cq[grepl("synonymous", consequence)] <- "syn"
  cq[grepl("missense", consequence)]   <- "mis"
  cq[grepl("frameshift|stop_gained|splice_donor|splice_acceptor|start_lost|stop_lost", consequence)] <- "LGD"
  cq[grepl("synonymous", consequence) & !is.na(splice_ai) & splice_ai >= 0.8] <- "LGD"
  cq
}
# Three D-mis tiers by MisFit1.5_S, with the lowest tier requiring AlphaMissense support.
MISSENSE_BINS <- list(
  mis1 = function(df) df$cq == "mis" & !is.na(unfactor(df$MisFit1.5_S)) & unfactor(df$MisFit1.5_S) >= 0.1,
  mis2 = function(df) df$cq == "mis" & !is.na(unfactor(df$MisFit1.5_S)) &
                       unfactor(df$MisFit1.5_S) < 0.1 & unfactor(df$MisFit1.5_S) >= 0.03,
  mis3 = function(df) df$cq == "mis" & !is.na(unfactor(df$MisFit1.5_S)) & unfactor(df$MisFit1.5_S) < 0.03 &
                       !is.na(unfactor(df$AlphaMissense)) & unfactor(df$AlphaMissense) >= 0.7
)

dnv <- read.delim(file.path(IN_DIR, "autism_DNVs.txt"))
dnv_aff <- dnv[dnv$Pheno == "Affected", ]; dnv_unaff <- dnv[dnv$Pheno == "Unaffected", ]
base <- read.delim(file.path(IN_DIR, "background_mutation_rate.txt"), header = FALSE,col.names = c("Gene", "Transcript", "Chrom", "Mu_lgd", "Mu_mis1", "Mu_mis2", "Mu_mis3"))

inh <- read.delim(gzfile(file.path(IN_DIR, "inherited_variants_nonEUR_relatedEUR.txt.gz")))
inh$cq <- classify_cq(inh$consequence, inh$SpliceAI_max)
inh <- inh[inh$i_asd == 2, ]
inh_trans <- inh[inh$iGT == 1, ]; inh_untrans <- inh[inh$iGT == 0, ]

case_var <- read.delim(gzfile(file.path(IN_DIR, "cases_rmDNVs_rare_variants.txt.gz")))
ctrl_var <- read.delim(gzfile(file.path(IN_DIR, "UKBB_controls_rare_variants.txt.gz")))
case_var$cq <- classify_cq(case_var$canon_Consequence, case_var$SpliceAI_max)
ctrl_var$cq <- classify_cq(ctrl_var$canon_Consequence, ctrl_var$SpliceAI_max)

genes <- unique(base$Gene)
count_by_gene <- function(gene, mask, genes) as.integer(table(factor(gene[mask], levels = genes)))
source_counts <- function(df, gene_col, prefix) {
  g <- df[[gene_col]]
  out <- data.frame(
    lof  = count_by_gene(g, df$cq == "LGD", genes),
    mis1 = count_by_gene(g, MISSENSE_BINS$mis1(df), genes),
    mis2 = count_by_gene(g, MISSENSE_BINS$mis2(df), genes),
    mis3 = count_by_gene(g, MISSENSE_BINS$mis3(df), genes)
  )
  setNames(out, paste0(prefix, c("_lof", "_mis1", "_mis2", "_mis3")))
}

counts <- cbind(
  source_counts(dnv_aff,     "HGNC",       "dn_aff"),
  source_counts(dnv_unaff,   "HGNC",       "dn_unaff"),
  source_counts(inh_trans,   "gene",       "inh_trans"),
  source_counts(inh_untrans, "gene",       "inh_untrans"),
  source_counts(case_var,    "canon_gene", "cc_case"),
  source_counts(ctrl_var,    "canon_gene", "cc_ctrl")
)

x <- cbind(base, counts)
x$Mu_lgd[x$Mu_lgd == 0]   <- 1e-8
x$Mu_mis1[x$Mu_mis1 == 0] <- 1e-8
x$Mu_mis2[x$Mu_mis2 == 0] <- 1e-8

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
write.table(x, file.path(OUT_DIR, "bbTADA_gene_count.txt"),sep = "\t", quote = FALSE, row.names = FALSE)
