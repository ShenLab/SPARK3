# Transmission disequilibrium test (TDT) for inherited variants in autism trios.

INH_FILE <- "Data/autism_transmission_variants.txt.gz"

# Mutually exclusive bins stratifying missense variants by MisFit-S and LoF variants by Shet (gene-level selection coefficient).
MISFIT_TIERS <- list(
  "S>=0.03"       = function(s) s >= 0.03,
  "0.005<=S<0.03" = function(s) s >= 0.005 & s < 0.03,
  "1e-4<=S<0.005" = function(s) s >= 1e-4  & s < 0.005,
  "S<1e-4"        = function(s) s < 1e-4
)
SHET_TIERS <- list(
  "Shet>=0.1"       = function(s) s >= 0.1,
  "0.03<=Shet<0.1"  = function(s) s >= 0.03 & s < 0.1,
  "5e-4<=Shet<0.03" = function(s) s >= 5e-4 & s < 0.03,
  "Shet<5e-4"       = function(s) s < 5e-4
)

# Variant classification
variant_class_matrix <- function(inh) {
  splice_ai <- suppressWarnings(as.numeric(inh$SpliceAI_max))
  is_lof <- grepl("frameshift|stop_gained|splice_donor|splice_acceptor", inh$consequence) |
    (!is.na(splice_ai) & splice_ai >= 0.8)
  is_mis <- grepl("missense", inh$consequence) & !is_lof
  is_syn <- grepl("synonymous", inh$consequence) & !is_lof
  is_lof_hc <- is_lof & inh$LoF == "HC"

  # !is.na(x) guards go first so a tier predicate never returns NA for a missing score.
  misfit <- suppressWarnings(as.numeric(inh$MisFit1.5_S))
  mis_tiers <- vapply(MISFIT_TIERS, function(f) is_mis & !is.na(misfit) & f(misfit), logical(nrow(inh)))

  shet <- suppressWarnings(as.numeric(inh$Shet))
  lof_tiers <- vapply(SHET_TIERS, function(f) is_lof_hc & !is.na(shet) & f(shet), logical(nrow(inh)))

  cbind(Synonymous = is_syn, Missense = is_mis, mis_tiers, LoF = is_lof_hc, lof_tiers)
}

# Transmission test for one variant class: binom.test of transmitted vs. non-transmitted parental alleles against the null of 50% transmission.
tdt_test <- function(in_class, iGT) {
  n_trans    <- sum(in_class & iGT == 1)
  n_nontrans <- sum(in_class & iGT == 0)
  test <- binom.test(n_trans, n_trans + n_nontrans, p = 0.5)
  data.frame(
    n_trans = n_trans, n_nontrans = n_nontrans,
    ratio = n_trans / (n_trans + n_nontrans),
    LowCI = test$conf.int[1], UpCI = test$conf.int[2], P = test$p.value
  )
}
tdt_by_class <- function(class_mat, iGT, phenotype) {
  rows <- lapply(colnames(class_mat), function(cls) {
    cbind(VariantClass = cls, tdt_test(class_mat[, cls], iGT))
  })
  cbind(Phenotype = phenotype, do.call(rbind, rows))
}

inh <- read.delim(INH_FILE, stringsAsFactors = FALSE)
class_mat <- variant_class_matrix(inh)

tdt <- rbind(
  tdt_by_class(class_mat[inh$i_asd == 2, ], inh$iGT[inh$i_asd == 2], "Affected"),
  tdt_by_class(class_mat[inh$i_asd == 1, ], inh$iGT[inh$i_asd == 1], "Unaffected")
)
tdt$Significance <- ifelse(tdt$P < 0.01, "**", ifelse(tdt$P < 0.05, "*", ""))

print(tdt, row.names = FALSE)
