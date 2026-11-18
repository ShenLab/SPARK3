# Fold enrichment of de novo variants (DNVs) in autistic probands vs. unaffected siblings.

DNV_FILE   <- "Data/autism_DNVs.txt"
TRIO_FILE  <- "Data/autism_trio_samples.txt"

# MisFit-S cutoffs defining nested damaging-missense (D-mis) tiers.
MISFIT_THRESHOLDS <- c(0.001, 0.003, 0.005, 0.01, 0.03, 0.1)

# Sample strata, each defined by a predicate on the trio-offspring table.
# `%in%` rather than `==` throughout: Sex and FamHist carry missing values, and `==` would yield NA, which silently injects NA rows when used to subset.
STRATA <- list(
  All       = function(x) rep(TRUE, nrow(x)),
  Female    = function(x) x$Sex %in% "Female",
  Male      = function(x) x$Sex %in% "Male",
  Simplex   = function(x) x$FamHist %in% "Simplex",
  Multiplex = function(x) x$FamHist %in% "Multiplex"
)

variant_class_matrix <- function(dnv, thresholds = MISFIT_THRESHOLDS) {
  is_mis <- dnv$cq %in% "mis"
  misfit <- suppressWarnings(as.numeric(dnv$`MisFit1.5_S`))

  dmis <- vapply(
    thresholds,
    function(cutoff) is_mis & !is.na(misfit) & misfit >= cutoff,
    logical(nrow(dnv))
  )
  colnames(dmis) <- sprintf("Dmis_MisFit_S%g", thresholds)

  cbind(Syn = dnv$cq %in% "syn", Mis = is_mis, dmis, LGD = dnv$cq %in% "LGD")
}

count_by_stratum <- function(dnv, ped, class_mat, strata = STRATA) {
  variants <- t(vapply(
    strata,
    function(keep) colSums(class_mat[keep(dnv), , drop = FALSE]),
    numeric(ncol(class_mat))
  ))
  trios <- vapply(strata, function(keep) sum(keep(ped)), numeric(1))
  list(variants = variants, trios = trios)
}

# ---------------------------------------------------------------------------
# Enrichment test
# ---------------------------------------------------------------------------
poisson_enrichment <- function(affected, unaffected) {
  strata  <- rownames(affected$variants)
  classes <- colnames(affected$variants)

  rows <- lapply(strata, function(stratum) {
    n_aff   <- affected$trios[[stratum]]
    n_unaff <- unaffected$trios[[stratum]]

    do.call(rbind, lapply(classes, function(cls) {
      k_aff   <- as.integer(affected$variants[stratum, cls])
      k_unaff <- as.integer(unaffected$variants[stratum, cls])
      test    <- poisson.test(c(k_aff, k_unaff), c(n_aff, n_unaff))

      data.frame(
        Stratum        = stratum,
        VariantClass   = cls,
        n_dnv_aff      = k_aff,
        n_dnv_unaff    = k_unaff,
        n_trios_aff    = n_aff,
        n_trios_unaff  = n_unaff,
        rate_aff       = k_aff / n_aff,
        rate_unaff     = k_unaff / n_unaff,
        RateDiff       = k_aff / n_aff - k_unaff / n_unaff,
        FoldEnrichment = unname(test$estimate),
        LowCI          = test$conf.int[1],
        UpCI           = test$conf.int[2],
        P              = test$p.value,
        stringsAsFactors = FALSE
      )
    }))
  })

  do.call(rbind, rows)
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
ped <- read.delim(TRIO_FILE, stringsAsFactors = FALSE)
dnv <- read.delim(DNV_FILE,  stringsAsFactors = FALSE)

# Sex and family history live on the sample table; carry them onto the variants.
dnv <- merge(dnv, ped[, c("IID", "Sex", "FamHist")], by = "IID", all.x = TRUE)
class_mat <- variant_class_matrix(dnv)
affected <- count_by_stratum(
  dnv[dnv$Pheno %in% "Affected", ],
  ped[ped$Pheno %in% "Affected", ],
  class_mat[dnv$Pheno %in% "Affected", , drop = FALSE]
)
unaffected <- count_by_stratum(
  dnv[dnv$Pheno %in% "Unaffected", ],
  ped[ped$Pheno %in% "Unaffected", ],
  class_mat[dnv$Pheno %in% "Unaffected", , drop = FALSE]
)

burden <- poisson_enrichment(affected, unaffected)

print(burden[burden$Stratum == "All", ], row.names = FALSE)
