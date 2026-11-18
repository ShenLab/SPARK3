# Overrepresentation of Allen Brain Atlas cell-type and GO Cellular Component terms among LoF- vs missense-predominant, and known vs novel, ASD risk genes.
#
# For each term in a library, compares how many genes from each group overlap that term via a 2x2 Fisher's exact test (group1 vs group2, not against genome background), then cross-references
# each group's own Enrichr overrepresentation result (background = all genes in the library) for the same term.

suppressPackageStartupMessages(library(httr))

GENE_TABLE   <- "Data/FDR0.01_bbTADA_genes.txt"
OUT_DIR      <- "Data"
LIBRARIES    <- c("Allen_Brain_Atlas_up", "GO_Cellular_Component_2026")
ENRICHR_BASE <- "https://maayanlab.cloud/Enrichr"

# ---------------------------------------------------------------------------
# Gene groups
# ---------------------------------------------------------------------------
y <- read.csv(GENE_TABLE, header = TRUE, sep = "\t")
y$BF_mis   <- y$BF_mis1 * y$BF_mis2 * y$BF_mis3
y$all_BF   <- y$BF_lof * y$BF_mis
y$perc_lof <- log10(y$BF_lof) / log10(y$all_BF)
y$perc_mis <- log10(y$BF_mis) / log10(y$all_BF)

groups <- list(
  LoF      = y$gene[y$perc_lof >= 0.7 & y$gene != "SORCS3"],
  Missense = union(y$gene[y$perc_mis >= 0.7 | y$gene == "SORCS3"], "SORCS2"),
  Known    = y$gene[(y$gene.score %in% c("1", "1S")) | y$ASD_studies != "."],
  Novel    = union(y$gene[!(y$gene.score %in% c("1", "1S")) & y$ASD_studies == "."], "SORCS2")
)

# ---------------------------------------------------------------------------
# Enrichr API (https://maayanlab.cloud/Enrichr/help#api)
# ---------------------------------------------------------------------------

# Download a full gene-set library as a named list: term -> character vector of gene symbols. Allen Brain Atlas-style libraries encode a per-gene score as "GENE,score"; stripping at the
# first comma works for both that format and libraries with plain gene symbols.
download_library <- function(library_name) {
  lines <- readLines(paste0(ENRICHR_BASE, "/geneSetLibrary?mode=text&libraryName=", library_name))
  split_lines <- strsplit(lines, "\t", fixed = TRUE)
  terms <- vapply(split_lines, `[`, character(1), 1)
  genes <- lapply(split_lines, function(x) sub(",.*", "", x[-(1:2)]))
  setNames(genes, terms)
}

# Submit a gene list to Enrichr and return its overrepresentation table for one library (background = all genes in that library, i.e. the standard Enrichr enrichment result).
enrichr_enrich <- function(genes, description, library_name) {
  added <- POST(paste0(ENRICHR_BASE, "/addList"),
                body = list(list = paste(genes, collapse = "\n"), description = description),
                encode = "multipart")
  user_list_id <- content(added, as = "parsed", type = "application/json")$userListId
  stopifnot(!is.null(user_list_id))
  Sys.sleep(1)  # be polite to a free public API across the repeated calls below
  exported <- GET(paste0(ENRICHR_BASE, "/export"),
                   query = list(userListId = user_list_id, filename = description, backgroundType = library_name))
  read.delim(text = content(exported, as = "text", encoding = "UTF-8"))
}

# Attach one group's Enrichr result columns to a term table, prefixed by `label` so the merged table stays self-documenting (e.g. Novel_P.value, Known_P.value).
add_enrichr_columns <- function(term_table, genes, label, library_name) {
  res <- enrichr_enrich(genes, label, library_name)
  non_term <- colnames(res) != "Term"
  colnames(res)[non_term] <- paste0(label, "_", colnames(res)[non_term])
  merge(term_table, res, by = "Term", all.x = TRUE)
}

# ---------------------------------------------------------------------------
# Two-group Fisher's exact test, term by term
# ---------------------------------------------------------------------------
#' For every term in `library`, a 2x2 Fisher's exact test of gene1 vs gene2 membership in that
#' term (group1 vs group2, not against genome background).
fisher_overrep <- function(library, genes1, genes2) {
  do.call(rbind, lapply(names(library), function(term) {
    term_genes <- library[[term]]
    match1 <- intersect(genes1, term_genes); match2 <- intersect(genes2, term_genes)
    ct1 <- length(match1); ct2 <- length(match2)
    p <- fisher.test(matrix(c(ct1, length(genes1) - ct1, ct2, length(genes2) - ct2), nrow = 2, byrow = TRUE))$p.value
    data.frame(Term = term, ct1 = ct1, frac1 = ct1 / length(genes1), ct2 = ct2, frac2 = ct2 / length(genes2),
               pval = p, gene1 = paste(match1, collapse = ","), gene2 = paste(match2, collapse = ","),
               stringsAsFactors = FALSE)
  }))
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

for (library_name in LIBRARIES) {
  lib <- download_library(library_name)

  novel_vs_known <- fisher_overrep(lib, groups$Novel, groups$Known)
  novel_vs_known <- add_enrichr_columns(novel_vs_known, groups$Novel, "Novel", library_name)
  novel_vs_known <- add_enrichr_columns(novel_vs_known, groups$Known, "Known", library_name)

  mis_vs_lof <- fisher_overrep(lib, groups$Missense, groups$LoF)
  mis_vs_lof <- add_enrichr_columns(mis_vs_lof, groups$Missense, "Missense", library_name)
  mis_vs_lof <- add_enrichr_columns(mis_vs_lof, groups$LoF, "LoF", library_name)

  tag <- if (library_name == "Allen_Brain_Atlas_up") "cells_up" else "go"
  write.table(novel_vs_known, file.path(OUT_DIR, paste0("novel_known_", tag, ".txt")), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(mis_vs_lof,     file.path(OUT_DIR, paste0("mis_lof_", tag, ".txt")),     sep = "\t", quote = FALSE, row.names = FALSE)
}
