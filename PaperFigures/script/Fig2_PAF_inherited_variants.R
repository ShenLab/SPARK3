# Population attributable fraction (PAF) of de novo vs. inherited D-mis/LoF burden, overall and by sex / family history.
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(ggpattern); library(openxlsx)})

FILE    <- "Data/variant_rate_DNVs_CaseCtrl.xlsx"
OUT_DIR <- "PaperFigures/output"

# Sample strata: which case-control and DNV columns to pull, the plot's y-axis range, and the row label used in Supplementary Table 17, section (3).
GROUPS <- list(
  list(label = "All",       sheet_label = "All autistic", case_col = "Case",     ctrl_col = "UKBB_controls",
       aff_col = "aff",     unaff_col = "unaff",     y_lim = c(0, 27), y_breaks = seq(0, 25, 5)),
  list(label = "Male",      sheet_label = "Males",        case_col = "Male",     ctrl_col = "Male_controls",
       aff_col = "M_aff",   unaff_col = "M_unaff",   y_lim = c(0, 28), y_breaks = seq(0, 28, 5)),
  list(label = "Female",    sheet_label = "Females",      case_col = "Female",   ctrl_col = "Female_controls",
       aff_col = "F_aff",   unaff_col = "F_unaff",   y_lim = c(0, 28), y_breaks = seq(0, 28, 5)),
  list(label = "Simplex",   sheet_label = "Simplex",      case_col = "Simplex",  ctrl_col = "UKBB_controls",
       aff_col = "SP_aff",  unaff_col = "SP_unaff",  y_lim = c(0, 31), y_breaks = seq(0, 30, 5)),
  list(label = "Multiplex", sheet_label = "Multiplex",    case_col = "Multiplex", ctrl_col = "UKBB_controls",
       aff_col = "MP_aff",  unaff_col = "MP_unaff",  y_lim = c(0, 31), y_breaks = seq(0, 30, 5))
)

unfactor <- function(s) as.numeric(as.character(s))
to_paf   <- function(burden) 1 - exp(-1 * burden)

CI <- function(case_n, ctrl_n, case_rate, ctrl_rate) {
  delta <- case_rate - ctrl_rate
  se <- sqrt(case_rate / case_n + ctrl_rate / ctrl_n)
  c(delta - 1.96 * se, delta + 1.96 * se)
}
SE_rate_diff <- function(n1, n2, rate1, rate2) sqrt(rate1 / n1 + rate2 / n2)

# Case-control, de novo, and inherited rate differences (+ 95% CIs) underlying the PAF.
# Inherited burden = total case-control burden minus de novo burden; its CI combines the case-control and DNV standard errors in quadrature (independent samples), not by
# subtracting the two CIs directly.
paf_components <- function(case_n, ctrl_n, case_rate, ctrl_rate, aff_n, unaff_n, aff_rate, unaff_rate) {
  cc_d  <- case_rate - ctrl_rate
  cc_ci <- CI(case_n, ctrl_n, case_rate, ctrl_rate)
  se_cc <- SE_rate_diff(case_n, ctrl_n, case_rate, ctrl_rate)

  dn_d  <- aff_rate - unaff_rate
  dn_ci <- CI(aff_n, unaff_n, aff_rate, unaff_rate)
  se_dn <- SE_rate_diff(aff_n, unaff_n, aff_rate, unaff_rate)

  in_d  <- cc_d - dn_d
  se_in <- sqrt(se_cc^2 + se_dn^2)
  in_ci <- c(in_d - 1.96 * se_in, in_d + 1.96 * se_in)

  list(cc_d = cc_d, cc_ci = cc_ci, dn_d = dn_d, dn_ci = dn_ci, in_d = in_d, in_ci = in_ci)
}

# De novo and inherited PAF for one variant class (D-mis or LoF), with 95% CIs. CIs are mapped through PAF = 1 - exp(-Delta) at their endpoints (an exact, monotonic
# transform), rather than linearized with the delta method.
paf_data <- function(case_n, ctrl_n, case_rate, ctrl_rate, aff_n, unaff_n, aff_rate, unaff_rate, variant_type) {
  comp <- paf_components(case_n, ctrl_n, case_rate, ctrl_rate, aff_n, unaff_n, aff_rate, unaff_rate)
  rbind(
    data.frame(Burden = comp$dn_d, CI1 = comp$dn_ci[1], CI2 = comp$dn_ci[2],
               PAF = to_paf(comp$dn_d), PAF_CI1 = to_paf(comp$dn_ci[1]), PAF_CI2 = to_paf(comp$dn_ci[2]),
               Variant_type = variant_type, Mode = "DNVs"),
    data.frame(Burden = comp$in_d, CI1 = comp$in_ci[1], CI2 = comp$in_ci[2],
               PAF = to_paf(comp$in_d), PAF_CI1 = to_paf(comp$in_ci[1]), PAF_CI2 = to_paf(comp$in_ci[2]),
               Variant_type = variant_type, Mode = "Inherited")
  )
}

# A table summarize case-control / DNV / inherited excess and PAF for one variant class in one stratum.
paf_table_row <- function(case_n, ctrl_n, case_rate, ctrl_rate, aff_n, unaff_n, aff_rate, unaff_rate) {
  comp <- paf_components(case_n, ctrl_n, case_rate, ctrl_rate, aff_n, unaff_n, aff_rate, unaff_rate)
  data.frame(
    `Case N` = case_n, `Control N` = ctrl_n,
    `Rate in cases` = case_rate, `Rate in controls` = ctrl_rate,
    `Case-ctrl excess` = comp$cc_d, `Case-ctrl excess CI1` = comp$cc_ci[1], `Case-ctrl excess CI2` = comp$cc_ci[2],
    `DNV rate in autistic individuals` = aff_rate, `DNV rate in siblings without autism` = unaff_rate,
    `DNV excess` = comp$dn_d, `DNV excess CI1` = comp$dn_ci[1], `DNV excess CI2` = comp$dn_ci[2],
    `Inherited excess` = comp$in_d, `Inherited excess CI1` = comp$in_ci[1], `Inherited excess CI2` = comp$in_ci[2],
    `DNV PAF` = to_paf(comp$dn_d), `DNV PAF CI1` = to_paf(comp$dn_ci[1]), `DNV PAF CI2` = to_paf(comp$dn_ci[2]),
    `Inherited PAF` = to_paf(comp$in_d), `Inherited PAF CI1` = to_paf(comp$in_ci[1]), `Inherited PAF CI2` = to_paf(comp$in_ci[2]),
    check.names = FALSE
  )
}

paf_plot <- function(df, y_lim, y_breaks, text_size) {
  df <- df %>% mutate(
    Variant_type = factor(Variant_type, levels = c("LoF", "D-mis")),
    Mode = factor(Mode, levels = c("DNVs", "Inherited")),
    PAF_percent = PAF * 100
  )
  ggplot(df, aes(x = Gene_groups, y = PAF_percent, fill = Variant_type, pattern = Mode)) +
    geom_col_pattern(position = position_dodge(width = 0.9), color = "black", pattern_fill = "black",
                      pattern_angle = 45, pattern_density = 0.05, pattern_spacing = 0.05,
                      pattern_key_scale_factor = 0.6) +
    geom_errorbar(aes(ymin = PAF_CI1 * 100, ymax = PAF_CI2 * 100),
                  position = position_dodge(width = 0.9), width = 0.2, size = 1, color = "gray20") +
    scale_fill_manual(values = c("LoF" = "sienna2", "D-mis" = "mediumpurple4")) +
    scale_pattern_manual(values = c("DNVs" = "stripe", "Inherited" = "none")) +
    guides(fill = guide_legend(override.aes = list(pattern = "none")),
           pattern = guide_legend(override.aes = list(
             fill = c("gray90", "gray90"), pattern_fill = c("black", NA),
             pattern_density = c(0.05, 0), pattern = c("stripe", "none")))) +
    geom_text(aes(y = PAF_CI2 * 100, label = sprintf("%.1f%%", PAF_percent)),
              position = position_dodge(width = 0.9), vjust = -1, size = text_size) +
    geom_text(aes(y = PAF_CI2 * 100, label = "***"),
              position = position_dodge(width = 0.9), vjust = 0.5, size = 12) +
    scale_y_continuous(limits = y_lim, breaks = y_breaks, expand = expansion(mult = c(0.015, 0.05))) +
    labs(y = "Population attributable fraction (%)", fill = "Variant class", pattern = "Inheritance") +
    theme_classic(base_size = 14) +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.title.y = element_text(size = 25), axis.text.y = element_text(size = 25),
          legend.title = element_text(size = 25), legend.text = element_text(size = 20))
}

#' Build the DNV/Inherited x D-mis/LoF PAF table and summary plot for one stratum.
run_group <- function(cc1, dnv1, group) {
  case_n <- unfactor(cc1[[group$case_col]][1]);  ctrl_n <- unfactor(cc1[[group$ctrl_col]][1])
  case_Dmis <- unfactor(cc1[[group$case_col]][4]); ctrl_Dmis <- unfactor(cc1[[group$ctrl_col]][4])
  case_Dlof <- unfactor(cc1[[group$case_col]][5]); ctrl_Dlof <- unfactor(cc1[[group$ctrl_col]][5])

  aff_n <- unfactor(dnv1[[group$aff_col]][1]);     unaff_n <- unfactor(dnv1[[group$unaff_col]][1])
  aff_Dmis <- unfactor(dnv1[[group$aff_col]][4]);  unaff_Dmis <- unfactor(dnv1[[group$unaff_col]][4])
  aff_Dlof <- unfactor(dnv1[[group$aff_col]][5]);  unaff_Dlof <- unfactor(dnv1[[group$unaff_col]][5])

  df <- rbind(
    paf_data(case_n, ctrl_n, case_Dmis, ctrl_Dmis, aff_n, unaff_n, aff_Dmis, unaff_Dmis, "D-mis"),
    paf_data(case_n, ctrl_n, case_Dlof, ctrl_Dlof, aff_n, unaff_n, aff_Dlof, unaff_Dlof, "LoF")
  )
  df$Gene_groups <- "All"

  plot <- paf_plot(df, group$y_lim, group$y_breaks, text_size = 10)

  dn_d <- sum(df$Burden[df$Mode == "DNVs"]); in_d <- sum(df$Burden[df$Mode == "Inherited"])
  summary <- c(dn_paf = sprintf("%.1f%%", to_paf(dn_d) * 100),
               in_paf = sprintf("%.1f%%", to_paf(in_d) * 100))

  table <- rbind(
    cbind(Individuals = group$sheet_label, `Variant class` = "D-mis",
          paf_table_row(case_n, ctrl_n, case_Dmis, ctrl_Dmis, aff_n, unaff_n, aff_Dmis, unaff_Dmis)),
    cbind(Individuals = group$sheet_label, `Variant class` = "LoF",
          paf_table_row(case_n, ctrl_n, case_Dlof, ctrl_Dlof, aff_n, unaff_n, aff_Dlof, unaff_Dlof))
  )

  list(plot = plot, summary = summary, table = table)
}

cc  <- read.xlsx(FILE, sheet = 1); cc1 <- cc[c(1:3, 13, 24), 2:ncol(cc)]
colnames(cc1) <- colnames(cc)[2:ncol(cc)]

dnv <- read.xlsx(FILE, sheet = 2); dnv1 <- dnv[c(2:4, 14, 25), 2:ncol(dnv)]
colnames(dnv1) <- c("aff", "unaff", "delta", "SP_aff", "SP_unaff", "SP_delta", "MP_aff", "MP_unaff", "MP_delta", "M_aff", "M_unaff", "M_delta", "F_aff", "F_unaff", "F_delta")

results <- setNames(lapply(GROUPS, run_group, cc1 = cc1, dnv1 = dnv1), sapply(GROUPS, `[[`, "label"))
paf_summary <- sapply(results, `[[`, "summary")
print(paf_summary)

paf_table <- do.call(rbind, lapply(results, `[[`, "table"))
rownames(paf_table) <- NULL
print(paf_table)

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
pdf(file.path(OUT_DIR, "PAF_overall_sample_set.pdf"), width = 8, height = 10)
for (label in c("Simplex", "Multiplex", "Male", "Female")) print(results[[label]]$plot)
dev.off()
pdf(file.path(OUT_DIR, "PAF_overall.pdf"), width = 12, height = 10)
print(results[["All"]]$plot)
dev.off()
