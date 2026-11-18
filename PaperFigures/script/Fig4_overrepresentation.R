unfactor<-function(s){as.numeric(as.character(s))}
library(ggplot2);library(data.table);library(dplyr);library(minpack.lm);library(patchwork);library(epitools);library(openxlsx);library(ggrepel)
file<-'Data/BrainRegions_Cellular_analysis.xlsx'
region1<-read.xlsx(file,sheet=1);cell1<-read.xlsx(file,sheet=2)
colnames(region1)<-region1[1,];region1<-region1[-1,];colnames(cell1)<-cell1[1,];cell1<-cell1[-1,]
region2<-read.xlsx(file,sheet=3);cell2<-read.xlsx(file,sheet=4)
colnames(region2)<-region2[1,];region2<-region2[-1,];colnames(cell2)<-cell2[1,];cell2<-cell2[-1,]

pseudo_or<-function(overlap,term_genes,gene_set,total_genes){
  a<-overlap;b<-term_genes-a
  c<-gene_set-a;d<-total_genes-term_genes-gene_set+a
  cont_table <- matrix(c(a, b, c, d), nrow = 2)
  per_test<-fisher.test(cont_table)
  or<-per_test$estimate
  p<-per_test$p.value
  return(c(or,p))
}
enrichr_plot<-function(data,gene1,gene2,p1_thrd,p2_thrd,n_gene1,n_gene2,dot_range,dot_breaks,dot_limits){
  data<-data[!(data$Count1==0&data$Count2==0),];data<-data[!is.na(data$Category),]
  data$Count1 <- as.numeric(data$Count1);data$Count2 <- as.numeric(data$Count2);data$Fisher_P <- as.numeric(data$Fisher_P)
  data$Odds_Ratio1 <- as.numeric(data$Odds_Ratio1);data$Odds_Ratio2 <- as.numeric(data$Odds_Ratio2)
  data$Enrich_P1 <- as.numeric(data$Enrich_P1);data$Enrich_P2 <- as.numeric(data$Enrich_P2)
  for(i in 1:nrow(data)){
    if(data$Count1[i]==0){term_genes<-unfactor(c(unlist(strsplit(data$Overlap2[i],'/'))))[2];new<-pseudo_or(1,term_genes+1,n_gene1+1,19000);data$Odds_Ratio1[i]<-new[1];data$Enrich_P1[i]<-new[2]}
    if(data$Count2[i]==0){term_genes<-unfactor(c(unlist(strsplit(data$Overlap1[i],'/'))))[2];new<-pseudo_or(1,term_genes+1,n_gene2+1,19000);data$Odds_Ratio2[i]<-new[1];data$Enrich_P2[i]<-new[2]}
  }
  data$log10_OR1 <- log10(data$Odds_Ratio1);data$log10_OR2 <- log10(data$Odds_Ratio2);data$logFisher <- -log10(data$Fisher_P)
  min_limit <- floor(min(c(data$log10_OR1, data$log10_OR2), na.rm = TRUE));max_limit <- ceiling(max(c(data$log10_OR1, data$log10_OR2), na.rm = TRUE))
  # terms meeting this now get a bold black point border instead of a text label. NA-safe: a plain
  # `Enrich_P1 < p1_thrd | Enrich_P2 < p2_thrd` yields NA (not FALSE) whenever one p-value is NA and
  # the other doesn't clear its threshold, and scale_color_manual below has no mapping for NA, so
  # those rows were silently dropped from the plot entirely instead of just drawn without a border.
  data$sig <- (!is.na(data$Enrich_P1) & data$Enrich_P1 < p1_thrd) | (!is.na(data$Enrich_P2) & data$Enrich_P2 < p2_thrd)
  
  if(length(unique(data$Category))>5){
    # 13 Allen Brain Atlas Categorys (finer than the 7-way Category grouping); region plots use
    # the gene's own Category column (added to sheets 1 & 4 of the xlsx) instead of Category.
    category_colors <- c("Isocortex" = "#D62728", "Olfactory areas" = "#9467BD", "Hippocampus" = "#FF9F40",
                         "Cortical subplate" = "#17BECF", "Striatum" = "#1F77B4", "Pallidum" = "#AEC7E8",
                         "Amygdala" = "#E377C2", "Thalamus" = "#2CA02C", "Hypothalamus" = "#98DF8A",
                         "Midbrain" = "#8C564B", "Pons" = "#C49C94", "Medulla" = "#7F7F7F", "Cerebellum" = "#E8C9A0")
    legend_title<-'Brain region'
    # only rows with a Category are plottable (Category was only assigned to significant,
    # already-Category-labelled rows, same as Category itself)
    data <- data[!is.na(data$Category), ]
    data$PlotGroup <- data$Category
    # draw order (z-order): geom_point renders rows in the order they appear in `data`, so later
    # rows sit on top when points overlap. Striatum (blue, #1F77B4) and Isocortex (red, #D62728)
    # are the direct successors of the old Subcortical/Cortical colors, so keep the same priority:
    # every other Category first, Striatum next, Isocortex last.
    draw_order <- c(setdiff(unique(data$PlotGroup), c('Striatum','Isocortex')), 'Striatum', 'Isocortex')
    data <- data[order(match(data$PlotGroup, draw_order)), ]
    # border colour: black for significant terms; for the rest, SAME colour as the dot's own fill
    # (not "transparent"/stroke=0) -- PDF/PostScript treats a 0-width line as a "hairline" that's
    # always drawn at ~1px regardless of width, so a truly invisible border needs same-colour-as-
    # fill, not zero width, to avoid a visible ring surviving into the saved PDF.
    data$border_col <- ifelse(data$sig, 'black', unname(category_colors[as.character(data$PlotGroup)]))
    x_breaks <- c(0.1, 0.5, 1, 2, 5, 10); x_limits <- NULL  # data-driven, as before
    y_breaks <- x_breaks; y_limits <- x_limits
  }
  if(length(unique(data$Category))<=5){
    category_colors <- c("Postsynaptic density-related" = "#1f77b4","Nucleus" = "#d62728","Others" = "#2ca02c")
    legend_title<-'Cellular component'
    data$Category<-factor(data$Category,levels=c('Nucleus','Postsynaptic density-related','Others'))
    data$PlotGroup <- data$Category
    # draw order (z-order): Others first, then Postsynaptic density-related (blue), then Nucleus
    # (red) last, so blue sits above Others and red sits above both.
    draw_order <- c('Others', 'Postsynaptic density-related', 'Nucleus')
    data <- data[order(match(data$PlotGroup, draw_order)), ]
    data$border_col <- ifelse(data$sig, 'black', unname(category_colors[as.character(data$PlotGroup)]))
    x_breaks <- c(0.5, 1, 2, 5, 10, 15, 20, 25); x_limits <- c(0.5, 25)
    y_breaks <- x_breaks; y_limits <- x_limits
  }
  
  p <- ggplot(data, aes(x = Odds_Ratio1, y = Odds_Ratio2)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50", linewidth = 0.6) +
    geom_vline(xintercept = 1, color = "black", linewidth = 0.6) +
    geom_hline(yintercept = 1, color = "black", linewidth = 0.6) +
    geom_point(aes(size = Count1 + Count2, fill = PlotGroup, color = I(border_col)),
               shape = 21, stroke = 1.3,
               position = position_jitter(width = 0.01, height = 0.01), alpha = 0.85) +
    coord_equal(
      xlim = if (is.null(x_limits)) c(min(data$Odds_Ratio1, data$Odds_Ratio2, na.rm = TRUE),
                                      max(data$Odds_Ratio1, data$Odds_Ratio2, na.rm = TRUE)) else x_limits,
      ylim = if (is.null(y_limits)) c(min(data$Odds_Ratio1, data$Odds_Ratio2, na.rm = TRUE),
                                      max(data$Odds_Ratio1, data$Odds_Ratio2, na.rm = TRUE)) else y_limits
    ) +
    scale_x_log10(
      name = bquote("Odds ratio of "*.(gene1)),
      breaks = x_breaks,
      labels = x_breaks,
      limits = x_limits
    ) +
    scale_y_log10(
      name = bquote("Odds ratio of "*.(gene2)),
      breaks = y_breaks,
      labels = y_breaks,
      limits = y_limits
    ) +
    scale_fill_manual(
      values = category_colors,
      name = legend_title
    ) +
    scale_size_continuous(name = "Number of genes") +
    guides(
      # size = 6 makes legend dots larger; color = NA removes the border from the legend swatches
      # specifically (legend keys otherwise pick up whichever row's border_col ggplot samples) --
      # this only affects the legend glyphs, not the actual plotted points
      fill = guide_legend(override.aes = list(size = 6, color = NA))
    ) +
    theme_classic(base_size = 20) +
    scale_size_continuous(
      name = "Number of genes",
      range = dot_range, 
      breaks = dot_breaks,
      limits = dot_limits
    )+
    theme(
      axis.line = element_blank(),
      axis.ticks = element_line(size = 0.6),
      axis.text = element_text(size = 18),
      axis.title = element_text(size = 22, face = "bold"),
      legend.title = element_text(size = 18),
      legend.text = element_text(size = 16),
      plot.title = element_blank(),
      legend.position = "right"
    )
  return(p)
}

novel_known_region<-enrichr_plot(region1,'newly identified candidate genes','known genes',5e-5,5e-5,137,97,c(2,6),c(5, 10, 15),c(1,max(unfactor(region1$Count1)+unfactor(region1$Count2))));
novel_known_cell<-enrichr_plot(cell1,'newly identified candidate genes','known genes',1.6e-4,1.6e-4,137,97,c(2,6),c(25,50,75,100),c(1,max(unfactor(cell1$Count1)+unfactor(cell1$Count2))))
mis_lgd_region<-enrichr_plot(region2,'missense genes','LoF genes',5e-5,5e-5,100,94,c(3,6),c(1, 5, 10),c(1,max(unfactor(region2$Count1)+unfactor(region2$Count2))));
mis_lgd_cell<-enrichr_plot(cell2,'missense genes','LoF genes',1.6e-4,1.6e-4,100,94,c(2,6),c(25,50,75,100),c(1,max(unfactor(cell2$Count1)+unfactor(cell2$Count2))))
pdf('PaperFigures/output/Figure4.pdf',width = 12,height = 10)
novel_known_region
novel_known_cell
mis_lgd_region;mis_lgd_cell
dev.off()

