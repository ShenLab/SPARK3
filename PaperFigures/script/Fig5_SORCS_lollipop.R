unfactor<-function(s){as.numeric(as.character(s))}
library(Gviz);library(rtracklayer);library(trackViewer);library(RColorBrewer);library(readr);library(GenomicFeatures);library(openxlsx)
library(ggplot2);library(dplyr);library(gggenes);library(GenomicRanges);library(gridExtra);library(tidyr);library(stringr);library(data.table)

# Canonical protein length (aa), from UniProt's REST API rather than Ensembl's BioMart: BioMart's
# live web service (useEnsembl()/getBM()) is frequently unreachable (main site and all mirrors),
# and this is the only thing it was used for here -- UniProt's reviewed (Swiss-Prot) entry for
# each gene gives the same canonical sequence length directly.
get_protein_length <- function(gene) {
  r <- httr::GET("https://rest.uniprot.org/uniprotkb/search",
                 query = list(query = paste0("gene:", gene, " AND organism_id:9606 AND reviewed:true"),
                              fields = "accession,length", format = "json"))
  res <- jsonlite::fromJSON(httr::content(r, as = "text", encoding = "UTF-8"))$results
  res$sequence$length[1]
}
protein_lengths <- sapply(c("SORCS1","SORCS2","SORCS3"), get_protein_length)
print(protein_lengths)
reformat_variant1 <- function(variant, GeneEff) {
  if (GeneEff %in% c("frameshift", "frameshift_variant")) {
    match <- regmatches(variant, regexec("^(\\d+)([A-Z]+).*>\\d+\\*?", variant))
    if (length(match[[1]]) > 2) {pos <- match[[1]][2];ref <- substr(match[[1]][3], 1, 1);variant <- paste0(ref, pos, "fs*")}} else {
      variant <- gsub("(\\d+)([A-Z])>\\d+([A-Z])", "\\2\\1\\3", variant)
      variant <- gsub("(\\d+)([A-Z])>\\d+\\*", "\\2\\1*", variant)
    }
  return(variant)
}
reformat_variant2 <- function(AAchange) {
  sub(".*p\\.", "", AAchange)
}

dnv_aff<-read.csv('Data/SORCS_genes_dnv_aff.txt',header = T,sep = '\t');
inh_trans<-read.csv('Data/SORCS_genes_inherited.txt',header = T,sep='\t')
case_var<-read.csv('Data/SORCS_genes_case_var.txt',header = T,sep='\t')
all_ctrl<-read.csv('Data/all_SORCS_ctrl_var.txt',header = T,sep='\t')
all_case<-read.csv('Data/all_SORCS_case_var.txt',header = T,sep='\t')
case_n<-33319;ctrl_n<-142251
domains_file<-'Data/SORCS_domain_regions.xlsx'
all_domains<-read.xlsx(domains_file, sheet=1);plot_domains<-read.xlsx(domains_file,sheet = 1)

pdf('PaperFigures/output/Figure5.pdf',width = 25,height = 8)
for(gene in c('SORCS1','SORCS2','SORCS3')){
  protein_length<-protein_lengths[[gene]]
  domains<-all_domains[all_domains$gene==gene,]
  
  dnv_aff1<-dnv_aff[dnv_aff$HGNC==gene,]
  dnv_aff2<-dnv_aff1[,c('HGNC','canon_amino_acid_change','GeneEff','cq')];colnames(dnv_aff2)<-c('gene','AAchange','GeneEff','cq')
  if(nrow(dnv_aff2)>0){dnv_aff2$aachange<-do.call(c,lapply(1:nrow(dnv_aff2),function(s){reformat_variant2(as.character(dnv_aff2$AAchange[s]))}));dnv_aff2$position<-as.integer((sub("^[A-Za-z]+([0-9]+).*$", "\\1",dnv_aff2$aachange)))}
  inh_trans1<-inh_trans[inh_trans$gene==gene,]
  inh_trans2<-inh_trans1[,c('gene','canon_amino_acid_change','consequence','cq')];colnames(inh_trans2)<-c('gene','AAchange','GeneEff','cq')
  if(nrow(inh_trans2)>0){inh_trans2$aachange<-do.call(c,lapply(1:nrow(inh_trans2),function(s){reformat_variant1(as.character(inh_trans2[s,'AAchange']),inh_trans2[s,'GeneEff'])}));inh_trans2$position<- as.numeric(gsub("^([0-9]+).*", "\\1", inh_trans2$AAchange))}
  case_var1<-case_var[case_var$canon_gene==gene,];if(gene=='SORCS3'){case_var1<-case_var1[case_var1$POS!='105252822',]}
  case_var2<-case_var1[,c('canon_gene','canon_amino_acid_change','canon_Consequence','cq')];colnames(case_var2)<-c('gene','AAchange','GeneEff','cq')
  if(nrow(case_var2)>0){case_var2$aachange<-do.call(c,lapply(1:nrow(case_var2),function(s){reformat_variant1(as.character(case_var2[s,'AAchange']),case_var2[s,'GeneEff'])}));case_var2$position<- as.numeric(gsub("^([0-9]+).*", "\\1", case_var2$AAchange))}
  datasets <- list(dnv_aff2 = dnv_aff2, inh_trans2 = inh_trans2, case_var2 = case_var2)
  datasets <- lapply(names(datasets), function(name) {df <- datasets[[name]];if (nrow(df) > 0) {df$dataset <- name;return(df)} else {return(NULL)}})
  all_variants <- do.call(rbind, Filter(Negate(is.null), datasets))
  
  variant_counts <- aggregate(count ~ position + dataset + cq, data = transform(all_variants, count = 1), FUN = length)
  dataset_colors <- c("dnv_aff2" = "#E41A1C","inh_trans2" = "#377EB8","case_var2" = "#4DAF4A")
  variant_counts <- transform(variant_counts,color = dataset_colors[dataset],shape = ifelse(cq == "mis", "circle", "diamond"))
  variant_counts <- variant_counts %>%mutate(pos = ifelse(dataset %in% c("dnv_aff2", "inh_trans2", "case_var2"), "top", "bottom"))
  z<-variant_counts
  
  sample.gr <- GRanges(gene, IRanges(unfactor(z$position), width=1, names=NULL))
  sample.gr$shape<-z$shape;sample.gr$color<-z$color
  sample.gr$legend<-z$cq
  sample.gr$border<-'gray30';label.parameter.gp.black <- gpar(col="white")     #label.parameter.gp.black <- gpar(col="gray30")
  sample.gr$label.parameter.rot<-60
  sample.gr$score<-unfactor(z$count)
  sample.gr$SNPsideID<-z$pos
  
  features<-GRanges(gene,ranges = IRanges(start = domains$Domain_start, end = domains$Domain_end,names=domains$Domain_name),domain = domains$Domain_name)
  unique_domains <- unique(features$domain)
  domain_colors <- setNames(c('#FFFEE0', '#F4E9FF', '#D5F6FB', '#F1BEB5'), unique_domains);features$fill <- domain_colors[features$domain]
  
  features$height<-rep(0.05,length(features$fill))
  xaxis<-c(seq(0,floor(protein_length/100) * 100,100),protein_length);names(xaxis)<-xaxis;names(xaxis)[length(xaxis)]<-paste0(protein_length,"aa")
  legend<-list(labels=c("DNVs_affected","Transmitted","Variant_cases"),col=c("#E41A1C","#377EB8","#4DAF4A"),fill=c("#E41A1C","#377EB8","#4DAF4A"))
  lolliplot(sample.gr, features,ranges=GRanges(gene, IRanges(0, protein_length)),yaxis=c(0,max(z$count)),xaxis=xaxis,legend=legend,ylab=paste0("# ",gene," mutations"))
}
dev.off()

