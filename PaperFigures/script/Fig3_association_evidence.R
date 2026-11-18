library(ggplot2);library(dplyr);library(tidyr);library(patchwork);library(ggrepel)
unfactor<-function(s){as.numeric(as.character(s))}
y<-read.csv('Data/FDR0.01_bbTADA_genes.txt',header = T,sep='\t')
y<-y[order(y$BF,decreasing = T),]
y$BF_dnv <- y$BF_dn_lof * y$BF_dn_mis1 * y$BF_dn_mis2 * y$BF_dn_mis3
y$all_BF<-y$BF_lof*y$BF_mis1*y$BF_mis2*y$BF_mis3
y$perc_dn<-log10(y$BF_dnv)/log10(y$all_BF);y[which(y$perc_dn>1),'perc_dn']<-1

know<-y[(y$gene.score=='1'|y$gene.score=='1S')|y$ASD_studies!='.',];unknow<-y[y$gene.score!='1' & y$gene.score!='1S' & y$ASD_studies=='.',]
know1<-know[,c('gene','BF_lof','BF_mis1','BF_mis2','BF_mis3','perc_dn')];unknow1<-unknow[,c('gene','BF_lof','BF_mis1','BF_mis2','BF_mis3','perc_dn')]

###Separeate known and novel genes into scatter plot
know2<-cbind(know1[,c('gene','BF_lof')],know1$BF_mis1*know1$BF_mis2*know1$BF_mis3,know1$perc_dn);colnames(know2)<-c('gene','BF_lof','BF_mis',"ratio")
unknow2<-cbind(unknow1[,c('gene','BF_lof')],unknow1$BF_mis1*unknow1$BF_mis2*unknow1$BF_mis3,unknow1$perc_dn);colnames(unknow2)<-c('gene','BF_lof','BF_mis',"ratio")
pdf('PaperFigures/output/Figure3.pdf',width = 10,height = 7)
ggplot(know2, aes(x = log10(BF_lof), y = log10(BF_mis))) +geom_abline(slope = 1, intercept = 0, color = "gray60", linetype = "dashed") +geom_point(aes(fill = ratio), shape = 21, color = "black", size = 4) +
  geom_text_repel(aes(label = ifelse(log10(BF_lof) > 15 | log10(BF_mis) > 8, gene, "")),size = 4.5,max.overlaps = Inf,segment.color = "gray80",segment.linetype = "dashed",segment.size = 0.3) +
  scale_fill_gradientn(colors = c("blue", "yellow", "red"),name = "Proportion of DNV BF in total BF") +
  scale_x_continuous(name = "log10(BF of LoF)",breaks = c(0, 5, 10, 15, 20, 30, 40, 50),limits = c(0, max(log10(know2$BF_lof))),expand = c(0.05, 0.05)) +scale_y_continuous(name = "log10(BF of D-mis)",breaks = c(0, 5, 10, 15, 20, 30),limits = c(0, max(log10(know2$BF_mis))),expand = c(0.05, 0.05)) +
  theme_minimal() +theme(axis.text = element_text(size = 12),axis.title = element_text(size = 14),panel.grid.major = element_line(color = "gray90"),panel.grid.minor = element_blank())
ggplot(unknow2, aes(x = log10(BF_lof), y = log10(BF_mis))) +geom_abline(slope = 1, intercept = 0, color = "gray60", linetype = "dashed") +geom_point(aes(fill = ratio), shape = 21, color = "black", size = 4) +
  geom_text_repel(aes(label = ifelse(log10(BF_lof) > 5 | log10(BF_mis) > 5, gene, "")),size = 4.5,max.overlaps = Inf,segment.color = "gray80",segment.linetype = "dashed",segment.size = 0.3) +
  scale_fill_gradientn(colors = c("blue", "yellow", "red"),name = "Proportion of DNV BF in total BF") +
  scale_x_continuous(name = "log10(BF of LoF)",breaks = c(0, 5, 10, 15, 20, 30, 40, 50),limits = c(0, max(log10(unknow2$BF_lof))),expand = c(0.05, 0.05)) +scale_y_continuous(name = "log10(BF of D-mis)",breaks = c(0, 5, 10, 15, 20, 30),limits = c(0, max(log10(unknow2$BF_mis))),expand = c(0.05, 0.05)) +
  theme_minimal() +theme(axis.text = element_text(size = 12),axis.title = element_text(size = 14),panel.grid.major = element_line(color = "gray90"),panel.grid.minor = element_blank())
dev.off()


