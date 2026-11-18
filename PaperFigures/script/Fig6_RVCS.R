###Figure 6a: overall distribution: only EUR in affected vs unaffected siblings vs Population controls
library(ggplot2);library(dplyr);library(tidyr);library(patchwork);library(ggrepel);library(openxlsx);library(RColorBrewer);library(ggpattern);library(Cairo)
unfactor<-function(s){as.numeric(as.character(s))}
aff<-read.csv('Data/RVCS_EUR_autistic_individuals.txt',header = T,sep='\t')
unaff<-read.csv('Data/RVCS_EUR_siblings.txt',header = T,sep='\t')
ctrl<-read.csv('Data/RVCS_EUR_controls.txt',header = T,sep='\t')

df<-rbind(aff,unaff,ctrl)
df$group <- factor(df$group, levels = c("Autistic individuals", "Siblings without autism", "Population controls"))
colors <- c("Autistic individuals" = "#E63032","Siblings without autism" = "#377EB8","Population controls" = "#4DAF4A")
group_means <- df %>%group_by(group) %>%summarise(mean_S = mean(S));known_group_mean<-group_means
all_overall <- ggplot(df, aes(x = S, color = group)) +geom_density(size = 1.2) +scale_x_log10() + scale_color_manual(values = colors) +labs(x = "Rare variant composite score", y = "Density", title = NULL) +theme_minimal() +
  theme(text = element_text(size = 20),legend.title = element_blank(),legend.position = 'right',legend.background = element_rect(fill = "white", color = NA),legend.text = element_text(size = 20),axis.text.y = element_text(size =15),axis.title.y = element_text(size = 20)) +
  geom_vline(data = group_means, aes(xintercept = mean_S, color = group), linetype = "dashed", size = 1)
print(all_overall)
CairoPDF("PaperFigures/output/Fig6a_RVCS_distribution.pdf", width = 12, height = 3)
print(all_overall)
dev.off()

###Figure 6b: association between RVCS and clinical phenotypes
library(ggplot2);library(dplyr);library(tidyr);library(patchwork);library(ggrepel);library(openxlsx);library(RColorBrewer)
unfactor<-function(s){as.numeric(as.character(s))}
x<-read.csv('Data/RVCS_all_autistic_individuals.txt',header = T,sep='\t')
imp<-read.csv('Data/Cognitive_impairment.txt',header = T,sep='\t')
basic<-read.csv('Data/Basic_screen.txt',header = T,sep='\t')

##binary derived cognitive impairment, ADHD>5yr, Seizures:any age, depression>10yr, bipolar>15yr, schizophrenia>15yr,anxiety>10 
decile<-function(phe,phe_name,phe_label,geneset){
  group_names<-10:1
  if(phe_label=='Cognitive impairment'){
    phe<-phe[!is.na(phe[,phe_name]),];thrd <- rev(quantile(phe$S, probs = seq(0, 1, length.out = 11), na.rm = TRUE)[1:10]);thrd8<-thrd[8]
    y<-phe;j<-phe_name
    middle_y<-y[y$S<thrd8,];n_y_y<-nrow(middle_y[middle_y[,j]==1,]);n_y_n<-nrow(middle_y[middle_y[,j]==0,])
    y1<-y[y$S>=thrd[1],];n_y1_y<-nrow(y1[y1[,j]==1,]);n_y1_n<-nrow(y1[y1[,j]==0,]);y1_n<-rbind(c(n_y1_y,n_y1_n),c(n_y_y,n_y_n));y1_p<-fisher.test(y1_n,alternative = 'two.sided')$p.value;rr1<-(y1_n[1,1]/y1_n[1,2])/(y1_n[2,1]/y1_n[2,2]);y1_ci<-fisher.test(y1_n)$conf.int[1:2];data1<-c(rr1,y1_ci,y1_p)
    y2<-y[y$S>=thrd[2] & y$S<thrd[1],];n_y2_y<-nrow(y2[y2[,j]==1,]);n_y2_n<-nrow(y2[y2[,j]==0,]);y2_n<-rbind(c(n_y2_y,n_y2_n),c(n_y_y,n_y_n));y2_p<-fisher.test(y2_n,alternative = 'two.sided')$p.value;rr2<-(y2_n[1,1]/y2_n[1,2])/(y2_n[2,1]/y2_n[2,2]);y2_ci<-fisher.test(y2_n)$conf.int[1:2];data2<-c(rr2,y2_ci,y2_p)
    y3<-y[y$S>=thrd[3] & y$S<thrd[2],];n_y3_y<-nrow(y3[y3[,j]==1,]);n_y3_n<-nrow(y3[y3[,j]==0,]);y3_n<-rbind(c(n_y3_y,n_y3_n),c(n_y_y,n_y_n));y3_p<-fisher.test(y3_n,alternative = 'two.sided')$p.value;rr3<-(y3_n[1,1]/y3_n[1,2])/(y3_n[2,1]/y3_n[2,2]);y3_ci<-fisher.test(y3_n)$conf.int[1:2];data3<-c(rr3,y3_ci,y3_p)
    y4<-y[y$S>=thrd[4] & y$S<thrd[3],];n_y4_y<-nrow(y4[y4[,j]==1,]);n_y4_n<-nrow(y4[y4[,j]==0,]);y4_n<-rbind(c(n_y4_y,n_y4_n),c(n_y_y,n_y_n));y4_p<-fisher.test(y4_n,alternative = 'two.sided')$p.value;rr4<-(y4_n[1,1]/y4_n[1,2])/(y4_n[2,1]/y4_n[2,2]);y4_ci<-fisher.test(y4_n)$conf.int[1:2];data4<-c(rr4,y4_ci,y4_p)
    y5<-y[y$S>=thrd[5] & y$S<thrd[4],];n_y5_y<-nrow(y5[y5[,j]==1,]);n_y5_n<-nrow(y5[y5[,j]==0,]);y5_n<-rbind(c(n_y5_y,n_y5_n),c(n_y_y,n_y_n));y5_p<-fisher.test(y5_n,alternative = 'two.sided')$p.value;rr5<-(y5_n[1,1]/y5_n[1,2])/(y5_n[2,1]/y5_n[2,2]);y5_ci<-fisher.test(y5_n)$conf.int[1:2];data5<-c(rr5,y5_ci,y5_p)
    y6<-y[y$S>=thrd[6] & y$S<thrd[5],];n_y6_y<-nrow(y6[y6[,j]==1,]);n_y6_n<-nrow(y6[y6[,j]==0,]);y6_n<-rbind(c(n_y6_y,n_y6_n),c(n_y_y,n_y_n));y6_p<-fisher.test(y6_n,alternative = 'two.sided')$p.value;rr6<-(y6_n[1,1]/y6_n[1,2])/(y6_n[2,1]/y6_n[2,2]);y6_ci<-fisher.test(y6_n)$conf.int[1:2];data6<-c(rr6,y6_ci,y6_p)
    y7<-y[y$S>=thrd[7] & y$S<thrd[6],];n_y7_y<-nrow(y7[y7[,j]==1,]);n_y7_n<-nrow(y7[y7[,j]==0,]);y7_n<-rbind(c(n_y7_y,n_y7_n),c(n_y_y,n_y_n));y7_p<-fisher.test(y7_n,alternative = 'two.sided')$p.value;rr7<-(y7_n[1,1]/y7_n[1,2])/(y7_n[2,1]/y7_n[2,2]);y7_ci<-fisher.test(y7_n)$conf.int[1:2];data7<-c(rr7,y7_ci,y7_p)
    y8<-y[y$S>=thrd[8] & y$S<thrd[7],];n_y8_y<-nrow(y8[y8[,j]==1,]);n_y8_n<-nrow(y8[y8[,j]==0,]);y8_n<-rbind(c(n_y8_y,n_y8_n),c(n_y_y,n_y_n));y8_p<-fisher.test(y8_n,alternative = 'two.sided')$p.value;rr8<-(y8_n[1,1]/y8_n[1,2])/(y8_n[2,1]/y8_n[2,2]);y8_ci<-fisher.test(y8_n)$conf.int[1:2];data8<-c(rr8,y8_ci,y8_p)
    y9<-y[y$S>=thrd[9] & y$S<thrd[8],];n_y9_y<-nrow(y9[y9[,j]==1,]);n_y9_n<-nrow(y9[y9[,j]==0,]);y9_n<-rbind(c(n_y9_y,n_y9_n),c(n_y_y,n_y_n));y9_p<-fisher.test(y9_n,alternative = 'two.sided')$p.value;rr9<-(y9_n[1,1]/y9_n[1,2])/(y9_n[2,1]/y9_n[2,2]);y9_ci<-fisher.test(y9_n)$conf.int[1:2];data9<-c(rr9,y9_ci,y9_p)
    y10<-y[y$S<thrd[9],];n_y10_y<-nrow(y10[y10[,j]==1,]);n_y10_n<-nrow(y10[y10[,j]==0,]);y10_n<-rbind(c(n_y10_y,n_y10_n),c(n_y_y,n_y_n));y10_p<-fisher.test(y10_n,alternative = 'two.sided')$p.value;rr10<-(y10_n[1,1]/y10_n[1,2])/(y10_n[2,1]/y10_n[2,2]);y10_ci<-fisher.test(y10_n)$conf.int[1:2];data10<-c(rr10,y10_ci,y10_p)
  }else{
    y<-phe;j<-phe_name
    if(phe_label=='ADHD'){y<-y[unfactor(y$age_at_eval_months)>=60,]};if(phe_label=='Bipolar'){y<-y[unfactor(y$age_at_eval_months)>=180,]};if(phe_label=='Schizophrenia'){y<-y[unfactor(y$age_at_eval_months)<=180,]}
    if(phe_label=='Depression'){y<-y[unfactor(y$age_at_eval_months)>=120,]};if(phe_label=='Anxiety'){y<-y[unfactor(y$age_at_eval_months)>=120,]}
    thrd <- rev(quantile(y$S, probs = seq(0, 1, length.out = 11), na.rm = TRUE)[1:10]);thrd8<-thrd[8]
    middle_y<-y[y$S<thrd8,];n_y_y<-nrow(middle_y[!is.na(middle_y[,j]) & middle_y[,j]==1,]);n_y_n<-nrow(middle_y[is.na(middle_y[,j]),])
    y1<-y[y$S>=thrd[1],];n_y1_y<-nrow(y1[!is.na(y1[,j]) & y1[,j]==1,]);n_y1_n<-nrow(y1[is.na(y1[,j]),]);y1_n<-rbind(c(n_y1_y,n_y1_n),c(n_y_y,n_y_n));y1_p<-fisher.test(y1_n,alternative = 'two.sided')$p.value;rr1<-(y1_n[1,1]/y1_n[1,2])/(y1_n[2,1]/y1_n[2,2]);y1_ci<-fisher.test(y1_n)$conf.int[1:2];data1<-c(rr1,y1_ci,y1_p)
    y2<-y[y$S>=thrd[2] & y$S<thrd[1],];n_y2_y<-nrow(y2[!is.na(y2[,j]) & y2[,j]==1,]);n_y2_n<-nrow(y2[is.na(y2[,j]),]);y2_n<-rbind(c(n_y2_y,n_y2_n),c(n_y_y,n_y_n));y2_p<-fisher.test(y2_n,alternative = 'two.sided')$p.value;rr2<-(y2_n[1,1]/y2_n[1,2])/(y2_n[2,1]/y2_n[2,2]);y2_ci<-fisher.test(y2_n)$conf.int[1:2];data2<-c(rr2,y2_ci,y2_p)
    y3<-y[y$S>=thrd[3] & y$S<thrd[2],];n_y3_y<-nrow(y3[!is.na(y3[,j]) & y3[,j]==1,]);n_y3_n<-nrow(y3[is.na(y3[,j]),]);y3_n<-rbind(c(n_y3_y,n_y3_n),c(n_y_y,n_y_n));y3_p<-fisher.test(y3_n,alternative = 'two.sided')$p.value;rr3<-(y3_n[1,1]/y3_n[1,2])/(y3_n[2,1]/y3_n[2,2]);y3_ci<-fisher.test(y3_n)$conf.int[1:2];data3<-c(rr3,y3_ci,y3_p)
    y4<-y[y$S>=thrd[4] & y$S<thrd[3],];n_y4_y<-nrow(y4[!is.na(y4[,j]) & y4[,j]==1,]);n_y4_n<-nrow(y4[is.na(y4[,j]),]);y4_n<-rbind(c(n_y4_y,n_y4_n),c(n_y_y,n_y_n));y4_p<-fisher.test(y4_n,alternative = 'two.sided')$p.value;rr4<-(y4_n[1,1]/y4_n[1,2])/(y4_n[2,1]/y4_n[2,2]);y4_ci<-fisher.test(y4_n)$conf.int[1:2];data4<-c(rr4,y4_ci,y4_p)
    y5<-y[y$S>=thrd[5] & y$S<thrd[4],];n_y5_y<-nrow(y5[!is.na(y5[,j]) & y5[,j]==1,]);n_y5_n<-nrow(y5[is.na(y5[,j]),]);y5_n<-rbind(c(n_y5_y,n_y5_n),c(n_y_y,n_y_n));y5_p<-fisher.test(y5_n,alternative = 'two.sided')$p.value;rr5<-(y5_n[1,1]/y5_n[1,2])/(y5_n[2,1]/y5_n[2,2]);y5_ci<-fisher.test(y5_n)$conf.int[1:2];data5<-c(rr5,y5_ci,y5_p)
    y6<-y[y$S>=thrd[6] & y$S<thrd[5],];n_y6_y<-nrow(y6[!is.na(y6[,j]) & y6[,j]==1,]);n_y6_n<-nrow(y6[is.na(y6[,j]),]);y6_n<-rbind(c(n_y6_y,n_y6_n),c(n_y_y,n_y_n));y6_p<-fisher.test(y6_n,alternative = 'two.sided')$p.value;rr6<-(y6_n[1,1]/y6_n[1,2])/(y6_n[2,1]/y6_n[2,2]);y6_ci<-fisher.test(y6_n)$conf.int[1:2];data6<-c(rr6,y6_ci,y6_p)
    y7<-y[y$S>=thrd[7] & y$S<thrd[6],];n_y7_y<-nrow(y7[!is.na(y7[,j]) & y7[,j]==1,]);n_y7_n<-nrow(y7[is.na(y7[,j]),]);y7_n<-rbind(c(n_y7_y,n_y7_n),c(n_y_y,n_y_n));y7_p<-fisher.test(y7_n,alternative = 'two.sided')$p.value;rr7<-(y7_n[1,1]/y7_n[1,2])/(y7_n[2,1]/y7_n[2,2]);y7_ci<-fisher.test(y7_n)$conf.int[1:2];data7<-c(rr7,y7_ci,y7_p)
    y8<-y[y$S>=thrd[8] & y$S<thrd[7],];n_y8_y<-nrow(y8[!is.na(y8[,j]) & y8[,j]==1,]);n_y8_n<-nrow(y8[is.na(y8[,j]),]);y8_n<-rbind(c(n_y8_y,n_y8_n),c(n_y_y,n_y_n));y8_p<-fisher.test(y8_n,alternative = 'two.sided')$p.value;rr8<-(y8_n[1,1]/y8_n[1,2])/(y8_n[2,1]/y8_n[2,2]);y8_ci<-fisher.test(y8_n)$conf.int[1:2];data8<-c(rr8,y8_ci,y8_p)
    y9<-y[y$S>=thrd[9] & y$S<thrd[8],];n_y9_y<-nrow(y9[!is.na(y9[,j]) & y9[,j]==1,]);n_y9_n<-nrow(y9[is.na(y9[,j]),]);y9_n<-rbind(c(n_y9_y,n_y9_n),c(n_y_y,n_y_n));y9_p<-fisher.test(y9_n,alternative = 'two.sided')$p.value;rr9<-(y9_n[1,1]/y9_n[1,2])/(y9_n[2,1]/y9_n[2,2]);y9_ci<-fisher.test(y9_n)$conf.int[1:2];data9<-c(rr9,y9_ci,y9_p)
    y10<-y[y$S<thrd[9],];n_y10_y<-nrow(y10[!is.na(y10[,j]) & y10[,j]==1,]);n_y10_n<-nrow(y10[is.na(y10[,j]),]);y10_n<-rbind(c(n_y10_y,n_y10_n),c(n_y_y,n_y_n));y10_p<-fisher.test(y10_n,alternative = 'two.sided')$p.value;rr10<-(y10_n[1,1]/y10_n[1,2])/(y10_n[2,1]/y10_n[2,2]);y10_ci<-fisher.test(y10_n)$conf.int[1:2];data10<-c(rr10,y10_ci,y10_p)
  }
  phe1<-rbind(data1,data2,data3,data4,data5,data6,data7,data8,data9,data10)
  phe2<-cbind.data.frame(group_names,phe1);phe2$Pheno<-phe_label;colnames(phe2)<-c('Decile','OR','CI1','CI2','P','Pheno')
  phe3<-phe2[nrow(phe2):1,];phe3$Geneset<-geneset
  phe3<-phe3 %>% mutate(sig_label = case_when(P < 0.001 ~ "***",P < 0.01  ~ "**",P < 0.05  ~ "*",TRUE ~ ""))
  return(phe3)
}

all_phe_list<-list()
all_colors<-c('#7A52A6','#F49600','#0431FE','#FF8AD7','#29A15C','#DD201C','#7AB0DE')
imp_all<-merge(x,imp,by.x='ID',by.y='subject_sp_id')
imp_all1<-decile(imp_all,'derived_cog_impair','Cognitive impairment','All genes')
imp_all1$CI1 <- as.numeric(imp_all1$CI1);imp_all1$CI2 <- as.numeric(imp_all1$CI2);imp_all1$OR  <- as.numeric(imp_all1$OR)
imp_all1$Decile <- factor(imp_all1$Decile, levels = as.character(1:10));imp_all1$DecileNum <- as.numeric(as.character(imp_all1$Decile))
p <- ggplot(imp_all1, aes(x = DecileNum, y = OR)) +geom_ribbon(aes(ymin = CI1, ymax = CI2, group = 1),fill = all_colors[1], alpha = 0.25, color = NA) +
  geom_line(aes(group = 1), color = all_colors[1], size = 1.2) +geom_point(color = all_colors[1], size = 4) +
  geom_text(aes(label = sig_label), vjust = -0.05, color = all_colors[1], size = 6) +geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
  scale_x_continuous(breaks = 1:10, labels = 1:10) +scale_y_continuous(breaks = c(0.8, 1, 1.2, 1.4, 1.6, 1.8), limits = c(0.6, 2.0)) +
  labs(x = "Decile", y = "Odds ratio (OR)", title = unique(imp_all1$Pheno)) +theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, size = 15),axis.title = element_text(size = 15),axis.text = element_text(size = 15))
all_phe_list[[1]]<-p

phe_all<-merge(x,basic,by.x='ID',by.y='subject_sp_id')
term<-c('neuro_sz','behav_adhd','mood_bipol','mood_dep','mood_anx','schiz');term_names<-c('Seizure','ADHD','Bipolar','Depression','Anxiety','Schizophrenia')
known_colors<-c('#F49600','#0431FE','#FF8AD7','#29A15C','#DD201C','#7AB0DE');other_colors<-c('#F9CB80','#8098FD','#FEC6EC','#9CCBA7','#FEBEBD','#BDD7F1')
for(i in 1:length(term)){
  phe_all1<-decile(phe_all,term[i],term_names[i],'All genes')
  phe_all1$CI1 <- as.numeric(phe_all1$CI1);phe_all1$CI2 <- as.numeric(phe_all1$CI2);phe_all1$OR  <- as.numeric(phe_all1$OR)
  phe_all1[which(phe_all1$CI1<0.6),'CI1']<-0.6;phe_all1[which(phe_all1$CI2>2.0),'CI2']<-2
  phe_all1$Decile <- factor(phe_all1$Decile, levels = as.character(1:10));phe_all1$DecileNum <- as.numeric(as.character(phe_all1$Decile))
  p <- ggplot(phe_all1, aes(x = DecileNum, y = OR)) +geom_ribbon(aes(ymin = CI1, ymax = CI2, group = 1),fill = known_colors[i], alpha = 0.25, color = NA) +
    geom_line(aes(group = 1), color = known_colors[i], size = 1.2) +geom_point(color = known_colors[i], size = 4) +
    geom_text(aes(label = sig_label), vjust = -0.05, color = known_colors[i], size = 6) +geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
    scale_x_continuous(breaks = 1:10, labels = 1:10) +scale_y_continuous(breaks = c(0.8, 1, 1.2, 1.4, 1.6, 1.8), limits = c(0.6,2.0)) +
    labs(x = "Decile", y = "Odds ratio (OR)", title = unique(phe_all1$Pheno)) +theme_minimal(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5, size = 15),axis.title = element_text(size = 15),axis.text = element_text(size = 15))
  all_phe_list[[i+1]]<-p
}
all_combined_plot <- wrap_plots(all_phe_list, ncol = 3)
pdf('PaperFigures/output/Fig6b_RVCS_phenotypes.pdf',width = 12,height = 10)
print(all_combined_plot)
dev.off()


