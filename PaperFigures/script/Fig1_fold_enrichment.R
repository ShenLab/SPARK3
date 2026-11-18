unfactor<-function(s){as.numeric(as.character(s))}
library(ggplot2);library(data.table);library(dplyr);library(minpack.lm);library(patchwork);library(epitools)
dnv<-read.csv('Data/autism_DNVs.txt',header = T,sep='\t')
trios<-read.csv('Data/autism_trio_samples.txt',header = T,sep='\t')
aff_dnv_n<-nrow(trios[trios$Pheno=='Affected',]);unaff_dnv_n<-nrow(trios[trios$Pheno=='Unaffected',])
dnv_aff<-dnv[dnv$Pheno=='Affected',];dnv_unaff<-dnv[dnv$Pheno=='Unaffected',]

inh<-read.csv('Data/autism_Inherited_variants.txt',header = T,sep='\t')
aff_inh_n<-46590

##########   Plot de novo and inherited separate   #########
### group the LoF and missense into Shet or MisFit-S bins
effect <- function(data1, data2, variant_class,inheritance, n1, n2, step = NULL) {
  if (variant_class == 'LGD') s <- 'Shet';if (variant_class == 'mis') s <- 'MisFit1.5_S'
  if(variant_class=='LGD'&inheritance=='DeNovo'){n_windows=100};if(variant_class=='mis'&inheritance=='DeNovo'){n_windows=120};if(variant_class=='LGD'&inheritance=='Inherited'){n_windows=100};if(variant_class=='mis'&inheritance=='Inherited'){n_windows=120};
  n <- nrow(data1);window_size <- floor(n / n_windows)
  selection <- data1[order(data1[, s]), s]
  if (is.null(step)) step <- floor(window_size / 2)   # sliding step (default = half window overlap)
  starts <- seq(1, n - window_size + 1, by = step)
  win <- data.frame(group = seq_along(starts),min_val = NA, max_val = NA,rate1 = NA, rate2 = NA,burden = NA, ci1 = NA, ci2 = NA, p = NA)
  
  for (i in seq_along(starts)) {
    idx <- starts[i]:(starts[i] + window_size - 1)
    win$min_val[i] <- min(selection[idx], na.rm = TRUE);win$max_val[i] <- max(selection[idx], na.rm = TRUE)
    count1 <- nrow(data1[unfactor(data1[, s]) >= win[i, 'min_val'] & unfactor(data1[, s]) <  win[i, 'max_val'], ])
    count2 <- nrow(data2[unfactor(data2[, s]) >= win[i, 'min_val'] & unfactor(data2[, s]) <  win[i, 'max_val'], ])
    if(count1>0 & count2>0){
      if(inheritance=='DeNovo'){a <- poisson.test(c(count1, count2), c(n1, n2));win$burden[i]<-a$estimate;win$ci1[i]<-a$conf.int[1];win$ci2[i]<-a$conf.int[2];win$p[i]<-a$p.value}
      if(inheritance=='Inherited'){a<-binom.test(count1,count1+count2,p=0.5);win$burden[i]<-a$estimate;win$ci1[i]<-a$conf.int[1];win$ci2[i]<-a$conf.int[2];win$p[i]<-a$p.value}
    }
    win$rate1[i]<-count1/n1;win$rate2[i]<-count2/n2
  }
  win<-win[win$rate1!=0 & win$rate2!=0,];win<-win[win$min_val!='Inf',];if(win[nrow(win),'max_val']<selection[length(selection)]){win[nrow(win),'max_val']<-1}
  return(win)
}
inh_effect<-function(count1,count2){
  a<-binom.test(count1,count1+count2,p=0.5);overtr<-a$estimate;ci<-a$conf.int;ci1<-ci[1];ci2<-ci[2];p<-a$p.value
  return(c(overtr,ci1,ci2,p))
}

lof_aff<-dnv_aff[dnv_aff$cq=='LGD' & dnv_aff$LoF=='HC' & !is.na(dnv_aff$Shet),];lof_unaff<-dnv_unaff[dnv_unaff$cq=='LGD' & dnv_unaff$LoF=='HC' & !is.na(dnv_unaff$Shet),];mis_aff<-dnv_aff[dnv_aff$cq=='mis'& !is.na(dnv_aff$MisFit1.5_S),];mis_unaff<-dnv_unaff[dnv_unaff$cq=='mis'& !is.na(dnv_unaff$MisFit1.5_S),]
dnv_lof_eff<-effect(lof_aff,lof_unaff,'LGD','DeNovo',aff_dnv_n,unaff_dnv_n);dnv_mis_eff<-effect(mis_aff,mis_unaff,'mis','DeNovo',aff_dnv_n,unaff_dnv_n)
dnv_lof_eff$Variant <- "LoF";dnv_mis_eff$Variant <- "Missense"
dnv_data <- bind_rows(dnv_lof_eff, dnv_mis_eff)
lof_aff<-dnv_aff[dnv_aff$cq=='LGD' & dnv_aff$LoF=='HC',];lof_unaff<-dnv_unaff[dnv_unaff$cq=='LGD' & dnv_unaff$LoF=='HC',];mis_aff<-dnv_aff[dnv_aff$cq=='mis'& !is.na(dnv_aff$MisFit1.5_S),];mis_unaff<-dnv_unaff[dnv_unaff$cq=='mis'& !is.na(dnv_unaff$MisFit1.5_S),]
lof_aff_ct1<-nrow(lof_aff[unfactor(lof_aff$Shet)>=0.1,]);lof_unaff_ct1<-nrow(lof_unaff[unfactor(lof_unaff$Shet)>=0.1,]);a<-poisson.test(c(lof_aff_ct1,lof_unaff_ct1),c(aff_dnv_n,unaff_dnv_n));lof_dnv1_burden<-a$estimate;lof_dnv1_ci1<-a$conf.int[1];lof_dnv1_ci2<-a$conf.int[2];lof_dnv1_p<-a$p.value
lof_aff_ct2<-nrow(lof_aff[unfactor(lof_aff$Shet)<0.1&unfactor(lof_aff$Shet)>=0.01,]);lof_unaff_ct2<-nrow(lof_unaff[unfactor(lof_unaff$Shet)<0.1&unfactor(lof_unaff$Shet)>=0.01,]);a<-poisson.test(c(lof_aff_ct2,lof_unaff_ct2),c(aff_dnv_n,unaff_dnv_n));lof_dnv2_burden<-a$estimate;lof_dnv2_ci1<-a$conf.int[1];lof_dnv2_ci2<-a$conf.int[2];lof_dnv2_p<-a$p.value
lof_aff_ct3<-nrow(lof_aff[unfactor(lof_aff$Shet)<0.01&unfactor(lof_aff$Shet)>=1e-3,]);lof_unaff_ct3<-nrow(lof_unaff[unfactor(lof_unaff$Shet)<0.01&unfactor(lof_unaff$Shet)>=1e-3,]);a<-poisson.test(c(lof_aff_ct3,lof_unaff_ct3),c(aff_dnv_n,unaff_dnv_n));lof_dnv3_burden<-a$estimate;lof_dnv3_ci1<-a$conf.int[1];lof_dnv3_ci2<-a$conf.int[2];lof_dnv3_p<-a$p.value
lof_aff_ct4<-nrow(lof_aff[unfactor(lof_aff$Shet)<1e-3&unfactor(lof_aff$Shet)>=1e-4,]);lof_unaff_ct4<-nrow(lof_unaff[unfactor(lof_unaff$Shet)<1e-3&unfactor(lof_unaff$Shet)>=1e-4,]);a<-poisson.test(c(lof_aff_ct4,lof_unaff_ct4),c(aff_dnv_n,unaff_dnv_n));lof_dnv4_burden<-a$estimate;lof_dnv4_ci1<-a$conf.int[1];lof_dnv4_ci2<-a$conf.int[2];lof_dnv4_p<-a$p.value
lof_aff_ct5<-nrow(lof_aff[unfactor(lof_aff$Shet)<1e-4,]);lof_unaff_ct5<-nrow(lof_unaff[unfactor(lof_unaff$Shet)<1e-4,]);a<-poisson.test(c(lof_aff_ct5,lof_unaff_ct5),c(aff_dnv_n,unaff_dnv_n));lof_dnv5_burden<-a$estimate;lof_dnv5_ci1<-a$conf.int[1];lof_dnv5_ci2<-a$conf.int[2];lof_dnv5_p<-a$p.value
mis_aff_ct1<-nrow(dnv_aff[dnv_aff$cq=='mis' & unfactor(dnv_aff$MisFit1.5_S)>=0.1,]);mis_unaff_ct1<-nrow(dnv_unaff[dnv_unaff$cq=='mis' & unfactor(dnv_unaff$MisFit1.5_S)>=0.1,]);a<-poisson.test(c(mis_aff_ct1,mis_unaff_ct1),c(aff_dnv_n,unaff_dnv_n));mis_dnv1_burden<-a$estimate;mis_dnv1_ci1<-a$conf.int[1];mis_dnv1_ci2<-a$conf.int[2];mis_dnv1_p<-a$p.value
mis_aff_ct2<-nrow(dnv_aff[dnv_aff$cq=='mis'&unfactor(dnv_aff$MisFit1.5_S)<0.1&unfactor(dnv_aff$MisFit1.5_S)>=0.01,]);mis_unaff_ct2<-nrow(dnv_unaff[dnv_unaff$cq=='mis'&unfactor(dnv_unaff$MisFit1.5_S)<0.1&unfactor(dnv_unaff$MisFit1.5_S)>=0.01,]);a<-poisson.test(c(mis_aff_ct2,mis_unaff_ct2),c(aff_dnv_n,unaff_dnv_n));mis_dnv2_burden<-a$estimate;mis_dnv2_ci1<-a$conf.int[1];mis_dnv2_ci2<-a$conf.int[2];mis_dnv2_p<-a$p.value
mis_aff_ct3<-nrow(dnv_aff[dnv_aff$cq=='mis'&unfactor(dnv_aff$MisFit1.5_S)<0.01&unfactor(dnv_aff$MisFit1.5_S)>=1e-3,]);mis_unaff_ct3<-nrow(dnv_unaff[dnv_unaff$cq=='mis'&unfactor(dnv_unaff$MisFit1.5_S)<0.01&unfactor(dnv_unaff$MisFit1.5_S)>=1e-3,]);a<-poisson.test(c(mis_aff_ct3,mis_unaff_ct3),c(aff_dnv_n,unaff_dnv_n));mis_dnv3_burden<-a$estimate;mis_dnv3_ci1<-a$conf.int[1];mis_dnv3_ci2<-a$conf.int[2];mis_dnv3_p<-a$p.value
mis_aff_ct4<-nrow(dnv_aff[dnv_aff$cq=='mis'&unfactor(dnv_aff$MisFit1.5_S)<1e-3&unfactor(dnv_aff$MisFit1.5_S)>=1e-4,]);mis_unaff_ct4<-nrow(dnv_unaff[dnv_unaff$cq=='mis'&unfactor(dnv_unaff$MisFit1.5_S)<1e-3&unfactor(dnv_unaff$MisFit1.5_S)>=1e-4,]);a<-poisson.test(c(mis_aff_ct4,mis_unaff_ct4),c(aff_dnv_n,unaff_dnv_n));mis_dnv4_burden<-a$estimate;mis_dnv4_ci1<-a$conf.int[1];mis_dnv4_ci2<-a$conf.int[2];mis_dnv4_p<-a$p.value
mis_aff_ct5<-nrow(dnv_aff[dnv_aff$cq=='mis'&unfactor(dnv_aff$MisFit1.5_S)<1e-4,]);mis_unaff_ct5<-nrow(dnv_unaff[dnv_unaff$cq=='mis'&unfactor(dnv_unaff$MisFit1.5_S)<1e-4,]);a<-poisson.test(c(mis_aff_ct5,mis_unaff_ct5),c(aff_dnv_n,unaff_dnv_n));mis_dnv5_burden<-a$estimate;mis_dnv5_ci1<-a$conf.int[1];mis_dnv5_ci2<-a$conf.int[2];mis_dnv5_p<-a$p.value
dnv_bin<-rbind(c('>=0.1',0.32,'LoF',lof_aff_ct1/aff_dnv_n,lof_unaff_ct1/unaff_dnv_n,lof_dnv1_burden,lof_dnv1_ci1,lof_dnv1_ci2,lof_dnv1_p),c('>=0.1',0.32,'Missense',mis_aff_ct1/aff_dnv_n,mis_unaff_ct1/unaff_dnv_n,mis_dnv1_burden,mis_dnv1_ci1,mis_dnv1_ci2,mis_dnv1_p),
               c('[0.01, 0.1)',0.032,'LoF',lof_aff_ct2/aff_dnv_n,lof_unaff_ct2/unaff_dnv_n,lof_dnv2_burden,lof_dnv2_ci1,lof_dnv2_ci2,lof_dnv2_p),c('[0.01, 0.1)',0.032,'Missense',mis_aff_ct2/aff_dnv_n,mis_unaff_ct2/unaff_dnv_n,mis_dnv2_burden,mis_dnv2_ci1,mis_dnv2_ci2,mis_dnv2_p),
               c('[1e-3, 0.01)',3.2e-3,'LoF',lof_aff_ct3/aff_dnv_n,lof_unaff_ct3/unaff_dnv_n,lof_dnv3_burden,lof_dnv3_ci1,lof_dnv3_ci2,lof_dnv3_p),c('[1e-3, 0.01)',3.2e-3,'Missense',mis_aff_ct3/aff_dnv_n,mis_unaff_ct3/unaff_dnv_n,mis_dnv3_burden,mis_dnv3_ci1,mis_dnv3_ci2,mis_dnv3_p),
               c('[1e-4, 1e-3)',3.2e-4,'LoF',lof_aff_ct4/aff_dnv_n,lof_unaff_ct4/unaff_dnv_n,lof_dnv4_burden,lof_dnv4_ci1,lof_dnv4_ci2,lof_dnv4_p),c('[1e-4, 1e-3)',3.2e-4,'Missense',mis_aff_ct4/aff_dnv_n,mis_unaff_ct4/unaff_dnv_n,mis_dnv4_burden,mis_dnv4_ci1,mis_dnv4_ci2,mis_dnv4_p),
               c('<1e-4',5e-5,'LoF',lof_aff_ct5/aff_dnv_n,lof_unaff_ct5/unaff_dnv_n,lof_dnv5_burden,lof_dnv5_ci1,lof_dnv5_ci2,lof_dnv5_p),c('<1e-4',5e-5,'Missense',mis_aff_ct5/aff_dnv_n,mis_unaff_ct5/unaff_dnv_n,mis_dnv5_burden,mis_dnv5_ci1,mis_dnv5_ci2,mis_dnv5_p))
dnv_bin<-as.data.frame(dnv_bin);colnames(dnv_bin)<-c('bin','selection','Variant','rate1','rate2','Burden','CI1','CI2','p');dnv_bin$p<-unfactor(dnv_bin$p);dnv_bin$Burden<-unfactor(dnv_bin$Burden);dnv_bin$CI1<-unfactor(dnv_bin$CI1);dnv_bin$CI2<-unfactor(dnv_bin$CI2);dnv_bin$selection<-unfactor(dnv_bin$selection)
dnv_bin <- dnv_bin %>%mutate(sig_label = case_when(p < 0.001 ~ "***",p < 0.005 ~ "**",p < 0.01 ~ "*",TRUE ~ ""))
dnv_plot<-ggplot() +geom_point(data = dnv_data,aes(x = (max_val + min_val)/2, y = burden, color = Variant),alpha = 0.2, size = 2.5) +
  geom_point(data = dnv_bin,aes(x = selection, y = Burden, color = Variant),position = position_dodge(0.3), size = 5) +
  geom_errorbar(data = dnv_bin,aes(x = selection, ymin = CI1, ymax = CI2, color = Variant),position = position_dodge(0.3), width = 0.3,size=1) +
  geom_text(data = dnv_bin,aes(x = selection, y = CI2, label = sig_label, color = Variant),vjust = 0.1, position = position_dodge(0.3),size = 10, show.legend = FALSE) +
  geom_text(data = dnv_bin,aes(x = selection, y = 7, label = bin),inherit.aes = FALSE, vjust = 1.5, size = 8) +
  scale_x_log10(breaks = c(1e-4, 1e-3, 1e-2, 1e-1, 1),labels = c("1e-4", "1e-3", "0.01", "0.1", "1")) +scale_y_log10(breaks = c(1, 1.5, 2, 3, 5, 8, 10, 16),labels = scales::label_number()) +
  coord_cartesian(ylim = c(0.7, 16)) +scale_color_manual(values = c("LoF" = "#FC840D",  "Missense" = "#257ab6")) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray30") +geom_vline(xintercept = c(1e-4, 1e-3, 1e-2, 1e-1, 1),linetype = "dotted", color = "gray30") +theme_minimal(base_size = 14) +labs(x = "Selection coefficient", y = "Enrichment vs without autism",color = "Variant") +
  theme(legend.position = c(0.15, 0.9),axis.title = element_text(size = 30),axis.text.x = element_text(size = 25),axis.text.y = element_text(size = 25),legend.text = element_text(size = 30),legend.title = element_text(size = 30),panel.grid.major.x = element_blank())

lof_tr<-inh[inh$cq=='LGD' & inh$LoF=='HC' & inh$iGT==1 & !is.na(inh$Shet),];lof_nontr<-inh[inh$cq=='LGD' & inh$LoF=='HC' & inh$iGT==0 & !is.na(inh$Shet),];mis_tr<-inh[inh$cq=='mis'& inh$iGT==1 & !is.na(inh$MisFit1.5_S),];mis_nontr<-inh[inh$cq=='mis' & inh$iGT==0 & !is.na(inh$MisFit1.5_S),]
inh_lof_eff<-effect(lof_tr,lof_nontr,'LGD','Inherited',aff_inh_n,aff_inh_n);inh_mis_eff<-effect(mis_tr,mis_nontr,'mis','Inherited',aff_inh_n,aff_inh_n)
inh_lof_eff$Variant <- "LoF";inh_mis_eff$Variant <- "Missense"
inh_data <- bind_rows(inh_lof_eff, inh_mis_eff)
lof_tr<-inh[inh$cq=='LGD' & inh$LoF=='HC' & inh$iGT==1,];lof_nontr<-inh[inh$cq=='LGD' & inh$LoF=='HC' & inh$iGT==0,];mis_tr<-inh[inh$cq=='mis'& inh$iGT==1,];mis_nontr<-inh[inh$cq=='mis' & inh$iGT==0,]
lof_tr_ct1<-nrow(lof_tr[unfactor(lof_tr$Shet)>=0.1,]);lof_nontr_ct1<-nrow(lof_nontr[unfactor(lof_nontr$Shet)>=0.1,]);lof_inh_p1<-inh_effect(lof_tr_ct1,lof_nontr_ct1)
lof_tr_ct2<-nrow(lof_tr[unfactor(lof_tr$Shet)<0.1&unfactor(lof_tr$Shet)>=0.01,]);lof_nontr_ct2<-nrow(lof_nontr[unfactor(lof_nontr$Shet)<0.1&unfactor(lof_nontr$Shet)>=0.01,]);lof_inh_p2<-inh_effect(lof_tr_ct2,lof_nontr_ct2)
lof_tr_ct3<-nrow(lof_tr[unfactor(lof_tr$Shet)<0.01&unfactor(lof_tr$Shet)>=1e-3,]);lof_nontr_ct3<-nrow(lof_nontr[unfactor(lof_nontr$Shet)<0.01&unfactor(lof_nontr$Shet)>=1e-3,]);lof_inh_p3<-inh_effect(lof_tr_ct3,lof_nontr_ct3)
lof_tr_ct4<-nrow(lof_tr[unfactor(lof_tr$Shet)<1e-3&unfactor(lof_tr$Shet)>=1e-4,]);lof_nontr_ct4<-nrow(lof_nontr[unfactor(lof_nontr$Shet)<1e-3&unfactor(lof_nontr$Shet)>=1e-4,]);lof_inh_p4<-inh_effect(lof_tr_ct4,lof_nontr_ct4)
lof_tr_ct5<-nrow(lof_tr[unfactor(lof_tr$Shet)<1e-4,]);lof_nontr_ct5<-nrow(lof_nontr[unfactor(lof_nontr$Shet)<1e-4,]);a<-binom.test(lof_tr_ct5,lof_tr_ct5+lof_nontr_ct5,p=0.5);lof_inh_p5<-inh_effect(lof_tr_ct5,lof_nontr_ct5)
mis_tr_ct1<-nrow(mis_tr[unfactor(mis_tr$MisFit1.5_S)>=0.1,]);mis_nontr_ct1<-nrow(mis_nontr[unfactor(mis_nontr$MisFit1.5_S)>=0.1,]);mis_inh_p1<-inh_effect(mis_tr_ct1,mis_nontr_ct1)
mis_tr_ct2<-nrow(mis_tr[unfactor(mis_tr$MisFit1.5_S)>=0.01&unfactor(mis_tr$MisFit1.5_S)<0.1,]);mis_nontr_ct2<-nrow(mis_nontr[unfactor(mis_nontr$MisFit1.5_S)>=0.01&unfactor(mis_nontr$MisFit1.5_S)<0.1,]);mis_inh_p2<-inh_effect(mis_tr_ct2,mis_nontr_ct2)
mis_tr_ct3<-nrow(mis_tr[unfactor(mis_tr$MisFit1.5_S)>=1e-3&unfactor(mis_tr$MisFit1.5_S)<0.01,]);mis_nontr_ct3<-nrow(mis_nontr[unfactor(mis_nontr$MisFit1.5_S)>=1e-3&unfactor(mis_nontr$MisFit1.5_S)<0.01,]);mis_inh_p3<-inh_effect(mis_tr_ct3,mis_nontr_ct3)
mis_tr_ct4<-nrow(mis_tr[unfactor(mis_tr$MisFit1.5_S)>=1e-4&unfactor(mis_tr$MisFit1.5_S)<1e-3,]);mis_nontr_ct4<-nrow(mis_nontr[unfactor(mis_nontr$MisFit1.5_S)>=1e-4&unfactor(mis_nontr$MisFit1.5_S)<1e-3,]);mis_inh_p4<-inh_effect(mis_tr_ct4,mis_nontr_ct4)
mis_tr_ct5<-nrow(mis_tr[unfactor(mis_tr$MisFit1.5_S)<1e-4,]);mis_nontr_ct5<-nrow(mis_nontr[unfactor(mis_nontr$MisFit1.5_S)<1e-4,]);a<-binom.test(mis_tr_ct5,mis_tr_ct5+mis_nontr_ct5,p=0.5);mis_inh_p5<-inh_effect(mis_tr_ct5,mis_nontr_ct5)
inh_bin<-rbind(c('>=0.1',0.32,'LoF',lof_tr_ct1/aff_inh_n,lof_nontr_ct1/aff_inh_n,lof_inh_p1),c('>=0.1',0.32,'Missense',mis_tr_ct1/aff_inh_n,mis_nontr_ct1/aff_inh_n,mis_inh_p1),
               c('[0.01, 0.1)',0.032,'LoF',lof_tr_ct2/aff_inh_n,lof_nontr_ct2/aff_inh_n,lof_inh_p2),c('[0.01, 0.1)',0.032,'Missense',mis_tr_ct2/aff_inh_n,mis_nontr_ct2/aff_inh_n,mis_inh_p2),
               c('[1e-3, 0.01)',3.2e-3,'LoF',lof_tr_ct3/aff_inh_n,lof_nontr_ct3/aff_inh_n,lof_inh_p3),c('[1e-3, 0.01)',3.2e-3,'Missense',mis_tr_ct3/aff_inh_n,mis_nontr_ct3/aff_inh_n,mis_inh_p3),
               c('[1e-4, 1e-3)',3.2e-4,'LoF',lof_tr_ct4/aff_inh_n,lof_nontr_ct4/aff_inh_n,lof_inh_p4),c('[1e-4, 1e-3)',3.2e-4,'Missense',mis_tr_ct4/aff_inh_n,mis_nontr_ct4/aff_inh_n,mis_inh_p4),
               c('<1e-4',5e-5,'LoF',lof_tr_ct5/aff_inh_n,lof_nontr_ct5/aff_inh_n,lof_inh_p5),c('<1e-4',5e-5,'Missense',mis_tr_ct5/aff_inh_n,mis_nontr_ct5/aff_inh_n,mis_inh_p5))
inh_bin<-as.data.frame(inh_bin);colnames(inh_bin)<-c('bin','selection','Variant','rate1','rate2','Burden','CI1','CI2','p');inh_bin$p<-unfactor(inh_bin$p);inh_bin$Burden<-unfactor(inh_bin$Burden);inh_bin$CI1<-unfactor(inh_bin$CI1);inh_bin$CI2<-unfactor(inh_bin$CI2);inh_bin$selection<-unfactor(inh_bin$selection)
inh_bin <- inh_bin %>%mutate(sig_label = case_when(p < 0.001 ~ "***",p < 0.005 ~ "**",p < 0.01 ~ "*",TRUE ~ ""))
inh_plot<-ggplot() +geom_point(data = inh_data,aes(x = (max_val + min_val)/2, y = burden, color = Variant),alpha = 0.2, size = 2.5) +
  geom_point(data = inh_bin,aes(x = selection, y = Burden, color = Variant),position = position_dodge(0.3), size = 5) +
  geom_errorbar(data = inh_bin,aes(x = selection, ymin = CI1, ymax = CI2, color = Variant),position = position_dodge(0.3), width = 0.3,size=1) +
  geom_text(data = inh_bin,aes(x = selection, y = CI2, label = sig_label, color = Variant),vjust = 0.1, position = position_dodge(0.3),size = 10, show.legend = FALSE) +
  geom_text(data = inh_bin,aes(x = selection, y = 0.55, label = bin),inherit.aes = FALSE, vjust = 1.5, size = 8) +
  scale_x_log10(breaks = c(1e-4, 1e-3, 1e-2, 1e-1, 1),labels = c("1e-4", "1e-3", "0.01", "0.1", "1")) +scale_y_continuous(breaks = seq(0.5, 0.57, 0.01),labels = sprintf("%.2f", seq(0.5, 0.57, 0.01)))+
  coord_cartesian(ylim = c(0.49,0.57)) +scale_color_manual(values = c("LoF" = "#E84445",  "Missense" = "#1999B2")) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray30") +geom_vline(xintercept = c(1e-4, 1e-3, 1e-2, 1e-1, 1),linetype = "dotted", color = "gray30") +theme_minimal(base_size = 14) +labs(x = "Selection coefficient", y = "% transmission",color = "Variant") +
  theme(legend.position = c(0.15, 0.9),axis.title = element_text(size = 30),axis.text.x = element_text(size = 25),axis.text.y = element_text(size = 25),legend.text = element_text(size = 30),legend.title = element_text(size = 30),panel.grid.major.x = element_blank())

pdf('PaperFigures/output/Figure1.pdf',width = 11,height = 10)
dnv_plot;inh_plot
dev.off()


