library(GenomicRanges); library(stringr); library(openxlsx); library(readxl)
unfactor<-function(s){as.numeric(as.character(s))}
### Change to working directory where the folder was cloned
source('beta_binomial_TADA/script/bbTADA_function.R')
### File name for supplementary table
info<-read.csv('beta_binomial_TADA/output/bbTADA_prior_LoF_Dmis.txt',header=T,sep='\t')  ### Information on genes: constraint, priors
count<-read.csv('beta_binomial_TADA/output/bbTADA_gene_count.txt',header=T,sep='\t')     ### Loading gene count tables
count_gene<-cbind(count$Gene,1:nrow(count));colnames(count_gene)<-c('gene','num');count_gene<-as.data.frame(count_gene)
info1<-merge(info,count_gene,by='gene',all.y=T)
info2<-info1[order(unfactor(info1$num)),1:(ncol(info1)-1)]
info<-info2

### Setting sample sizes
n_aff_dnv<-34247;n_unaff_dnv<-10666
n_aff_inh<-18898
n_case<-33319;n_ctrl<-142251

### Bayes Factor calculations
beta.dn <- 0.2

### DN LoF
BF_dn_lof <- BF_DN_SNV(count_case=count$dn_aff_lof, count_con=count$dn_unaff_lof, n_case=n_aff_dnv, n_con=n_unaff_dnv, mut=count$Mu_lgd, gamma.dn=info$prior_dn_lof, beta.dn=beta.dn)
### DN mis1
BF_dn_mis1 <- BF_DN_SNV(count_case=count$dn_aff_mis1, count_con=count$dn_unaff_mis1, n_case=n_aff_dnv, n_con=n_unaff_dnv, mut=count$Mu_mis1, gamma.dn=info$prior_dn_mis1, beta.dn=beta.dn)
### DN mis2
BF_dn_mis2 <- BF_DN_SNV(count_case=count$dn_aff_mis2, count_con=count$dn_unaff_mis2, n_case=n_aff_dnv, n_con=n_unaff_dnv, mut=count$Mu_mis2, gamma.dn=info$prior_dn_mis2, beta.dn=beta.dn)
### DN mis3
BF_dn_mis3 <- BF_DN_SNV(count_case=count$dn_aff_mis3, count_con=count$dn_unaff_mis3, n_case=n_aff_dnv, n_con=n_unaff_dnv, mut=count$Mu_mis3, gamma.dn=info$prior_dn_mis3, beta.dn=beta.dn)

########################################################################
### IN LoF
BF_in_lof <- BF_CC_SNV(count_case=count$inh_trans_lof, count_con=count$inh_untrans_lof, n_case= n_aff_inh, n_con=n_aff_inh, mut=count$Mu_lgd, gamma.cc=info$prior_inh_lof)
### IN mis1
BF_in_mis1 <- BF_CC_SNV(count_case=count$inh_trans_mis1, count_con=count$inh_untrans_mis1, n_case= n_aff_inh, n_con=n_aff_inh, mut=count$Mu_mis1, gamma.cc=info$prior_inh_mis1)
### IN mis2
BF_in_mis2 <- BF_CC_SNV(count_case=count$inh_trans_mis2, count_con=count$inh_untrans_mis2, n_case= n_aff_inh, n_con=n_aff_inh, mut=count$Mu_mis2, gamma.cc=info$prior_inh_mis2)
### IN mis3
BF_in_mis3 <- BF_CC_SNV(count_case=count$inh_trans_mis3, count_con=count$inh_untrans_mis3, n_case= n_aff_inh, n_con=n_aff_inh, mut=count$Mu_mis3, gamma.cc=info$prior_inh_mis3)

########################################################################
### Case-control: n-dependent beta-binomial for the case/control split rho(m) = a/(a + b*m^c), estimated from synonymous. 
### rho is indexed by the gene's SYNONYMOUS case + control count (m_syn),
### the same value for LoF and all mis classes of that gene
bb<-read.csv('beta_binomial_TADA/output/ndep_rho_abc_by_class.txt',header = T,sep='\t')
bb_a <- bb$a[1]; bb_b <- bb$b[1]; bb_c <- bb$c[1]

syn <- read.csv('beta_binomial_TADA/input/case_ctrl_VarCount_perGene.txt', header = TRUE, sep = '\t')
syn<-syn[,c('gene','syn_count_case','syn_count_ctrl')];colnames(syn) <- c('gene', 'syn_case', 'syn_ctrl')
m_syn <- unfactor(syn$syn_case)[match(count$gene, syn$gene)] + unfactor(syn$syn_ctrl)[match(count$gene, syn$gene)]
cat(sum(is.na(m_syn) | m_syn == 0), 'genes have no synonymous count; their BF_cc stays binomial\n')
### CC LoF
BF_cc_lof <- BF_CC_SNV_BB(count_case=count$cc_case_lof, count_con=count$cc_ctrl_lof, n_case=n_case, n_con=n_ctrl, mut=count$Mu_lgd, gamma.cc=info$prior_cc_lof, m_syn=m_syn, a=bb_a, b=bb_b, c=bb_c)
### CC mis1
BF_cc_mis1 <- BF_CC_SNV_BB(count_case=count$cc_case_mis1, count_con=count$cc_ctrl_mis1, n_case=n_case, n_con=n_ctrl, mut=count$Mu_mis1, gamma.cc=info$prior_cc_mis1, m_syn=m_syn, a=bb_a, b=bb_b, c=bb_c)
### CC mis2
BF_cc_mis2 <- BF_CC_SNV_BB(count_case=count$cc_case_mis2, count_con=count$cc_ctrl_mis2, n_case=n_case, n_con=n_ctrl, mut=count$Mu_mis2, gamma.cc=info$prior_cc_mis2, m_syn=m_syn, a=bb_a, b=bb_b, c=bb_c)
### CC mis3
BF_cc_mis3 <- BF_CC_SNV_BB(count_case=count$cc_case_mis3, count_con=count$cc_ctrl_mis3, n_case=n_case, n_con=n_ctrl, mut=count$Mu_mis3, gamma.cc=info$prior_cc_mis3, m_syn=m_syn, a=bb_a, b=bb_b, c=bb_c)

########################################################################
### TADA outcome calculation
########################################################################

### Aggregate BF
BF_asd <- cbind(pmax(BF_dn_lof*BF_in_lof*BF_cc_lof, 1),
                pmax(BF_dn_mis1*BF_in_mis1*BF_cc_mis1, 1),
                pmax(BF_dn_mis2*BF_in_mis2*BF_cc_mis2, 1),
                pmax(BF_dn_mis3*BF_in_mis3*BF_cc_mis3, 1))
colnames(BF_asd)<-c('BF_lof','BF_mis1','BF_mis2','BF_mis3')
### FDR calculation
BF<-apply(BF_asd, 1, prod);BF[is.na(BF)]<-1
qval_asd <- Bayesian.FDR(BF, pi0 = 1 - 0.05);colnames(qval_asd)<-c('PPA','FDR');qval_asd<-as.data.frame(qval_asd)
x<-cbind(count[,c('Gene','Chrom')],info[,c('pLI','Shet')],BF_asd,BF,qval_asd,info[,4:ncol(info)],count[,4:ncol(count)])
y<-cbind(count[,c('Gene','Chrom')],info[,c('pLI','Shet')],BF_dn_lof,BF_in_lof,BF_cc_lof,pmax(BF_dn_lof*BF_in_lof*BF_cc_lof, 1),BF_dn_mis1,BF_in_mis1,BF_cc_mis1,pmax(BF_dn_mis1*BF_in_mis1*BF_cc_mis1, 1),BF_dn_mis2,BF_in_mis2,BF_cc_mis2,
         pmax(BF_dn_mis2*BF_in_mis2*BF_cc_mis2, 1),BF_dn_mis3,BF_in_mis3,BF_cc_mis3,pmax(BF_dn_mis3*BF_in_mis3*BF_cc_mis3, 1),BF,qval_asd,info[,4:ncol(info)],count[,4:ncol(count)])
colnames(y)<-c("gene","chr", "pLI","Shet","BF_dn_lof","BF_in_lof","BF_cc_lof","BF_lof","BF_dn_mis1","BF_in_mis1","BF_cc_mis1","BF_mis1","BF_dn_mis2","BF_in_mis2","BF_cc_mis2","BF_mis2","BF_dn_mis3","BF_in_mis3","BF_cc_mis3","BF_mis3","BF", 
               "PPA","FDR","prior_dn_lof","prior_dn_mis1","prior_dn_mis2","prior_dn_mis3","prior_inh_lof","prior_inh_mis1","prior_inh_mis2","prior_inh_mis3","prior_cc_lof","prior_cc_mis1","prior_cc_mis2","prior_cc_mis3","Mu_lgd","Mu_mis1","Mu_mis2","Mu_mis3",
               "dn_aff_lof","dn_aff_mis1","dn_aff_mis2","dn_aff_mis3","dn_unaff_lof","dn_unaff_mis1","dn_unaff_mis2","dn_unaff_mis3","inh_trans_lof","inh_trans_mis1","inh_trans_mis2","inh_trans_mis3","inh_untrans_lof","inh_untrans_mis1","inh_untrans_mis2",
               "inh_untrans_mis3","cc_case_lof","cc_case_mis1","cc_case_mis2","cc_case_mis3","cc_ctrl_lof","cc_ctrl_mis1","cc_ctrl_mis2","cc_ctrl_mis3")

##to get exome-wide p values
N <- 18698
y$pseudo_p <- NA
ord <- order(y$FDR)
y$pseudo_p[ord] <- y$FDR[ord] * seq_len(nrow(y)) / N
write.table(y,'beta_binomial_TADA/output/bbTADA_BF_FDR.txt',col.names = T,row.names = F,sep='\t',quote = F)


