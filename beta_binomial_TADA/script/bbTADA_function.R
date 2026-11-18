## Original code and conception: Xin He (He et al.,2013), Jiebiao Wang, Lambertus Klei, and Jack Fu (Fu et al., 2022)
## Modifications and extensions: Lu Qiao
## Functions used to apply the TADA framework

################################################################################################
### log Bayes factor calculation for de novo variant contribution to a gene
## Input:
# x.dn: the de novo variant count 
# n.dn: the sample size (number of families)
# mu: the mutation rate (of this type of mutational events in this gene)
# Prior distribution of risk: gamma ~ Gamma(gamma.dn*beta.dn, beta.dn)
## Output:
# lBF: log Bayes factor evidence for de novo variants of this type in this gene
################################################################################################
log.bayes.factor.dn <- function(x.dn, n.dn, mu,gamma.dn, beta.dn){
  marg.lik0 <- dpois(x.dn, 2*n.dn*mu, log=TRUE)
  marg.lik1 <- dnbinom(x.dn, gamma.dn*beta.dn, beta.dn/(beta.dn+2*n.dn*mu), log=TRUE)
  lBF <- marg.lik1-marg.lik0
  return (lBF=lBF)
}

################################################################################################
### Bayes factor calculation for case-control/inherited contribution to a gene
## Input: 
# x.cc is a vector of case/transmitted and control/untransmitted variant counts in this gene
# n.cc is the number of case and control samples (or number of probands when evaluating inherited variants)
# gamma.cc is a vector of the prior mean on risk
# Prior distribution of q|H1: Gamma(rho1, nu1)
# Prior distribution of q|H0: Gamma(rho0, nu0)
## Output:
# BF: Bayes factor evidence for case-control/inherited contribution to a gene
################################################################################################
bayes.factor.cc <- function(x.cc, n.cc, gamma.cc, rho1, nu1, rho0, nu0){
  marglik0.cc <- evidence.null.cc(x.cc, n.cc, rho0, nu0)
  marglik1.cc <- evidence.alt.cc(x.cc, n.cc, gamma.cc, rho1, nu1)
  BF.cn <- marglik1.cc$cn / marglik0.cc$cn
  ## Assuming null and alt of the control model is identical
  if(is.na(BF.cn)){BF.cn<-1} 
  BF.ca <- marglik1.cc$ca / marglik0.cc$ca
  BF <- BF.cn * BF.ca
  return(BF=BF)
}

################################################################################################
### Helper function to bayes.factor.cc
### Calculates the evidence of case-control/inherited counts under null model
## Input: 
# x.cc is a vector of case/transmitted and control/untransmitted variant counts in this gene
# n.cc is the number of case and control samples (or number of probands when evaluating inherited variants)
# Prior distribution of q|H0: Gamma(rho0, nu0)
## Output:
# cn: marginal likelihood of control/untransmitted variant counts under null model
# ca: marginal likelihood of case/transmitted variant counts under null model
################################################################################################
evidence.null.cc <- function(x.cc, n.cc, rho0, nu0) {
  marglik0.ctrl.log <- log(dnbinom(x.cc$cn, rho0, nu0/(nu0+n.cc$cn)))
  marglik0.case.log <- log(dnbinom(x.cc$ca, rho0+x.cc$cn, (nu0+n.cc$cn)/(nu0+n.cc$cn+n.cc$ca)))
  marglik0.log <- marglik0.ctrl.log + marglik0.case.log
  
  return (list(cn=exp(marglik0.ctrl.log), ca=exp(marglik0.case.log)))#, total=exp(marglik0.log)))
}

################################################################################################
### Helper function to bayes.factor.cc
### Calculates the evidence of case-control/inherited counts under alternate model
## Input: 
# x.cc is a vector of case/transmitted and control/untransmitted variant counts in this gene
# n.cc is the number of case and control samples (or number of probands when evaluating inherited variants)
# gamma.cc is a vector of the prior mean on risk
## Output:
# cn: marginal likelihood of control/untransmitted variant counts under alternate model
# ca: marginal likelihood of case/transmitted variant counts under alternate model
################################################################################################
evidence.alt.cc <- function(x.cc, n.cc, gamma.cc, rho1, nu1){
  marglik1.ctrl <- dnbinom(x.cc$cn, rho1, nu1/(nu1+n.cc$cn))
  marglik1.case <- dnbinom(x.cc$ca, rho1+x.cc$cn, (nu1+n.cc$cn)/(nu1+n.cc$ca*gamma.cc+n.cc$cn))
  marglik1 <- marglik1.ctrl * marglik1.case
  return (list(cn=marglik1.ctrl, ca=marglik1.case))#, total=marglik1))
}

################################################################################################
### Wrapper function to obtain de novo Bayes Factors for SNV/indels, over a list of genes for a given variant type
### Makes calls to log.bayes.factor.dn for each gene
## Input: 
# count_case: a vector of de novo variant counts over each gene in cases
# count_con: a vector of de novo variant counts over each gene in controls
# n_case: the number of cases
# n_con: the number of controls
# mut: a vector of mutation rates for each gene
# Prior distribution of RR: gamma ~ Gamma(gamma.dn*beta.dn, beta.dn)
## Output:
# BF: a vector of Bayes factors of de novo evidence from SNV/indel for each gene for a given variant type
################################################################################################
BF_DN_SNV <- function(count_case, count_con, n_case, n_con, mut, gamma.dn, beta.dn=.2){
  ### If there are siblings, evaluate BF evidence for siblings to offset evidence for probands
  if(n_con>0){
    BF_con <- exp(sapply(1:length(mut), function(i) log.bayes.factor.dn(x.dn=count_con[i], n.dn=n_con, mu=mut[i], beta.dn=beta.dn, gamma.dn=gamma.dn[i])))
  }else{
    BF_con <- rep(1, length(mut))
  }
  BF_con[BF_con<1] <- 1
  
  ### BF evidence calculation for probands
  BF_case <- exp(sapply(1:length(mut), function(i) log.bayes.factor.dn(x.dn=count_case[i], n.dn=n_case, mu=mut[i], beta.dn=beta.dn, gamma.dn=gamma.dn[i])))
  BF_case[BF_case<1] <- 1
  
  ### Combine proband and sibling BF evidence
  ### Enforce BF back to 1 if no observations in probands or siblings
  BF <- pmax(BF_case/BF_con, 1) 
  BF[which(count_case==0 & count_con==0)] <- 1
  BF[is.na(BF)] <- 1
  return(BF)
}

################################################################################################
### Wrapper function to obtain case-control/inherited Bayes Factors for SNV/indels, over a list of genes for a given variant type
### Makes calls to bayes.factor.cc for each gene
## Input: 
# count_case: a vector of case/transmitted variant counts over each gene in cases
# count_con: a vector of control/untransmitted variant counts over each gene in controls
# n_case: the number of cases
# n_con: the number of controls (case-control analysis) or number of cases (inherited analysis)
# mut: a vector of mutation rates for each gene
# gamma.cc: a vector of the prior mean on risk
# nu: a hyperparameter that controls nu0 and nu1 in determining prior distributions, passed to bayes.factor.cc
## Output:
# BF: a vector of Bayes factors of case-control/inherited evidence from SNV/indel for each gene for a given variant type
################################################################################################
BF_CC_SNV <- function(count_case, count_con, n_case, n_con, mut, gamma.cc, nu=5000){
  rho.in = nu * sum(count_con, na.rm=T) / (2*length(mut)*n_con)
  rho.in = as.numeric(rho.in)*mut/mean(mut, na.rm=TRUE)
  BF = sapply((1:length(mut)), function(i){
    bayes.factor.cc(x.cc = data.frame(ca=count_case[i], cn=count_con[i]), n.cc=data.frame(ca=n_case, cn=n_con), gamma.cc=gamma.cc[i], rho1=rho.in[i], nu1=nu, rho0=rho.in[i], nu0=nu)
  })
  
  ### Enforce BF back to 1 if no mutation rate or no observed variant counts
  BF[mut==0] <- 1
  BF[which((count_case+count_con)==0)] <- 1
  return(BF) 
}

################################################################################################
### Beta-binomial extension of the TADA case-control Bayes factor Source, with BF_CC_SNV_BB(), a drop-in replacement for BF_CC_SNV() in the case-control part.
###
### The original TADA model is x_ctrl ~ Pois(q * n_ctrl),  x_case ~ Pois(q * gamma * n_case),  q ~ Gamma(rho, nu)
### For any q, this joint likelihood splits exactly into two pieces:
###   P(x_case, x_ctrl) = P(m) * P(x_case | m),   m = x_case + x_ctrl
###   (1) P(m)          : NB(rho, nu / (nu + n_ctrl + gamma * n_case))
###   (2) P(x_case | m) : Binomial(m, p), p = gamma*n_case / (gamma*n_case + n_ctrl)
### So the original BF_cc = BF_total(m) * BF_split.
### Piece (2), the case/control split, is the binomial that the synonymous QQ plot showed to be over-dispersed. We replace only that piece with the
### n-dependent beta-binomial used in beta_binomial_TADA/script/bbTADA_overdispersion_parameters.R:
###   x_case | m ~ BetaBin(m, p*M, (1-p)*M),  M = (1-rho)/rho, rho = a / (a + b * m_syn^c)
### m_syn is the gene's rare synonymous case + control count. Over-dispersion is treated as a property of the gene, calibrated where we expect no signal
### (synonymous, lambda ~ 1), and the same rho is used for LoF and all mis classes.
### Genes with no synonymous count (m_syn NA or 0) get the binomial (no correction). Same mean as the binomial, larger variance. As b -> Inf, rho_m -> 0 and
### BF_CC_SNV_BB() returns exactly the original BF_CC_SNV().
################################################################################################

## log beta-binomial pmf
ldbetabinom <- function(x, m, al, be) {
  lchoose(m, x) + lbeta(x + al, m - x + be) - lbeta(al, be)
}

## log BF for one gene
log.bayes.factor.cc.bb <- function(x.ca, x.cn, n.ca, n.cn, gamma.cc, rho, nu, m_syn,a,b,c) {
  m <- x.ca + x.cn
  if (is.na(m) || m == 0) return(0)
  
  ## (1) total count m, as in the original TADA (rho1 = rho0, nu1 = nu0)
  lm1 <- dnbinom(m, rho, nu / (nu + n.cn + gamma.cc * n.ca), log = TRUE)
  lm0 <- dnbinom(m, rho, nu / (nu + n.cn + n.ca), log = TRUE)
  
  ## (2) case/control split, beta-binomial with rho indexed by the gene's synonymous m
  p1 <- gamma.cc * n.ca / (gamma.cc * n.ca + n.cn)
  p0 <- n.ca / (n.ca + n.cn)
  if (is.na(m_syn) || m_syn <= 0) {            # no synonymous data: binomial
    ls1 <- dbinom(x.ca, m, p1, log = TRUE)
    ls0 <- dbinom(x.ca, m, p0, log = TRUE)
  } else {
    rho_g <- a / (a + b * m_syn^c)
    M <- (1 - rho_g) / rho_g
    ls1 <- ldbetabinom(x.ca, m, p1 * M, (1 - p1) * M)
    ls0 <- ldbetabinom(x.ca, m, p0 * M, (1 - p0) * M)
  }
  
  (lm1 - lm0) + (ls1 - ls0)
}

## Wrapper with the same arguments as BF_CC_SNV(), plus m_syn (one value per gene, same gene order as count_case) and a, b, c
BF_CC_SNV_BB <- function(count_case, count_con, n_case, n_con, mut, gamma.cc, m_syn, nu = 5000,a,b,c) {
  rho.in = nu * sum(count_con, na.rm = T) / (2 * length(mut) * n_con)
  rho.in = as.numeric(rho.in) * mut / mean(mut, na.rm = TRUE)
  BF = exp(sapply(seq_along(mut), function(i)
    log.bayes.factor.cc.bb(count_case[i], count_con[i], n_case, n_con,
                           gamma.cc[i], rho.in[i], nu, m_syn[i], a = a, b = b, c = c)))
  
  ### Enforce BF back to 1 if no mutation rate or no observed variant counts
  BF[mut == 0] <- 1
  BF[which((count_case + count_con) == 0)] <- 1
  return(BF)
}

################################################################################################
### Function to transform total BF evidence to FDR values
### Bayesian FDR control (PMID:19822692, Section2.3)
## Input: 
# BF: total Bayes factor matrix (each row is a gene, multiple columns allowed)
# pi0: estimated proportion of non-risk genes
## Output:
# FDR: a vector of FDR values corresponding to the rows of the input matrix
################################################################################################
Bayesian.FDR <- function(BF, pi0) {
  # order the BF in decreasing order, need to retain order to get results back in proper order 
  i.order=order(BF, decreasing = T)
  BF=BF[i.order];BF<-BF[!is.na(BF)]
  # convert BFs to PPA (posterior probability of alternative model)
  pi <- 1-pi0
  q <- pi*BF/(1-pi+pi*BF) # PPA
  q0 <- 1 - q # posterior probability of null model
  
  # the FDR at each PPA cutoff
  FDR=cumsum(q0)/(1:length(BF))
  
  # reorder to the original order
  FDR[i.order]=FDR
  q[i.order]=q
  return (cbind(q,FDR))
}


################################################################################################
### Helper function to distribute BF evidence inversely to LOEUF score for CNV evidence
## Input: 
# x: a vector of log Bayes factors evaluated for CNV events, with names labeled by a "," concatenated string of constituent genes
# info: matrix that contains gene level information (Supplementary Table 8: gene names, constraint)
## Output:
# evidence: a matrix containing the log Bayes (lBF) factor contributed to each gene(gene)from evaluated CNVs
################################################################################################
split_weight <- function(x, info){
  genes <- unlist(str_split(names(x), ","))
  oes <- info$LOEUF[match(genes, info$gene_gencodeV33)]
  weights <- sum(oes)/oes
  weights <- weights/sum(weights)
  return(evidence=data.frame(gene=genes, lBF=as.numeric(x)*weights))
}

################################################################################################
### Function to backtransform FDR to be on p-value scale
## Input: 
# x: a vector of FDRs
## Output:
# p: a vector of p-values
################################################################################################
q2p <- function(x){
  ec <- x
  ec <- sort(ec)
  ec <- ec*(1:length(ec))/(length(ec))/max(ec)    
  return(p=ec[names(x)])
}
