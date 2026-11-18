### rare variant composite score (RVCS) per gene
args=commandArgs(trailingOnly=TRUE)
library(data.table);library(dplyr);library(openxlsx)
unfactor<-function(s){as.numeric(as.character(s))}
known<-read.xlsx('Supplementary Table.xlsx',sheet = 3)   ##509 autism/NDD known genes
known<-known[-1,];colnames(known)<-known[1,];known<-known[-1,]
pli<-read.table('beta_binomial_TADA/input/gnomAD4_pLI_GeneBayes_Shet.txt',header=T,sep='\t')

id=as.character(args[1])   ###gene symbols
inDir=as.character(args[2]) ###Directory storing variants in autistic individuals, siblings without autism, UKBB controls one by one 
file=paste0(inDir,"/",id,".txt")  ### variant in each individuals
x<-read.csv(file,header=T,sep='\t')
x$consequence='.'
x[which(x$canon_Consequence %like% "synonymous"),"consequence"]<-'syn'
x[which(x$canon_Consequence %like% "missense"),"consequence"]<-'mis'
x[grep('frameshift|splice_acceptor|splice_donor|stop_gained',x$canon_Consequence),'consequence']<-'LGD'
x[which(unfactor(x$SpliceAI_max)>=0.8 & x$canon_Consequence == "synonymous"),"consequence"]<-"LGD"
y<-merge(x,pli,by.x='canon_gene',by.y='gene',all.x=T);x<-y
z<-x[unfactor(x$topmed_AF)<=1e-4,]

x<-z  ### all genes
if(nrow(x)>0){
	x1<-x[(x$consequence=='mis' &!is.na(unfactor(x$MisFit1.5_S))& unfactor(x$MisFit1.5_S)>=0.01)|(x$consequence=='LGD' & !is.na(x$LoF) & x$LoF=='HC' & !is.na(x$Shet) & unfactor(x$Shet)>=0.01),]
	x2<-x1[!is.na(x1$ID),];x1<-x2
	if(nrow(x1)>0){
        	for(i in 1:nrow(x1)){
                	if(i==1){S=1}
                	if(x1[i,'consequence']=='mis'){s<-1-unfactor(x1[i,'MisFit1.5_S'])}
                	if(x1[i,'consequence']=='LGD'){s<-1-unfactor(x1[i,'Shet'])}
       			S<-S*s
        	}
        	S1<-1-S
        	all<-matrix(c(id,S1,nrow(x1)),nrow=1)
	}else{all<-c(id,0,0)}
}else{all<-c(id,0,0)}

x2<-merge(z,known,by.x='canon_gene',by.y='gene.symbol')
x<-x2  ### known genes
if(nrow(x)>0){
	x1<-x[(x$consequence=='mis' & unfactor(x$MisFit1.5_S)>=0.01)|(x$consequence=='LGD' & !is.na(x$LoF) & x$LoF=='HC' & !is.na(x$Shet) & unfactor(x$Shet)>=0.01),]
	x2<-x1[!is.na(x1$ID),];x1<-x2
	if(nrow(x1)>0){
        	for(i in 1:nrow(x1)){
                	if(i==1){S=1}
                	if(x1[i,'consequence']=='mis'){s<-1-unfactor(x1[i,'MisFit1.5_S'])}
                	if(x1[i,'consequence']=='LGD'){s<-1-unfactor(x1[i,'Shet'])}
       		S<-S*s
        	}
        	S1<-1-S
        	know<-matrix(c(S1,nrow(x1)),nrow=1)
	}else{know<-c(0,0)}
}else{know<-c(0,0)}

x3 <- z[!(z$canon_gene %in% known$gene.symbol), ]
x<-x3   ### unknown genes
if(nrow(x)>0){
	x1<-x[(x$consequence=='mis' & unfactor(x$MisFit1.5_S)>=0.01)|(x$consequence=='LGD' & !is.na(x$LoF) & x$LoF=='HC' & !is.na(x$Shet) & unfactor(x$Shet)>=0.01),]
	x2<-x1[!is.na(x1$ID),];x1<-x2
	if(nrow(x1)>0){
		for(i in 1:nrow(x1)){
			if(i==1){S=1}
			if(x1[i,'consequence']=='mis'){s<-1-unfactor(x1[i,'MisFit1.5_S'])}
			if(x1[i,'consequence']=='LGD'){s<-1-unfactor(x1[i,'Shet'])}
			S<-S*s
		}
		S1<-1-S
		unknown<-matrix(c(S1,nrow(x1)),nrow=1)
	}else{unknown<-c(0,0)}
}else{unknown<-c(0,0)}
final<-matrix(c(all,know,unknown),nrow=1)
write.table(final,paste0(inDir,'/',id,'_S.txt'),col.names=F,row.names=F,sep='\t',quote=F)

