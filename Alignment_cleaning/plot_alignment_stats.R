
setwd("/Users/fa8/volumes/fa8_network/117/fede/Orthofinder/Proteomes_2nd_round_OneIsoform/OrthoFinder/Results_Jul25_1/ALI_STATS/")
sps = c("dmel","phyg","bcop-grc","bcop","bimp-grc","bimp","ling-grc","ling","aaphi","orobi","contarinia")
k = list()
for(sp in sps) {
	k[[sp]] = read.table(paste(sp,".STATS",sep=""),sep="\t",header=F)
}

covered_frac = list()
unique_frac = list()
ident_frac = list()
for(sp in sps) {
	covered_frac[[sp]] = k[[sp]][,3]
	unique_frac[[sp]]  = k[[sp]][,4]
	ident_frac[[sp]]  = k[[sp]][,5]
}

# library(ggplot2)
# rr$name = paste(rep("het",nrow(rr)),rr$het,sep="")
# p <- ggplot(covered_frac, aes(x=name, y=ref_frac, fill=name)) + # fill=name allow to automatically dedicate a color for each group
#   geom_violin()

# devtools::install_github("TomKellyGenetics/vioplot", ref = "dev")
# library("vioplot")
# 
# vioplot(x)

# Make a data frame of same size (sampling data):
# Max 4000 points:
  rr = data.frame(matrix(nrow=0,ncol=3))
  colnames(rr) = c("sp","cov","uniq")
  for(sp in sps) {
  	kk = sample(1:length(covered_frac[[sp]]),4000)
  	jj = cbind(rep(sp,4000),covered_frac[[sp]][kk],unique_frac[[sp]][kk])
  	rr = rbind(rr,jj)
  }
  colnames(rr) = c("sp","cov_frac","uniq_frac")
  rr$cov_frac = as.numeric(rr$cov_frac)
  rr$uniq_frac = as.numeric(rr$uniq_frac)
  library(ggplot2)
  p <- ggplot(rr, aes(x=sp, y=cov_frac, fill=sp)) + # fill=name allow to automatically dedicate a color for each group
       geom_violin()
  p <- ggplot(rr, aes(x=sp, y=uniq_frac, fill=sp)) + # fill=name allow to automatically dedicate a color for each group
       geom_violin()

par(mar=c(10,4,4,4))
par(mfrow=c(1,3))
boxplot(covered_frac,las=2,ylab="Covered fraction",outline=F)
boxplot(unique_frac,las=2,ylab="Unique fraction",outline=F)
boxplot(ident_frac,las=2,ylab="Identity fraction",outline=F)

# Now stratify ident fractions by completeness...
par(mfrow=c(3,4))
cor.tests = list()
for(sp in sps) {
	cor_ = cor(ident_frac[[sp]],covered_frac[[sp]])
	plot(ident_frac[[sp]],covered_frac[[sp]], main=paste(sp,"cor=",cor_),xlab="Ident frac",ylab="Covered frac",pch=20,col=rgb(.2,.2,.2,.3))
	
}





