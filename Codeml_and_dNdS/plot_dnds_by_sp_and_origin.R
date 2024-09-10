
setwd("/Users/fa8/Desktop/GRCs")

# scp -p farm22-head1:/lustre/scratch126/tol/teams/jaron/users/fede/CDS_seqs/CODEML_Results/DNDS_ALL.tsv .
# 

origins = read.table("ALL_CLASSIFICATIONS.ok",sep="\t",header=FALSE)
dnds = read.table("DNDS_ALL.tsv",sep="\t",header=FALSE, quote = "\"", skipNul = TRUE)

colnames(origins) = c("tree","sp","gene","origin","blen")
origins$gene = paste(origins$sp,"_grc_",origins$gene,sep="")
colnames(dnds) = c("gene","wbg","wfg","blen","lnl0","lnl1","ali","tree")

par(mfrow=c(3,2))
for(sp in c("bcop","bimp","ling")) {
	ce = origins[which(origins$sp==sp&origins$origin=="cecidomyiidae"),]
	sc = origins[which(origins$sp==sp&origins$origin=="sciaridae"),]
	
	dnds_ce = dnds[which(dnds$gene %in% ce$gene),]
	dnds_sc = dnds[which(dnds$gene %in% sc$gene),]
	
	plot(dnds_ce$wbg,dnds_ce$wfg,xlab="Tree dN/dS",ylab="Tip dN/dS",main=paste(sp,", origin=cecidomyiidae",sep=""),ylim=c(0,1),xlim=c(0,.2))
	abline(coef=c(0,1),col="orange",lty=2)
	plot(dnds_sc$wbg,dnds_sc$wfg,xlab="Tree dN/dS",ylab="Tip dN/dS",main=paste(sp,", origin=sciaridae",sep=""),ylim=c(0,1),xlim=c(0,.2))
	abline(coef=c(0,1),col="orange",lty=2)
}

par(mfrow=c(4,2))
dnds_ = dnds[grep("grc",dnds$gene,invert=T),]
for(sp in c("dmel","orobi","contarinia","aaphi","phyg","bcop","bimp","ling")) {
	
	dnds___ = dnds_[grep(sp,dnds_$gene),]
	
	plot(dnds___$wbg,dnds___$wfg,xlab="Tree dN/dS",ylab="Tip dN/dS",main=sp,ylim=c(0,1),xlim=c(0,.2))
	abline(coef=c(0,1),col="orange",lty=2)
}








