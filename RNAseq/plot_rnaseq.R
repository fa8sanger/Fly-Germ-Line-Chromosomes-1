setwd("/Users/fa8/Desktop/GRCs")

ling = read.table("Ling_RNAseq_by_Chr.tsv",row.names=1,sep="\t",header=T)

totals = apply(ling,1,sum)

ling_props = ling/totals
grcs = apply(ling_props[,grep("GRC",colnames(ling_props))],1,sum)

library("viridis")           
par(mfrow=c(1,2))
par(mar=c(10,4,4,4))
barplot(t(as.matrix(ling_props)),las=2,cex.names=.6,ylab="fraction of total expression",ylim=c(0,1.2),col=viridis(n=ncol(ling_props)))
legend("topleft",legend=colnames(ling_props),col=viridis(n=ncol(ling_props)),cex=.6,pch=15,ncol=2)
barplot(grcs,las=2,cex.names=.6,ylab="fraction of total expression",main="GRCs",col="darkcyan")


setwd("/Users/fa8/Desktop/GRCs")

bcop = read.table("bcop_RNAseq_by_Chr.tsv",row.names=1,sep="\t",header=T)

totals = apply(bcop,1,sum)

bcop_props = bcop/totals
grcs = apply(bcop_props[,grep("GRC",colnames(bcop_props))],1,sum)

library("viridis")           
par(mfrow=c(1,2))
par(mar=c(10,4,4,4))
barplot(t(as.matrix(bcop_props)),las=2,cex.names=.6,ylab="fraction of total expression",ylim=c(0,1.2),col=viridis(n=ncol(bcop_props)))
legend("topleft",legend=colnames(bcop_props),col=viridis(n=ncol(bcop_props)),cex=.6,pch=15,ncol=2)
barplot(grcs,las=2,cex.names=.6,ylab="fraction of total expression",main="GRCs",col="darkcyan")

##########################################################################################
# With new version: proper read counts: proper pairs, uniquily mapped, etc:
setwd("/Users/fa8/Desktop/GRCs")

ling = read.table("Ling.readcounts-Unique-f2-F256.tsv",row.names=1,sep="\t",header=T)

totals = apply(ling,1,sum)

ling_props = ling/totals
grcs = apply(ling_props[,grep("GRC",colnames(ling_props))],1,sum)

library("viridis")           
par(mfrow=c(1,2))
par(mar=c(11,4,4,4))
barplot(t(as.matrix(ling_props)),las=2,cex.names=.5,ylab="fraction of total expression",ylim=c(0,1.2),col=viridis(n=ncol(ling_props)))
legend("topleft",legend=colnames(ling_props),col=viridis(n=ncol(ling_props)),cex=.6,pch=15,ncol=2)
barplot(grcs,las=2,cex.names=.5,ylab="fraction of total expression",main="GRCs",col="darkcyan")


setwd("/Users/fa8/Desktop/GRCs")

bcop = read.table("Bcop.readcounts-Unique-f2-F256.tsv",row.names=1,sep="\t",header=T)

totals = apply(bcop,1,sum)

bcop = bcop[names(totals[which(totals>1e6)]),]
totals = totals[rownames(bcop)]

bcop_props = bcop/totals
grcs = apply(bcop_props[,grep("GRC",colnames(bcop_props))],1,sum)

library("viridis")           
par(mfrow=c(1,2))
par(mar=c(10,4,4,4))
barplot(t(as.matrix(bcop_props)),las=2,cex.names=.6,ylab="fraction of total expression",ylim=c(0,1.2),col=viridis(n=ncol(bcop_props)))
legend("topleft",legend=colnames(bcop_props),col=viridis(n=ncol(bcop_props)),cex=.6,pch=15,ncol=2)
barplot(grcs,las=2,cex.names=.6,ylab="fraction of total expression",main="GRCs",col="darkcyan")

