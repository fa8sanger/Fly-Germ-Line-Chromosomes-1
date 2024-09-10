setwd("~/Desktop/GRCs/") 

tt = read.table("OG_BASED_SYNTENY_TABLE.tsv",sep="\t",header=FALSE)
colnames(tt) = c("og","genome","gene","chr","start","end","genome2","gene2","chr2","start2","end2")
qtt = tt[,1:6]
qtt$id = paste(qtt$genome,":",qtt$gene,sep="")
qtt = qtt[!duplicated(qtt),]
ttt = tt[,7:11]
ttt = ttt[!duplicated(ttt),]
ttt$id = paste(ttt$genome2,":",ttt$gene2,sep="")

div = read.table("ALL_DIVERGENCES.dS_version.tsv",sep="\t",header=TRUE)
div$id = paste(div$sp,":",div$gene,sep="")

# Add chrom info to div:
div2 = merge(div, qtt, by.x="id", by.y="id")
div2$tid = paste(div2$tsp,":",div2$tgene,sep="")
div3 = merge(div2, ttt, by.x="tid", by.y="id")


# For Bcop & origin=sciaridae, plot:
#   GRC1-SUPER_X
#   GRC1-SUPER_1
#   GRC1-SUPER_2
#   GRC1-SUPER_3
#   GRC2-SUPER_X
#   GRC2-SUPER_1
#   GRC2-SUPER_2
#   GRC2-SUPER_3

par(mfrow=c(1,3))
qsp="bcop_grc"
tsp="bcop_core"
origin="sciaridae"
qchrs=c("SUPER_GRC1","SUPER_GRC2")
tchrs=c("SUPER_X","SUPER_1","SUPER_2","SUPER_3")
div_ = div3[which(div3$sp==qsp & div3$tsp==tsp & div3$origin==origin & div3$chr %in% qchrs & div3$chr2 %in% tchrs),]
div_$factor = paste(gsub("SUPER_","",div_$chr)," vs ",gsub("SUPER_","",div_$chr2),sep="")
par(mar=c(10,4,4,4))
boxplot(dS~factor,data=div_,main=paste(qsp,"vs",tsp,"by origin=",origin),xlab="",ylab="dS",las=2,cex.names=.6,cex.main=.6,varwidth=TRUE,ylim=c(0,2),notch=T)

j=list()
for(qchr in qchrs) {
	div_ = div3[which(div3$sp==qsp & div3$tsp==tsp & div3$origin==origin & div3$chr==qchr & div3$chr2 %in% tchrs),]
	div_=div_[order(div_$start),]
	div_=div_[order(div_$chr),]
	j[[qchr]] = rle(div_$chr2)
	j[[qchr]] = as.data.frame(cbind(j[[qchr]]$values, j[[qchr]]$lengths))
	j[[qchr]][,2] = as.numeric(j[[qchr]][,2])
	cat("max for",qchr,"=",max(j[[qchr]][,2]),"\n")
}


qsp="ling_grc"
tsp="ling_core"
origin="sciaridae"
qchrs=c("SUPER_GRC1","SUPER_GRC2")
tchrs=c("SUPER_X","SUPER_1","SUPER_2","SUPER_3")
div_ = div3[which(div3$sp==qsp & div3$tsp==tsp & div3$origin==origin & div3$chr %in% qchrs & div3$chr2 %in% tchrs),]
div_$factor = paste(gsub("SUPER_","",div_$chr)," vs ",gsub("SUPER_","",div_$chr2),sep="")
par(mar=c(10,4,4,4))
boxplot(dS~factor,data=div_,main=paste(qsp,"vs",tsp,"by origin=",origin),xlab="",ylab="dS",las=2,cex.names=.6,cex.main=.6,varwidth=TRUE,ylim=c(0,2),notch=T)

j=list()
for(qchr in qchrs) {
	div_ = div3[which(div3$sp==qsp & div3$tsp==tsp & div3$origin==origin & div3$chr==qchr & div3$chr2 %in% tchrs),]
	div_=div_[order(div_$start),]
	div_=div_[order(div_$chr),]
	j[[qchr]] = rle(div_$chr2)
	j[[qchr]] = as.data.frame(cbind(j[[qchr]]$values, j[[qchr]]$lengths))
	j[[qchr]][,2] = as.numeric(j[[qchr]][,2])
	cat("max for",qchr,"=",max(j[[qchr]][,2]),"\n")
}

qsp="bimp_grc"
tsp="bimp_core"
origin="sciaridae"
qchrs=c("SUPER_GRC")
tchrs=c("SUPER_X","SUPER_1","SUPER_2","SUPER_3")
div_ = div3[which(div3$sp==qsp & div3$tsp==tsp & div3$origin==origin & div3$chr %in% qchrs & div3$chr2 %in% tchrs),]
div_$factor = paste(gsub("SUPER_","",div_$chr)," vs ",gsub("SUPER_","",div_$chr2),sep="")
par(mar=c(10,4,4,4))
boxplot(dS~factor,data=div_,main=paste(qsp,"vs",tsp,"by origin=",origin),xlab="",ylab="dS",las=2,cex.names=.6,cex.main=.6,varwidth=TRUE,ylim=c(0,2),notch=T)

j=list()
for(qchr in qchrs) {
	div_ = div3[which(div3$sp==qsp & div3$tsp==tsp & div3$origin==origin & div3$chr==qchr & div3$chr2 %in% tchrs),]
	div_=div_[order(div_$start),]
	div_=div_[order(div_$chr),]
	j[[qchr]] = rle(div_$chr2)
	j[[qchr]] = as.data.frame(cbind(j[[qchr]]$values, j[[qchr]]$lengths))
	j[[qchr]][,2] = as.numeric(j[[qchr]][,2])
	cat("max for",qchr,"=",max(j[[qchr]][,2]),"\n")
}



# Now cecido:
par(mfrow=c(1,3))
qsp="bcop_grc"
tsp="aaphi"
origin="cecidomyiidae"
qchrs=c("SUPER_GRC1","SUPER_GRC2")
tchrs=c("CM059996.1","CM059997.1","CM059998.1","CM059999.1")
div_ = div3[which(div3$sp==qsp & div3$tsp==tsp & div3$origin==origin & div3$chr %in% qchrs & div3$chr2 %in% tchrs),]
div_$factor = paste(gsub("SUPER_","",div_$chr)," vs ",gsub("SUPER_","",div_$chr2),sep="")
par(mar=c(10,4,4,4))
boxplot(dS~factor,data=div_,main=paste(qsp,"vs",tsp,"by origin=",origin),xlab="",ylab="dS",las=2,cex.names=.6,cex.main=.6,varwidth=TRUE,ylim=c(0,3),notch=TRUE)
zz = as.data.frame(matrix(nrow=0,ncol=3))
j=list()
for(qchr in qchrs) {
	div_ = div3[which(div3$sp==qsp & div3$tsp==tsp & div3$origin==origin & div3$chr==qchr & div3$chr2 %in% tchrs),]
	div_=div_[order(div_$start),]
	div_=div_[order(div_$chr),]
	j[[qchr]] = rle(div_$chr2)
	j[[qchr]] = as.data.frame(cbind(j[[qchr]]$values, j[[qchr]]$lengths))
	j[[qchr]][,2] = as.numeric(j[[qchr]][,2])
	cat(qsp,tsp,"max for",qchr,"=",max(j[[qchr]][,2]),"\n")
	zz = rbind(zz,cbind(rep(paste(qsp,qchr,sep=":"),length(j[[qchr]][,1])),j[[qchr]]))
}

qsp="ling_grc"
tsp="aaphi"
origin="cecidomyiidae"
qchrs=c("SUPER_GRC1","SUPER_GRC2")
tchrs=c("CM059996.1","CM059997.1","CM059998.1","CM059999.1")
div_ = div3[which(div3$sp==qsp & div3$tsp==tsp & div3$origin==origin & div3$chr %in% qchrs & div3$chr2 %in% tchrs),]
div_$factor = paste(gsub("SUPER_","",div_$chr)," vs ",gsub("SUPER_","",div_$chr2),sep="")
par(mar=c(10,4,4,4))
boxplot(dS~factor,data=div_,main=paste(qsp,"vs",tsp,"by origin=",origin),xlab="",ylab="dS",las=2,cex.names=.6,cex.main=.6,varwidth=TRUE,ylim=c(0,3),notch=TRUE)
j=list()
for(qchr in qchrs) {
	div_ = div3[which(div3$sp==qsp & div3$tsp==tsp & div3$origin==origin & div3$chr==qchr & div3$chr2 %in% tchrs),]
	div_=div_[order(div_$start),]
	div_=div_[order(div_$chr),]
	j[[qchr]] = rle(div_$chr2)
	j[[qchr]] = as.data.frame(cbind(j[[qchr]]$values, j[[qchr]]$lengths))
	j[[qchr]][,2] = as.numeric(j[[qchr]][,2])
	cat(qsp,tsp,"max for",qchr,"=",max(j[[qchr]][,2]),"\n")
	zz = rbind(zz,cbind(rep(paste(qsp,qchr,sep=":"),length(j[[qchr]][,1])),j[[qchr]]))
}

qsp="bimp_grc"
tsp="aaphi"
origin="cecidomyiidae"
qchrs=c("SUPER_GRC")
tchrs=c("CM059996.1","CM059997.1","CM059998.1","CM059999.1")
div_ = div3[which(div3$sp==qsp & div3$tsp==tsp & div3$origin==origin & div3$chr %in% qchrs & div3$chr2 %in% tchrs),]
div_$factor = paste(gsub("SUPER_","",div_$chr)," vs ",gsub("SUPER_","",div_$chr2),sep="")
par(mar=c(10,4,4,4))
boxplot(dS~factor,data=div_,main=paste(qsp,"vs",tsp,"by origin=",origin),xlab="",ylab="dS",las=2,cex.names=.6,cex.main=.6,varwidth=TRUE,ylim=c(0,3),notch=TRUE)
j=list()
for(qchr in qchrs) {
	div_ = div3[which(div3$sp==qsp & div3$tsp==tsp & div3$origin==origin & div3$chr==qchr & div3$chr2 %in% tchrs),]
	div_=div_[order(div_$start),]
	div_=div_[order(div_$chr),]
	j[[qchr]] = rle(div_$chr2)
	j[[qchr]] = as.data.frame(cbind(j[[qchr]]$values, j[[qchr]]$lengths))
	j[[qchr]][,2] = as.numeric(j[[qchr]][,2])
	cat(qsp,tsp,"max for",qchr,"=",max(j[[qchr]][,2]),"\n")
	zz = rbind(zz,cbind(rep(paste(qsp,qchr,sep=":"),length(j[[qchr]][,1])),j[[qchr]]))
}

colnames(zz) = c("sp_grc","target_chr","seq_len_rep")
par(mar=c(15,4,4,4))
par(mfrow=c(1,1))
boxplot(seq_len_rep~sp_grc,data=zz,las=2,xlab="",ylab="Sequential repetitions",cex.names=.6,ylim=c(0,70))

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
  p <- ggplot(zz, aes(x=sp_grc, y=seq_len_rep, fill=sp_grc)) + # fill=name allow to automatically dedicate a color for each group
       geom_violin()


# Now with orobi:
par(mfrow=c(1,3))
qsp="bcop_grc"
tsp="orobi"
origin="cecidomyiidae"
qchrs=c("SUPER_GRC1","SUPER_GRC2")
tchrs = c("CM052014.1","CM052015.1","CM052016.1","CM052017.1")
div_ = div3[which(div3$sp==qsp & div3$tsp==tsp & div3$origin==origin & div3$chr %in% qchrs & div3$chr2 %in% tchrs),]
div_$factor = paste(gsub("SUPER_","",div_$chr)," vs ",gsub("SUPER_","",div_$chr2),sep="")
par(mar=c(10,4,4,4))
boxplot(dS~factor,data=div_,main=paste(qsp,"vs",tsp,"by origin=",origin),xlab="",ylab="dS",las=2,cex.names=.6,cex.main=.6,varwidth=TRUE)

qsp="ling_grc"
tsp="orobi"
origin="cecidomyiidae"
qchrs=c("SUPER_GRC1","SUPER_GRC2")
tchrs = c("CM052014.1","CM052015.1","CM052016.1","CM052017.1")
div_ = div3[which(div3$sp==qsp & div3$tsp==tsp & div3$origin==origin & div3$chr %in% qchrs & div3$chr2 %in% tchrs),]
div_$factor = paste(gsub("SUPER_","",div_$chr)," vs ",gsub("SUPER_","",div_$chr2),sep="")
par(mar=c(10,4,4,4))
boxplot(dS~factor,data=div_,main=paste(qsp,"vs",tsp,"by origin=",origin),xlab="",ylab="dS",las=2,cex.names=.6,cex.main=.6,varwidth=TRUE)

qsp="bimp_grc"
tsp="orobi"
origin="cecidomyiidae"
qchrs=c("SUPER_GRC")
tchrs = c("CM052014.1","CM052015.1","CM052016.1","CM052017.1")
div_ = div3[which(div3$sp==qsp & div3$tsp==tsp & div3$origin==origin & div3$chr %in% qchrs & div3$chr2 %in% tchrs),]
div_$factor = paste(gsub("SUPER_","",div_$chr)," vs ",gsub("SUPER_","",div_$chr2),sep="")
par(mar=c(10,4,4,4))
boxplot(dS~factor,data=div_,main=paste(qsp,"vs",tsp,"by origin=",origin),xlab="",ylab="dS",las=2,cex.names=.6,cex.main=.6,varwidth=TRUE)



# Compare divergences between bcop grc and all the other cecidomyiidae:
# Then stratify by chromosome perhaps
par(mfrow=c(2,3))
qsp = "bcop_grc"
origin = "cecidomyiidae"
tsps = c("orobi","aaphi","contarinia")
qchr = "SUPER_GRC1"
div_ = div2[which(div2$sp==qsp & div2$tsp %in% tsps & div2$origin==origin & div2$chr==qchr),]
boxplot(dS~tsp,data=div_,ylim=c(0,3),main=paste(qsp,qchr),notch=TRUE)
qsp = "bcop_grc"
origin = "cecidomyiidae"
tsps = c("orobi","aaphi","contarinia")
qchr = "SUPER_GRC2"
div_ = div2[which(div2$sp==qsp & div2$tsp %in% tsps & div2$origin==origin & div2$chr==qchr),]
boxplot(dS~tsp,data=div_,ylim=c(0,3),main=paste(qsp,qchr),notch=TRUE)
qsp = "ling_grc"
origin = "cecidomyiidae"
tsps = c("orobi","aaphi","contarinia")
qchr = "SUPER_GRC1"
div_ = div2[which(div2$sp==qsp & div2$tsp %in% tsps & div2$origin==origin & div2$chr==qchr),]
boxplot(dS~tsp,data=div_,ylim=c(0,3),main=paste(qsp,qchr),notch=TRUE)
qsp = "ling_grc"
origin = "cecidomyiidae"
tsps = c("orobi","aaphi","contarinia")
qchr = "SUPER_GRC2"
div_ = div2[which(div2$sp==qsp & div2$tsp %in% tsps & div2$origin==origin & div2$chr==qchr),]
boxplot(dS~tsp,data=div_,ylim=c(0,3),main=paste(qsp,qchr),notch=TRUE)
qsp = "bimp_grc"
origin = "cecidomyiidae"
tsps = c("orobi","aaphi","contarinia")
qchr = "SUPER_GRC"
div_ = div2[which(div2$sp==qsp & div2$tsp %in% tsps & div2$origin==origin & div2$chr==qchr),]
boxplot(dS~tsp,data=div_,ylim=c(0,3),main=paste(qsp,qchr),notch=TRUE)





qsp = "bcop_grc"
origin = "cecidomyiidae"
tsps = c("orobi","aaphi","contarinia")
qchr = "SUPER_GRC1"
div_ = div2[which(div2$sp==qsp & div2$tsp %in% tsps & div2$origin==origin & div2$chr==qchr),]
counts = vector()
counts[tsps] = 0
for(gene in unique(div_$gene.x)) {
	mm = div_[which(div_$gene.x==gene),]
	mm = mm[order(mm$dS,decreasing=F),]
	best = mm[1,"tsp"]
	counts[best] = counts[best]+1
}
cat("--------------------------------------\n")
cat(qsp,qchr,"\n")
counts

qsp = "bcop_grc"
origin = "cecidomyiidae"
tsps = c("orobi","aaphi","contarinia")
qchr = "SUPER_GRC2"
div_ = div2[which(div2$sp==qsp & div2$tsp %in% tsps & div2$origin==origin & div2$chr==qchr),]
counts = vector()
counts[tsps] = 0
for(gene in unique(div_$gene.x)) {
	mm = div_[which(div_$gene.x==gene),]
	mm = mm[order(mm$dS,decreasing=F),]
	best = mm[1,"tsp"]
	counts[best] = counts[best]+1
}
cat("--------------------------------------\n")
cat(qsp,qchr,"\n")
counts

qsp = "ling_grc"
origin = "cecidomyiidae"
tsps = c("orobi","aaphi","contarinia")
qchr = "SUPER_GRC1"
div_ = div2[which(div2$sp==qsp & div2$tsp %in% tsps & div2$origin==origin & div2$chr==qchr),]
counts = vector()
counts[tsps] = 0
for(gene in unique(div_$gene.x)) {
	mm = div_[which(div_$gene.x==gene),]
	mm = mm[order(mm$dS,decreasing=F),]
	best = mm[1,"tsp"]
	counts[best] = counts[best]+1
}
cat("--------------------------------------\n")
cat(qsp,qchr,"\n")
counts

qsp = "ling_grc"
origin = "cecidomyiidae"
tsps = c("orobi","aaphi","contarinia")
qchr = "SUPER_GRC2"
div_ = div2[which(div2$sp==qsp & div2$tsp %in% tsps & div2$origin==origin & div2$chr==qchr),]
counts = vector()
counts[tsps] = 0
for(gene in unique(div_$gene.x)) {
	mm = div_[which(div_$gene.x==gene),]
	mm = mm[order(mm$dS,decreasing=F),]
	best = mm[1,"tsp"]
	counts[best] = counts[best]+1
}
cat("--------------------------------------\n")
cat(qsp,qchr,"\n")
counts

qsp = "bimp_grc"
origin = "cecidomyiidae"
tsps = c("orobi","aaphi","contarinia")
qchr = "SUPER_GRC"
div_ = div2[which(div2$sp==qsp & div2$tsp %in% tsps & div2$origin==origin & div2$chr==qchr),]
counts = vector()
counts[tsps] = 0
for(gene in unique(div_$gene.x)) {
	mm = div_[which(div_$gene.x==gene),]
	mm = mm[order(mm$dS,decreasing=F),]
	best = mm[1,"tsp"]
	counts[best] = counts[best]+1
}
cat("--------------------------------------\n")
cat(qsp,qchr,"\n")
counts



bcop_grc SUPER_GRC1 
> counts
     orobi      aaphi contarinia 
        87        897        107 
bcop_grc SUPER_GRC2 
> counts
     orobi      aaphi contarinia 
        87        537        112 
ling_grc SUPER_GRC1 
> counts
     orobi      aaphi contarinia 
        56        208         65 
ling_grc SUPER_GRC2 
> counts
     orobi      aaphi contarinia 
        27        608         44 
bimp_grc SUPER_GRC 
> counts
     orobi      aaphi contarinia 
        37        129         54 







