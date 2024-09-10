setwd("~/Desktop/GRCs")

origins = read.table("ORIGINS_ONTO_CHRS.tsv",sep="\t",header=FALSE)
colnames(origins) = c("sp","chr","start","end","gene","origin")


origins$spchr = paste(origins$sp,origins$chr,sep="_")
origins = origins[grep("GRC",origins$chr),]
origins$id = paste(origins$sp,origins$gene,sep=":")

for(sp in c("bcop","bimp","ling")) {
	file = paste(sp,"_grc.gene.bed",sep="")
	kk = read.table(file,sep="\t",header=FALSE)
	colnames(kk) = c("chr","start","end","gene")
	kk$id = paste(sp,":",kk$gene,".t1",sep="")
	kk = kk[-which(kk$id %in% origins$id),]
	kk_ = data.frame(sp=sp,chr=kk$chr,start=kk$start,end=kk$end,gene=paste(kk$gene,".t1",sep=""),origin="Uncertain",spchr=paste(sp,kk$chr,sep="_"),id=kk$id)
	origins = rbind(origins,kk_)
}

origins = origins[grep("unloc",origins$chr,invert=T),]

origins$colour = ""
origins[which(origins$origin == "cecidomyiidae"),"colour"] = "darkcyan"
origins[which(origins$origin == "sciaridae"),"colour"] = "firebrick"
origins[which(origins$origin == "Uncertain"),"colour"] = "gray"

par(mfrow=c(1,1))
par(mar=c(4,4,4,4))
origins = origins[sample(1:nrow(origins)),]
#origins = origins[order(origins$origin,decreasing=T),]

#for(spchr in sort(unique(origins$spchr))) {
#	kk = origins[which(origins$spchr==spchr),]
#	plot(kk$pos,rep(2,nrow(kk)),main=spchr,xlab="",ylab="",pch=15,col=kk$colour,yaxt="n",cex=2,bty="n")
#	legend("topleft",legend=c("cecidomyiidae origin","sciaridae origin"),pch=15,col=c("darkcyan","firebrick"))
#}
plot(NULL,NULL,yaxt="n",bty="n",xlim=c(0,max(origins$end)+500000),ylim=c(0,7),ylab="",xlab="Chrom. position",main="Gene origins")
yat = 1
for(spchr in sort(unique(origins$spchr))) {
	kk = origins[which(origins$spchr==spchr),]
	#plot(NULL,NULL,main=spchr,xlab="",ylab="",pch=15,col=kk$colour,yaxt="n",cex=2,bty="n",xlim=c(0,max(kk$end)),ylim=c(0,2))
	#segments(kk$start,1.2,kk$end,0.8,col=kk$colour)
	rect(kk$start, yat-0.2, kk$end, yat + 0.2, density = NULL, angle = 45, col = kk$colour, border = kk$colour,lwd=.2)
	text(0,yat+.5,spchr,pos=4)
	yat = yat + 1
}
legend("topleft",legend=c("cecidomyiidae origin","sciaridae origin","unknown"),pch=15,col=c("darkcyan","firebrick","gray"))





