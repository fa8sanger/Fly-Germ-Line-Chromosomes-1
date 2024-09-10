setwd("~/Desktop/GRCs/")

tt = read.table("OG_BASED_SYNTENY_TABLE.tsv",sep="\t",header=FALSE)
colnames(tt) = c("og","genome1","gene1","chr1","start1","end1","genome2","gene2","chr2","start2","end2")

#sort(table(tt$og),decreasing=T)[1:20]
og_counts = table(tt$og)
keep_these = names(og_counts[which(og_counts < 1000)])

tt = tt[which(tt$og %in% keep_these),]

origins = read.table("ORIGINS_ONTO_CHRS.tsv",sep="\t",header=FALSE)
colnames(origins) = c("sp","chr","start","end","gene","origin")


origins$spchr = paste(origins$sp,origins$chr,sep="_")
origins = origins[grep("GRC",origins$chr),]
origins$id = paste(paste(origins$sp,"_grc",sep=""),origins$gene,sep=":")
cecido_ids = origins[which(origins$origin=="cecidomyiidae"),"id"]
sciara_ids = origins[which(origins$origin=="sciaridae"),"id"]

cecido_ids = gsub("\\.t1","",cecido_ids)
sciara_ids = gsub("\\.t1","",sciara_ids)


par(mfrow=c(1,2))
from_genome = "bimp_grc"
from_chr    = "SUPER_GRC"
to_genome   = "bimp_core"
to_chr      = "SUPER_1"

jj = tt[which(tt$genome1==from_genome & tt$chr1==from_chr & tt$genome2==to_genome & tt$chr2==to_chr),]
plot(NULL,NULL,xlim=c(0,max(jj$end1)),ylim=c(0,max(jj$end2)),xlab=paste(from_genome,from_chr),ylab=paste(to_genome,to_chr))
segments(jj$start1,jj$start2,jj$end1,jj$end2,col="gray")
jj$id1 = paste(jj$genome1,jj$gene1,sep=":")
jj$id2 = paste(jj$genome2,jj$gene2,sep=":")
jj_s = jj[which(jj$id1 %in% sciara_ids | jj$id2 %in% sciara_ids),]
segments(jj_s$start1,jj_s$start2,jj_s$end1,jj_s$end2,col="firebrick",lwd=2)
jj_c = jj[which(jj$id1 %in% cecido_ids | jj$id2 %in% cecido_ids),]
segments(jj_c$start1,jj_c$start2,jj_c$end1,jj_c$end2,col="darkcyan",lwd=2)

from_genome = "bimp_grc"
from_chr    = "SUPER_GRC"
to_genome   = "bimp_core"
to_chr      = "SUPER_X"

jj = tt[which(tt$genome1==from_genome & tt$chr1==from_chr & tt$genome2==to_genome & tt$chr2==to_chr),]
plot(NULL,NULL,xlim=c(0,max(jj$end1)),ylim=c(0,max(jj$end2)),xlab=paste(from_genome,from_chr),ylab=paste(to_genome,to_chr))
segments(jj$start1,jj$start2,jj$end1,jj$end2,col="orange")
jj$id1 = paste(jj$genome1,jj$gene1,sep=":")
jj$id2 = paste(jj$genome2,jj$gene2,sep=":")
jj_s = jj[which(jj$id1 %in% sciara_ids | jj$id2 %in% sciara_ids),]
segments(jj_s$start1,jj_s$start2,jj_s$end1,jj_s$end2,col="firebrick",lwd=2)
jj_c = jj[which(jj$id1 %in% cecido_ids | jj$id2 %in% cecido_ids),]
segments(jj_c$start1,jj_c$start2,jj_c$end1,jj_c$end2,col="darkcyan",lwd=2)

# GRC1 vs GRC2:
par(mfrow=c(1,3))
from_genome = "bcop_grc"
from_chr    = "SUPER_GRC1"
to_genome   = "ling_grc"
to_chr      = "SUPER_GRC2"

jj = tt[which(tt$genome1==from_genome & tt$chr1==from_chr & tt$genome2==to_genome & tt$chr2==to_chr),]
plot(NULL,NULL,xlim=c(0,max(jj$end1)),ylim=c(0,max(jj$end2)),xlab=paste(from_genome,from_chr),ylab=paste(to_genome,to_chr))
segments(jj$start1,jj$start2,jj$end1,jj$end2,col="gray")
jj$id1 = paste(jj$genome1,jj$gene1,sep=":")
jj$id2 = paste(jj$genome2,jj$gene2,sep=":")
jj_s = jj[which(jj$id1 %in% sciara_ids | jj$id2 %in% sciara_ids),]
segments(jj_s$start1,jj_s$start2,jj_s$end1,jj_s$end2,col="firebrick",lwd=2)
jj_c = jj[which(jj$id1 %in% cecido_ids | jj$id2 %in% cecido_ids),]
segments(jj_c$start1,jj_c$start2,jj_c$end1,jj_c$end2,col="darkcyan",lwd=2)

from_genome = "bcop_grc"
from_chr    = "SUPER_GRC1"
to_genome   = "bcop_core"
to_chr      = "SUPER_X"

jj = tt[which(tt$genome1==from_genome & tt$chr1==from_chr & tt$genome2==to_genome & tt$chr2==to_chr),]
plot(NULL,NULL,xlim=c(0,max(jj$end1)),ylim=c(0,max(jj$end2)),xlab=paste(from_genome,from_chr),ylab=paste(to_genome,to_chr))
segments(jj$start1,jj$start2,jj$end1,jj$end2,col="orange")
jj$id1 = paste(jj$genome1,jj$gene1,sep=":")
jj$id2 = paste(jj$genome2,jj$gene2,sep=":")
jj_s = jj[which(jj$id1 %in% sciara_ids | jj$id2 %in% sciara_ids),]
segments(jj_s$start1,jj_s$start2,jj_s$end1,jj_s$end2,col="firebrick",lwd=2)
jj_c = jj[which(jj$id1 %in% cecido_ids | jj$id2 %in% cecido_ids),]
segments(jj_c$start1,jj_c$start2,jj_c$end1,jj_c$end2,col="darkcyan",lwd=2)

from_genome = "bcop_grc"
from_chr    = "SUPER_GRC1"
to_genome   = "bcop_core"
to_chr      = "SUPER_1"

jj = tt[which(tt$genome1==from_genome & tt$chr1==from_chr & tt$genome2==to_genome & tt$chr2==to_chr),]
plot(NULL,NULL,xlim=c(0,max(jj$end1)),ylim=c(0,max(jj$end2)),xlab=paste(from_genome,from_chr),ylab=paste(to_genome,to_chr))
segments(jj$start1,jj$start2,jj$end1,jj$end2,col="orange")
jj$id1 = paste(jj$genome1,jj$gene1,sep=":")
jj$id2 = paste(jj$genome2,jj$gene2,sep=":")
jj_s = jj[which(jj$id1 %in% sciara_ids | jj$id2 %in% sciara_ids),]
segments(jj_s$start1,jj_s$start2,jj_s$end1,jj_s$end2,col="firebrick",lwd=2)
jj_c = jj[which(jj$id1 %in% cecido_ids | jj$id2 %in% cecido_ids),]
segments(jj_c$start1,jj_c$start2,jj_c$end1,jj_c$end2,col="darkcyan",lwd=2)



##########################################################################################
# Microsyntenty and gene loss plot
##########################################################################################
from_genome = "bimp_grc"
from_chr    = "SUPER_GRC"
from_start  = 2.8e6
from_end    = 3e6
to_genome   = "bimp_core"
to_chr      = "SUPER_1"
to_start    = 2.8e6
to_end      = 3e6

from_genes = tt[which(tt$chr1==from_chr & tt$start1>=from_start & tt$end1<=from_end),c(4:6 )]
to_genes   = tt[which(tt$chr2==to_chr   & tt$start2>=to_start   & tt$end2<=to_end  ),c(9:11)]
colnames(from_genes) = colnames(to_genes) = c("chr","start","end")

plot(NULL,NULL,xlim=c(from_start-2000,from_end),ylim=c(to_start-2000,to_end),
     xlab=paste(from_genome,": ",from_chr,":",from_start,"-",from_end,sep=""),
     ylab=paste(to_genome,  ": ",to_chr,  ":",to_start,  "-",to_end,sep=""))

rect(from_genes[,"start"],to_start-1500, from_genes[,"end"], to_start+1500, col="black")
rect(from_start-1500,to_genes[,"start"], from_start+1500, to_genes[,"end"], col="black")

jj = tt[which(tt$genome1==from_genome & tt$chr1==from_chr & tt$genome2==to_genome & 
              tt$chr2==to_chr & tt$start1>=from_start & tt$end1<=from_end & tt$start2>=to_start & 
              tt$end2<=to_end),]

# Find genes not in an OG with the other species and flag them
from_genes_with_ortology = tt[which(tt$chr1==from_chr & tt$start1>=from_start & tt$end1<=from_end & tt$genome2==to_genome  ),c(4:6 )]
to_genes_with_ortology   = tt[which(tt$chr2==to_chr   & tt$start2>=to_start   & tt$end2<=to_end   & tt$genome1==from_genome),c(4:6 )]
colnames(from_genes_with_ortology) = colnames(to_genes_with_ortology) = c("chr","start","end")

rect(from_genes_with_ortology[,"start"],to_start-1500, from_genes_with_ortology[,"end"], to_start+1500, col="red")
rect(from_start-1500,to_genes_with_ortology[,"start"], from_start+1500, to_genes_with_ortology[,"end"], col="red")



segments(jj$start1,jj$start2,jj$end1,jj$end2,col="gray")
jj$id1 = paste(jj$genome1,jj$gene1,sep=":")
jj$id2 = paste(jj$genome2,jj$gene2,sep=":")
jj_s = jj[which(jj$id1 %in% sciara_ids | jj$id2 %in% sciara_ids),]
segments(jj_s$start1,jj_s$start2,jj_s$end1,jj_s$end2,col="firebrick",lwd=2)
jj_c = jj[which(jj$id1 %in% cecido_ids | jj$id2 %in% cecido_ids),]
segments(jj_c$start1,jj_c$start2,jj_c$end1,jj_c$end2,col="darkcyan",lwd=2)

legend("topright",legend=c("with orthology","without_orthology"),col=c("red","black"),pch=15)



##########################################################################################
# Another segment:
from_genome = "bimp_grc"
from_chr    = "SUPER_GRC"
from_start  = 1.08e7
from_end    = 1.1e7
to_genome   = "bimp_core"
to_chr      = "SUPER_1"
to_start    = 1.6e7
to_end      = 1.63e7

from_genes = tt[which(tt$chr1==from_chr & tt$start1>=from_start & tt$end1<=from_end),c(4:6 )]
to_genes   = tt[which(tt$chr2==to_chr   & tt$start2>=to_start   & tt$end2<=to_end  ),c(9:11)]
colnames(from_genes) = colnames(to_genes) = c("chr","start","end")

plot(NULL,NULL,xlim=c(from_start-2000,from_end),ylim=c(to_start-2000,to_end),
     xlab=paste(from_genome,": ",from_chr,":",from_start,"-",from_end,sep=""),
     ylab=paste(to_genome,  ": ",to_chr,  ":",to_start,  "-",to_end,sep=""))

rect(from_genes[,"start"],to_start-1500, from_genes[,"end"], to_start+1500, col="black")
rect(from_start-1500,to_genes[,"start"], from_start+1500, to_genes[,"end"], col="black")

jj = tt[which(tt$genome1==from_genome & tt$chr1==from_chr & tt$genome2==to_genome & 
              tt$chr2==to_chr & tt$start1>=from_start & tt$end1<=from_end & tt$start2>=to_start & 
              tt$end2<=to_end),]

# Find genes not in an OG with the other species and flag them
from_genes_with_ortology = tt[which(tt$chr1==from_chr & tt$start1>=from_start & tt$end1<=from_end & tt$genome2==to_genome  ),c(4:6 )]
to_genes_with_ortology   = tt[which(tt$chr2==to_chr   & tt$start2>=to_start   & tt$end2<=to_end   & tt$genome1==from_genome),c(9:11)]
colnames(from_genes_with_ortology) = colnames(to_genes_with_ortology) = c("chr","start","end")

# from_genes_with_ortology = tt[which(tt$genome1==from_genome & tt$chr1==from_chr & tt$genome2==to_genome & 
#               tt$chr2==to_chr & tt$start1>=from_start & tt$end1<=from_end & tt$start2>=to_start & 
#               tt$end2<=to_end),c(4:6 )]
# to_genes_with_ortology   = tt[which(tt$genome1==from_genome & tt$chr1==from_chr & tt$genome2==to_genome & 
#               tt$chr2==to_chr & tt$start1>=from_start & tt$end1<=from_end & tt$start2>=to_start & 
#               tt$end2<=to_end),c(9:11)]
# colnames(from_genes_with_ortology) = colnames(to_genes_with_ortology) = c("chr","start","end")

rect(from_genes_with_ortology[,"start"],to_start-1500, from_genes_with_ortology[,"end"], to_start+1500, col="red",border="red")
rect(from_start-1500,to_genes_with_ortology[,"start"], from_start+1500, to_genes_with_ortology[,"end"], col="red",border="red")

segments(jj$start1,jj$start2,jj$end1,jj$end2,col="gray")
jj$id1 = paste(jj$genome1,jj$gene1,sep=":")
jj$id2 = paste(jj$genome2,jj$gene2,sep=":")
jj_s = jj[which(jj$id1 %in% sciara_ids | jj$id2 %in% sciara_ids),]
segments(jj_s$start1,jj_s$start2,jj_s$end1,jj_s$end2,col="firebrick",lwd=2)
jj_c = jj[which(jj$id1 %in% cecido_ids | jj$id2 %in% cecido_ids),]
segments(jj_c$start1,jj_c$start2,jj_c$end1,jj_c$end2,col="darkcyan",lwd=2)

legend("topright",legend=c("with orthology","without_orthology"),col=c("red","black"),pch=15)



##########################################################################################
##########################################################################################

par(mfrow=c(2,2))
from_genome = "bcop_grc"
from_chr    = "SUPER_GRC1"
to_genome   = "ling_grc"
to_chr      = "SUPER_GRC2"

jj = tt[which(tt$genome1==from_genome & tt$chr1==from_chr & tt$genome2==to_genome & tt$chr2==to_chr),]
plot(NULL,NULL,xlim=c(0,max(jj$end1)),ylim=c(0,max(jj$end2)),xlab=paste(from_genome,from_chr),ylab=paste(to_genome,to_chr))
segments(jj$start1,jj$start2,jj$end1,jj$end2,col="orange")

from_genome = "bcop_grc"
from_chr    = "SUPER_GRC1"
to_genome   = "ling_grc"
to_chr      = "SUPER_GRC1"

jj = tt[which(tt$genome1==from_genome & tt$chr1==from_chr & tt$genome2==to_genome & tt$chr2==to_chr),]
plot(NULL,NULL,xlim=c(0,max(jj$end1)),ylim=c(0,max(jj$end2)),xlab=paste(from_genome,from_chr),ylab=paste(to_genome,to_chr))
segments(jj$start1,jj$start2,jj$end1,jj$end2,col="orange")

from_genome = "bcop_grc"
from_chr    = "SUPER_GRC2"
to_genome   = "ling_grc"
to_chr      = "SUPER_GRC2"

jj = tt[which(tt$genome1==from_genome & tt$chr1==from_chr & tt$genome2==to_genome & tt$chr2==to_chr),]
plot(NULL,NULL,xlim=c(0,max(jj$end1)),ylim=c(0,max(jj$end2)),xlab=paste(from_genome,from_chr),ylab=paste(to_genome,to_chr))
segments(jj$start1,jj$start2,jj$end1,jj$end2,col="orange")

from_genome = "bcop_grc"
from_chr    = "SUPER_GRC2"
to_genome   = "ling_grc"
to_chr      = "SUPER_GRC1"

jj = tt[which(tt$genome1==from_genome & tt$chr1==from_chr & tt$genome2==to_genome & tt$chr2==to_chr),]
plot(NULL,NULL,xlim=c(0,max(jj$end1)),ylim=c(0,max(jj$end2)),xlab=paste(from_genome,from_chr),ylab=paste(to_genome,to_chr))
segments(jj$start1,jj$start2,jj$end1,jj$end2,col="orange")




# Now ling GRC2 vs aaphi
par(mfrow=c(2,2))
from_genome = "ling_grc"
from_chr    = "SUPER_GRC2"
to_genome   = "aaphi"
to_chr      = "CM059996.1"

jj = tt[which(tt$genome1==from_genome & tt$chr1==from_chr & tt$genome2==to_genome & tt$chr2==to_chr),]
plot(NULL,NULL,xlim=c(0,max(jj$end1)),ylim=c(0,max(jj$end2)),xlab=paste(from_genome,from_chr),ylab=paste(to_genome,to_chr))
segments(jj$start1,jj$start2,jj$end1,jj$end2,col="orange")

from_genome = "ling_grc"
from_chr    = "SUPER_GRC2"
to_genome   = "aaphi"
to_chr      = "CM059997.1"

jj = tt[which(tt$genome1==from_genome & tt$chr1==from_chr & tt$genome2==to_genome & tt$chr2==to_chr),]
plot(NULL,NULL,xlim=c(0,max(jj$end1)),ylim=c(0,max(jj$end2)),xlab=paste(from_genome,from_chr),ylab=paste(to_genome,to_chr))
segments(jj$start1,jj$start2,jj$end1,jj$end2,col="orange")

from_genome = "ling_grc"
from_chr    = "SUPER_GRC2"
to_genome   = "aaphi"
to_chr      = "CM059998.1"

jj = tt[which(tt$genome1==from_genome & tt$chr1==from_chr & tt$genome2==to_genome & tt$chr2==to_chr),]
plot(NULL,NULL,xlim=c(0,max(jj$end1)),ylim=c(0,max(jj$end2)),xlab=paste(from_genome,from_chr),ylab=paste(to_genome,to_chr))
segments(jj$start1,jj$start2,jj$end1,jj$end2,col="orange")

from_genome = "ling_grc"
from_chr    = "SUPER_GRC2"
to_genome   = "aaphi"
to_chr      = "CM059999.1"

jj = tt[which(tt$genome1==from_genome & tt$chr1==from_chr & tt$genome2==to_genome & tt$chr2==to_chr),]
plot(NULL,NULL,xlim=c(0,max(jj$end1)),ylim=c(0,max(jj$end2)),xlab=paste(from_genome,from_chr),ylab=paste(to_genome,to_chr))
segments(jj$start1,jj$start2,jj$end1,jj$end2,col="orange")



quartz()
par(mfrow=c(2,2))
from_genome = "bcop_grc"
from_chr    = "SUPER_GRC1"
to_genome   = "aaphi"
to_chr      = "CM059996.1"

jj = tt[which(tt$genome1==from_genome & tt$chr1==from_chr & tt$genome2==to_genome & tt$chr2==to_chr),]
plot(NULL,NULL,xlim=c(0,max(jj$end1)),ylim=c(0,max(jj$end2)),xlab=paste(from_genome,from_chr),ylab=paste(to_genome,to_chr))
segments(jj$start1,jj$start2,jj$end1,jj$end2,col="orange")

from_genome = "bcop_grc"
from_chr    = "SUPER_GRC1"
to_genome   = "aaphi"
to_chr      = "CM059997.1"

jj = tt[which(tt$genome1==from_genome & tt$chr1==from_chr & tt$genome2==to_genome & tt$chr2==to_chr),]
plot(NULL,NULL,xlim=c(0,max(jj$end1)),ylim=c(0,max(jj$end2)),xlab=paste(from_genome,from_chr),ylab=paste(to_genome,to_chr))
segments(jj$start1,jj$start2,jj$end1,jj$end2,col="orange")

from_genome = "bcop_grc"
from_chr    = "SUPER_GRC1"
to_genome   = "aaphi"
to_chr      = "CM059998.1"

jj = tt[which(tt$genome1==from_genome & tt$chr1==from_chr & tt$genome2==to_genome & tt$chr2==to_chr),]
plot(NULL,NULL,xlim=c(0,max(jj$end1)),ylim=c(0,max(jj$end2)),xlab=paste(from_genome,from_chr),ylab=paste(to_genome,to_chr))
segments(jj$start1,jj$start2,jj$end1,jj$end2,col="orange")

from_genome = "bcop_grc"
from_chr    = "SUPER_GRC1"
to_genome   = "aaphi"
to_chr      = "CM059999.1"

jj = tt[which(tt$genome1==from_genome & tt$chr1==from_chr & tt$genome2==to_genome & tt$chr2==to_chr),]
plot(NULL,NULL,xlim=c(0,max(jj$end1)),ylim=c(0,max(jj$end2)),xlab=paste(from_genome,from_chr),ylab=paste(to_genome,to_chr))
segments(jj$start1,jj$start2,jj$end1,jj$end2,col="orange")



quartz()
par(mfrow=c(2,2))
from_genome = "bcop_grc"
from_chr    = "SUPER_GRC2"
to_genome   = "aaphi"
to_chr      = "CM059996.1"

jj = tt[which(tt$genome1==from_genome & tt$chr1==from_chr & tt$genome2==to_genome & tt$chr2==to_chr),]
plot(NULL,NULL,xlim=c(0,max(jj$end1)),ylim=c(0,max(jj$end2)),xlab=paste(from_genome,from_chr),ylab=paste(to_genome,to_chr))
segments(jj$start1,jj$start2,jj$end1,jj$end2,col="orange")

from_genome = "bcop_grc"
from_chr    = "SUPER_GRC2"
to_genome   = "aaphi"
to_chr      = "CM059997.1"

jj = tt[which(tt$genome1==from_genome & tt$chr1==from_chr & tt$genome2==to_genome & tt$chr2==to_chr),]
plot(NULL,NULL,xlim=c(0,max(jj$end1)),ylim=c(0,max(jj$end2)),xlab=paste(from_genome,from_chr),ylab=paste(to_genome,to_chr))
segments(jj$start1,jj$start2,jj$end1,jj$end2,col="orange")

from_genome = "bcop_grc"
from_chr    = "SUPER_GRC2"
to_genome   = "aaphi"
to_chr      = "CM059998.1"

jj = tt[which(tt$genome1==from_genome & tt$chr1==from_chr & tt$genome2==to_genome & tt$chr2==to_chr),]
plot(NULL,NULL,xlim=c(0,max(jj$end1)),ylim=c(0,max(jj$end2)),xlab=paste(from_genome,from_chr),ylab=paste(to_genome,to_chr))
segments(jj$start1,jj$start2,jj$end1,jj$end2,col="orange")

from_genome = "bcop_grc"
from_chr    = "SUPER_GRC2"
to_genome   = "aaphi"
to_chr      = "CM059999.1"

jj = tt[which(tt$genome1==from_genome & tt$chr1==from_chr & tt$genome2==to_genome & tt$chr2==to_chr),]
plot(NULL,NULL,xlim=c(0,max(jj$end1)),ylim=c(0,max(jj$end2)),xlab=paste(from_genome,from_chr),ylab=paste(to_genome,to_chr))
segments(jj$start1,jj$start2,jj$end1,jj$end2,col="orange")




# Now ling GRC2 vs aaphi
par(mfrow=c(2,2))
from_genome = "ling_grc"
from_chr    = "SUPER_GRC1"
to_genome   = "aaphi"
to_chr      = "CM059996.1"

jj = tt[which(tt$genome1==from_genome & tt$chr1==from_chr & tt$genome2==to_genome & tt$chr2==to_chr),]
plot(NULL,NULL,xlim=c(0,max(jj$end1)),ylim=c(0,max(jj$end2)),xlab=paste(from_genome,from_chr),ylab=paste(to_genome,to_chr))
segments(jj$start1,jj$start2,jj$end1,jj$end2,col="orange")

from_genome = "ling_grc"
from_chr    = "SUPER_GRC1"
to_genome   = "aaphi"
to_chr      = "CM059997.1"

jj = tt[which(tt$genome1==from_genome & tt$chr1==from_chr & tt$genome2==to_genome & tt$chr2==to_chr),]
plot(NULL,NULL,xlim=c(0,max(jj$end1)),ylim=c(0,max(jj$end2)),xlab=paste(from_genome,from_chr),ylab=paste(to_genome,to_chr))
segments(jj$start1,jj$start2,jj$end1,jj$end2,col="orange")

from_genome = "ling_grc"
from_chr    = "SUPER_GRC1"
to_genome   = "aaphi"
to_chr      = "CM059998.1"

jj = tt[which(tt$genome1==from_genome & tt$chr1==from_chr & tt$genome2==to_genome & tt$chr2==to_chr),]
plot(NULL,NULL,xlim=c(0,max(jj$end1)),ylim=c(0,max(jj$end2)),xlab=paste(from_genome,from_chr),ylab=paste(to_genome,to_chr))
segments(jj$start1,jj$start2,jj$end1,jj$end2,col="orange")

from_genome = "ling_grc"
from_chr    = "SUPER_GRC1"
to_genome   = "aaphi"
to_chr      = "CM059999.1"

jj = tt[which(tt$genome1==from_genome & tt$chr1==from_chr & tt$genome2==to_genome & tt$chr2==to_chr),]
plot(NULL,NULL,xlim=c(0,max(jj$end1)),ylim=c(0,max(jj$end2)),xlab=paste(from_genome,from_chr),ylab=paste(to_genome,to_chr))
segments(jj$start1,jj$start2,jj$end1,jj$end2,col="orange")


# ling_grc vs orobi
par(mfrow=c(2,2))
from_genome = "ling_grc"
from_chr    = "SUPER_GRC2"
to_genome   = "orobi"
to_chr      = "CM052014.1"

jj = tt[which(tt$genome1==from_genome & tt$chr1==from_chr & tt$genome2==to_genome & tt$chr2==to_chr),]
plot(NULL,NULL,xlim=c(0,max(jj$end1)),ylim=c(0,max(jj$end2)),xlab=paste(from_genome,from_chr),ylab=paste(to_genome,to_chr))
segments(jj$start1,jj$start2,jj$end1,jj$end2,col="orange")
jj$id1 = paste(jj$genome1,jj$gene1,sep=":")
jj$id2 = paste(jj$genome2,jj$gene2,sep=":")
jj_s = jj[which(jj$id1 %in% sciara_ids | jj$id2 %in% sciara_ids),]
segments(jj_s$start1,jj_s$start2,jj_s$end1,jj_s$end2,col="firebrick",lwd=2)
jj_c = jj[which(jj$id1 %in% cecido_ids | jj$id2 %in% cecido_ids),]
segments(jj_c$start1,jj_c$start2,jj_c$end1,jj_c$end2,col="darkcyan",lwd=2)

from_genome = "ling_grc"
from_chr    = "SUPER_GRC2"
to_genome   = "orobi"
to_chr      = "CM052015.1"

jj = tt[which(tt$genome1==from_genome & tt$chr1==from_chr & tt$genome2==to_genome & tt$chr2==to_chr),]
plot(NULL,NULL,xlim=c(0,max(jj$end1)),ylim=c(0,max(jj$end2)),xlab=paste(from_genome,from_chr),ylab=paste(to_genome,to_chr))
segments(jj$start1,jj$start2,jj$end1,jj$end2,col="orange")
jj$id1 = paste(jj$genome1,jj$gene1,sep=":")
jj$id2 = paste(jj$genome2,jj$gene2,sep=":")
jj_s = jj[which(jj$id1 %in% sciara_ids | jj$id2 %in% sciara_ids),]
segments(jj_s$start1,jj_s$start2,jj_s$end1,jj_s$end2,col="firebrick",lwd=2)
jj_c = jj[which(jj$id1 %in% cecido_ids | jj$id2 %in% cecido_ids),]
segments(jj_c$start1,jj_c$start2,jj_c$end1,jj_c$end2,col="darkcyan",lwd=2)

from_genome = "ling_grc"
from_chr    = "SUPER_GRC2"
to_genome   = "orobi"
to_chr      = "CM052016.1"

jj = tt[which(tt$genome1==from_genome & tt$chr1==from_chr & tt$genome2==to_genome & tt$chr2==to_chr),]
plot(NULL,NULL,xlim=c(0,max(jj$end1)),ylim=c(0,max(jj$end2)),xlab=paste(from_genome,from_chr),ylab=paste(to_genome,to_chr))
segments(jj$start1,jj$start2,jj$end1,jj$end2,col="orange")
jj$id1 = paste(jj$genome1,jj$gene1,sep=":")
jj$id2 = paste(jj$genome2,jj$gene2,sep=":")
jj_s = jj[which(jj$id1 %in% sciara_ids | jj$id2 %in% sciara_ids),]
segments(jj_s$start1,jj_s$start2,jj_s$end1,jj_s$end2,col="firebrick",lwd=2)
jj_c = jj[which(jj$id1 %in% cecido_ids | jj$id2 %in% cecido_ids),]
segments(jj_c$start1,jj_c$start2,jj_c$end1,jj_c$end2,col="darkcyan",lwd=2)

from_genome = "ling_grc"
from_chr    = "SUPER_GRC2"
to_genome   = "orobi"
to_chr      = "CM052017.1"

jj = tt[which(tt$genome1==from_genome & tt$chr1==from_chr & tt$genome2==to_genome & tt$chr2==to_chr),]
plot(NULL,NULL,xlim=c(0,max(jj$end1)),ylim=c(0,max(jj$end2)),xlab=paste(from_genome,from_chr),ylab=paste(to_genome,to_chr))
segments(jj$start1,jj$start2,jj$end1,jj$end2,col="orange")
jj$id1 = paste(jj$genome1,jj$gene1,sep=":")
jj$id2 = paste(jj$genome2,jj$gene2,sep=":")
jj_s = jj[which(jj$id1 %in% sciara_ids | jj$id2 %in% sciara_ids),]
segments(jj_s$start1,jj_s$start2,jj_s$end1,jj_s$end2,col="firebrick",lwd=2)
jj_c = jj[which(jj$id1 %in% cecido_ids | jj$id2 %in% cecido_ids),]
segments(jj_c$start1,jj_c$start2,jj_c$end1,jj_c$end2,col="darkcyan",lwd=2)





quartz()
par(mfrow=c(2,2))
from_genome = "bimp_grc"
from_chr    = "SUPER_GRC"
to_genome   = "aaphi"
to_chr      = "CM059996.1"

jj = tt[which(tt$genome1==from_genome & tt$chr1==from_chr & tt$genome2==to_genome & tt$chr2==to_chr),]
plot(NULL,NULL,xlim=c(0,max(jj$end1)),ylim=c(0,max(jj$end2)),xlab=paste(from_genome,from_chr),ylab=paste(to_genome,to_chr))
segments(jj$start1,jj$start2,jj$end1,jj$end2,col="orange")

from_genome = "bimp_grc"
from_chr    = "SUPER_GRC"
to_genome   = "aaphi"
to_chr      = "CM059997.1"

jj = tt[which(tt$genome1==from_genome & tt$chr1==from_chr & tt$genome2==to_genome & tt$chr2==to_chr),]
plot(NULL,NULL,xlim=c(0,max(jj$end1)),ylim=c(0,max(jj$end2)),xlab=paste(from_genome,from_chr),ylab=paste(to_genome,to_chr))
segments(jj$start1,jj$start2,jj$end1,jj$end2,col="orange")

from_genome = "bimp_grc"
from_chr    = "SUPER_GRC"
to_genome   = "aaphi"
to_chr      = "CM059998.1"

jj = tt[which(tt$genome1==from_genome & tt$chr1==from_chr & tt$genome2==to_genome & tt$chr2==to_chr),]
plot(NULL,NULL,xlim=c(0,max(jj$end1)),ylim=c(0,max(jj$end2)),xlab=paste(from_genome,from_chr),ylab=paste(to_genome,to_chr))
segments(jj$start1,jj$start2,jj$end1,jj$end2,col="orange")

from_genome = "bimp_grc"
from_chr    = "SUPER_GRC"
to_genome   = "aaphi"
to_chr      = "CM059999.1"

jj = tt[which(tt$genome1==from_genome & tt$chr1==from_chr & tt$genome2==to_genome & tt$chr2==to_chr),]
plot(NULL,NULL,xlim=c(0,max(jj$end1)),ylim=c(0,max(jj$end2)),xlab=paste(from_genome,from_chr),ylab=paste(to_genome,to_chr))
segments(jj$start1,jj$start2,jj$end1,jj$end2,col="orange")


