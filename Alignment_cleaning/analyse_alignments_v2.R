# Usage:
# cd /lustre/scratch126/tol/teams/jaron/users/fede/Orthofinder/Proteomes_2nd_round_OneIsoform/OrthoFinder/Results_Jul25_1/
# mkdir ALI_STATS
# Make a perl script that reads the orthogroup table and decides which alis to process
# and calls the analyse_alignments thingy
# Or I can do everything on a single R script, without a for loop. Won't take too long I think

# DO THIS FIRST
# cd /lustre/scratch126/tol/teams/jaron/users/fede/Orthofinder/Proteomes_2nd_round_OneIsoform/OrthoFinder/Results_Jul25_1/MultipleSequenceAlignments
# for a in OG*fa
#  do perl -pi -e 's/^ +//g' $a
# done
# 
# Then:
# cd /lustre/scratch126/tol/teams/jaron/users/fede/Orthofinder/Proteomes_2nd_round_OneIsoform/OrthoFinder/Results_Jul25_1/
# bsub5000 -e kk -o kk Rscript analyse_alignments_v2.R

# WHEN FINISHED:
# rm *STATS
# for sp in dmel phyg bcop-grc bcop bimp-grc bimp ling-grc ling aaphi orobi contarinia
#  do grep -P "$sp\t" *stats.tsv >> $sp.STATS
# done

# args = commandArgs(trailingOnly=TRUE)
# alignment = args[1]

library(seqinr)
setwd("/lustre/scratch126/tol/teams/jaron/users/fede/Orthofinder/Proteomes_2nd_round_OneIsoform/OrthoFinder/Results_Jul25_1")
ortg = read.table("Orthogroups/Orthogroups.GeneCount.tsv",sep="\t",header=T,row.names=1)
ortg = ortg[,1:(ncol(ortg)-1)]
presence = ortg
presence[] = 0
presence[ortg>0] = 1

# Requirements:
# * Max of 4 ort per species
# * Mean < 2
# * At least 4 species > 0
max_per_ortg  = apply(ortg,1,max)
mean_per_ortg = apply(ortg,1,mean)
sps_per_ortg  = apply(ortg,1,function(x) sum(x>0))
stats = as.data.frame(cbind(max_per_ortg,mean_per_ortg,sps_per_ortg))
good_ortg = stats[which(stats$max_per_ortg <= 4 & stats$sps_per_ortg >= 4 & stats$mean_per_ortg <= 2),]
good_ortg_names = rownames(good_ortg)
apply(ortg[good_ortg_names,],2,sum)
apply(presence[good_ortg_names,],2,sum)


#mkdir("ALISTATS")

cat(length(good_ortg_names), " good orthologous groups (max_per_sp <= 4 & mean_per_group <= 2 & sps_per_group >= 4)\n")
count = 0
for(og_name in good_ortg_names) {
	count = count + 1
	cat("Processing:",og_name,count,"/",length(good_ortg_names),"...\n")
	kk = read.alignment(paste("MultipleSequenceAlignments/",og_name,".fa",sep=""),format="fasta")
	kk = as.matrix.alignment(kk)
	kk[grep("[a-z]",kk)] = "1"
	kk[grep("-",kk)] = "0"
	kk[grep("\\*",kk)] = "0"
	#kk = as.numeric(kk)
	storage.mode(kk) <- "numeric"
	
	
	rownames(kk) = gsub("_grc","-grc",rownames(kk))
	species = sapply(rownames(kk),function(x) unlist(strsplit(x,"_"))[1])
	# table(species)
	
	genes = rownames(kk)
	genes = gsub("dmel_","",genes)
	genes = gsub("bimp_core_","",genes)
	genes = gsub("bimp_grc_","",genes)
	genes = gsub("bcop_core_","",genes)
	genes = gsub("bcop_grc_","",genes)
	genes = gsub("ling_core_","",genes)
	genes = gsub("ling_grc_","",genes)
	genes = gsub("contarinia_","",genes)
	genes = gsub("aaphi_","",genes)
	genes = gsub("orobi_","",genes)
	genes = gsub("phyg_","",genes)
	
	# Later I'll see in which cases the two genes are consecutive in the genome
	
	# Columns with at least 30% aas: (and in at least 3 seqs): 
	perc30 = vector(length=ncol(kk))
	perc30[] = 0
	perc30[which(apply(kk,2,sum)/nrow(kk) > 0.3 & apply(kk,2,sum) >= 3)] = 1
	
	per_gene = apply(kk,1,function(x) sum(x == perc30 & x == 1)) / sum(perc30)
	per_gene_unique = apply(kk,1,function(x) sum(x != perc30 & x == 1)) / sum(perc30)
	
	# It should be easy to calculate complementarity between consecutive genes
	
	results = cbind(species,genes,per_gene, per_gene_unique)
	colnames(results) = c("species","gene","perc_covered","perc_unique")
	write.table(results,file=paste("ALI_STATS/",og_name,".stats.tsv",sep=""),sep="\t",row.names=F,col.names=T,quote=F)
	#cat(warnings())
}


# In the end, do a violin plot!
# One per species!!!!!



