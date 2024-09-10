# Usage:
# bsub5000 -e LOG -o LOG Rscript analyse_alignments_with_identity_v2.R 
# 
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
# Load CDS sequences in an appropriate format
CDS_dir="/lustre/scratch126/tol/teams/jaron/users/fede/CDS_seqs/"
CDS4alis_dir="/lustre/scratch126/tol/teams/jaron/users/fede/CDS_seqs/CDS4alis/"
AA4alis_dir="/lustre/scratch126/tol/teams/jaron/users/fede/CDS_seqs/AA4alis/"
sps = c("aaphi","bcop_core","bcop_grc","bimp_core","bimp_grc","contarinia","dmel","ling_core","ling_grc","orobi","phyg")
cds_list = list()
for(sp in sps) {
	fasta=paste(CDS_dir,sp,".fa",sep="")
	cat("Loading",fasta,"...\n")
	jj = read.fasta(fasta,seqtype="DNA",as.string=TRUE)
	jj = unlist(sapply(jj,function(x) as.character(x)))
	cds_list[[sp]] = jj
}


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
	#kk = read.alignment(ali,format="fasta")
	kk = as.matrix.alignment(kk)
	aa = kk
	most_freq = apply(aa,2,function(x) names(sort(table(x),decreasing=T)[1]))
	most_freq[which(most_freq=="-")] = "@" # to avoid counting them
	kk[grep("[a-z]",kk)] = "1"
	kk[grep("-",kk)] = "0"
	kk[grep("\\*",kk)] = "0"
	#kk = as.numeric(kk)
	storage.mode(kk) <- "numeric"
	
	rownames(kk) = gsub("_grc","-grc",rownames(kk))
	rownames(kk) = gsub("_core","-core",rownames(kk))
	species = sapply(rownames(kk),function(x) unlist(strsplit(x,"_"))[1])
	# table(species)
	
	genes = rownames(kk)
	genes = gsub("dmel_","",genes)
	genes = gsub("bimp-core_","",genes)
	genes = gsub("bimp-grc_","",genes)
	genes = gsub("bcop-core_","",genes)
	genes = gsub("bcop-grc_","",genes)
	genes = gsub("ling-core_","",genes)
	genes = gsub("ling-grc_","",genes)
	genes = gsub("contarinia_","",genes)
	genes = gsub("aaphi_","",genes)
	genes = gsub("orobi_","",genes)
	genes = gsub("phyg_","",genes)
	
	# Later I'll see in which cases the two genes are consecutive in the genome
	
	# Columns with at least 30% aas: (and in at least 3 seqs): 
	perc30 = vector(length=ncol(kk))
	perc30[] = 0
	perc30[which(apply(kk,2,sum)/nrow(kk) > 0.3 & apply(kk,2,sum) >= 3)] = 1
	
	# This is for the entire ali, but for each gene I have to focus on what it got aligned
	per_gene_eq_consensus = apply(aa,1,function(x) sum(x==most_freq & perc30 == 1 & x!="-"))
	
	per_gene = apply(kk,1,function(x) sum(x == perc30 & x == 1)) 
	per_gene_unique = apply(kk,1,function(x) sum(x != perc30 & x == 1)) / sum(perc30)
	
	ident_per_gene = per_gene_eq_consensus/per_gene

	per_gene_frac = apply(kk,1,function(x) sum(x == perc30 & x == 1)) / sum(perc30)

	# It should be easy to calculate complementarity between consecutive genes
	
	results = as.data.frame(cbind(species,genes,per_gene_frac, per_gene_unique,ident_per_gene))
	colnames(results) = c("species","gene","perc_covered","perc_unique","ident_per_gene")
	write.table(results,file=paste("ALI_STATS/",og_name,".stats.tsv",sep=""),sep="\t",row.names=F,col.names=T,quote=F)
	#cat(warnings())
	
	######################################################################################
	# Now decide which ones to use for subsequent analysis! CDS alis and AA alis         #
	######################################################################################
	good_ones = which(results$perc_covered > 0.80 & results$perc_unique < 0.10 & results$ident_per_gene > 0.40)
	cat("   removing",nrow(results)-length(good_ones),"/",nrow(results),"sps=",length(species),"genes...\n");
	cat("      ",length(which(results$perc_covered <= 0.8)),"<0.8\n")
	cat("      ",length(which(results$perc_unique  >= 0.1)),">0.1\n")
	cat("      ",length(which(results$ident_per_gene <= 0.4)),"<0.4\n")
	species = gsub("-","_",species)
	species = species[good_ones]
	genes   = genes[good_ones]
	if(length(species) == 0) {
		next
	}
	aa_ali = aa[good_ones,]
	cds_seqs = vector()
	for(i in c(1:length(species))) {
		sp = species[i]
		gene = genes[i]
		if(sp == "phyg") {
			gene = unlist(strsplit(gene,"_"))[3]
		}
		#cat(sp,gene,"\n")
		zz = grep(gene,names(cds_list[[sp]]),value=T)
		if(length(zz) != 1) {
			cat("  ERROR?", sp, gene, length(zz),"\n")
			cat("  Choosing", zz[1],"\n")
			zz = zz[1]
		}
		# get the CDS sequences and save them!
		cds_seqs[paste(sp,gene,sep="_")] = cds_list[[sp]][zz]
		
		# and save the aa alignments as well		
	}
	if(length(cds_seqs) >= 3) {
		write.fasta(sapply(cds_seqs,function(x) s2c(x)),names(cds_seqs),as.string=FALSE,file.out=paste(CDS4alis_dir,og_name,".fa",sep=""))
		aa_seqs = sapply(1:nrow(aa_ali), function(x) paste(aa_ali[x,],collapse=""))
		lineas = vector()
		for(i in c(1:length(aa_seqs))) {
			lineas = c(lineas,paste(">",names(cds_seqs)[i],sep=""),aa_seqs[i])
		}
		fileConn<-file(paste(AA4alis_dir,og_name,".fa",sep=""))
		writeLines(lineas,fileConn)
		close(fileConn)

		#jorl = as.alignment(nb = nrow(aa_ali), nam = names(cds_seqs), seq = aa_seqs)
		#write.fasta(jorl,names=names(cds_seqs),file.out=paste(AA4alis_dir,og_name,".fa",sep=""))
		#write.fasta(sequences=sapply(aa_seqs, function(x) s2c(x)),names=names(cds_seqs),as.string=FALSE,file.out=paste(AA4alis_dir,og_name,".fa",sep=""))
	} 
	######################################################################################
	# Keep track of kept and modified orthogroups: original and modified OGs
	######################################################################################

	######################################################################################
}



# In the end, do a violin plot!
# One per species!!!!!



