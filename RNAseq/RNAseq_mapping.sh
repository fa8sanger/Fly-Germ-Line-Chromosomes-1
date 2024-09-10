
export MODULEPATH=/software/treeoflife/shpc/current/views/tol:/software/treeoflife/custom-installs/modules:/software/modules
module add fastp/0.23.4--hadf994f_3
# paso de fastqc

# For Ling I don't need to trim the reads

# Fastp for Bcop:
cd /lustre/scratch125/casm/team268im/fa8/117/fede/RNAseq/Bcop
export MODULEPATH=/software/treeoflife/shpc/current/views/tol:/software/treeoflife/custom-installs/modules:/software/modules
module add fastp/0.23.4--hadf994f_3
mkdir fastp
for dir in 13F_1 13F_2 13F_3 13M_1 13M_2 13M_3 13M_4 15F_1 15F_2 15F_3 15M_1 15M_2 15M_3 21F_1 21F_2 21F_3 21M_1 21M_2 21M_3 21M_4 4F_1 4F_2 4F_3 4M_1 4M_2 4M_3 Fbody1 Fbody2 Fbody3 Fgerm1 Fgerm2 Fgerm3 Mbody1 Mbody2 Mbody3 Mgerm1 Mgerm2 Mgerm3 
 do echo $dir
 cd $dir
 bsub5000 -e /lustre/scratch125/casm/team268im/fa8/117/fede/RNAseq/Bcop/$dir.err -o /lustre/scratch125/casm/team268im/fa8/117/fede/RNAseq/Bcop/$dir.err "fastp -i *_1.fastq.gz -I *_2.fastq.gz -o /lustre/scratch125/casm/team268im/fa8/117/fede/RNAseq/Bcop/fastp/$dir.1.fastp.fastq.gz -O /lustre/scratch125/casm/team268im/fa8/117/fede/RNAseq/Bcop/fastp/$dir.2.fastp.fastq.gz --length_required 33 --cut_front --cut_tail --cut_mean_quality 20 --thread 1"
 cd ..
 done

# Index reference genomes with hisat2:
cd /lustre/scratch125/casm/team268im/fa8/117/fede/Genomes
mkdir hisat2
bsub10000 -e kk  -o kk  hisat2-build idLycInge5.1.primary.fa              hisat2/Ling
bsub10000 -e kk2 -o kk2 hisat2-build idBraCopr2_1_HAP1.primary.curated.fa hisat2/Bcop

# And run hisat2:
# Ling:
cd /lustre/scratch125/casm/team268im/fa8/117/fede/RNAseq/Ling
mkdir hisat2
for input in L10_combined L11_combined L12_EKRN230017329-1A_HYCMYDSX5_L3 L13_EKRN230017330-1A_HYCMYDSX5_L3 L14_EKRN230017331-1A_HYCMYDSX5_L3 L15_combined L16_EKRN230017333-1A_HYCMYDSX5_L3 L17_EKRN230017334-1A_HYCMYDSX5_L3 L18_EKRN230017335-1A_HYCMYDSX5_L3 L1_EKRN230017318-1A_HY5YHDSX5_L2 L2_EKRN230017319-1A_HY5YHDSX5_L2 L3_EKRN230017320-1A_HY5YHDSX5_L2 L4_EKRN230017321-1A_HYCJ7DSX5_L4 L5_EKRN230017322-1A_HY5YHDSX5_L2 L6_EKRN230017323-1A_HY5YHDSX5_L2 L7_EKRN230017324-1A_HYCMYDSX5_L3 L8_EKRN230017325-1A_HYCMYDSX5_L3 L9_EKRN230017326-1A_HYCMYDSX5_L3 
 do echo $input
 bsub20000 -e $input.err -o $input.err "hisat2 --dta -x /lustre/scratch125/casm/team268im/fa8/117/fede/Genomes/hisat2/Ling -1 \"$input\"_1.trimmed.fq.gz -2 \"$input\"_2.trimmed.fq.gz | samtools view -bS | samtools sort -n -o hisat2/$input.query_sorted.bam"
 done
# Bcop:
cd /lustre/scratch125/casm/team268im/fa8/117/fede/RNAseq/Bcop/fastp
mkdir hisat2
for input in 13F_1 13F_2 13F_3 13M_1 13M_2 13M_3 13M_4 15F_1 15F_2 15F_3 15M_1 15M_2 15M_3 21F_1 21F_2 21F_3 21M_1 21M_2 21M_3 21M_4 4F_1 4F_2 4F_3 4M_1 4M_2 4M_3 Fbody1 Fbody2 Fbody3 Fgerm1 Fgerm2 Fgerm3 Mbody1 Mbody2 Mbody3 Mgerm1 Mgerm2 Mgerm3 
 do echo $input
 bsub20000 -e $input.err -o $input.err "hisat2 --dta -x /lustre/scratch125/casm/team268im/fa8/117/fede/Genomes/hisat2/Bcop -1 $input.1.fastp.fastq.gz -2 $input.2.fastp.fastq.gz | samtools view -bS | samtools sort -n -o hisat2/$input.query_sorted.bam"
 done


# Coordinate sort (and index) all the resulting bams:
cd /lustre/scratch125/casm/team268im/fa8/117/fede/RNAseq/Ling/hisat2
for input in L10_combined L11_combined L12_EKRN230017329-1A_HYCMYDSX5_L3 L13_EKRN230017330-1A_HYCMYDSX5_L3 L14_EKRN230017331-1A_HYCMYDSX5_L3 L15_combined L16_EKRN230017333-1A_HYCMYDSX5_L3 L17_EKRN230017334-1A_HYCMYDSX5_L3 L18_EKRN230017335-1A_HYCMYDSX5_L3 L1_EKRN230017318-1A_HY5YHDSX5_L2 L2_EKRN230017319-1A_HY5YHDSX5_L2 L3_EKRN230017320-1A_HY5YHDSX5_L2 L4_EKRN230017321-1A_HYCJ7DSX5_L4 L5_EKRN230017322-1A_HY5YHDSX5_L2 L6_EKRN230017323-1A_HY5YHDSX5_L2 L7_EKRN230017324-1A_HYCMYDSX5_L3 L8_EKRN230017325-1A_HYCMYDSX5_L3 L9_EKRN230017326-1A_HYCMYDSX5_L3 
 do echo $input
  bsub5000 -e $input.err -o $input.err "samtools sort -o $input.coord_sorted.bam $input.query_sorted.bam; samtools index $input.coord_sorted.bam "
 done
cd /lustre/scratch125/casm/team268im/fa8/117/fede/RNAseq/Bcop/fastp/hisat2
for input in 13F_1 13F_2 13F_3 13M_1 13M_2 13M_3 13M_4 15F_1 15F_2 15F_3 15M_1 15M_2 15M_3 21F_1 21F_2 21F_3 21M_1 21M_2 21M_3 21M_4 4F_1 4F_2 4F_3 4M_1 4M_2 4M_3 Fbody1 Fbody2 Fbody3 Fgerm1 Fgerm2 Fgerm3 Mbody1 Mbody2 Mbody3 Mgerm1 Mgerm2 Mgerm3 
 do echo $input
  bsub5000 -e $input.err -o $input.err "samtools sort -o $input.coord_sorted.bam $input.query_sorted.bam; samtools index $input.coord_sorted.bam "
 done


# Merge (and index) all the bams into one:
cd /lustre/scratch125/casm/team268im/fa8/117/fede/RNAseq/Ling/hisat2
bsub5000 -e kk.err -o kk.err "samtools merge All_rnaseq.bam *coord_sorted.bam; samtools index All_rnaseq.bam"
cd /lustre/scratch125/casm/team268im/fa8/117/fede/RNAseq/Bcop/fastp/hisat2
bsub5000 -e kk.err -o kk.err "samtools merge All_rnaseq.bam *coord_sorted.bam; samtools index All_rnaseq.bam"


# cd /lustre/scratch125/casm/team268im/fa8/117/fede/RNAseq/Ling/hisat2
# for a in *coord_sorted.bam; do echo $a `samtools view -c $a SUPER_GRC1` `samtools view -c $a SUPER_1`; done
# cd /lustre/scratch125/casm/team268im/fa8/117/fede/RNAseq/Bcop/fastp/hisat2
# for a in *coord_sorted.bam; do echo $a `samtools view -c $a SUPER_GRC1` `samtools view -c $a SUPER_1`; done


# Ling:
cd /lustre/scratch125/casm/team268im/fa8/117/fede/RNAseq/Ling/hisat2
echo "bam SUPER_1 SUPER_X SUPER_2 SUPER_3 SUPER_GRC1 SUPER_GRC1_unloc_1 SUPER_GRC2" | tr ' ' '\t' > Ling.readcounts.tsv
for a in *coord_sorted.bam;   
 do echo $a `samtools view -c $a SUPER_1` `samtools view -c $a SUPER_X` `samtools view -c $a SUPER_2` `samtools view -c $a SUPER_3` `samtools view -c $a SUPER_GRC1` `samtools view -c $a SUPER_GRC1_unloc_1` `samtools view -c $a SUPER_GRC2`| tr ' ' '\t'  >> Ling.readcounts.tsv ;  
 done 
# Bcop:
cd /lustre/scratch125/casm/team268im/fa8/117/fede/RNAseq/Bcop/fastp/hisat2
echo "bam SUPER_1 SUPER_X SUPER_2 SUPER_3 SUPER_GRC1 SUPER_GRC2" | tr ' ' '\t' > Bcop.readcounts.tsv
for a in *coord_sorted.bam; 
 do echo $a `samtools view -c $a SUPER_1` `samtools view -c $a SUPER_X` `samtools view -c $a SUPER_2` `samtools view -c $a SUPER_3` `samtools view -c $a SUPER_GRC1`  `samtools view -c $a SUPER_GRC2` | tr ' ' '\t' >> Bcop.readcounts.tsv; 
done 

# Ignoring alignments with NH != 1 (only want those mapping uniquely for this estimation), only proper pairs (-f 2) , not secondary alignments (-F 256):
# Ling:
cd /lustre/scratch125/casm/team268im/fa8/117/fede/RNAseq/Ling/hisat2
echo "bam SUPER_1 SUPER_X SUPER_2 SUPER_3 SUPER_GRC1 SUPER_GRC1_unloc_1 SUPER_GRC2" | tr ' ' '\t' > Ling.readcounts-Unique-f2-F256.tsv
for a in *coord_sorted.bam;   
 do echo $a `samtools view -f 2 -F 256 $a SUPER_1 | grep -P "NH:i:1$" | wc -l` `samtools view -f 2 -F 256 $a SUPER_X | grep -P "NH:i:1$" | wc -l` `samtools view -f 2 -F 256 $a SUPER_2 | grep -P "NH:i:1$" | wc -l` `samtools view -f 2 -F 256 $a SUPER_3 | grep -P "NH:i:1$" | wc -l` `samtools view -f 2 -F 256 $a SUPER_GRC1 | grep -P "NH:i:1$" | wc -l` `samtools view -f 2 -F 256 $a SUPER_GRC1_unloc_1 | grep -P "NH:i:1$" | wc -l` `samtools view -f 2 -F 256 $a SUPER_GRC2 | grep -P "NH:i:1$" | wc -l`| tr ' ' '\t'  >> Ling.readcounts-Unique-f2-F256.tsv ;  
 done 
# Bcop:
cd /lustre/scratch125/casm/team268im/fa8/117/fede/RNAseq/Bcop/fastp/hisat2
echo "bam SUPER_1 SUPER_X SUPER_2 SUPER_3 SUPER_GRC1 SUPER_GRC2" | tr ' ' '\t' > Bcop.readcounts-Unique-f2-F256.tsv
for a in *coord_sorted.bam; 
 do echo $a `samtools view -f 2 -F 256 $a SUPER_1 | grep -P "NH:i:1$" | wc -l` `samtools view -f 2 -F 256 $a SUPER_X | grep -P "NH:i:1$" | wc -l` `samtools view -f 2 -F 256 $a SUPER_2 | grep -P "NH:i:1$" | wc -l` `samtools view -f 2 -F 256 $a SUPER_3 | grep -P "NH:i:1$" | wc -l` `samtools view -f 2 -F 256 $a SUPER_GRC1 | grep -P "NH:i:1$" | wc -l`  `samtools view -f 2 -F 256 $a SUPER_GRC2 | grep -P "NH:i:1$" | wc -l` | tr ' ' '\t' >> Bcop.readcounts-Unique-f2-F256.tsv; 
done 

@SQ     SN:SUPER_1      LN:76508704
@SQ     SN:SUPER_X      LN:67901076
@SQ     SN:SUPER_2      LN:56451848
@SQ     SN:SUPER_3      LN:52603560
@SQ     SN:SUPER_GRC1   LN:36229871
@SQ     SN:SUPER_GRC1_unloc_1   LN:1196638
@SQ     SN:SUPER_GRC2   LN:34896639

@SQ     SN:SUPER_1      LN:76508704
@SQ     SN:SUPER_X      LN:67901076
@SQ     SN:SUPER_2      LN:56451848
@SQ     SN:SUPER_3      LN:52603560
@SQ     SN:SUPER_GRC1   LN:36229871
@SQ     SN:SUPER_GRC1_unloc_1   LN:1196638
@SQ     SN:SUPER_GRC2   LN:34896639


# Ling:
L10_combined.coord_sorted.bam 17085 23856725
L11_combined.coord_sorted.bam 25356 23360093
L12_EKRN230017329-1A_HYCMYDSX5_L3.coord_sorted.bam 20011 23704263
L13_EKRN230017330-1A_HYCMYDSX5_L3.coord_sorted.bam 56177 23565762
L14_EKRN230017331-1A_HYCMYDSX5_L3.coord_sorted.bam 39609 19435227
L15_combined.coord_sorted.bam 52107 22000782
L16_EKRN230017333-1A_HYCMYDSX5_L3.coord_sorted.bam 35156 23226206
L17_EKRN230017334-1A_HYCMYDSX5_L3.coord_sorted.bam 29779 21891476
L18_EKRN230017335-1A_HYCMYDSX5_L3.coord_sorted.bam 47208 20786861
L1_EKRN230017318-1A_HY5YHDSX5_L2.coord_sorted.bam 81991 16825003
L2_EKRN230017319-1A_HY5YHDSX5_L2.coord_sorted.bam 43051 10365216
L3_EKRN230017320-1A_HY5YHDSX5_L2.coord_sorted.bam 60912 15805729
L4_EKRN230017321-1A_HYCJ7DSX5_L4.coord_sorted.bam 69735 31700377
L5_EKRN230017322-1A_HY5YHDSX5_L2.coord_sorted.bam 59666 35310035
L6_EKRN230017323-1A_HY5YHDSX5_L2.coord_sorted.bam 101171 29907120
L7_EKRN230017324-1A_HYCMYDSX5_L3.coord_sorted.bam 109845 20880356
L8_EKRN230017325-1A_HYCMYDSX5_L3.coord_sorted.bam 105062 20023390
L9_EKRN230017326-1A_HYCMYDSX5_L3.coord_sorted.bam 77253 17115671

# Bcop:
13F_1.coord_sorted.bam 25609 56787
13F_2.coord_sorted.bam 35312 93607
13F_3.coord_sorted.bam 59205 187337
13M_1.coord_sorted.bam 52428 182481
13M_2.coord_sorted.bam 52178 180430
13M_3.coord_sorted.bam 50256 178511
13M_4.coord_sorted.bam 38088 129172
15F_1.coord_sorted.bam 82866 289311
15F_2.coord_sorted.bam 53023 170026
15F_3.coord_sorted.bam 54576 167524
15M_1.coord_sorted.bam 24816 82149
15M_2.coord_sorted.bam 20830 65800
15M_3.coord_sorted.bam 22972 75222
21F_1.coord_sorted.bam 23393 65652
21F_2.coord_sorted.bam 30114 80481
21F_3.coord_sorted.bam 33275 88999
21M_1.coord_sorted.bam 19010 60472
21M_2.coord_sorted.bam 23506 76343
21M_3.coord_sorted.bam 38394 137493
21M_4.coord_sorted.bam 40129 144997
4F_1.coord_sorted.bam 31101 61235
4F_2.coord_sorted.bam 25313 65108
4F_3.coord_sorted.bam 27732 70094
4M_1.coord_sorted.bam 24837 84184
4M_2.coord_sorted.bam 22884 76800
4M_3.coord_sorted.bam 16548 52725
Fbody1.coord_sorted.bam 970219 2430153
Fbody2.coord_sorted.bam 805784 1990553
Fbody3.coord_sorted.bam 897977 2085473
Fgerm1.coord_sorted.bam 934596 2631525
Fgerm2.coord_sorted.bam 745834 1996187
Fgerm3.coord_sorted.bam 709482 1646940
Mbody1.coord_sorted.bam 983857 2510645
Mbody2.coord_sorted.bam 951818 2430691
Mbody3.coord_sorted.bam 848903 2094680
Mgerm1.coord_sorted.bam 697174 2096900
Mgerm2.coord_sorted.bam 553347 1676364
Mgerm3.coord_sorted.bam 591072 1569524


process trimReads {
        publishDir params.outdir, mode:'copy'

        input:
        tuple val(sample_id), path(reads)

        output:
        tuple val(sample_id), path('fastp/*.fastp.fastq.gz') , optional:true, emit: reads

        script:
        """
        mkdir fastp
        fastp -i ${reads[0]} -I ${reads[1]} -o fastp/${sample_id}.1.fastp.fastq.gz -O fastp/${sample_id}.2.fastp.fastq.gz --length_required 33 --cut_front --cut_tail --cut_mean_quality 20 --thread ${task.cpus}
        """
}

process fastqc {
        tag "FASTQC on $sample_id"
        
        input:
        tuple val(sample_id), path(reads)

        output:
        path("fastqc/${sample_id}")

        script:
        """
        mkdir -p fastqc/${sample_id}
        fastqc -o fastqc/${sample_id} --nogroup -f fastq -q ${reads} 
        """
}

process multiqc {
        publishDir params.outdir, mode:'move'

        input:
        val(prefix)
        path('*')

        output:
        path("${prefix}_multiqc")

        script:
        """
        multiqc -o ${prefix}_multiqc . 
        """
}

process indexGenomeHisat2 {

        input:
        path genome_f

        output:
        path "hisat2", emit: index

        script:
        """
        mkdir hisat2
        hisat2-build ${genome_f} hisat2/genome
        """
}

process mapToGenomeHisat2 {
        publishDir params.outdir, mode:'copy'

        input:
        path(index)
        tuple val(sample_id), path(reads)

        output:
        path("hisat2/bams/${sample_id}.query_sorted.bam*")

        script:
        """
        mkdir -p hisat2/bams
        INDEX=`find -L ./ -name "*.1.ht2" | sed 's/\\.1.ht2\$//'`
        hisat2 --dta -x \$INDEX -1 ${reads[0]} -2 ${reads[1]} | samtools view -bS | samtools sort -n -o hisat2/bams/${sample_id}.query_sorted.bam
        """
}