cd /lustre/scratch126/tol/teams/jaron/users/fede
mkdir Blast_synteny
install_dir=/lustre/scratch125/casm/team268im/fa8/119/bin/ncbi-blast-2.15.0+/bin
cd /lustre/scratch126/tol/teams/jaron/users/fede/Blast_synteny
aaphi_genome=/lustre/scratch126/tol/teams/jaron/users/fede/Trash/GCA_030463065.1_ASM3046306v1_genomic.fna
lycorella_genome=/lustre/scratch126/tol/teams/jaron/data/assemblies_Sanger/insects/Lycoriella_ingenua/assembly/curated/idLycInge5.1/idLycInge5.1.primary.curated.fa
bimp_genome=/lustre/scratch126/tol/teams/jaron/data/assemblies_Sanger/insects/Bradysia_impatiens/assembly/curated/idBraImpa2.1/idBraImpa2.1.primary.curated.fa

# modify chr/scaffold/contig names:
cp -p $aaphi_genome ./aaphi_genome.fa
perl -pi -e 's/^>/>aaphi_/' aaphi_genome.fa
cp -p $bimp_genome ./bimp_genome.fa
perl -pi -e 's/^>/>bimp_/' bimp_genome.fa
# Remove GRCs from bimp:
perl remove_grcs.pl bimp_genome.fa > bimp_genome.Core.fa
cat aaphi_genome.fa bimp_genome.Core.fa > aaphi+bimp.fa

$install_dir/makeblastdb -in aaphi+bimp.fa -parse_seqids -dbtype nucl


cat /lustre/scratch126/tol/teams/jaron/users/sasha/diptera/Bradysia_impatiens/buscos/final/grcs_only/results_insecta_odb10/run_insecta_odb10/busco_sequences/single_copy_busco_sequences/*faa > bimp_grc_single_busco_proteins.fa


$install_dir/tblastn -db aaphi+bimp.fa -query bimp_grc_single_busco_proteins.fa -out GRC_vs_aaphi+bimpCore.out -outfmt 6 -evalue 0.0001 -num_descriptions 1 -num_alignments 1

$install_dir/tblastn -db aaphi+bimp.fa -query bimp_grc_single_busco_proteins.fa -out GRC_vs_aaphi+bimpCore.mts1.out -outfmt 6 -evalue 0.0001 -max_target_seqs 1

max_target_seqs


cd /lustre/scratch126/tol/teams/jaron/users/fede/Blast
install_dir=/lustre/scratch125/casm/team268im/fa8/119/bin/ncbi-blast-2.15.0+/bin
for a in *.fa *.fasta; 
 do echo $a
 $install_dir/makeblastdb -in $a -parse_seqids -dbtype prot
 done

install_dir=/lustre/scratch125/casm/team268im/fa8/119/bin/ncbi-blast-2.15.0+/bin
bsub5000 -e kk1 -o kk1 $install_dir/blastp -db contarinia.fa -query ling_braker2.grc.aa -out ling_vs_contarinia.grc.bls -evalue 0.0001 -num_descriptions 5 -num_alignments 5
bsub5000 -e kk2 -o kk3 $install_dir/blastp -db uniprotkb_drosophila_melanogaster_AND_m_2024_07_16.fasta -query ling_braker2.grc.aa -out ling_vs_dmel_uniprot.grc.bls -evalue 0.0001 -num_descriptions 5 -num_alignments 5
bsub5000 -e kk3 -o kk3 $install_dir/blastp -db dmel-all-translation-r6.58.fasta -query ling_braker2.grc.aa -out ling_vs_dmel_flybase.grc.bls -evalue 0.0001 -num_descriptions 5 -num_alignments 5

install_dir=/lustre/scratch125/casm/team268im/fa8/119/bin/ncbi-blast-2.15.0+/bin
bsub5000 -e kk1 -o kk1 $install_dir/blastp -db contarinia.fa -query ling_braker3.core.aa -out ling_vs_contarinia.core.bls -evalue 0.0001 -num_descriptions 5 -num_alignments 5
bsub5000 -e kk2 -o kk3 $install_dir/blastp -db uniprotkb_drosophila_melanogaster_AND_m_2024_07_16.fasta -query ling_braker3.core.aa -out ling_vs_dmel_uniprot.core.bls -evalue 0.0001 -num_descriptions 5 -num_alignments 5
bsub5000 -e kk3 -o kk3 $install_dir/blastp -db dmel-all-translation-r6.58.fasta -query ling_braker3.core.aa -out ling_vs_dmel_flybase.core.bls -evalue 0.0001 -num_descriptions 5 -num_alignments 5


Braker gene predictions for Lycoriella
* 21909 core genome, protein based
* 6845 grc genome, protein based
* 17536 core genome, protein+rnaseq based
I find it strange that braker3 predicts less than braker2 for the core genome

Out of 7641 GRC proteins (some of the 6845 genes have more than one isoform predicted)
* 6388 have homologs in Contarinia (Blast based)
* 5687 in Dmel (flybase, n=30802)
* 5934 in Dmel (uniprot, n=42676)

A quick look to the Blast alignments show many hits for the GRC with good conservation

Although drosophila and Lycoriella are distantly related, I was thinking on establishing
orthology relationships between them to take advantage of the available data on drosophila,
to for example look for expression profiles (tissue, development) of the genes in 
Lycoriellas GRCs.

Orthology relationships with Contarinia would be useful to estimate simple pairwise dN/dS
(or tree-based later)



# Probar Orthofinder:
https://github.com/davidemms/OrthoFinder

# dN/dS pairwise alignments:
# https://github.com/veg/hyphy/issues/1273

Hyphy in:
/lustre/scratch126/tol/teams/jaron/users/fede/software/bin

https://drostlab.github.io/orthologr/articles/dNdS_estimation.html

THIS SOUNDS GREAT! AND IN R
https://rdrr.io/github/HajkD/orthologr/man/compute_dnds.html


Dear @bioinfowheat,

If you have an alignment with two sequences, the simplest way to get pairwise estimates would be using the FitMG94.bf analysis from https://github.com/veg/hyphy-analyses/tree/master/FitMG94

hyphy /path/to/FitMG94.bf --alignment /path/to/2.fas --tree neighbor-joining --output /path/to/2.json
/path/to/2.json will have a lot of details on the fit result and is a JSON file. In particular you will see something like

"Standard MG94":{
     "AIC-c":3069.060099885265,
     "Confidence Intervals":{
       "non-synonymous/synonymous rate ratio":{
         "LB":0.05932492302024387,
         "UB":0.09858175408945982
        }
      },
     "Log Likelihood":-1517.976078455871,
     "Rate Distributions":{
       "Substitution rate from nucleotide A to nucleotide C":1.654785163501859,
       "Substitution rate from nucleotide A to nucleotide G":1,
       "Substitution rate from nucleotide A to nucleotide T":0.4123593781163388,
       "Substitution rate from nucleotide C to nucleotide G":0.920591030862722,
       "Substitution rate from nucleotide C to nucleotide T":1.083971553866486,
       "Substitution rate from nucleotide G to nucleotide T":0.4226548876413209,
       "non-synonymous/synonymous rate ratio":0.07708666449963816
      },
     "display order":1,
     "estimated parameters":16
    }
  }
The point estimate of dN/dS is "non-synonymous/synonymous rate ratio" : 0.07708666449963816 and the confidence interval is reported in

"Confidence Intervals":{
       "non-synonymous/synonymous rate ratio":{
         "LB":0.05932492302024387,
         "UB":0.09858175408945982
        }
      }
You can read JSON with Python/R etc or you can use jq (a very powerful command line tools to extract values). For example if you all you care about is the dN/dS estimate for a bunch of files you can do something like (with pathnames adjusted of course)

hyphy FitMG94.bf --alignment /Users/sergei/Desktop/2.txt --tree neighbor-joining --output pairwise.json >/dev/null; jq '.["fits"]["Standard MG94"]["Rate Distributions"]["non-synonymous/synonymous rate ratio"]' pairwise.json
Best,
Sergei