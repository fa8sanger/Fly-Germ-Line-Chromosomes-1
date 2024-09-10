# TRASH_run.sh assembly.fa --o output.path
cd /lustre/scratch125/casm/team268im/fa8/117/fede/Trash
mkdir Bimp_output
mkdir Lyco_output
mkdir Bcop_output
mkdir Aaphi_output # Aphid killer
mkdir Orobi_output # Obolodiplosis robiniae, Locust gall midge


# Trying now with absolute path... it failed before :-( 
cp -p /lustre/scratch126/tol/teams/jaron/data/assemblies_Sanger/insects/Bradysia_impatiens/assembly/curated/idBraImpa2.1/idBraImpa2.1.primary.curated.fa .
cp -p /lustre/scratch126/tol/teams/jaron/data/assemblies_Sanger/insects/Bradysia_coprophila/assembly/draft/idBraCopr2.20240426/idBraCopr2.20240426.primary.fa.gz .
cp -p /lustre/scratch126/tol/teams/jaron/data/assemblies_Sanger/insects/Lycoriella_ingenua/assembly/curated/idLycInge5.1/idLycInge5.1.primary.fa.gz .
cp -p /lustre/scratch123/tol/teams/grit/tom_mathers/curations/idBraCopr2_1/idBraCopr2_1_HAP1.primary.curated.fa .
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/030/463/065/GCA_030463065.1_ASM3046306v1/GCA_030463065.1_ASM3046306v1_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/028/476/595/GCA_028476595.1_CAF_Oro_1.0/GCA_028476595.1_CAF_Oro_1.0_genomic.fna.gz
gunzip *.gz


# Run only for the large chromosomes:
perl filter_contigs.pl idBraImpa2.1.primary.curated.fa           1e6 > bimp.fa
perl filter_contigs.pl idBraCopr2_1_HAP1.primary.curated.fa      1e6 > bcop.fa
perl filter_contigs.pl idLycInge5.1.primary.fa                   1e6 > lyco.fa
perl filter_contigs.pl GCA_030463065.1_ASM3046306v1_genomic.fna  1e6 > aaphi.fa
perl filter_contigs.pl GCA_028476595.1_CAF_Oro_1.0_genomic.fna   1e6 > orobi.fa


bsub5000 -q long -e bimp_trash.err  -o bimp_trash.err   ./TRASH_run.sh /lustre/scratch125/casm/team268im/fa8/117/fede/Trash/bimp.fa     --o /lustre/scratch125/casm/team268im/fa8/117/fede/Trash/Bimp_output 
bsub5000 -q long -e bcop_trash.err  -o bcop_trash.err   ./TRASH_run.sh /lustre/scratch125/casm/team268im/fa8/117/fede/Trash/bcop.fa     --o /lustre/scratch125/casm/team268im/fa8/117/fede/Trash/Bcop_output 
bsub5000 -q long -e lyco_trash.err  -o lyco_trash.err   ./TRASH_run.sh /lustre/scratch125/casm/team268im/fa8/117/fede/Trash/lyco.fa     --o /lustre/scratch125/casm/team268im/fa8/117/fede/Trash/Lyco_output 
bsub5000 -q long -e aaphi_trash.err -o aaphi_trash.err  ./TRASH_run.sh /lustre/scratch125/casm/team268im/fa8/117/fede/Trash/aaphi.fa    --o /lustre/scratch125/casm/team268im/fa8/117/fede/Trash/Aaphi_output 
bsub5000 -q long -e orobi_trash.err -o orobi_trash.err  ./TRASH_run.sh /lustre/scratch125/casm/team268im/fa8/117/fede/Trash/orobi.fa    --o /lustre/scratch125/casm/team268im/fa8/117/fede/Trash/Orobi_output 


