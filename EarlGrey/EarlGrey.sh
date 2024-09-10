# Modified by fa8 based on Sam's (se13) script 

export MODULEPATH=/software/treeoflife/shpc/current/views/tol:/software/treeoflife/custom-installs/modules:/software/modules

module load braker3/3.0.8--hdfd78af_0
module load earlgrey/3.0-c1
module load bedtools/2.31.1--hf5e1c6e_1 

cd /lustre/scratch125/casm/team268im/fa8/117/fede/EarlGrey
ln -s /lustre/scratch125/casm/team268im/fa8/117/fede/Genomes/GCA_028476595.1_CAF_Oro_1.0_genomic.fna Orobi/genome.fa
ln -s /lustre/scratch125/casm/team268im/fa8/117/fede/Genomes/GCA_030463065.1_ASM3046306v1_genomic.fna Aaphi/genome.fa
ln -s /lustre/scratch125/casm/team268im/fa8/117/fede/Genomes/idBraCopr2_1_HAP1.primary.curated.fa Bcop/genome.fa
ln -s /lustre/scratch125/casm/team268im/fa8/117/fede/Genomes/idBraImpa2.1.primary.curated.fa Bimp/genome.fa
cp -p /lustre/scratch125/casm/team268im/fa8/117/fede/Genomes/idLycInge5.1.primary.fa.gz Ling/genome.fa.gz
gunzip Ling/genome.fa.gz


# Simplify headers:
cd /lustre/scratch125/casm/team268im/fa8/117/fede/EarlGrey
for species in Orobi Aaphi Bcop Bimp Ling
 do echo $species
 cd $species
 awk '{print $1}' genome.fa > genome.simple_header.fa
 cd ..
 done

# Remove softmasking:
cd /lustre/scratch125/casm/team268im/fa8/117/fede/EarlGrey
for species in Orobi Aaphi Bcop Bimp Ling
 do echo $species
 cd $species
 awk '/^>/ {print($0)}; /^[^>]/ {print(toupper($0))}' genome.simple_header.fa > genome.simple_header.unmasked.fa
 cd ..
 done


# Run EarlGrey
export MODULEPATH=/software/treeoflife/shpc/current/views/tol:/software/treeoflife/custom-installs/modules:/software/modules
module load braker3/3.0.8--hdfd78af_0
module load earlgrey/3.0-c1
module load bedtools/2.31.1--hf5e1c6e_1 
cd /lustre/scratch125/casm/team268im/fa8/117/fede/EarlGrey
for species in Orobi Aaphi Bcop Bimp Ling
 do echo $species
 cd $species
 mkdir results_dir
 cpus=8
 bsub60000 -q basement -e $species.earl.err -o $species.earl.err -n $cpus -R span[ptile=$cpus] earlGrey -g genome.simple_header.fa -s $species -o results_dir/ -t $cpus
 cd ..
 done


# Separate GRCs from core genome:
cd /lustre/scratch125/casm/team268im/fa8/117/fede/EarlGrey
for species in Bcop Bimp Ling
 do echo $species
 cd $species
 perl /lustre/scratch125/casm/team268im/fa8/117/fede/EarlGrey/separate_genomes.pl genome.simple_header.unmasked.fa
 cd ..
 done


# Running EarlGrey separately for core and GRCs
export MODULEPATH=/software/treeoflife/shpc/current/views/tol:/software/treeoflife/custom-installs/modules:/software/modules
module load braker3/3.0.8--hdfd78af_0
module load earlgrey/3.0-c1
module load bedtools/2.31.1--hf5e1c6e_1 
cd /lustre/scratch125/casm/team268im/fa8/117/fede/EarlGrey
for species in Bcop Bimp Ling
 do echo $species
 cd $species
 mkdir core_results_dir
 mkdir grc_results_dir
 rm $species.core.err
 rm $species.grc.err
 cpus=12
 bsub60000 -q basement -e $species.grc.err  -o $species.grc.err  -n $cpus -R span[ptile=$cpus] earlGrey -g grc.fa  -s $species -o grc_results_dir/  -t $cpus
 bsub60000 -q basement -e $species.core.err -o $species.core.err -n $cpus -R span[ptile=$cpus] earlGrey -g core.fa -s $species -o core_results_dir/ -t $cpus
 cd ..
 done





export MODULEPATH=/software/treeoflife/shpc/current/views/tol:/software/treeoflife/custom-installs/modules:/software/modules
module load braker3/3.0.8--hdfd78af_0
module load earlgrey/3.0-c1
module load bedtools/2.31.1--hf5e1c6e_1 
cd /lustre/scratch125/casm/team268im/fa8/117/fede/EarlGrey
 species=Orobi
 cd $species
 cpus=8
 bsub -M 100000 -R "select[mem>100000] rusage[mem=100000]" -q basement -e $species.earl.new.err -o $species.earl.new.err -n $cpus -R span[ptile=$cpus] earlGrey -g genome.simple_header.fa -s $species -o results_dir/ -t $cpus


# And mask fasta
# I remask the unmasked genome (not original softmasked genome) using earlGrey repeat bed.
bedtools maskfasta -soft -fi genome.simple_header.unmasked.fa -bed earlGrey_repeats.bed -fo genome.simple_header.earlGrey_masked.fa














 
 