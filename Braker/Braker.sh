
# Braker tutorial:
# https://www.youtube.com/watch?v=UXTkJ4mUkyg

# Braker
# I run braker as follows as braker 3 (using RNAseq data) or braker 2 (without RNAseq data). 
# With braker 3 I use the best available prot db. For species where I need to use braker2 
# I added the results of my species with a braker3 run to the prot db. I'm not sure what 
# is ideal, but it is worth thinking about what you expect to be annotating when choosing 
# your prot db. I haven't actually run braker3 with the module only braker2, but this 
# should work fine. There are other ways of providing rnaseq data to braker3 detailed on 
# the github. I have a nextflow pipeline for mapping rnaseq reads too but it is the first 
# one I made so let me know if you want it -- I could improve it a little.
# 
# Resources
# memory '40G'
# cpus 32
# queue 'long'
# 
# Braker recommended prot dbs : https://bioinf.uni-greifswald.de/bioinf/partitioned_odb11/
# 
# '--useexisting: Use the present config and parameter files if they exist for 'species'; 
# will overwrite original parameters if BRAKER performs an AUGUSTUS training.'

wget https://bioinf.uni-greifswald.de/bioinf/partitioned_odb11/Arthropoda.fa.gz 
ln -s  /lustre/scratch126/tol/teams/jaron/users/sasha/diptera/annotations/Contarinia_nasturtii/GCF_009176525.2 .
cat GCF_009176525.2/protein.faa Arthropoda.fa > Arthropoda+contarina.fa


module add augustus/3.5.0--pl5321h700735d_3
module add braker3/3.0.8--hdfd78af_0


# Mask genomes for EarlGrey repeats:
cd /lustre/scratch125/casm/team268im/fa8/117/fede/EarlGrey/Aaphi
bedtools maskfasta -soft -fi genome.simple_header.unmasked.fa -bed results_dir/Aaphi_EarlGrey/Aaphi_summaryFiles/Aaphi.filteredRepeats.bed -fo genome.simple_header.earlGrey_masked.fa
cd /lustre/scratch125/casm/team268im/fa8/117/fede/EarlGrey/Orobi
bedtools maskfasta -soft -fi genome.simple_header.unmasked.fa -bed results_dir/Orobi_EarlGrey/Orobi_summaryFiles/Orobi.filteredRepeats.bed -fo genome.simple_header.earlGrey_masked.fa


# Install Augustus:
cd /lustre/scratch125/casm/team268im/fa8/119/bin
git clone https://github.com/Gaius-Augustus/Augustus.git
cd Augustus
# switch off MySQL usage by setting MYSQL = false in common.mk --> I had an error before
make augustus
module add bamtools/2.5.2--hdcf5f25_3
make auxprogs # not all of them compiled successfully
# export PATH=~/augustus/bin:~/augustus/scripts:$PATH
# export AUGUSTUS_CONFIG_PATH=/my_path_to_AUGUSTUS/augustus/config/

# Install Prothint:
cd /lustre/scratch125/casm/team268im/fa8/119/bin
git clone https://github.com/gatech-genemark/ProtHint.git

# Install GeneMark:
cd /lustre/scratch125/casm/team268im/fa8/119/bin
git clone https://github.com/gatech-genemark/GeneMark-ETP.git
cd GeneMark-ETP
# export PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/ProtHint/bin:$PATH
./check_install.pl


# Make protein db headers simpler:
cd /lustre/scratch125/casm/team268im/fa8/117/fede/Braker
awk '{print $1}' Arthropoda+contarina.fa > Arthropoda+contarina.simple_header.fa


## module add singularity
## singularity build braker3 docker://teambraker/braker3:latest
## singularity exec braker3 braker.pl [OPTIONS]
## 
## 
## speciesName=Aaphi
## module purge
## module load singularity
## cd /lustre/scratch125/casm/team268im/fa8/117/fede/Braker/$speciesName
## masked_genome=/lustre/scratch125/casm/team268im/fa8/117/fede/EarlGrey/$speciesName/genome.simple_header.earlGrey_masked.fa
## cpus=32
## #prot_db=/lustre/scratch125/casm/team268im/fa8/117/fede/Braker/Arthropoda+contarina.fa
## prot_db=/lustre/scratch125/casm/team268im/fa8/117/fede/Braker/Arthropoda.fa
## singularity exec -B ${PWD}:${PWD} braker3 braker.pl --genome=$masked_genome workingdir=/lustre/scratch125/casm/team268im/fa8/117/fede/Braker/$speciesName --GENEMARK_PATH=${ETP}/gmes --softmasking --workingdir=. --threads $cpus --species=$speciesName --gff3 --prot_seq=$prot_db  --useexisting

# Trying braker3 again:
#  module add singularity
#  singularity build braker3 docker://teambraker/braker3:latest
#  
#  masked_genome=/lustre/scratch125/casm/team268im/fa8/117/fede/Braker/Ling/core.fa
#  prot_db=/lustre/scratch125/casm/team268im/fa8/117/fede/Braker/Arthropoda+contarina.fa
#  rnaseq=/lustre/scratch125/casm/team268im/fa8/117/fede/RNAseq/Ling/hisat2/All_rnaseq.bam
#  export AUGUSTUS_CONFIG_PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/Augustus/config
#  #export AUGUSTUS_BIN_PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/Augustus/bin
#  export AUGUSTUS_SCRIPTS_PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/Augustus/scripts
#  export GENEMARK_PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/GeneMark-ETP/bin
#  export PROTHINT_PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/ProtHint/bin
#
#  singularity exec -B ${PWD}:${PWD}  braker3 braker.pl --genome=$masked_genome workingdir=/lustre/scratch125/casm/team268im/fa8/117/fede/Braker/$speciesName --GENEMARK_PATH=${ETP}/gmes --bam=$rnaseq --softmasking --workingdir=. --threads $cpus --species=$speciesName --gff3 --prot_seq=$prot_db  --useexisting

# Braker2 (with prot db)
# Resources
# memory '40G'
# cpus 32
# queue 'long'
# 
# export PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/Augustus/bin:/lustre/scratch125/casm/team268im/fa8/119/bin/Augustus/scripts:$PATH
# export PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/GeneMark-ETP/bin:$PATH
# export PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/ProtHint/bin:$PATH
# module add braker2/2.1.6--hdfd78af_5
module add bamtools/2.5.2--hdcf5f25_3
module add augustus/3.5.0--pl5321h700735d_3
module add braker3/3.0.8--hdfd78af_0
speciesName=Aaphi
cd /lustre/scratch125/casm/team268im/fa8/117/fede/Braker/$speciesName
masked_genome=/lustre/scratch125/casm/team268im/fa8/117/fede/EarlGrey/$speciesName/genome.simple_header.earlGrey_masked.fa
cpus=16
prot_db=/lustre/scratch125/casm/team268im/fa8/117/fede/Braker/Arthropoda+contarina.simple_header.fa
export AUGUSTUS_CONFIG_PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/Augustus/config
#export AUGUSTUS_BIN_PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/Augustus/bin
export AUGUSTUS_SCRIPTS_PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/Augustus/scripts
export GENEMARK_PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/GeneMark-ETP/bin
export PROTHINT_PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/ProtHint/bin
#braker.pl --genome=$masked_genome --softmasking --workingdir=. --threads $cpus --species=speciesName --gff3 --prot_seq=$prot_db  --useexisting
bsub40000 -q long -e braker.err -o braker.err -n $cpus -R span[ptile=$cpus]  braker.pl --genome=$masked_genome --softmasking --workingdir=. --threads $cpus --species=$speciesName --gff3 --prot_seq=$prot_db  --useexisting


module add bamtools/2.5.2--hdcf5f25_3
module add augustus/3.5.0--pl5321h700735d_3
module add braker3/3.0.8--hdfd78af_0
speciesName=Orobi
cd /lustre/scratch125/casm/team268im/fa8/117/fede/Braker/$speciesName
masked_genome=/lustre/scratch125/casm/team268im/fa8/117/fede/EarlGrey/$speciesName/genome.simple_header.earlGrey_masked.fa
cpus=16
prot_db=/lustre/scratch125/casm/team268im/fa8/117/fede/Braker/Arthropoda+contarina.simple_header.fa
export AUGUSTUS_CONFIG_PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/Augustus/config
#export AUGUSTUS_BIN_PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/Augustus/bin
export AUGUSTUS_SCRIPTS_PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/Augustus/scripts
export GENEMARK_PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/GeneMark-ETP/bin
export PROTHINT_PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/ProtHint/bin
#braker.pl --genome=$masked_genome --softmasking --workingdir=. --threads $cpus --species=speciesName --gff3 --prot_seq=$prot_db  --useexisting
bsub40000 -q long -e braker.err -o braker.err -n $cpus -R span[ptile=$cpus]  braker.pl --genome=$masked_genome --softmasking --workingdir=. --threads $cpus --species=$speciesName --gff3 --prot_seq=$prot_db  --useexisting


Run Braker2 (prot only) for the GRC of Lycoriella
Run Braker3 (prot+RNAseq) for the core genome of Lycoriella
Look for orthologs in Drosophila... to transfer annotations, do studies, etc (ORTHOFINDER?)
I could also use annotations from Contarinia

##########################################################################################
# Lycoriella:
  cd /lustre/scratch125/casm/team268im/fa8/117/fede/Braker/Ling
  repeats_bed=/lustre/scratch126/tol/teams/jaron/users/sasha/diptera/Lycoriella_ingenua/TEs/results/idLyncInge5_EarlGrey/idLyncInge5_summaryFiles/idLyncInge5.filteredRepeats.bed 
  genome_unmasked=/lustre/scratch126/tol/teams/jaron/users/sasha/diptera/Lycoriella_ingenua/TEs/idLycInge5.1.primary.fa
  rnaseq=/lustre/scratch125/casm/team268im/fa8/117/fede/RNAseq/Ling/hisat2/All_rnaseq.bam
  bedtools maskfasta -soft -fi $genome_unmasked -bed $repeats_bed -fo /lustre/scratch125/casm/team268im/fa8/117/fede/Braker/Ling/genome.simple_header.earlGrey_masked.fa

# Separate GRC and core genome
  perl /lustre/scratch125/casm/team268im/fa8/117/fede/EarlGrey/separate_genomes.pl genome.simple_header.earlGrey_masked.fa

# Run Braker3 for core genome:
  #  cd /lustre/scratch126/tol/teams/jaron/users/fede/Braker/Ling
  #  rnaseq=/lustre/scratch126/tol/teams/jaron/users/fede/RNAseq/Ling/hisat2/All_rnaseq.bam
  #  module add bamtools/2.5.2--hdcf5f25_3
  #  module add augustus/3.5.0--pl5321h700735d_3
  #  module add braker3/3.0.8--hdfd78af_0
  #  masked_genome=core.fa
  #  cpus=16
  #  prot_db=/lustre/scratch126/tol/teams/jaron/users/fede/Braker/Arthropoda+contarina.simple_header.fa
  #  export AUGUSTUS_CONFIG_PATH=/lustre/scratch126/tol/teams/jaron/users/fede/software/Augustus/config
  #  #export AUGUSTUS_BIN_PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/Augustus/bin
  #  export AUGUSTUS_SCRIPTS_PATH=/lustre/scratch126/tol/teams/jaron/users/fede/software/Augustus/scripts
  #  export GENEMARK_PATH=/lustre/scratch126/tol/teams/jaron/users/fede/software/GeneMark-ETP/bin
  #  export PROTHINT_PATH=/lustre/scratch126/tol/teams/jaron/users/fede/software/ProtHint/bin
  #  mkdir Braker3_core
  #  bsub40000 -q long -e braker3core.err -o braker3core.err -n $cpus -R span[ptile=$cpus]  braker.pl --bam $rnaseq--genome=$masked_genome --softmasking --workingdir=Braker3_core --threads $cpus --species=Lycoriella.core.rnaseq --gff3 --prot_seq=$prot_db  --useexisting
  Now using bsub < braker_job.sh!


  # cd /lustre/scratch125/casm/team268im/fa8/117/fede/Braker/Ling
  # rnaseq=/lustre/scratch125/casm/team268im/fa8/117/fede/RNAseq/Ling/hisat2/All_rnaseq.bam
  # module add bamtools/2.5.2--hdcf5f25_3
  # module add augustus/3.5.0--pl5321h700735d_3
  # module add braker3/3.0.8--hdfd78af_0
  # masked_genome=core.fa
  # cpus=16
  # prot_db=/lustre/scratch125/casm/team268im/fa8/117/fede/Braker/Arthropoda+contarina.simple_header.fa
  # export AUGUSTUS_CONFIG_PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/Augustus/config
  # #export AUGUSTUS_BIN_PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/Augustus/bin
  # export AUGUSTUS_SCRIPTS_PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/Augustus/scripts
  # export GENEMARK_PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/GeneMark-ETP/bin
  # export PROTHINT_PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/ProtHint/bin
  # mkdir Braker3_core
  # bsub40000 -q long -e braker3core.err -o braker3core.err -n $cpus -R span[ptile=$cpus]  braker.pl --bam $rnaseq--genome=$masked_genome --softmasking --workingdir=Braker3_core --threads $cpus --species=Lycoriella.core.rnaseq --gff3 --prot_seq=$prot_db  --useexisting

# Run Braker2 for core genome:
  cd /lustre/scratch125/casm/team268im/fa8/117/fede/Braker/Ling
  module add bamtools/2.5.2--hdcf5f25_3
  module add augustus/3.5.0--pl5321h700735d_3
  module add braker3/3.0.8--hdfd78af_0
  masked_genome=core.fa
  cpus=16
  prot_db=/lustre/scratch125/casm/team268im/fa8/117/fede/Braker/Arthropoda+contarina.simple_header.fa
  export AUGUSTUS_CONFIG_PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/Augustus/config
  #export AUGUSTUS_BIN_PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/Augustus/bin
  export AUGUSTUS_SCRIPTS_PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/Augustus/scripts
  export GENEMARK_PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/GeneMark-ETP/bin
  export PROTHINT_PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/ProtHint/bin
  #braker.pl --genome=$masked_genome --softmasking --workingdir=. --threads $cpus --species=speciesName --gff3 --prot_seq=$prot_db  --useexisting
  mkdir Braker2_core
  bsub40000 -q long -e braker2core.err -o braker2core.err -n $cpus -R span[ptile=$cpus]  braker.pl --genome=$masked_genome --softmasking --workingdir=Braker2_core --threads $cpus --species=Lycoriella.core --gff3 --prot_seq=$prot_db  --useexisting

# Run Braker2 for GRC:
  cd /lustre/scratch125/casm/team268im/fa8/117/fede/Braker/Ling
  module add bamtools/2.5.2--hdcf5f25_3
  module add augustus/3.5.0--pl5321h700735d_3
  module add braker3/3.0.8--hdfd78af_0
  masked_genome=grc.fa
  cpus=16
  prot_db=/lustre/scratch125/casm/team268im/fa8/117/fede/Braker/Arthropoda+contarina.simple_header.fa
  export AUGUSTUS_CONFIG_PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/Augustus/config
  #export AUGUSTUS_BIN_PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/Augustus/bin
  export AUGUSTUS_SCRIPTS_PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/Augustus/scripts
  export GENEMARK_PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/GeneMark-ETP/bin
  export PROTHINT_PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/ProtHint/bin
  #braker.pl --genome=$masked_genome --softmasking --workingdir=. --threads $cpus --species=speciesName --gff3 --prot_seq=$prot_db  --useexisting
  mkdir Braker2_grc
  bsub40000 -q long -e braker2grc.err -o braker2grc.err -n $cpus -R span[ptile=$cpus]  braker.pl --genome=$masked_genome --softmasking --workingdir=Braker2_grc --threads $cpus --species=Lycoriella.grc --gff3 --prot_seq=$prot_db  --useexisting


# Compare annotations:
  /lustre/scratch125/casm/team268im/fa8/117/fede/software/gffcompare/gffcompare -r ./data/lycoriella_braker3_core/braker3/lycoriella/braker.gff3 -o braker2_vs_3.3asref.out ./Braker2_core/braker.gff3
  /lustre/scratch125/casm/team268im/fa8/117/fede/software/gffcompare/gffcompare -r ./Braker2_core/braker.gff3 -o braker2_vs_3.2asref.out ./data/lycoriella_braker3_core/braker3/lycoriella/braker.gff3
  

##########################################################################################
# Bradysia coprophila:
  cd /lustre/scratch125/casm/team268im/fa8/117/fede/Braker/Bcop
  cp -p /lustre/scratch126/tol/teams/jaron/users/sasha/diptera/Bradysia_coprophila/TEs/idBraCopr2.1.primary.masked.fa .
  rnaseq=/lustre/scratch125/casm/team268im/fa8/117/fede/RNAseq/Bcop/fastp/hisat2/All_rnaseq.bam

# Separate GRC and core genome
  cd /lustre/scratch125/casm/team268im/fa8/117/fede/Braker/Bcop
  perl /lustre/scratch125/casm/team268im/fa8/117/fede/EarlGrey/separate_genomes.pl idBraCopr2.1.primary.masked.fa

# Run Braker3 for core:
  cd /lustre/scratch125/casm/team268im/fa8/117/fede/Braker/Bcop
  mkdir logs
  # Copy everything onto data:
  mkdir data
  mkdir data/bcop_braker3_core
  mkdir data/bcop_braker3_core/braker3
  mkdir data/bcop_braker3_core/braker3/bcop
  bsub < braker_job.bcop.sh

# Run Braker2 for GRC:
  cd /lustre/scratch125/casm/team268im/fa8/117/fede/Braker/Bcop
  module add bamtools/2.5.2--hdcf5f25_3
  module add augustus/3.5.0--pl5321h700735d_3
  module add braker3/3.0.8--hdfd78af_0
  masked_genome=grc.fa
  cpus=16
  prot_db=/lustre/scratch125/casm/team268im/fa8/117/fede/Braker/Arthropoda+contarina.simple_header.fa
  export AUGUSTUS_CONFIG_PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/Augustus/config
  #export AUGUSTUS_BIN_PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/Augustus/bin
  export AUGUSTUS_SCRIPTS_PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/Augustus/scripts
  export GENEMARK_PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/GeneMark-ETP/bin
  export PROTHINT_PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/ProtHint/bin
  mkdir Braker2_grc
  bsub40000 -q long -e braker2grc.err -o braker2grc.err -n $cpus -R span[ptile=$cpus]  braker.pl --genome=$masked_genome --softmasking --workingdir=Braker2_grc --threads $cpus --species=Bcop.grc --gff3 --prot_seq=$prot_db  --useexisting



##########################################################################################
# Bradysia impatiens:
  cd /lustre/scratch125/casm/team268im/fa8/117/fede/Braker/Bimp
  # reps=/lustre/scratch126/tol/teams/jaron/projects/springtails_haploid_selection/data/results/bradysia_earlgrey/bradysiaImpatiens_EarlGrey/bradysiaImpatiens_summaryFiles/bradysiaImpatiens.filteredRepeats.bed
  # unmasked_genome=/lustre/scratch126/tol/teams/jaron/data/assemblies_Sanger/insects/Bradysia_impatiens/assembly/curated/idBraImpa2.1/idBraImpa2.1.primary.curated.fa
  # masked_genome=idBraImpa2.1.primary.curated.masked.fa
  # 
  # cd /lustre/scratch125/casm/team268im/fa8/117/fede/Braker/Bimp
  # bedtools maskfasta -soft -fi $unmasked_genome -bed $reps -fo $masked_genome

  masked_genome=/lustre/scratch126/tol/teams/jaron/data/diptera/Bradysia_impatiens/genome/idBraImpa2.1.primary.masked.fa
  perl /lustre/scratch125/casm/team268im/fa8/117/fede/EarlGrey/separate_genomes.pl $masked_genome
  
# Run Braker2 for GRC:
  cd /lustre/scratch125/casm/team268im/fa8/117/fede/Braker/Bimp
  module add bamtools/2.5.2--hdcf5f25_3
  module add augustus/3.5.0--pl5321h700735d_3
  module add braker3/3.0.8--hdfd78af_0
  masked_genome=grc.fa
  cpus=16
  prot_db=/lustre/scratch125/casm/team268im/fa8/117/fede/Braker/Arthropoda+contarina.simple_header.fa
  export AUGUSTUS_CONFIG_PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/Augustus/config
  #export AUGUSTUS_BIN_PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/Augustus/bin
  export AUGUSTUS_SCRIPTS_PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/Augustus/scripts
  export GENEMARK_PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/GeneMark-ETP/bin
  export PROTHINT_PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/ProtHint/bin
  mkdir Braker2_grc
  bsub40000 -q long -e braker2grc.err -o braker2grc.err -n $cpus -R span[ptile=$cpus]  braker.pl --genome=$masked_genome --softmasking --workingdir=Braker2_grc --threads $cpus --species=Bimp.grc --gff3 --prot_seq=$prot_db  --useexisting

# Run Braker2 for core:
  cd /lustre/scratch125/casm/team268im/fa8/117/fede/Braker/Bimp
  module add bamtools/2.5.2--hdcf5f25_3
  module add augustus/3.5.0--pl5321h700735d_3
  module add braker3/3.0.8--hdfd78af_0
  masked_genome=core.fa
  cpus=16
  prot_db=/lustre/scratch125/casm/team268im/fa8/117/fede/Braker/Arthropoda+contarina.simple_header.fa
  export AUGUSTUS_CONFIG_PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/Augustus/config
  #export AUGUSTUS_BIN_PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/Augustus/bin
  export AUGUSTUS_SCRIPTS_PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/Augustus/scripts
  export GENEMARK_PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/GeneMark-ETP/bin
  export PROTHINT_PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/ProtHint/bin
  mkdir Braker2_core
  bsub40000 -q long -e braker2core.err -o braker2core.err -n $cpus -R span[ptile=$cpus]  braker.pl --genome=$masked_genome --softmasking --workingdir=Braker2_core --threads $cpus --species=Bimp.core --gff3 --prot_seq=$prot_db  --useexisting

  cd /lustre/scratch125/casm/team268im/fa8/117/fede/Braker/Bimp
  module add bamtools/2.5.2--hdcf5f25_3
  module add augustus/3.5.0--pl5321h700735d_3
  module add braker3/3.0.8--hdfd78af_0
  masked_genome=grcSuper.fa
  cpus=16
  prot_db=/lustre/scratch125/casm/team268im/fa8/117/fede/Braker_onlyArthropoda/Arthropoda+contarina.simple_header.fa
  export AUGUSTUS_CONFIG_PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/Augustus/config
  #export AUGUSTUS_BIN_PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/Augustus/bin
  export AUGUSTUS_SCRIPTS_PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/Augustus/scripts
  export GENEMARK_PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/GeneMark-ETP/bin
  export PROTHINT_PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/ProtHint/bin
  mkdir Braker2_grcChrLevel
  #bsub40000 -q long -e braker2grc.err -o braker2grc.err -n $cpus -R span[ptile=$cpus]  braker.pl --genome=$masked_genome --softmasking --workingdir=Braker2_grc --threads $cpus --species=Bimp.grc --gff3 --prot_seq=$prot_db  --useexisting
  bsub40000 -q long -e braker2grcChrLevel.err -o braker2grcChrLevel.err -n $cpus -R span[ptile=$cpus]  braker.pl --genome=$masked_genome --softmasking --workingdir=Braker2_grcChrLevel --threads $cpus --species=Bimp.grc --gff3 --prot_seq=$prot_db  --useexisting

# Try braker2 with singularity:
  cd /lustre/scratch125/casm/team268im/fa8/117/fede/Braker/Bimp
  mkdir logs
  # Copy everything onto data:
  mkdir data
  mkdir data/bimp_braker3_core
  mkdir data/bimp_braker3_core/braker3
  mkdir data/bimp_braker3_core/braker3/bimp
  bsub < braker_job.bimp.sh




##########################################################################################
# NOTES

https://github.com/samebdon/springtail_haploid_selection/blob/main/scripts/describe_gtf.py
module add ISG/python/3.12.3
python3.12 -m pip install pandas
python3.12 -m pip install docopt
python3.12 describe_gtf.py -f braker.gff3


# AGAT (I usually run agat on braker.gtf's to tidy them up)
https://agat.readthedocs.io/en/latest/
 module load ISG/singularity/3.11.4
 singularity pull docker://quay.io/biocontainers/agat:1.0.0--pl5321hdfd78af_0
 singularity run agat_1.0.0--pl5321hdfd78af_0.sif



## export PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/Augustus/bin:/lustre/scratch125/casm/team268im/fa8/119/bin/Augustus/scripts:$PATH
## module add braker2/2.1.6--hdfd78af_5
## module add bamtools/2.5.2--hdcf5f25_3
## # module add augustus/3.5.0--pl5321h700735d_3
## # module add braker3/3.0.8--hdfd78af_0
## speciesName=Orobi
## cd /lustre/scratch125/casm/team268im/fa8/117/fede/Braker/$speciesName
## masked_genome=/lustre/scratch125/casm/team268im/fa8/117/fede/EarlGrey/$speciesName/genome.simple_header.earlGrey_masked.fa
## cpus=32
## prot_db=/lustre/scratch125/casm/team268im/fa8/117/fede/Braker/Arthropoda+contarina.fa
## bsub40000 -q long -e braker.err -o braker.err -n $cpus -R span[ptile=$cpus]  braker.pl --genome=$masked_genome --softmasking --workingdir=. --threads $cpus --species=speciesName --gff3 --prot_seq=$prot_db --useexisting --AUGUSTUS_CONFIG_PATH=.



# OR
# Braker3 (with prot db + RNAseq)
braker.pl --genome=genome.simple_header.earlGrey_masked.fa --softmasking --workingdir=. --threads cpus --species=speciesName --gff3 --prot_seq=prot_db.fa --useexisting --bam=sample1.rnaseq.query_sorted.bam,sample2.rnaseq.query_sorted.bam 



# AGAT (I usually run agat on braker.gtf's to tidy them up)
# https://agat.readthedocs.io/en/latest/ 
# Bioinformatics Web Server - University of Greifswald
# 




