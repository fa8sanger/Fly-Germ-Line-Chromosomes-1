# Braker second round, including: 
# Contarinia, Bcop core, Ling core, pseudolycoriella hygida

cdt
mkdir Braker_2nd_round
cd Braker_2nd_round
bcop_prots=/lustre/scratch126/tol/teams/jaron/users/fede/Braker/Bcop/data/bcop_braker3_core/braker3/bcop/braker.aa
ling_prots=/lustre/scratch126/tol/teams/jaron/users/fede/Braker/Ling/data/lycoriella_braker3_core/braker3/lycoriella/braker.aa
# pseudolycoriella hygida
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/029/228/625/GCA_029228625.1_FCFRP_Bhyg_1.0/GCA_029228625.1_FCFRP_Bhyg_1.0_genomic.gff.gz .
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/029/228/625/GCA_029228625.1_FCFRP_Bhyg_1.0/GCA_029228625.1_FCFRP_Bhyg_1.0_cds_from_genomic.fna.gz 
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/029/228/625/GCA_029228625.1_FCFRP_Bhyg_1.0/GCA_029228625.1_FCFRP_Bhyg_1.0_translated_cds.faa.gz
gunzip GCA_029228625.1_FCFRP_Bhyg_1.0_translated_cds.faa.gz
awk '{print $1}'  GCA_029228625.1_FCFRP_Bhyg_1.0_translated_cds.faa >  GCA_029228625.1_FCFRP_Bhyg_1.0_translated_cds.simple_header.faa
phyg_prots=GCA_029228625.1_FCFRP_Bhyg_1.0_translated_cds.simple_header.faa

perl -pi -e 's/^>lcl\|/>phyg-/' GCA_029228625.1_FCFRP_Bhyg_1.0_translated_cds.simple_header.faa
cp $bcop_prots ./bcop_prots.fa
perl -pi -e 's/^>/>bcop-/' bcop_prots.fa
cp $ling_prots ./ling_prots.fa
perl -pi -e 's/^>/>ling-/' ling_prots.fa

cp -p ../Braker/Arthropoda+contarina.simple_header.fa .
cat Arthropoda+contarina.simple_header.fa > Arth+cont+bcop+ling+phyg.fa
echo >> Arth+cont+bcop+ling+phyg.fa
cat bcop_prots.fa >> Arth+cont+bcop+ling+phyg.fa
echo >> Arth+cont+bcop+ling+phyg.fa
cat ling_prots.fa >> Arth+cont+bcop+ling+phyg.fa
echo >> Arth+cont+bcop+ling+phyg.fa
cat $phyg_prots >> Arth+cont+bcop+ling+phyg.fa

cp -p Arth+cont+bcop+ling+phyg.fa Arth+cont+bcop+ling+phyg.nospecialchars.fa
# perl -pi -e 's/:/_/g' Arth+cont+bcop+ling+phyg.nospecialchars.fa
# perl -pi -e 's/\|/_/g' Arth+cont+bcop+ling+phyg.nospecialchars.fa
perl -pi -e 's/[\t ]//g' Arth+cont+bcop+ling+phyg.nospecialchars.fa
perl -pi -e 's/[^A-Za-z0-9\.\r\n\>\*]/_/g' Arth+cont+bcop+ling+phyg.nospecialchars.fa


##########################################################################################
# Aaphi
  module add bamtools/2.5.2--hdcf5f25_3
  module add augustus/3.5.0--pl5321h700735d_3
  module add braker3/3.0.8--hdfd78af_0
  speciesName=Aaphi
  cd /lustre/scratch125/casm/team268im/fa8/117/fede/Braker_2nd_round/
  mkdir $speciesName
  cd $speciesName
  masked_genome=/lustre/scratch125/casm/team268im/fa8/117/fede/EarlGrey/$speciesName/genome.simple_header.earlGrey_masked.fa
  cpus=16
  prot_db=/lustre/scratch125/casm/team268im/fa8/117/fede/Braker_2nd_round/Arth+cont+bcop+ling+phyg.nospecialchars.fa
  export AUGUSTUS_CONFIG_PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/Augustus/config
  #export AUGUSTUS_BIN_PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/Augustus/bin
  export AUGUSTUS_SCRIPTS_PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/Augustus/scripts
  export GENEMARK_PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/GeneMark-ETP/bin
  export PROTHINT_PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/ProtHint/bin
  #braker.pl --genome=$masked_genome --softmasking --workingdir=. --threads $cpus --species=speciesName --gff3 --prot_seq=$prot_db  --useexisting
  bsub40000 -q normal -e braker.err -o braker.err -n $cpus -R span[ptile=$cpus]  braker.pl --genome=$masked_genome --softmasking --workingdir=. --threads $cpus --species=$speciesName --gff3 --prot_seq=$prot_db  --useexisting

# Orobi
  module add bamtools/2.5.2--hdcf5f25_3
  module add augustus/3.5.0--pl5321h700735d_3
  module add braker3/3.0.8--hdfd78af_0
  speciesName=Orobi
  cd /lustre/scratch125/casm/team268im/fa8/117/fede/Braker_2nd_round/
  mkdir $speciesName
  cd $speciesName
  masked_genome=/lustre/scratch125/casm/team268im/fa8/117/fede/EarlGrey/$speciesName/genome.simple_header.earlGrey_masked.fa
  cpus=16
  prot_db=/lustre/scratch125/casm/team268im/fa8/117/fede/Braker_2nd_round/Arth+cont+bcop+ling+phyg.nospecialchars.fa
  export AUGUSTUS_CONFIG_PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/Augustus/config
  #export AUGUSTUS_BIN_PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/Augustus/bin
  export AUGUSTUS_SCRIPTS_PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/Augustus/scripts
  export GENEMARK_PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/GeneMark-ETP/bin
  export PROTHINT_PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/ProtHint/bin
  #braker.pl --genome=$masked_genome --softmasking --workingdir=. --threads $cpus --species=speciesName --gff3 --prot_seq=$prot_db  --useexisting
  bsub40000 -q long -e braker.err -o braker.err -n $cpus -R span[ptile=$cpus]  braker.pl --genome=$masked_genome --softmasking --workingdir=. --threads $cpus --species=$speciesName --gff3 --prot_seq=$prot_db  --useexisting


##########################################################################################
# Ling
  cd /lustre/scratch125/casm/team268im/fa8/117/fede/Braker_2nd_round
  mkdir Ling
  cd Ling
  cp -p /lustre/scratch125/casm/team268im/fa8/117/fede/Braker/Ling/*fa .

# Ling core genome
  cd /lustre/scratch125/casm/team268im/fa8/117/fede/Braker_2nd_round/Ling
  mkdir logs
  # Create data DIR structure and copy files within... 
  bsub < braker_job.sh

# Braker2 for GRC:
  cd /lustre/scratch125/casm/team268im/fa8/117/fede/Braker_2nd_round/Ling
  module add bamtools/2.5.2--hdcf5f25_3
  module add augustus/3.5.0--pl5321h700735d_3
  module add braker3/3.0.8--hdfd78af_0
  masked_genome=grc.fa
  cpus=16
  prot_db=/lustre/scratch125/casm/team268im/fa8/117/fede/Braker_2nd_round/Arth+cont+bcop+ling+phyg.nospecialchars.fa
  export AUGUSTUS_CONFIG_PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/Augustus/config
  #export AUGUSTUS_BIN_PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/Augustus/bin
  export AUGUSTUS_SCRIPTS_PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/Augustus/scripts
  export GENEMARK_PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/GeneMark-ETP/bin
  export PROTHINT_PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/ProtHint/bin
  #braker.pl --genome=$masked_genome --softmasking --workingdir=. --threads $cpus --species=speciesName --gff3 --prot_seq=$prot_db  --useexisting
  mkdir Braker2_grc
  bsub40000 -q long -e braker2grc.err -o braker2grc.err -n $cpus -R span[ptile=$cpus]  braker.pl --genome=$masked_genome --softmasking --workingdir=Braker2_grc --threads $cpus --species=Lycoriella.grc --gff3 --prot_seq=$prot_db  --useexisting


##########################################################################################
# Bcop
  cd /lustre/scratch125/casm/team268im/fa8/117/fede/Braker_2nd_round
  mkdir Bcop
  cd Bcop
  cp -p /lustre/scratch125/casm/team268im/fa8/117/fede/Braker/Bcop/*fa .

# Bcop core genome
  cd /lustre/scratch125/casm/team268im/fa8/117/fede/Braker_2nd_round/Bcop
  mkdir logs
  # Create data DIR structure and copy files within... 
  bsub < braker_job.bcop.sh

# Bcop for GRC:
  cd /lustre/scratch125/casm/team268im/fa8/117/fede/Braker_2nd_round/Bcop
  module add bamtools/2.5.2--hdcf5f25_3
  module add augustus/3.5.0--pl5321h700735d_3
  module add braker3/3.0.8--hdfd78af_0
  masked_genome=grc.fa
  cpus=16
  prot_db=/lustre/scratch125/casm/team268im/fa8/117/fede/Braker_2nd_round/Arth+cont+bcop+ling+phyg.nospecialchars.fa
  export AUGUSTUS_CONFIG_PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/Augustus/config
  #export AUGUSTUS_BIN_PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/Augustus/bin
  export AUGUSTUS_SCRIPTS_PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/Augustus/scripts
  export GENEMARK_PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/GeneMark-ETP/bin
  export PROTHINT_PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/ProtHint/bin
  #braker.pl --genome=$masked_genome --softmasking --workingdir=. --threads $cpus --species=speciesName --gff3 --prot_seq=$prot_db  --useexisting
  mkdir Braker2_grc
  bsub40000 -q long -e braker2grc.err -o braker2grc.err -n $cpus -R span[ptile=$cpus]  braker.pl --genome=$masked_genome --softmasking --workingdir=Braker2_grc --threads $cpus --species=Bcop.grc --gff3 --prot_seq=$prot_db  --useexisting


##########################################################################################
# Bimp
  cd /lustre/scratch125/casm/team268im/fa8/117/fede/Braker_2nd_round
  mkdir Bimp
  cd Bimp
  masked_genome=/lustre/scratch126/tol/teams/jaron/data/diptera/Bradysia_impatiens/genome/idBraImpa2.1.primary.masked.fa
  perl /lustre/scratch125/casm/team268im/fa8/117/fede/EarlGrey/separate_genomes.pl $masked_genome
  
# Run Braker2 for GRC:
  cd /lustre/scratch125/casm/team268im/fa8/117/fede/Braker_2nd_round/Bimp
  module add bamtools/2.5.2--hdcf5f25_3
  module add augustus/3.5.0--pl5321h700735d_3
  module add braker3/3.0.8--hdfd78af_0
  masked_genome=grc.fa
  cpus=16
  prot_db=/lustre/scratch125/casm/team268im/fa8/117/fede/Braker_2nd_round/Arth+cont+bcop+ling+phyg.nospecialchars.fa
  export AUGUSTUS_CONFIG_PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/Augustus/config
  #export AUGUSTUS_BIN_PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/Augustus/bin
  export AUGUSTUS_SCRIPTS_PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/Augustus/scripts
  export GENEMARK_PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/GeneMark-ETP/bin
  export PROTHINT_PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/ProtHint/bin
  mkdir Braker2_grc
  #bsub40000 -q long -e braker2grc.err -o braker2grc.err -n $cpus -R span[ptile=$cpus]  braker.pl --genome=$masked_genome --softmasking --workingdir=Braker2_grc --threads $cpus --species=Bimp.grc --gff3 --prot_seq=$prot_db  --useexisting
  bsub40000 -q long -e braker2grc.err -o braker2grc.err -n $cpus -R span[ptile=$cpus]  braker.pl --genome=$masked_genome --softmasking --workingdir=Braker2_grc --threads $cpus --species=Bimp.grc --prot_seq=$prot_db  --useexisting

# Run Braker2 for core:
  cd /lustre/scratch125/casm/team268im/fa8/117/fede/Braker_2nd_round/Bimp
  module add bamtools/2.5.2--hdcf5f25_3
  module add augustus/3.5.0--pl5321h700735d_3
  module add braker3/3.0.8--hdfd78af_0
  masked_genome=core.fa
  cpus=16
  prot_db=/lustre/scratch125/casm/team268im/fa8/117/fede/Braker_2nd_round/Arth+cont+bcop+ling+phyg.nospecialchars.fa
  export AUGUSTUS_CONFIG_PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/Augustus/config
  #export AUGUSTUS_BIN_PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/Augustus/bin
  export AUGUSTUS_SCRIPTS_PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/Augustus/scripts
  export GENEMARK_PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/GeneMark-ETP/bin
  export PROTHINT_PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/ProtHint/bin
  mkdir Braker2_core
  bsub40000 -q long -e braker2core.err -o braker2core.err -n $cpus -R span[ptile=$cpus]  braker.pl --genome=$masked_genome --softmasking --workingdir=Braker2_core --threads $cpus --species=Bimp.core --gff3 --prot_seq=$prot_db  --useexisting

   cd /lustre/scratch125/casm/team268im/fa8/117/fede/Braker_2nd_round/Bimp
  module add bamtools/2.5.2--hdcf5f25_3
  module add augustus/3.5.0--pl5321h700735d_3
  module add braker3/3.0.8--hdfd78af_0
  masked_genome=grcSuper.fa
  cpus=16
  prot_db=/lustre/scratch125/casm/team268im/fa8/117/fede/Braker_2nd_round/Arth+cont+bcop+ling+phyg.nospecialchars.fa
  export AUGUSTUS_CONFIG_PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/Augustus/config
  #export AUGUSTUS_BIN_PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/Augustus/bin
  export AUGUSTUS_SCRIPTS_PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/Augustus/scripts
  export GENEMARK_PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/GeneMark-ETP/bin
  export PROTHINT_PATH=/lustre/scratch125/casm/team268im/fa8/119/bin/ProtHint/bin
  mkdir Braker2_grcChrLevel
  #bsub40000 -q long -e braker2grc.err -o braker2grc.err -n $cpus -R span[ptile=$cpus]  braker.pl --genome=$masked_genome --softmasking --workingdir=Braker2_grc --threads $cpus --species=Bimp.grc --gff3 --prot_seq=$prot_db  --useexisting
  bsub40000 -q long -e braker2grcChrLevel.err -o braker2grcChrLevel.err -n $cpus -R span[ptile=$cpus]  braker.pl --genome=$masked_genome --softmasking --workingdir=Braker2_grcChrLevel --threads $cpus --species=Bimp.grc --gff3 --prot_seq=$prot_db  --useexisting

# Try braker2 with singularity:
  cd /lustre/scratch125/casm/team268im/fa8/117/fede/Braker_2nd_round/Bimp
  mkdir logs
  # Copy everything onto data:
  mkdir data
  mkdir data/bimp_braker3_core
  mkdir data/bimp_braker3_core/braker3
  mkdir data/bimp_braker3_core/braker3/bimp
  bsub < braker_job.bimp.sh




# Monitoring:
 cd /lustre/scratch125/casm/team268im/fa8/117/fede/Braker_2nd_round
 tail -f */*err */logs/*err*

cd /lustre/scratch126/tol/teams/jaron/users/fede/Braker_2nd_round
find . -name "braker.gtf"
python3.12 analyze_exons.py -f ./Bcop/Braker2_grc/braker.gtf
python3.12 analyze_exons.py -f ./Bcop/data/bcop_braker3_core/braker3/bcop/braker.gtf
python3.12 analyze_exons.py -f ./Orobi/braker.gtf
python3.12 analyze_exons.py -f ./Ling/Braker2_grc/braker.gtf
python3.12 analyze_exons.py -f ./Ling/data/lycoriella_braker3_core/braker3/lycoriella/braker.gtf
python3.12 analyze_exons.py -f ./Bimp/Braker2_core/braker.gtf
python3.12 analyze_exons.py -f ./Bimp/Braker2_grc/braker.gtf
python3.12 analyze_exons.py -f ./Aaphi/braker.gtf
cd /lustre/scratch126/tol/teams/jaron/users/fede/Braker
find . -name "braker.gtf"
echo ./Bcop/Braker2_grc/braker.gtf                                       ; python3.12 analyze_exons.py -f ./Bcop/Braker2_grc/braker.gtf
echo ./Bcop/data/bcop_braker3_core/braker3/bcop/braker.gtf               ; python3.12 analyze_exons.py -f ./Bcop/data/bcop_braker3_core/braker3/bcop/braker.gtf
echo ./Orobi/braker.gtf                                                  ; python3.12 analyze_exons.py -f ./Orobi/braker.gtf
echo ./Ling/Braker2_core/braker.gtf                                      ; python3.12 analyze_exons.py -f ./Ling/Braker2_core/braker.gtf
echo ./Ling/Braker2_grc/braker.gtf                                       ; python3.12 analyze_exons.py -f ./Ling/Braker2_grc/braker.gtf
echo ./Ling/data/lycoriella_braker3_core/braker3/lycoriella/braker.gtf   ; python3.12 analyze_exons.py -f ./Ling/data/lycoriella_braker3_core/braker3/lycoriella/braker.gtf
echo ./Bimp/Braker2_core/braker.gtf                                      ; python3.12 analyze_exons.py -f ./Bimp/Braker2_core/braker.gtf
echo ./Bimp/Braker2_grc/braker.gtf                                       ; python3.12 analyze_exons.py -f ./Bimp/Braker2_grc/braker.gtf
echo ./Aaphi/braker.gtf                                                  ; python3.12 analyze_exons.py -f ./Aaphi/braker.gtf


##########################################################################################








