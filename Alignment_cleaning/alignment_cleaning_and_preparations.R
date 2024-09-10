# This script / NOTES is to explain what to do and prepare the data

cdt
mkdir CDS_seqs
cd CDS_seqs
cp -p /lustre/scratch126/tol/teams/jaron/users/fede/Braker_2nd_round/Bcop/Braker2_grc/braker.codingseq bcop_grc.fa
cp -p /lustre/scratch126/tol/teams/jaron/users/fede/Braker_2nd_round/Bcop/data/bcop_braker3_core/braker3/bcop/braker.codingseq bcop_core.fa
cp -p /lustre/scratch126/tol/teams/jaron/users/fede/Braker_2nd_round/Orobi/braker.codingseq orobi.fa
cp -p /lustre/scratch126/tol/teams/jaron/users/fede/Braker_2nd_round/Ling/Braker2_grc/braker.codingseq ling_grc.fa
cp -p /lustre/scratch126/tol/teams/jaron/users/fede/Braker_2nd_round/Ling/data/lycoriella_braker3_core/braker3/lycoriella/braker.codingseq ling_core.fa
cp -p /lustre/scratch126/tol/teams/jaron/users/fede/Braker_2nd_round/Bimp/Braker2_core/braker.codingseq bimp_core.fa
cp -p /lustre/scratch126/tol/teams/jaron/users/fede/Braker_2nd_round/Bimp/data/bimp_braker3_core/braker3/bimp/braker.codingseq bimp_grc.fa
cp -p /lustre/scratch126/tol/teams/jaron/users/fede/Braker_2nd_round/Aaphi/braker.codingseq aaphi.fa
wget http://ftp.flybase.net/genomes/Drosophila_melanogaster/current/fasta/dmel-all-CDS-r6.58.fasta.gz
gunzip dmel-all-CDS-r6.58.fasta.gz
mv dmel-all-CDS-r6.58.fasta dmel.fa
cp -p /lustre/scratch126/tol/teams/jaron/users/fede/GCF_009176525.2/cds_from_genomic.fna contarinia.fa
cp -p /lustre/scratch126/tol/teams/jaron/users/fede/Braker_2nd_round/GCA_029228625.1_FCFRP_Bhyg_1.0_cds_from_genomic.fna.gz phyg.fa.gz
gunzip phyg.fa.gz
perl -pi -e 's/ //g' dmel.fa # trick for longer seq names
# Now modify the script analyse_alignments_with_identity_v2.R so that it loads all of these sequences



* Analizar los alineamientos y eliminar las secuencias parciales o con demasiado unique
* Guardar esos alineamientos
* Buscar las secuencias de codones correspondientes. Esto va a ser challenging
    * Cargar todos los fasta en memoria y luego hacer greps de esas matrices. En R
    * Esta tarde: crear un directorio con todas las secuencias para cada especie
* Correr Prequal o algún otro para limpiar. También BGE (funciona a nivel de codones?)
    * /lustre/scratch126/tol/teams/jaron/users/fede/software/BMGE/src/BMGE.jar
    * /lustre/scratch126/tol/teams/jaron/users/fede/software/prequal
    * perl ~/perl5/bin/HmmCleaner.pl  OG0000329.fa

Secuencias de codones:
* Fácil para las anotadas con Braker: braker.codingseq
* Dmel: más complicado: 
* Contarinia: 
* Phyg: 

All proteomes used in Orthofinder 2nd round:
- cp -p /lustre/scratch126/tol/teams/jaron/users/fede/Braker_2nd_round/Bcop/Braker2_grc/braker.aa bcop_grc.fa
- cp -p /lustre/scratch126/tol/teams/jaron/users/fede/Braker_2nd_round/Bcop/data/bcop_braker3_core/braker3/bcop/braker.aa bcop_core.fa
- cp -p /lustre/scratch126/tol/teams/jaron/users/fede/Braker_2nd_round/Orobi/braker.aa orobi.fa
- cp -p /lustre/scratch126/tol/teams/jaron/users/fede/Braker_2nd_round/Ling/Braker2_grc/braker.aa ling_grc.fa
- cp -p /lustre/scratch126/tol/teams/jaron/users/fede/Braker_2nd_round/Ling/data/lycoriella_braker3_core/braker3/lycoriella/braker.aa ling_core.fa
- cp -p /lustre/scratch126/tol/teams/jaron/users/fede/Braker_2nd_round/Bimp/Braker2_core/braker.aa bimp_core.fa
- # cp -p /lustre/scratch126/tol/teams/jaron/users/fede/Braker_2nd_round/Bimp/Braker2_grc/braker.aa # This failed
- cp -p /lustre/scratch126/tol/teams/jaron/users/fede/Braker_2nd_round/Bimp/data/bimp_braker3_core/braker3/bimp/braker.aa bimp_grc.fa
- cp -p /lustre/scratch126/tol/teams/jaron/users/fede/Braker_2nd_round/Aaphi/braker.aa aaphi.fa
- cp -p ../Proteomes/contarinia.fa .
- cp -p ../Proteomes/dmel.fa .
- cp -p /lustre/scratch126/tol/teams/jaron/users/fede/Braker_2nd_round/GCA_029228625.1_FCFRP_Bhyg_1.0_translated_cds.faa phyg.fa
- dmel comes from: ../Blast/dmel-all-translation-r6.58.fasta.
    - Codon sequences here: http://ftp.flybase.net/genomes/Drosophila_melanogaster/current/fasta/dmel-all-CDS-r6.58.fasta.gz
- contarinia.fa comes from: /lustre/scratch126/tol/teams/jaron/users/fede/GCF_009176525.2/cds_from_genomic.fna

Now I have all I need to get those cds alignments… How to align CDS without using translatorx??
I think MACSE v2 could be a good choice: https://academic.oup.com/mbe/article/35/10/2582/5079334
PRANK has a codon model, but they say it is very slow
A MACSE pipeline combining with hmmcleaner: https://github.com/ranwez/MACSE_V2_PIPELINES
Code here: https://www.agap-ge2pop.org/macsee-pipelines/
Warning This programs mainly aims at removing long non-homologous fragments. Smaller ones will be kept, to avoid removing fragments that are really homologous at this step. If your analyses are sensitive to such non-homologous fragments (e.g. dN/dS estimation), we strongly advice to also use a post filtering of your alignment at the amino acid level (e.g. using HMMCleaner, BMGE or trimAl) and to report this AA masking/filtering at the nucleotide level using reportMaskAA2NT. —> I’ve got a strategy! nice
* java -jar macse_v2.07.jar -prog alignSequences
* java -jar macse_v2.07.jar -prog reportMaskAA2NT

Strategy:
* Retrieve codon sequences (removing problematic “genes”, partial or containing rubbish)
* Clean with prequal
* Align with macse
* Clean with BMGE
* Fix with macse using BMGE info
* Convert to phy and tree, dnds, etc…