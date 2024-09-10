# Copy soft-masked genome files:
cd /lustre/scratch126/tol/teams/jaron/users/fede/Genomes_SoftMasked
cp -p ../Braker_2nd_round/Bcop/core.fa bcop.core.masked.fa
cp -p ../Braker_2nd_round/Bcop/grc.fa bcop.grc.masked.fa
cp -p ../Braker_2nd_round/Bimp/core.fa bimp.core.masked.fa
cp -p ../Braker_2nd_round/Bimp/grc.fa bimp.grc.masked.fa
cp -p ../Braker_2nd_round/Ling/core.fa ling.core.masked.fa
cp -p ../Braker_2nd_round/Ling/grc.fa ling.grc.masked.fa
cp -p ../EarlGrey/Orobi/genome.simple_header.earlGrey_masked.fa orobi.masked.fa
cp -p ../EarlGrey/Aaphi/genome.simple_header.earlGrey_masked.fa aaphi.masked.fa


# Format genomes to make blast databases:
install_dir=/lustre/scratch125/casm/team268im/fa8/119/bin/ncbi-blast-2.15.0+/bin
cd /lustre/scratch126/tol/teams/jaron/users/fede/Genomes_SoftMasked
for a in *fa 
do echo $a
# masking
perl -pi -e 's/[acgt]/n/g' $a
$install_dir/makeblastdb -in $a -parse_seqids -dbtype nucl
done

# Run blast using the drosophila proteins as input:
install_dir=/lustre/scratch125/casm/team268im/fa8/119/bin/ncbi-blast-2.15.0+/bin
cd /lustre/scratch126/tol/teams/jaron/users/fede/Dro_vs_Genomes
for a in /lustre/scratch126/tol/teams/jaron/users/fede/Genomes_SoftMasked/*fa 
 do echo $a
 cpus=8
 bsub10000 -q long -n $cpus -R span[ptile=$cpus]  -e $a.bls.err    -o $a.bls.err    $install_dir/tblastn -db $a -query dmel.fa -out $a.bls.out     -evalue 0.001           -num_threads $cpus -num_descriptions  10 -num_alignments  10 
 bsub10000 -q long -n $cpus -R span[ptile=$cpus]  -e $a.bls.f6.err -o $a.bls.f6.err $install_dir/tblastn -db $a -query dmel.fa -out $a.bls.f6.out  -evalue 0.001 -outfmt 6 -num_threads $cpus #-max_target_seqs  10 
done
# Use -db_hard_mask to avoid repeats


# Copy GFF3 files:
cd /lustre/scratch126/tol/teams/jaron/users/fede/Genomes_SoftMasked
cp -p /lustre/scratch126/tol/teams/jaron/users/fede/Braker_2nd_round/Bcop/Braker2_grc/braker.gff3 bcop_grc.gff3
cp -p /lustre/scratch126/tol/teams/jaron/users/fede/Braker_2nd_round/Bcop/data/bcop_braker3_core/braker3/bcop/braker.gff3 bcop_core.gff3
cp -p /lustre/scratch126/tol/teams/jaron/users/fede/Braker_2nd_round/Orobi/braker.gff3 orobi.gff3
cp -p /lustre/scratch126/tol/teams/jaron/users/fede/Braker_2nd_round/Ling/Braker2_grc/braker.gff3 ling_grc.gff3
cp -p /lustre/scratch126/tol/teams/jaron/users/fede/Braker_2nd_round/Ling/data/lycoriella_braker3_core/braker3/lycoriella/braker.gff3 ling_core.gff3
cp -p /lustre/scratch126/tol/teams/jaron/users/fede/Braker_2nd_round/Bimp/Braker2_core/braker.gff3 bimp_core.gff3
# cp -p /lustre/scratch126/tol/teams/jaron/users/fede/Braker_2nd_round/Bimp/Braker2_grc/braker.aa # This failed
cp -p /lustre/scratch126/tol/teams/jaron/users/fede/Braker_2nd_round/Bimp/data/bimp_braker3_core/braker3/bimp/braker.gff3 bimp_grc.gff3
cp -p /lustre/scratch126/tol/teams/jaron/users/fede/Braker_2nd_round/Aaphi/braker.gff3 aaphi.gff3

# Convert GFF3 files to gene beds:
cd /lustre/scratch126/tol/teams/jaron/users/fede/Genomes_SoftMasked
for sp in bcop_core bcop_grc bimp_core bimp_grc ling_core ling_grc orobi aaphi
do echo $sp
 gff2bed < $sp.gff3 | grep -w "gene" | cut -f1-4 > $sp.gene.bed
done

# Convert Blast f6 results to bed files:
for a in *f6.out
 do echo $a
 cat $a | awk '{ if($10>$9) { print $2,$9-1,$10,$1":"$11; } else { print $2,$10-1,$9,$1":"$11; }  }' | tr ' ' '\t' > $a.bed
 done




When finished: 
* Intersect bed of blast hits (incl. dmel gene name and evalue) with gff3 of gene predictions
* Subtract bed of hits with gff3 of gene predictions
* For each dmel gene, find regions with enough hit and evalue that are distant enough
  from regions with good intersection and evalue
  
Buscar hits de un gen que no caen en genes y para los cuales ese gen no tiene ningún hit solapando genes cercanos
Necesito hacer una especie de bed para cada gen cubriendo territorio génico. 
A ese bed le añado como 50Kb a cada lado
Ahora busco hits de ese gen que no caigan en territorio génico, tenga un evalue decente y 
no esté cerca de los otros hits de ese gen



