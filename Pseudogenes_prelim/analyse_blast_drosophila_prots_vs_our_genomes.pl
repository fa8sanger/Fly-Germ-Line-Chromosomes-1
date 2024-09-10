#!/usr/bin/perl -w
use strict;

# Example: 
# cd /lustre/scratch126/tol/teams/jaron/users/fede/Genomes_SoftMasked
# perl analyse_blast_drosophila_prots_vs_our_genomes.pl bimp_grc

# for sp in bcop_core bcop_grc bimp_core bimp_grc ling_core ling_grc orobi aaphi
#   do echo $sp
#      perl analyse_blast_drosophila_prots_vs_our_genomes.pl $sp > $sp.RESULTS
#   done

#*****************************************************************************************
#   For each non-gene hit of a query gene, find the closest gene-hit for the query gene  *    
#   If too far, calculate some metric...                                                 *   
#   And report in a bed file                                                             *  
#*****************************************************************************************

my $sp = $ARGV[0]; 
# aaphi.masked.fa.bls.f6.out.bed
# my $gen = $ARGV[1]; # aaphi.gene.bed
# my $out = $ARGV[2]; # aaphi
my $bls = "$sp.masked.fa.bls.f6.out.bed";
$bls =~ s/_core/.core/;
$bls =~ s/_grc/.grc/;
my $gen = "$sp.gene.bed";

# Find genic & non-genic blast hits:
print STDERR "Executing: intersectBed -wo -a $bls -b $gen  > $sp.gene_hits.bed\n";
print STDERR "Executing: intersectBed -v  -a $bls -b $gen  > $sp.nongene_hits.bed\n";

system("intersectBed -wo -a $bls -b $gen  > $sp.gene_hits.bed   ");
system("intersectBed -v  -a $bls -b $gen  > $sp.nongene_hits.bed");

# Merge the nongene_hits... That would make things easier to process
print STDERR "sortBed  -i $sp.nongene_hits.bed                         > $sp.nongene_hits.sorted.bed\n";
print STDERR "mergeBed -i $sp.nongene_hits.sorted.bed -c 4 -o collapse > $sp.nongene_hits.merged.bed\n";

system("sortBed  -i $sp.nongene_hits.bed                         > $sp.nongene_hits.sorted.bed");
system("mergeBed -i $sp.nongene_hits.sorted.bed -c 4 -o collapse > $sp.nongene_hits.merged.bed");

# Merge the hits... That would make things easier to process
print STDERR "sortBed  -i $sp.gene_hits.bed                         > $sp.gene_hits.sorted.bed\n";
print STDERR "mergeBed -i $sp.gene_hits.sorted.bed -c 4 -o collapse > $sp.gene_hits.merged.bed\n";

system("sortBed  -i $sp.gene_hits.bed                         > $sp.gene_hits.sorted.bed");
system("mergeBed -i $sp.gene_hits.sorted.bed -c 4 -o collapse > $sp.gene_hits.merged.bed");

# Cargar los candidatos:
my %candidates;
my %candidate_genes;
open(I,"$sp.nongene_hits.merged.bed") || die "ñaañajjkjakjakj\n";
while(<I>) {
	chomp;
	my($chr,$start,$end,$genes) = split;
	my $id = "$chr:$start-$end";
	my @tmp = split(/,/,$genes);
	$candidates{$id} = 1;
	foreach my $t ( @tmp ) {
		my($g,$e) = (split(/:/,$t))[0,1];
		$candidate_genes{$id}->{$g} = 1;
	}
}
close(I);


# Ahora, para cada non-coding segment, busco los coding segments mas cercanos
# Si alguno de los coding segments cercanos son hits de alguno de los genes en el non-coding: discard
# Mi flankBed: 
my $MIN_DISTANCE = 50000; # empirically estimated... 50000 es una exageración pero bueno. Mejor 10000?
open(I, "$sp.gene_hits.merged.bed") || die "Noarlallajaj\n";
open(O, ">$sp.gene_hits.merged+$MIN_DISTANCE.bed") || die "!ñjaññaañ\n";
while(<I>) {
	chomp;
	my($chr,$start,$end,$info) = split;
	$start = $start-$MIN_DISTANCE;
	$end = $end+$MIN_DISTANCE;
	if($start < 0) {
		$start = 0;
	}
	print O "$chr\t$start\t$end\t$info\n";
}
close(O);
close(I);

# De los substraidos, intersect con $sp.gene_hits.merged+$MIN_DISTANCE.bed
system("intersectBed -wo -a $sp.nongene_hits.merged.bed -b $sp.gene_hits.merged+$MIN_DISTANCE.bed > $sp.nongene_hits.vs.flankgenehits.bed");

# Explore:
open(I, "$sp.nongene_hits.vs.flankgenehits.bed") || die "ñaañññkjkkkkkk\n";
while(<I>) {
	chomp;
	my($chr,$start,$end,$genes1,$genes2) = (split)[0,1,2,3,7];
	my $candidate_id = "$chr:$start-$end";
	my $overlap = 0;
	$genes1 .= ",";
	$genes2 .= ",";
	$genes1 =~ s/:.+?,/\t/g;
	$genes2 =~ s/:.+?,/\t/g;
	#print STDERR "$_\n$genes1\n$genes2\n\n\n";
	foreach my $g1 ( split(/\t/,$genes1) ) {
		foreach my $g2 ( split(/\t/,$genes2) ) {
			if($g1 eq $g2) {
				$overlap++;
			}
		}
	}
	if($overlap > 0) {
		$candidates{$candidate_id} = 0;
	} 
}
close(I);

my %counts_per_gene;
foreach my $id ( keys %candidates ) {
	my @tt = keys %{$candidate_genes{$id}};
	@" = ",";
	print "$id\t$candidates{$id}\t@tt\n";
	foreach my $t ( @tt ) {
		$counts_per_gene{$candidates{$id}}->{$t}++;
	}
}

foreach my $type ( keys %counts_per_gene ) {
	foreach my $gene ( keys %{$counts_per_gene{$type}} ) {
		print "GENE_INFO\t$type\t$gene\t$counts_per_gene{$type}->{$gene}\n";
	}
}

#
#sp.gene.bed: gene prediction coordinates
#bcop.core.masked.fa.bls.f6.out: blast output

__END__

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





