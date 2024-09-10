#!/usr/bin/perl -w

use strict;
$|=1;

# Recibe un species name (bcop, bimp or ling), fichero fai (for chr names and lengths), 
# un chr name (e.g. SUPER_GRC_1), y un gff file con las localizaciones de los genes (e.g 
# bcop-grc.gff3)
# 
# El origen (cecidomyiidae or sciaridae) lo saca de ALL_CLASSIFICATIONS.ok
# 
# Puedo tratar de simplificar estas cosas
# Y podría generar esto para todos los cromosomas con un tamaño > 1Mb y que contengan
# la palabrita GRC
# Y las localizaciones de los genes ya sé dónde están

my $bcop_grc_gff3 = "/lustre/scratch126/tol/teams/jaron/users/fede/Braker_2nd_round/Bcop/Braker2_grc/braker.gff3";
my $bimp_grc_gff3 = "/lustre/scratch126/tol/teams/jaron/users/fede/Braker_2nd_round/Bimp/data/bimp_braker3_core/braker3/bimp/braker.gff3";
my $ling_grc_gff3 = "/lustre/scratch126/tol/teams/jaron/users/fede/Braker_2nd_round/Ling/Braker2_grc/braker.gff3";

# Load gene locations
my %bcop;
my %bcopgene2chr;
open(I, "gff2bed < $bcop_grc_gff3 | grep -w \"gene\" | cut -f1-4 |") || die "Noalallal\n";
while(<I>) {
	chomp;
	my($chr,$start,$end,$gene) = split;
	next if($chr =~ /unloc/);
	$gene = "$gene.t1";
	$bcop{$chr}->{$gene} = "$start:$end"; #($start+$end)/2;
	$bcopgene2chr{$gene} = $chr;
}
close(I);
my %bimp;
my %bimpgene2chr;
open(I, "gff2bed < $bimp_grc_gff3 | grep -w \"gene\" | cut -f1-4 |") || die "Noalallal\n";
while(<I>) {
	chomp;
	my($chr,$start,$end,$gene) = split;
	next if($chr =~ /unloc/);
	$gene = "$gene.t1";
	$bimp{$chr}->{$gene} = "$start:$end"; #($start+$end)/2;
	$bimpgene2chr{$gene} = $chr;
}
close(I);
my %ling;
my %linggene2chr;
open(I, "gff2bed < $ling_grc_gff3 | grep -w \"gene\" | cut -f1-4 |") || die "Noalallal\n";
while(<I>) {
	chomp;
	my($chr,$start,$end,$gene) = split;
	next if($chr =~ /unloc/);
	$gene = "$gene.t1";
	$ling{$chr}->{$gene} = "$start:$end"; #($start+$end)/2;
	$linggene2chr{$gene} = $chr;
}
close(I);

print STDERR scalar(keys(%bcopgene2chr)), " genes in bcop\n";
print STDERR scalar(keys(%bimpgene2chr)), " genes in bimp\n";
print STDERR scalar(keys(%linggene2chr)), " genes in ling\n";

print STDERR keys %bcop, " in bcop\n";
print STDERR keys %bimp, " in bimp\n";
print STDERR keys %ling, " in ling\n";

# Load origins:
my %origins;
open(I,"/lustre/scratch126/tol/teams/jaron/users/fede/TreeClassifications/ALL_CLASSIFICATIONS.ok") || die "ñajañjaña\n";
while(<I>) {
	my($sp,$gene,$origin) = (split(/\t/,$_))[1,2,3];
	my($chr,$location) = ("","");
	if($sp eq "bcop") {
		$chr = $bcopgene2chr{$gene};
		$location = $bcop{$chr}->{$gene};
	} elsif($sp eq "bimp") {
		$chr = $bimpgene2chr{$gene};
		$location = $bimp{$chr}->{$gene};
	} elsif($sp eq "ling") {
		$chr = $linggene2chr{$gene};
		$location = $ling{$chr}->{$gene};
	}
	my($start,$end) = split(/:/,$location);
	$start = $start + 1;
	print "$sp\t$chr\t$start\t$end\t$gene\t$origin\n";
	#$origins{$sp}->{$gene} = $origin;
}
close(I);


__END__



cd Proteomes_2nd_round
cp -p /lustre/scratch126/tol/teams/jaron/users/fede/Braker_2nd_round/Bcop/Braker2_grc/braker.aa bcop_grc.fa
cp -p /lustre/scratch126/tol/teams/jaron/users/fede/Braker_2nd_round/Bcop/data/bcop_braker3_core/braker3/bcop/braker.aa bcop_core.fa
cp -p /lustre/scratch126/tol/teams/jaron/users/fede/Braker_2nd_round/Orobi/braker.aa orobi.fa
cp -p /lustre/scratch126/tol/teams/jaron/users/fede/Braker_2nd_round/Ling/Braker2_grc/braker.aa ling_grc.fa
cp -p /lustre/scratch126/tol/teams/jaron/users/fede/Braker_2nd_round/Ling/data/lycoriella_braker3_core/braker3/lycoriella/braker.aa ling_core.fa
cp -p /lustre/scratch126/tol/teams/jaron/users/fede/Braker_2nd_round/Bimp/Braker2_core/braker.aa bimp_core.fa
# cp -p /lustre/scratch126/tol/teams/jaron/users/fede/Braker_2nd_round/Bimp/Braker2_grc/braker.aa # This failed
cp -p /lustre/scratch126/tol/teams/jaron/users/fede/Braker_2nd_round/Bimp/data/bimp_braker3_core/braker3/bimp/braker.aa bimp_grc.fa
cp -p /lustre/scratch126/tol/teams/jaron/users/fede/Braker_2nd_round/Aaphi/braker.aa aaphi.fa
cp -p ../Proteomes/contarinia.fa .
cp -p ../Proteomes/dmel.fa .
cp -p /lustre/scratch126/tol/teams/jaron/users/fede/Braker_2nd_round/GCA_029228625.1_FCFRP_Bhyg_1.0_translated_cds.faa phyg.fa


