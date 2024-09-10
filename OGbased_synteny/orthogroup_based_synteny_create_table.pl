#!/usr/bin/perl -w


use strict;
$|= 1;

# Example run:
# cdt
# perl orthogroup_based_synteny_create_table.pl > OG_BASED_SYNTENY_TABLE.tsv


# genomes: bcop_core bcop_grc bimp_core bimp_grc ling_core ling_grc orobi aaphi
# gff files in: /lustre/scratch126/tol/teams/jaron/users/fede/Genomes_SoftMasked

# * Puedo representar la sintenia usando la pertenencia a orthogroups… 
#     * así lo voy a hacer para las regiones X - GRC, por ejemplo. O para comparar los GRCs de uno y otro
#     * así también podemos estudiar la pérdida de genes
#     * Necesito la pertenencia a los OGs y las coordenadas de los genes, eso es todo. Easy peasy
#     * Podría obtener todas las coordenadas para todos los genomas, con un solo script
#         * 1. Leo coords de genes, guardo en genes{sp}->{gene} = coords. Y coords{sp}->{coords} = gene
#         * 2. Leo OGs y creo grafo de coords de genes: cc{sp1}->{sp2}->{‘coords1’} = coords2; cc{sp2}->{sp1}->{‘coords2’} = coords1
#         * 3. Guardo en un fichero tabulado todos los pares y ya podemos jugar en R

my $GENE_PATH = "/lustre/scratch126/tol/teams/jaron/users/fede/Genomes_SoftMasked";
my $OG_FILE   = "/lustre/scratch126/tol/teams/jaron/users/fede/Orthofinder/Proteomes_2nd_round_OneIsoform/OrthoFinder/Results_Jul25_1/Orthogroups/Orthogroups.tsv";
my @genomes = ("bcop_core","bcop_grc","bimp_core","bimp_grc","ling_core","ling_grc","orobi","aaphi");
my %genes;
my %coords;

# Load gene coords:
foreach my $genome ( @genomes ) {
	print STDERR "Reading $genome...\n";
	open(I,"$GENE_PATH/$genome.gene.bed") || die "ñkjasñkjañjjj $genome\n";
	while(<I>) {
		chomp;
		my($chr,$start,$end,$gene) = split;
		$genes{$genome}->{$gene} = "$chr:$start:$end";
		$coords{$genome}->{"$chr:$start:$end"} = $gene;
		#print STDERR "$genome $gene = $chr:$start:$end\n";
	}
	close(I);
}
print STDERR scalar(keys(%coords)), " in coords\n";

# Load OG info:
open(I, $OG_FILE) || die "aññañaiakakkakak\n";
my $header = <I>;
chomp($header);
my @columns = ("aaphi","bcop_core","bcop_grc","bimp_core","bimp_grc","contarinia","dmel","ling_core","ling_grc","orobi","phyg");
my %ogs;
while(<I>) {
	chomp;
	my @tmp = split(/\t/,$_);
	my $og = $tmp[0];
	for(my $i=1; $i<=$#tmp; $i++) {
		my $genome = $columns[$i-1];
		next unless(exists($genes{$genome}));
		my $genes_in_og = $tmp[$i];
		$genes_in_og =~ s/\.t1//g;
		if(!defined($genes_in_og) || $genes_in_og eq "") {
			next;
		}
		# print STDERR "Ready to add $genes_in_og for $genome\n";
		if($genes_in_og =~ /,/) {
			$genes_in_og =~ s/ +//g;
			my @tt = split(/,/,$genes_in_og);
			foreach my $t ( @tt ) {
				$ogs{$og}->{$genome}->{$t} = 1;			
			}
		} else {
			$ogs{$og}->{$genome}->{$genes_in_og} = 1;
		}
	}
}
close(I);
print STDERR scalar(keys(%ogs)), " in ogs\n";

# Print each OG:
foreach my $og ( keys %ogs ) {
	my @gs = keys(%{$ogs{$og}});
	foreach my $genome1 ( @gs ) {
		foreach my $genome2 ( @gs ) {
			next if($genome1 eq $genome2);
			foreach my $gene1 ( keys %{$ogs{$og}->{$genome1}} ) {
				foreach my $gene2 ( keys %{$ogs{$og}->{$genome2}} ) {
					# print STDERR "1: $genes{$genome1}->{$gene1} [$genome1 / $gene1]\n";
					my($chr1,$start1,$end1) = split(/:/,$genes{$genome1}->{$gene1});
					print "$og\t";
					print "$genome1\t$gene1\t$chr1\t$start1\t$end1\t";
					my($chr2,$start2,$end2) = split(/:/,$genes{$genome2}->{$gene2});
					print "$genome2\t$gene2\t$chr2\t$start2\t$end2\n";
				}
			}
		}
	}
}


__END__

