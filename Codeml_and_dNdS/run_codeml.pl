#!/usr/bin/perl -w

use strict;
# module add paml/4.10.7--h031d066_0

# Example: 
# perl ../run_codeml.pl ../../CDS_seqs/CleanAlignments/OG0005035/OG0005035.final.codon.fa ../../CDS_seqs/CleanAlignments/OG0005035/OG0005035.final.aa.fa.contree .

# Running all the cases:
# cd /lustre/scratch126/tol/teams/jaron/users/fede/CDS_seqs/CleanAlignments
# mkdir /lustre/scratch126/tol/teams/jaron/users/fede/CDS_seqs/CleanAlignments.kk
# module add paml/4.10.7--h031d066_0
# mkdir LOGS_CODEML
# # og=OG0012602
# for og in OG*
#  do echo $og
#  mkdir /lustre/scratch126/tol/teams/jaron/users/fede/CDS_seqs/CleanAlignments.kk/$og/
#  bsub5000 -e LOGS_CODEML/$og.err -o LOGS_CODEML/$og.err perl /lustre/scratch126/tol/teams/jaron/users/fede/CODEML/run_codeml.pl /lustre/scratch126/tol/teams/jaron/users/fede/CDS_seqs/CleanAlignments/$og/$og.final.codon.fa /lustre/scratch126/tol/teams/jaron/users/fede/CDS_seqs/CleanAlignments/$og/$og.final.aa.fa.contree /lustre/scratch126/tol/teams/jaron/users/fede/CDS_seqs/CleanAlignments.kk/$og/
#  done
##########################################################################################
# This script receives a CDS alignment and a tree and an output path
# It optimises branch lengths with codeml to make it codon-based rather than aa-based
# Then, extract the tree from the output
# Then, it reads the sequence names from the alignment and for each of them adds a "#1"
#   so that codeml runs the two omegas model and calculates the dN/dS of that tip
# Then reads the result and keeps the omegas
# Finally, write a TSV file with OG id, gene, omega... Later even classification
##########################################################################################

my $cds_ali        = $ARGV[0];
my $tree           = $ARGV[1];
my $output_path    = $ARGV[2];
my $INSTALL_DIR    = "/lustre/scratch126/tol/teams/jaron/users/fede/CODEML";
my $CODEML_BLEN    = "optim_blen.ctl";
my $CODEML_2OMEGAS = "two_omegas.ctl";

my %gene_names;
open(I, $cds_ali) || die "noooooo\n";
while(<I>) {
	if(/^>/) {
		chomp;
		s/^>//;
		$gene_names{$_} = 1;
	}
}
close(I);

##########################################################################################
# Optimise branch lengths:
chdir($output_path);
#`cp -p $INSTALL_DIR/$CODEML_BLEN .`;
#`perl -pi -e s/ALIGNMENT_ALIGNMENT/$cds_ali/` $CODEML_BLEN;
#`perl -pi -e s/TREE_TREE/$tree/` $CODEML_BLEN;
open(I, "$INSTALL_DIR/$CODEML_BLEN") || die "Norlrllalal\n";
open(O, ">$CODEML_BLEN") || die "anñdjñaj\n";
while(<I>) {
	if(/ALIGNMENT_ALIGNMENT/) {
		s/ALIGNMENT_ALIGNMENT/$cds_ali/;
	} elsif(/TREE_TREE/) {
		s/TREE_TREE/$tree/;
	}
	print O $_;
}
close(O);
close(I);

system("codeml $CODEML_BLEN");

my $tree_blen_optim = "";
open(I, "blen_optimization_results.txt") || die "Norrrrrlll\n";
open(O, ">tree.optim_blen.tree") || die "ñkajñaja\n";
my($lnl_h0, $lnl_h1);
while(<I>) {
	if(/^tree length = /) {
		$_ = <I>;
		$_ = <I>;
		$_ = <I>;
		$tree_blen_optim = <I>;
		chomp($tree_blen_optim);
		print O $tree_blen_optim, "\n";
	} elsif(/^lnL/) {
		chomp;
		s/ +/ /g;
		$lnl_h0 = (split)[4]; 
	}
}
close(O);
close(I);
##########################################################################################


##########################################################################################
# Run the two omega models for each gene:
print '-'x80,"\n";
print '-'x80,"\n";
open(RR, ">RESULTS.tsv") || die "ñaá´á+anañkaj\n";
foreach my $gene ( keys %gene_names ) {
	print STDERR "Working with gene: $gene\n";
	$lnl_h1 = "";
	my $OUT_PREFIX = "$gene"."_2w";
	open(I, "$INSTALL_DIR/$CODEML_2OMEGAS") || die "Norlrllalal\n";
	open(O, ">$CODEML_2OMEGAS") || die "anñdjñaj\n";
	while(<I>) {
		if(/ALIGNMENT_ALIGNMENT/) {
			s/ALIGNMENT_ALIGNMENT/$cds_ali/;
		} elsif(/TREE_TREE/) {
			s/TREE_TREE/$OUT_PREFIX.2w.tree/;
		} elsif(/TWO_OMEGAS/) {
			# change output name
			s/TWO_OMEGAS/$OUT_PREFIX/;
		}
		print O $_;
	}
	close(O);
	close(I);
	
	my $tmp_tree = $tree_blen_optim;
	$tmp_tree =~ s/$gene/$gene #1/;
	open(O, ">$OUT_PREFIX.2w.tree") || die "Ñañañjoiiioooo\n";
	print O $tmp_tree,"\n";
	close(O);
	
	system("codeml $CODEML_2OMEGAS");
	
	# Read resulting omegas:
	open(I, $OUT_PREFIX) || die "añdjañjaaj\n";
	my($w1,$w2,$tip_length) = ("","","");
	while(<I>) {
		if(/^w \(dN/) {
			chomp;
			s/ +/ /g;
			($w1,$w2) = (split)[4,5];
		} elsif(/^tree length =/) { # get TIP length
			$_ = <I>;
			$_ = <I>;
			$_ = <I>;
			$tmp_tree = <I>;
			chomp($tmp_tree);
			# e.g. (contarinia_XP_031616982.1: 0.820698, (((bcop_grc_g12323.t1: 0.502882, (bimp_core_g1016.t1: 0.234230, bcop_core_g3805.t1: 0.148999): 0.258871)>
			$tmp_tree =~ s/[,\(\) ]/ /g;
			$tmp_tree =~ s/ +/ /g;
			$tmp_tree =~ s/: /:/g;
			print STDERR "TMP_TREE:$tmp_tree\n";
			my @tmp = split(/ /,$tmp_tree);
			foreach my $t ( @tmp ) {
				print STDERR "GENE in $t\n";
				if($t =~ /$gene/) {
					$tip_length = (split(/:/,$t))[-1];
					print "TIP_LENGTH for gene=$gene: tip_length=$tip_length (w1=$w1, w2=$w2)\n";
				}	
			}
		} elsif(/^lnL/) {
			chomp;
			s/ +/ /g;
			$lnl_h1 = (split)[4]; 
		}
	}
	close(I);	
	print RR "$gene\t$w1\t$w2\t$tip_length\t$lnl_h0\t$lnl_h1\t$cds_ali\t$tree\n";
}
close(RR);

##########################################################################################



__END__




