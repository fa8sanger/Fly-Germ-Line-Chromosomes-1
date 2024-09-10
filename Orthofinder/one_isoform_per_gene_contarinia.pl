#!/usr/bin/perl -w

use strict;


# First load the gff:
my %genes;
my %cds_lengths;
my %cds2gene;
open(I, "/lustre/scratch126/tol/teams/jaron/users/fede/GCF_009176525.2/genomic.ii.gff") || die "nooalslala\n";
print STDERR "Reading...\n";
while(<I>) {
	chomp;
	next if(/^#/);
	my @tmp = split(/\t/,$_);
	if($tmp[2] eq "CDS") {
		my @jj = split(/;/,$tmp[8]);
		my $cds_id = $jj[0];
		$cds_id = (split(/-/,$cds_id))[1];
		$cds_lengths{$cds_id} += $tmp[4]-$tmp[3]+1;
		my $gene_id = $jj[2];
		$gene_id = (split(/,/,$gene_id))[0];
		$genes{$gene_id}->{$cds_id} = 1;
		$cds2gene{$cds_id} = $gene_id;
		#print STDERR "CDS_ID=$cds_id        GENE_ID=$gene_id     LEN=$cds_lengths{$cds_id} \n";
	}
}
close(I);

print STDERR scalar(keys(%genes)), " genes\n";
print STDERR scalar(keys(%cds_lengths)), " cds_lengths\n";
print STDERR scalar(keys(%cds2gene)), " cds2gene\n";

my %chosen_cds;
foreach my $gene ( keys %genes ) {
	foreach my $cds ( sort { $cds_lengths{$b} <=> $cds_lengths{$a} } keys %{$genes{$gene}} ) {
		$chosen_cds{$cds} = 1;
		last;
	}
}

#exit;
# One isoform per gene

my $cogelo = 0;
while(<>) {
	if(/^>/) {
		s/^>//;
		my $id = (split)[0];
		if(exists($chosen_cds{$id})) {
			$cogelo = 1;
		} else {
			$cogelo = 0;
		}
		print ">$_" if($cogelo == 1);
	} else {
		print if($cogelo == 1);
	}
}



__END__
NW_022198046.1  Gnomon  exon    3842723 3842832 .       +       .       ID=exon-XM_031767620.1-5;Parent=rna-XM_031767620.1;Dbxref=GeneID:116340886,Genbank:XM_031767620.1;gbkey=mRNA;gene=LOC116340886;product=E3 ubiquitin-prote>
NW_022198046.1  Gnomon  exon    3842918 3843009 .       +       .       ID=exon-XM_031767620.1-6;Parent=rna-XM_031767620.1;Dbxref=GeneID:116340886,Genbank:XM_031767620.1;gbkey=mRNA;gene=LOC116340886;product=E3 ubiquitin-prote>
NW_022198046.1  Gnomon  exon    3843081 3843306 .       +       .       ID=exon-XM_031767620.1-7;Parent=rna-XM_031767620.1;Dbxref=GeneID:116340886,Genbank:XM_031767620.1;gbkey=mRNA;gene=LOC116340886;product=E3 ubiquitin-prote>
NW_022198046.1  Gnomon  exon    3843392 3843455 .       +       .       ID=exon-XM_031767620.1-8;Parent=rna-XM_031767620.1;Dbxref=GeneID:116340886,Genbank:XM_031767620.1;gbkey=mRNA;gene=LOC116340886;product=E3 ubiquitin-prote>
NW_022198046.1  Gnomon  exon    3843576 3843738 .       +       .       ID=exon-XM_031767620.1-9;Parent=rna-XM_031767620.1;Dbxref=GeneID:116340886,Genbank:XM_031767620.1;gbkey=mRNA;gene=LOC116340886;product=E3 ubiquitin-prote>
NW_022198046.1  Gnomon  exon    3843857 3846779 .       +       .       ID=exon-XM_031767620.1-10;Parent=rna-XM_031767620.1;Dbxref=GeneID:116340886,Genbank:XM_031767620.1;gbkey=mRNA;gene=LOC116340886;product=E3 ubiquitin-prot>
NW_022198046.1  Gnomon  CDS     3842042 3842099 .       +       0       ID=cds-XP_031623480.1;Parent=rna-XM_031767620.1;Dbxref=GeneID:116340886,Genbank:XP_031623480.1;Name=XP_031623480.1;gbkey=CDS;gene=LOC116340886;product=E3>
NW_022198046.1  Gnomon  CDS     3842173 3842340 .       +       2       ID=cds-XP_031623480.1;Parent=rna-XM_031767620.1;Dbxref=GeneID:116340886,Genbank:XP_031623480.1;Name=XP_031623480.1;gbkey=CDS;gene=LOC116340886;product=E3>
NW_022198046.1  Gnomon  CDS     3842474 3842645 .       +       2       ID=cds-XP_031623480.1;Parent=rna-XM_031767620.1;Dbxref=GeneID:116340886,Genbank:XP_031623480.1;Name=XP_031623480.1;gbkey=CDS;gene=LOC116340886;product=E3>
