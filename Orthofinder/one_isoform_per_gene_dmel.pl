#!/usr/bin/perl -w

use strict;

# One isoform per gene

my $cogelo = 0;
my %ya;
while(<>) {
	if(/^>/) {
		chomp;
		my @tmp = split;
		foreach my $t ( @tmp) {
			if($t =~ /^parent/ ) {
				my $gene = (split(/[=,]/,$t))[1];
				if(!exists($ya{$gene})) {
					$cogelo = 1;
					$ya{$gene} = 1;
				} else {
					$cogelo = 0;
				}
			}
		}
		print "$_\n" if($cogelo == 1);
	} else {
		print if($cogelo == 1);
	}
}



