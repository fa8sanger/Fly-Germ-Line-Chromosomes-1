#!/usr/bin/perl -w

use strict;

# One isoform per gene

my $cogelo = 0;
while(<>) {
	if(/^>/) {
		chomp;
		my $iso = (split(/\./,$_))[1];
		if($iso eq "t1") {
			$cogelo = 1;
		} else {
			$cogelo = 0;
		}	
		print "$_\n" if($cogelo == 1);
	} else {
		print if($cogelo == 1);
	}
}



