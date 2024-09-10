#!/usr/bin/perl -w

use strict;

my $cogelo = 1;
my $grc = 1;
open(GRC, ">grc.fa");
open(CORE,">core.fa");
while(<>) {
	if(/^>/) {
		if(/SUPER/) {
			$cogelo = 1;
		} else {
			$cogelo = 0;
		}
		if(/GRC/) {
			$grc = 1;
		} else {
			$grc = 0;
		}
		if($cogelo == 1 && $grc == 1) {
			print GRC $_;
		}
		if($cogelo == 1 && $grc == 0) {
			print CORE $_;
		}
	} else {
		if($cogelo == 1) {
			if($grc == 1) {
				print GRC $_;
			} else {
				print CORE $_;
			}
		}
	}
}
close(GRC);
close(CORE);