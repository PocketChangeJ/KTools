#!/usr/bin/perl

use strict;
use warnings;

while(<>)
{
	chomp;
	if ($_ =~ m/^@/) {
                print $_."\n";
        }else{
                my @ta = split(/\t/, $_);
                $#ta > 10 or next;
		if (($ta[1] & 0x4) || ($ta[1] & 0x8)) {
			next;
		}

                print "$_\n";
        }
}

