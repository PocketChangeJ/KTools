#!/usr/bin/perl

=head

=cut

use strict;
use warnings;
use IO::File;

my $usage = qq'
USAGE: perl $0 input_sam (RNASeq align to cDNA) > output

* align RNA-seq read to cDNA sequence to generate sam file
* run this script and geneate plot from output using excel

';

my $sam_file = shift || die $usage;

# parse sam file
my %seq_len;
my %dist_sens;
my %dist_anti;
my $max_pos = 0;

open(FH, $sam_file) || die $!;
while(<FH>)
{
        chomp;
        my @a = split(/\t/, $_);
        if ($_ =~ m/^\@/)
        {
                # @SQ     SN:Csa1M000010.1        LN:327
                if ($_ =~ m/^\@SQ/) {
                        $a[1] =~ s/SN://;
                        $a[2] =~ s/LN://;
                        $seq_len{$a[1]} = $a[2];
                        #print "$a[1] $a[2]\n";
                }
        }
        else
        {
                # HWI-ST397:445:C56GJACXX:1:2316:19064:100583     0       Csa6M242200.1   757     255     101M    *       0       0       CGGATGCTGGGGAATAGTGTG

		if ($a[1] == 16 || $a[1] == 0 || $a[1] == 163 || $a[1] == 99)
		{
			my $pos = $a[3];
			die "[ERR] frg len $a[8]\n" if $a[8] < 0;
			$pos = $a[3] + length($a[9]) - 1 if $a[1] == 16;
			$pos = $a[3] + int($a[8]/2) if ($a[1] == 163 || $a[1] == 99);

			my $pos_norm = int(($pos / $seq_len{$a[2]}) * 100 );
			$max_pos = $pos_norm if $pos_norm > $max_pos;

			print $_."\n" if $pos_norm > 100;

			if ($a[1] == 16 || $a[1] == 163) {
	                        if ( defined $dist_anti{$pos_norm} ) {
	                                $dist_anti{$pos_norm}++;
	                        } else {
	                                $dist_anti{$pos_norm} = 1;
	                        }
	                } else {
				if ( defined $dist_sens{$pos_norm} ) {
					$dist_sens{$pos_norm}++;
				} else {
					$dist_sens{$pos_norm} = 1;
				}
			}
		}
		else
		{
			next;
		}
        }
}
close(FH);

for(my $i=0; $i<=$max_pos; $i++) {
	my ($pos, $sens, $anti) = ($i, 0, 0);
	$sens = $dist_sens{$pos} if defined $dist_sens{$pos};
	$anti = $dist_anti{$pos} if defined $dist_anti{$pos};
	print "$pos\t$sens\t$anti\n";
}




