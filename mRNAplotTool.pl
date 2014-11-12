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
my %dist;

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
                if ($a[1] == 0) 
                {
                        my $pos = $a[3] + 50;
                        my $pos_norm = int(($pos / $seq_len{$a[2]}) * 100 );
                        if ( defined $dist{$pos_norm} ) {
                                $dist{$pos_norm}++;
                        } else {
                                $dist{$pos_norm} = 1;
                        }
                }
        }
}
close(FH);

foreach my $pos (sort {$a<=>$b} keys %dist)
{
        print $pos."\t".$dist{$pos}."\n";
}




