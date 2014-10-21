#!/usr/bin/perl

package align;

use strict;
use warnings;
use Exporter;
use FindBin;
use IO::File;

our @ISA = qw (Exporter);
our @EXPORT = qw (levenshtein print_bin);

my $bin = ${FindBin::Bin}."../m";


sub print_bin
{
	print $bin."\n";
}

=head2

 multi_align_sRNA -- align sRNA to referece with multiply method (BWA)
 
=cut
=head
sub multi_align_sRNA
{
	my ($input_file, $reference, $output_file, $parameters, $mhit_num, $debug) = @_;
	my $align_program = ${FindBin::RealBin}."/../bin/bwa";

	# align input_file to reference using bwa

	my $sai = $temp_folder."/bwa.sai";
	my $log = $temp_folder."/bwa.log";
	my $bwa_mhit_param = "-n $mhit_num";
		if ($mhit_num > 1 ) { $bwa_mhit_param = "-n $mhit_num -s" }
		Util::process_cmd("$align_program index -p $reference -a bwtsw $reference 2> $log", $debug) unless -s "$reference.amb";
		Util::process_cmd("$align_program aln $parameters $reference $input_file 1> $sai 2>> $log", $debug);
		Util::process_cmd("$align_program samse $bwa_mhit_param $reference $sai $input_file 1> $output_file 2>> $log", $debug);
	}

	return 1;
}
=cut

=head2
 levenshtein -- compute levenshtein distance for two seq, by Shan Gao
 usage: $distance =  levenshtein($seq1, $seq2);
=cut
sub levenshtein
{
	my ($s1, $s2) = @_;
	my ($len1, $len2) = (length $s1, length $s2);
	return $len2 if ($len1 == 0);
	return $len1 if ($len2 == 0);

	my %mat;

	for (my $i = 0; $i <= $len1; ++$i)
    	{
		for (my $j = 0; $j <= $len2; ++$j)
		{
            		$mat{$i}{$j} = 0;
            		$mat{0}{$j} = $j;
        	}
        	$mat{$i}{0} = $i;
    	}

	my @ar1 = split(//, $s1);
	my @ar2 = split(//, $s2);

	for (my $i = 1; $i <= $len1; ++$i)
    	{   
        	for (my $j = 1; $j <= $len2; ++$j)
        	{   
            		my $cost = ($ar1[$i-1] eq $ar2[$j-1]) ? 0 : 1;
            		$mat{$i}{$j} = min($mat{$i-1}{$j} + 1,
					$mat{$i}{$j-1} + 1,
					$mat{$i-1}{$j-1} + $cost);
        	}
    	}
	return $mat{$len1}{$len2};
}

sub min
{
	my @list = @_;
    	my $min = $list[0];

	foreach my $i (@list) {   
        	$min = $i if ($i < $min);
	}
	return $min;
}


1;


