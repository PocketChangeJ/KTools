#!/usr/bin/perl

=head1
 samtools -- connect command of samtools
=cut

use strict;
use warnings;
use IO::File;
use Getopt::Std;


my $debug = 1;

while(<>)
{
	chomp;
	mergeBAM($_);
}

sub mergeBAM
{
	my $list = shift;
	
	my ($plus_bam, $minus_bam, $out_bam) = ($list."_plus.bam", $list."_minus.bam", $list."_all.bam");
	run_cmd("samtools merge $out_bam $plus_bam $minus_bam");
	sortBAM($out_bam);
}

=head2
 convertBAM2SAM
=cut
sub convertBAM2SAM
{
	my ($inFiles, $param) = @_;
	
	my @files = split(/,/ , $inFiles); 

}


=head2
 convertSAM2BAN
=cut
sub convertSAM2BAM
{
	my ($inFiles, $param) = @_;
	
	my @files = split(/,/, $inFiles);


}

=head2
 sortBAM
=cut
sub sortBAM
{
	my $input_bam = shift;
	
	die "[ERR]input file name\n" unless $input_bam =~ m/\.bam$/;
	my $sort_prefix = $input_bam; $sort_prefix =~ s/\.bam$/_sort/; 
	my $sort_bam    = $sort_prefix.".bam";
	run_cmd("samtools sort $input_bam $sort_prefix");
	#die "[ERR]no sorted bam\n" unless -s $sort_bam;
	run_cmd("mv $sort_bam $input_bam");
}

sub run_cmd
{
	my $cmd = shift;
	print $cmd."\n" and return(1) if $debug;
	system($cmd) && die "[ERR]cmd: $cmd\n";
}

