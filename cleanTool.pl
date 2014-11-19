#!/usr/bin/perl

=head1
 cleanTool.pl -- tools for clean NGS reads
=cut
use strict;
use warnings;
use IO::File;
use FindBin;
use Getopt::Std;

my $version = 0.1;
my $debug = 0;

my %options;
getopts('a:b:c:d:e:f:g:i:j:k:l:m:n:o:p:q:r:s:t:u:v:w:x:y:z:h', \%options);
unless (defined $options{'t'} ) { usage($version); }

# checking parameters

if	($options{'t'} eq 'trimo')	{ clean_trimo(\%options, \@ARGV); }	# parse multi dataset
elsif	($options{'t'} eq 'align')	{ clean_align(\%options, \@ARGV); }	# parse multi dataset
elsif   ($options{'t'} eq 'sRNA')	{ clean_sRNA(\%options, \@ARGV); }
else	{ usage($version); }
#################################################################
# kentnf: subroutine						#
#################################################################

=head2
 trimo -- clean reads using trimmomatic
=cut
sub clean_trimo
{
	my ($options, $files) = @_;
	my $usage = qq'
USAGE $0 -t trimo [options] sample1_R1 ... sampleN_R1 | sample1_R1,sample1_R2 ... sampleN_R1,sampleN_R2

'; 

	# checking trimmomatic and adapter sequence
	my ($trim_bin, $adp_SE, $adp_PE);
	$trim_bin = ${FindBin::RealBin}."/bin/trimmomatic0.32.jar";
	$adp_SE = ${FindBin::RealBin}."/bin/adapters/TruSeq3-SE.fa";
	$adp_PE = ${FindBin::RealBin}."/bin/adapters/TruSeq3-PE.fa";
	my @require_files = ($trim_bin, $adp_SE, $adp_PE);
	
	foreach my $f (@require_files) {
		print "[ERR]no file $f\n" and exit unless -s $f;
	}
	
	my $qual = 15; my $score = 7; my $seed_mismatch = 2; my $min_len = 40; my $thread = 24;
	$qual = $$options{'q'} if defined $$options{'q'}  && $$options{'q'} >= 10;
	$thread = $$options{'p'} if defined $$options{'p'} && $$options{'p'} > 0;
	
	# check input files
	my @input_files;
	foreach my $f ( @$files ) 
	{
		if ($f =~ m/,/) 
		{
			my @a = split(/,/, $f);
			print STDERR "[ERR]file $a[0] or $a[1] not exist\n" and next unless (-s $a[0] && -s $a[1]);
			push(@input_files, $f);
		}
		else
		{
			print STDERR "[ERR]file $f not exist\n" and next unless -s $f;
			push(@input_files, $f);
		}

	}
	die "$usage\n" if (scalar(@input_files) == 0);

	# trim adapter lowqual
	my ($ss, $clip, $inout, $outprefix, $cmd); 
	foreach my $f (@input_files) 
	{
		if ($f =~ m/,/) {
			my @a = split(/,/, $f);
			$outprefix = $a[0];
			$outprefix =~ s/\.gz//;
			$outprefix =~ s/\.fq//; $outprefix =~ s/\.fastq//;
			$outprefix =~ s/_R1//;  $outprefix =~ s/_1//;
			my ($p1, $p2, $s1, $s2) = ($outprefix."_paired_R1.fastq", $outprefix."_single_R1.fastq", $outprefix."_paired_R2.fastq", $outprefix."_single_R2.fastq");
			$inout = "$a[0] $a[1] $p1 $p2 $s1 $s2";
			$ss = 'PE';
			$clip = "ILLUMINACLIP:$adp_PE:$seed_mismatch:30:$score:1:true";
		} else {
			$outprefix = $f;
			$outprefix =~ s/\.gz//;
			$outprefix =~ s/\.fq//; $outprefix =~ s/\.fastq//;
			$inout = "$f $outprefix"."_clean.fastq";
			$ss = 'SE';
			$clip = "ILLUMINACLIP:$adp_SE:$seed_mismatch:30:$score";
		}
		$cmd = "java -jar $trim_bin $ss -threads $thread -phred33 $inout $clip LEADING:3 TRAILING:3 SLIDINGWINDOW:4:$qual MINLEN:$min_len";
		print $cmd."\n";
	}

}

=head2
 align -- clean reads through align reads
=cut
sub clean_align
{

}

=head2
 clean -- clean sRNA reads
=cut
sub clean_sRNA
{

}

=head2
 usage: print usage information
=cut
sub usage
{
	print qq'
USAGE: $0 -t [tool] [options] input file

	trimo	trim adapter, low quality, and short reads using trimmomatic.	

';
	exit;
}

=head2
 usage: print pipeline or command 
=cut
sub pipeline
{
	print qq'
>> pipeline for $0
>> 1. clean mRNA
      	1.1 clean adapter, low quality, short reads
	1.2 clean rRNA or other contanmination
	
>> 2. clean sRNA/degradome
	2.1 clean adapter, short reads,
	2.2 clean rRNA or other contanmination

>> 3. clean DNA
	3.1 clean adapter, low quality, short reads
';
	exit;
}
