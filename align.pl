#!/usr/bin/perl

=head

 Align list of reads to reference 
 1. remove rRNA from read
 2. align small RNA reads to genome for next analysis
 3. align RNASeq reads to genome for Quality Checking

 Author: Yi Zheng
 update: 02/09/2014 -- change it to align 
 update: 04/29/2014 -- update USAGE information
 update: 03/08/2014 -- accept fastq file
 update: 07/08/2013

=cut

use strict;
use warnings;
use FindBin;
use IO::File;
use Getopt::Std;

my $version = 0.1;
my $debug = 0;
my %options;
getopts('a:b:c:d:e:f:g:i:j:k:l:m:n:o:p:q:r:s:t:u:v:w:x:y:z:h', \%options);
unless (defined $options{'t'} ) { $options{'t'} = 'align'; } #

if	($options{'t'} eq 'align')	{ align(\%options, \@ARGV); }
elsif	($options{'t'} eq 'aport')	{ align_aport(\%options, \@ARGV); }
else	{ usage($version); }

#################################################################
# kentnf: subroutine						#
#################################################################

=head2
 generate report for align
=cut
sub align_aport
{
	my ($options, $files) = @_;

	my $usage = qq'
USAGE: $0 bowtie_report > output_report
 
';
	print $usage and exit unless defined $$files[0];

	print "#sample\ttotal\tmap\tunmap\n";
	my ($start, $sample_name, $cmd, $total, $align, $failed);
	my (@a, @b, @c, @d);

	my $fh = IO::File->new($$files[0]) || die $!;
	while(<$fh>)
	{
        	chomp;
        	#if ($_ =~ m/^Locate/ || $_ =~ m/^Read/) {  $start = 0; next; }
        	if ( $_ =~ m/bowtie/ )  { $start = 1; }
        	else { $start = 0; next; }

        	if ($start == 1)
        	{
                	$cmd = $_;                   @a = split(/\s+/, $cmd);
			$sample_name = $a[scalar(@a)-2];
	                $total = <$fh>;  chomp($total); @b = split(/\s+/, $total);
        	        $align = <$fh>;  chomp($align); @c = split(/\s+/, $align); 
	                $failed = <$fh>; chomp($failed);@d = split(/\s+/, $failed);
        	        <$fh>;
                	print $sample_name."\t".$b[3]."\t".$c[8]."\t".$d[6]."\n";
        	}
	}
	$fh->close;
}

=head2
 align -- align read to reference
=cut
sub align
{
	my ($options, $files) = @_;

	my $usage = qq'
USAGE: $0 -t align [options] input

	-p bowtie or bwa (default: bwa)
	-r rRNA_database_index (default = /home/database/rRNA_silva111)
	-o output_prefix
	-a number of CPU (default = 24)

	-u generate unmapped (only for bowtie)
	-m generate mapped (only for bowtie)
	-n edit_distance for bwa / mismatch for bowtie (default = 3) 
	-k good alignments per read (default = 1)
	-b report all/10000 hits for bowtie bwa

  * the alignment type has two values: 1 or 2
    1 for remove rRNA from input reads, the unaligned reads will be kept
    2 for align reads to reference, the aligned reads will be kept

';
	print $usage and exit unless defined $$files[0];
	foreach my $f (@$files) {
		my @r = split(/,/, $f); 
		die "[ERR]input files: $f\n" if (scalar(@r) > 2 || scalar(@r) < 1);
		foreach my $r (@r) {
			die "[ERR]file not exist: $r\n" unless -s $r;
		}
	}

	my $aligner = 'bowtie';
	$aligner = 'bwa' if (defined $$options{'p'} && $$options{'p'} eq 'bwa');

	my $edit_distance = 3;
	$edit_distance = $$options{'n'} if (defined $$options{'n'} && $$options{'n'} >= 0);
	
	my $reference = '/home/database/rRNA_silva111';
	$reference = $$options{'r'} if defined $$options{'r'};
	die "[ERR]ref file $reference\n" unless -s $reference;
	
	my $cpu = 24;
	$cpu = $$options{'a'} if (defined $$options{'a'} && $$options{'a'} > 0);

	my $out_prefix = '';
	$out_prefix = $$options{'o'} if (defined $$options{'o'});

	my $tophit = 1;
	$tophit = $$options{'k'} if (defined $$options{'k'} && $$options{'k'} > 0);

	# re-generate parameters according to aligner
	if ($aligner eq 'bowtie') {
		$cpu = "-p $cpu"; $edit_distance = "-v $edit_distance";	
		$tophit = "-k $tophit";
		$tophit = "--all --best" if defined $$options{'b'};
	} else {
		$cpu = "-t $cpu"; $edit_distance = "-n $edit_distance -o $edit_distance -e $edit_distance";
		$tophit = "-n $tophit";
		$tophit = "-n 100000" if defined $$options{'b'};
	}

	check_index($reference, $aligner);

	# main run command
	foreach my $f (@$files) 
	{
		my @r = split(/,/, $f);

		# check and set format of input reads
		my $format;
		foreach my $r (@r) {
			if (defined $format) {
				die "[ERR]format error PE: $f\n" if ($format ne check_reads($r));
			} else {
				$format = check_reads($r);
			}
		}		

		if ($aligner eq 'bowtie') {
			$format = '-f' if $format eq 'fa';
			$format = '-q' if $format eq 'fq';
		}
		
		# get input_prefix;
		my $input_prefix = $r[0]; $input_prefix =~ s/\.gz$//;
		$input_prefix =~ s/\.(fa|fq|fasta|fastq)$//;

		my ($input_files, $sam, $unmap, $mapped);  $unmap = ''; $mapped = '';
		$sam = $input_prefix.$out_prefix.".sam";
		$unmap = "--un ".$input_prefix.$out_prefix."_unmap" if defined $$options{'u'};
		$mapped = "--al ".$input_prefix.$out_prefix."_mapped" if defined $$options{'m'};

		if ($aligner eq 'bowtie') {
			$input_files = $r[0];
			$input_files = "-1 $r[0] -2 $r[1]" if defined $r[1];
			my $align_cmd = "$aligner $edit_distance $tophit $cpu $format $unmap $mapped -S $reference $input_files /dev/null";
			run_cmd($align_cmd);
		}
		else
		{
			my $sai1 = $r[0].".sai";
			if (defined $r[1] && -s $r[1]) {
				my $sai2 = $r[1].".sai";
				run_cmd("$aligner aln $cpu $edit_distance -i 1 -l 15 -k 1 -f $sai1 $reference $r[0]");
				run_cmd("$aligner aln $cpu $edit_distance -i 1 -l 15 -k 1 -f $sai2 $reference $r[1]") if -s $r[1];
				run_cmd("bwa sampe $tophit -f $sam $reference $sai1 $sai2 $r[0] $r[1]");
			} else {
				run_cmd("$aligner aln $cpu $edit_distance -i 1 -l 15 -k 1 -f $sai1 $reference $r[0]");
				run_cmd("bwa samse $tophit -f $sam $reference $sai1 $r[0]");
			}
		}
	}
}

#################################################################
# kentnf: subroutine			 			#
#################################################################

=head2
 check_reads -- check read format
=cut
sub check_reads
{
	my $file = shift;
	die "[ERR]file not exist $file\n" unless -s $file;
	my $head = `less $file | head -n 1 `;
	my $format = 'fa';
	$format = 'fq' if $head =~ m/^@/;
	return $format;
}

=head2
 check_index -- check and built index for diff aligner
=cut
sub check_index
{
	my ($r, $aligner) = @_;

	if ($aligner eq 'bowtie') {
		my $build = 0;
		foreach my $f (($r.".1.ebwt", $r.".2.ebwt", $r.".3.ebwt", $r.".4.ebwt", $r.".rev.1.ebwt", $r.".rev.2.ebwt")) {
			warn "[WARN]ref index missing: $f\n" and $build = 1 unless -s $f;
		}
		run_cmd("bowtie-build $r $r") if $build;
	} elsif ($aligner eq 'bwa') {
		my $build = 0;
		foreach my $f ( ($r.".amb", $r.".ann", $r.".bwt", $r.".pac", $r.".sa", $r.".fai") ) {
                	warn "[WARN]ref index missing: $f\n" and $build = 1 unless -s $f;
        	}
		run_cmd("bwa index $r") if $build;
		run_cmd("samtools faidx $r") if $build;
	}
}

sub run_cmd
{
	my $cmd = shift;
	print "[CMD]".$cmd."\n";
	system($cmd) && die "[ERR][CMD]$cmd\n";
}

