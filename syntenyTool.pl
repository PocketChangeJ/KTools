#!/usr/bin/perl

=head
 -- generate synteny table using MCScanX
 -- compare genome using lastz
 -- compare genome using .......
=cut

use strict;
use warnings;
use IO::File;
use FindBin;
use Getopt::Std;

my $debug = 1;
my $version = 0.1;

my %options;
getopts('a:b:c:d:e:f:g:i:j:k:l:m:n:o:p:q:r:s:t:u:v:w:x:y:z:h', \%options);

unless (defined $options{'t'} ) { usage($version); }

if	($options{'t'} eq 'mcscanx_blast')	{ mcscanx_blast(\%options, \@ARGV); }
elsif	($options{'t'} eq 'mcscanx_gff')	{ mcscanx_gff(\%options, \@ARGV); }
elsif	($options{'t'} eq 'mcscanx')		{ mcscanx(\%options, \@ARGV); }
else	{ usage($version); }

=head2
 mcscanx_blast: prepare blast result for mcscanx
=cut
sub mcscanx_blast
{
	my ($options, $files) = @_;
	my $usage = qq'
USAGE $0 -t mcscan_blast input_A input_B name_A name_B

* input_A and input_B are pep sequences
* name_A and name_B are output name
* must provide 4 input files even if A and B are same

';
	print $usage and exit unless (scalar(@$files) == 4);
	print "[ERR]input file not exist\n" unless (-s $$files[0] && -s $$files[1]);
	my ($input_A, $input_B, $name_A, $name_B) = @$files;
	my $output;
	if (($name_A eq $name_B) && ($input_A eq $input_B)) {
		$output = $name_A.".blast";
	} else {
		$output = $name_A."_".$name_B.".blast";
	}
	run_cmd("formatdb -i $input_B -p T");
	run_cmd("blastall -i $input_A -d $input_B -p blastp -e 1e-10 -b 5 -v 5 -a 24 -m 8 -o $output");
	exit;
}

=head2
 mcscanx_gff: convert annotation bed to gff
=cut
sub mcscanx_gff
{
	my ($options, $files) = @_;
	my $usage = qq'
USAGE $0 -t mcscan_gff gene_position.bed name

* gene_position.bed: used for RNASeq analysis
* name is prefix of output file

';
	print $usage and exit unless (scalar(@$files) == 2);
	print "[ERR]input file not exist\n" unless -s $$files[0];
	my ($bed_file, $name) = @$files;

	my ($name_chr, $name_gff) = ($name.".chr", $name.".gff");
	print "[ERR]output file exist" if (-s $name_chr || -s $name_gff);

	# convert bed to gff
	# input : Gm01    27642   27977   Glyma01g00210.1 234     -
	# output: Gm01    Glyma01g00210.1 27643   27977
	my $out1 = IO::File->new(">".$name_gff) || die $!;
	my $out2 = IO::File->new(">".$name_chr) || die $!;

	my %change_chr; # key: old_chr, value: new_chr
	my $chr_order = 0;

	my $fin = IO::File->new($bed_file) || die $!;
	while(<$fin>)
	{
		chomp;
		next if $_ =~ m/^#/;
		my @a = split(/\t/, $_);
		my $chr;
		if (defined $change_chr{$a[0]}) {
			$chr = $change_chr{$a[0]};
		} else {
			$chr_order++;
			$chr = $name.$chr_order;
			$change_chr{$a[0]} = $chr;
			print $out2 "$a[0]\t$chr\n"
		}
		print $out1 $chr,"\t",$a[3],"\t",$a[1]+1,"\t",$a[2],"\n";
	}
	$fin->close;
	$out1->close;
	$out2->close;	
}

=head2
 mcscanx: run mcscanx
=cut
sub mcscanx
{
	my ($options, $files) = @_;
	my $usage = qq'
USAGE: $0 name_A name_B

* must provide 4 input files even if A and B are same
';

	print $usage and exit unless (scalar(@$files) == 2);
	my ($name_A, $name_B) = @$files;

	my $input_prefix;
	if ($name_A eq $name_B) {
		$input_prefix = $name_A;
	} else {
		$input_prefix = $name_A."_".$name_B;
	}

	# check input files: .blast .gff
	print "[ERR]input files\n" and exit unless (-s $input_prefix.".blast");

	unless (-s $input_prefix.".gff") {
		my ($gff_A, $gff_B) = ($name_A.".gff", $name_B.".gff");
		print "[ERR]input files\n" and exit unless (-s $gff_A && -s $gff_B);
		system("cat $gff_A $gff_B > $input_prefix.gff")
	}

	# check program
	my $mcscanx_bin = $FindBin::RealBin."/bin/MCScanX/MCScanX";
	die "[ERR]no mcscanx bin\n" unless -s $mcscanx_bin;

	# run mcscanx
	run_cmd("$mcscanx_bin $input_prefix");
}

=head2
 run_cmd: run command
=cut
sub run_cmd
{
	my $cmd = shift;
	print $cmd."\n" and return(1) if $debug;
	system($cmd) && die "[ERR]cmd: $cmd\n";
}

=head2
 usage: print usage information and pipeline
=cut
sub usage
{
	my $version = shift;
	my $usage = qq'
USAGE: $0 -t ToolOption 
	
	mcscanx_blast	blast protein sequence
	mcscanx_gff	generate gff
	mcscanx

example of pipeline:
	\$ perl syntenyTool.pl -t mcscanx_blast at_rep_pep CM_protein_v3.5_rep_pep.fasta
	\$ perl syntenyTool.pl -t mcscanx_blast sl_rep_pep CM_protein_v3.5_rep_pep.fasta

	\$ formatdb -i CM_protein_v3.5_rep_pep.fasta -p T
	\$ blastall -i sl_rep_pep -d CM_protein_v3.5_rep_pep.fasta -p blastp -e 1e-10 -b 5 -v 5 -a 24 -m 8 -o sl_cm.blast
	\$ blastall -i at_rep_pep -d CM_protein_v3.5_rep_pep.fasta -p blastp -e 1e-10 -b 5 -v 5 -a 24 -m 8 -o at_cm.blast

	\$ perl syntenyTool.pl -t mcscanx_gff arabidopsis_rep_gene.bed at
	\$ perl syntenyTool.pl -t mcscanx_gff melon_rep_gene.bed cm
	\$ perl syntenyTool.pl -t mcscanx_gff tomato_rep_ITAG2.3.bed sl

	\$ perl syntenyTool.pl -t mcscanx at cm
	\$ perl syntenyTool.pl -t mcscanx sl cm

	# generate plot for gene family : 
	/home/kentnf/pipeline/iTAK/synteny/plant_synteny.pl -i vv_bHLH \
		-a vv_gene_position -b vv.collinearity -c vv_chrSize \
		-x at_gene_position -y at_vv.collinearity -z at_chrSize
	will genrate two picture for result

	# how to combine two result into one ?
 
';

	print $usage;
	exit;
}





