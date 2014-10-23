#!/usr/bin/perl

=head1
 mRNAtool.pl -- tools for mRNA analysis
=cut
use strict;
use warnings;
use IO::File;
use FindBin;
use Statistics::Basic qw(:all nofill);
use Getopt::Std;

my $version = 0.1;
my $debug = 0;

my %options;
getopts('t:s:l:e:d:p:n:f:u:h', \%options);

unless (defined $options{'t'} ) { usage($version); }

# checking parameters
my %ss = ('SE'=>1,'SS'=>1,'PE'=>1,'PS'=>1);
$options{'s'} = "SS" unless defined $options{'s'};
die "[ERR]sequencing-method: $options{'s'}\n" unless defined $ss{$options{'s'}};

if (defined $options{'l'} && $options{'l'} ne "fr-firststrand" && $options{'l'} ne "fr-secondstrand" && $options{'l'} ne "fr-unstranded" ) {
	die "[ERR]library-type: $options{'l'}\n";
} else {
	$options{'l'} = "fr-firststrand";
}

if	($options{'t'} eq 'align')	{ rnaseq_align(\%options, \@ARGV); }	# parse multi dataset
elsif	($options{'t'} eq 'tport')	{ rnaseq_tport(\%options, \@ARGV); }	# parse multi dataset
elsif	($options{'t'} eq 'count')	{ rnaseq_count(\%options, \@ARGV); }	# parse multi dataset
elsif	($options{'t'} eq 'norm')	{ rnaseq_norm(\%options, \@ARGV);  }	# parse single dataset
elsif	($options{'t'} eq 'corre')	{ rnaseq_corre(\%options, \@ARGV); }	# parse single dataset
elsif	($options{'t'} eq 'unmap')	{ rnaseq_unmap(\%options, \@ARGV); }	# parse single dataset
elsif	($options{'t'} eq 'pipeline')	{ pipeline(); }				# print pipelines
else	{ usage($version); } 

#################################################################
# kentnf: subroutine for each function				#
#################################################################
=head2
 rnaseq_align: align rnaseq_read to reference genome, generate report
=cut
sub rnaseq_align
{
	my ($options, $files) = @_;

	my $usage = qq'
USAGE: $0 [options] input1.fq input2.fq ...... | input1_r1_fq,input1_r2.fq input2_r1_fq,input2_r2.fq ...
	
 -s  [String]   Sequencing method for read files: (required)
        PE (paired-end);
        SE (single-end);
        PS (paired strand-specific);
        SS (single strand-specific);
 -d  [String]   Indexed database of genome sequences (required)
 -l  [String]   Library type, fr-unstranded, fr-firststrand, fr-secondstrand (required)
 -p  [Integer]  number of CPUs used for megablast clustering (default = 1)
 -r  [Integer]  mate-inner-dist (default = 100)
 -v  [Integer]  mate-std-dev (default = 20)
 -n  [Integer]  read mismatches (default = 0)
 -g  [Integer]  max-multihits (default : 20)
 -f  [String]	gene feature file (optional)

* input file could be fq, fa, or gzip file

';
	# check & correct parameters
	die "[ERR]no database\n$usage" unless defined $$options{'d'};
	my @index_files = ( 
		$$options{'d'}.".1.bt2", $$options{'d'}.".2.bt2", 
		$$options{'d'}.".3.bt2", $$options{'d'}.".4.bt2", 
		$$options{'d'}.".rev.1.bt2", $$options{'d'}.".rev.2.bt2" 
	);
	my $cpu = 20; $cpu = $$options{'p'} if defined $$options{'p'};
	my $mismatch = 1; $mismatch = $$options{'n'} if defined $$options{'n'};

	foreach my $index_file (@index_files) {
        	unless(-s $index_file) { die "[ERR] $index_file not exist\n"; }
	}
	
	# set the bin folder to ENV path
	$ENV{'PATH'} = ${FindBin::RealBin}."/bin".":".$ENV{'PATH'};

	my $tophat_bin;
	if (-s "${FindBin::RealBin}/bin/tophat-2.0.11.Linux_x86_64/tophat2" ) {
        	$tophat_bin = "${FindBin::RealBin}/bin/tophat-2.0.11.Linux_x86_64/tophat2";
	} else {
        	$tophat_bin = "tophat2";
	}
	print "Topaht2:\n".$tophat_bin."\n";	
	
	# check input files
	# key: output_suffix; value: input_file
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
	die "[ERR]no input files\n" if (scalar(@input_files) == 0);
	
	# main 
	# 1. create header of report file
	my $report_tophat = "";
	if ($$options{'s'} =~ m/^S/) {
		$report_tophat.="#sample\ttotal\tmapped\t%mapped\tmhit\t%mhit\n";
	} else {
		$report_tophat.="#sample\ttotal\tmapped\t%mapped\tmhit\t%mhit";
		$report_tophat.="\tLeftMap\t%LeftMap\tLeftMhit\t%LeftMhit\tRightMap\t%RightMap\tRightMhit\t%RightMhit\n";
	}

	# 2. tophat mapping
	foreach my $f (@input_files) 
	{
		my ($input, $output);
		if ($f =~ m/,/) {
			my @a = split(/,/, $f);
			$input = "$a[0] $a[1]";
			$output = $a[0]."_tophat";
		} else {
			$input = $f;
			$output = $f."_tophat";
		}

		# run tophat
		my $cmd_tophat = "$tophat_bin -o $output --library-type $$options{'l'} -p $cpu --segment-mismatches $mismatch --read-mismatches $mismatch --read-edit-dist $mismatch --read-gap-length $mismatch --max-multihits 20 --segment-length 25 $$options{'d'} $input";
		run_cmd($cmd_tophat);
		
		# parse report file
		my $report_tophat_file = $output."/align_summary.txt";
		my $rinfo = `cat $report_tophat_file`;
		chomp($rinfo); my @r = split(/\n/, $rinfo);	
		my ($total, $mapped, $mhit, $lefttotal, $leftmap, $leftmhit, $righttotal, $rightmap, $rightmhit);
		if ($f =~ m/,/)
		{
			if ( $r[1] =~ m/Input\s+:\s+(\d+)/ )	{ $lefttotal	= $1; }
			if ( $r[2] =~ m/Mapped\s+:\s+(\d+)/ )	{ $leftmap	= $1; }
			if ( $r[3] =~ m/of these:\s+(\d+)/ )	{ $leftmhit	= $1; }
			if ( $r[5] =~ m/Input\s+:\s+(\d+)/ )	{ $righttotal	= $1; }
			if ( $r[6] =~ m/Mapped\s+:\s+(\d+)/ )	{ $rightmap	= $1; }
			if ( $r[7] =~ m/of these:\s+(\d+)/ )	{ $rightmhit	= $1; }
			if ( $r[10] =~ m/Aligned pairs:\s+(\d+)/ ) { $mapped = $1; }
			if ( $r[11] =~ m/of these:\s+(\d+)/ )	{ $mhit 	= $1; }
			print "[WARN]Left ($lefttotal) is not same as right ($righttotal)\n" if $lefttotal ne $righttotal;
			$total = $lefttotal;
			$report_tophat.="$f\t$total\t$mapped\t".sprintf("%.2f", ($mapped/$total)*100);
			$report_tophat.="\t$mhit\t".sprintf("%.2f", ($mhit/$mapped)*100);
			$report_tophat.="\t$leftmap\t".sprintf("%.2f", ($leftmap/$lefttotal)*100);
			$report_tophat.="\t$leftmhit\t".sprintf("%.2f", ($leftmhit/$leftmap)*100);
			$report_tophat.="\t$rightmap\t".sprintf("%.2f", ($rightmap/$righttotal)*100);
			$report_tophat.="\t$rightmhit\t".sprintf("%.2f", ($rightmhit/$rightmap)*100)."\n";
		}
		else
		{
			if ( $r[1] =~ m/Input\s+:\s+(\d+)/ )	{ $total = $1; }
			if ( $r[2] =~ m/Mapped\s+:\s+(\d+)/ )	{ $mapped = $1;}
			if ( $r[3] =~ m/of these:\s+(\d+)/ )	{ $mhit = $1; }
			$report_tophat.="$f\t$total\t$mapped\t".sprintf("%.2f", ($mapped/$total)*100);
			$report_tophat.="\t$mhit\t".sprintf("%.2f", ($mhit/$mapped)*100)."\n";
		}
	}

	# output report information to file
	my $tophat_report_file = "report_tophat.txt";
	my $out = IO::File->new(">".$tophat_report_file) || die $!;
	print $out $report_tophat;
	$out->close;
}

=head2
 rnaseq_tport: parse tophat report file to generate result
=cut
sub rnaseq_tport
{
	my ($options, $files) = @_;
	my $usage = qq'
USAGE: $0 align_summary1.txt align_summary2.txt ... align_summary3.txt > report_tophat.txt

';
	print $usage and exit if scalar(@$files) == 0;

	# checking input files	
	foreach my $f (@$files) {
		print "[ERR]no file $f\n" and exit unless -s $f;
	}

	# parse report file
	my $report_tophat;
	my $report_line_num = 0;
	foreach my $f (@$files) {
		my $rinfo = `cat $f`;
		chomp($rinfo);
		my @r = split(/\n/, $rinfo);
		my ($total, $mapped, $mhit, $lefttotal, $leftmap, $leftmhit, $righttotal, $rightmap, $rightmhit);

		if (scalar @r == 14)
		{
			if ( $r[1] =~ m/Input\s+:\s+(\d+)/ )    { $lefttotal    = $1; }
			if ( $r[2] =~ m/Mapped\s+:\s+(\d+)/ )   { $leftmap      = $1; }
			if ( $r[3] =~ m/of these:\s+(\d+)/ )    { $leftmhit     = $1; }
			if ( $r[5] =~ m/Input\s+:\s+(\d+)/ )    { $righttotal   = $1; }
			if ( $r[6] =~ m/Mapped\s+:\s+(\d+)/ )   { $rightmap     = $1; }
			if ( $r[7] =~ m/of these:\s+(\d+)/ )    { $rightmhit    = $1; }
			if ( $r[10] =~ m/Aligned pairs:\s+(\d+)/ ) { $mapped = $1; }
			if ( $r[11] =~ m/of these:\s+(\d+)/ )   { $mhit         = $1; }
			print "[WARN]Left ($lefttotal) is not same as right ($righttotal)\n" if $lefttotal ne $righttotal;
			$total = $lefttotal;
			$report_tophat.="$f\t$total\t$mapped\t".sprintf("%.2f", ($mapped/$total)*100);
			$report_tophat.="\t$mhit\t".sprintf("%.2f", ($mhit/$mapped)*100);
			$report_tophat.="\t$leftmap\t".sprintf("%.2f", ($leftmap/$lefttotal)*100);
			$report_tophat.="\t$leftmhit\t".sprintf("%.2f", ($leftmhit/$leftmap)*100);
			$report_tophat.="\t$rightmap\t".sprintf("%.2f", ($rightmap/$righttotal)*100);
			$report_tophat.="\t$rightmhit\t".sprintf("%.2f", ($rightmhit/$rightmap)*100)."\n";
			$report_line_num = 14 if scalar(@r) > $report_line_num;
		}
		elsif (scalar @r == 5)
		{
			if ( $r[1] =~ m/Input\s+:\s+(\d+)/ )    { $total = $1; }
			if ( $r[2] =~ m/Mapped\s+:\s+(\d+)/ )   { $mapped = $1;}
			if ( $r[3] =~ m/of these:\s+(\d+)/ )    { $mhit = $1; }
			$report_tophat.="$f\t$total\t$mapped\t".sprintf("%.2f", ($mapped/$total)*100);
			$report_tophat.="\t$mhit\t".sprintf("%.2f", ($mhit/$mapped)*100)."\n";
			$report_line_num = 5 if scalar(@r) > $report_line_num;
                }
		else
		{
			print "[WARN]erro report format $f\n";
		}
	}

	my $report_head = '';
	if ($report_line_num == 14) {
		$report_head.="#sample\ttotal\tmapped\t%mapped\tmhit\t%mhit";
                $report_head.="\tLeftMap\t%LeftMap\tLeftMhit\t%LeftMhit\tRightMap\t%RightMap\tRightMhit\t%RightMhit\n";
	} else {
		$report_head.="#sample\ttotal\tmapped\t%mapped\tmhit\t%mhit\n";
	}

	$report_tophat = $report_head.$report_tophat;

	my $tophat_report_file = "report_tophat.txt";
        my $out = IO::File->new(">".$tophat_report_file) || die $!;
        print $out $report_tophat;
        $out->close;
}

=head2
 rnaseq_count: count and normalization for rnaseq analysis
 generate raw count and rpkm dataset for statistics analysis
=cut
sub rnaseq_count
{
	my ($options, $files) = @_;

	my $usage = qq'
USAGE: $0 [options] input1.bam input2.bam ......
        
 -s  [String]   Sequencing method for read files: (required)
        PE (paired-end);
        SE (single-end);
        PS (paired strand-specific);
        SS (single strand-specific);
 -f  [String]	gene feature file
 -l  [String]   Library type, fr-unstranded, fr-firststrand, fr-secondstrand (required)
 -o  [String]	prefix of output file (default : exp)

';
	# check parameters
	die "[ERR]no input file\n$usage" if (scalar(@$files) == 0);
	die "[ERR]no feature file\n$usage" unless defined $$options{'f'} && -s $$options{'f'};
	my $feature_bed = $$options{'f'};
	
	die "[ERR]sequencing method\n$usage" unless defined $$options{'s'};
	my $sequencing = $$options{'s'};
	die "[ERR]library type\n$usage" unless defined $$options{'l'};
	my $library = $$options{'l'}; 

	my $output_prefix = "exp";
	$output_prefix = $$options{'o'} if defined $$options{'o'};
	my $output_exp_raw  = $output_prefix."_raw_count.txt";
	my $output_exp_norm = $output_prefix."_rpkm.txt";
	die "[ERR]raw count exist\n" if -s $output_exp_raw;
	die "[ERR]RPKM exist\n" if -s $output_exp_norm;
	my $norm_method = "RPKM";

	# check alignment file
	foreach my $f ( @$files ) { 
		die "[ERR]no alignment file $f\n" unless -s $f;
		die "[ERR]no bam file $f\n" unless $f =~ m/\.bam$/;
	}

	# load feature into to hash
	# set bin size
	# change it base on genome size and density of feature
	# set it small will use more memory, set it big will take more time to locate read to feature
	my $bin_size = 1000;
	my %feature_struc = load_feature_bed($feature_bed, $bin_size);

	# count expression for alignment
	# key: file_name, feature_id, strand (sense or antisense)
	# value: raw count
	my %exp;

	my %flag_plus = ('147' => 1, '99' => 1, '145' => 1, '97' => 1,
			'355' => 1, '403'=> 1, '353' => 1, '401'=> 1);
	my %flag_minus= ('163' => 1, '83' => 1, '161' => 1, '81' => 1,
			'339' => 1, '419'=> 1, '337' => 1, '417'=> 1);

	foreach my $f ( @$files ) {
		my $sam = $f; $sam =~ s/\.bam/\.sam/;
		run_cmd("samtools view -h -o $sam $f") unless -s $sam;

		# uniq reads for Paired reads
		my %PE_uniq = ();

		# code for count expression
		my $fh = IO::File->new($sam) || die $!;
		while(<$fh>)
		{
			chomp;
			next if $_ =~ m/^@/;
			my @a = split(/\t/, $_);
			my $chr = $a[2];
			my $read_start = $a[3];				# get read start for single end
			my $mapped_len = parse_cigar($a[5]);
			my $read_end = $a[3] + $mapped_len - 1;		# get read end for single end
			my $flag = $a[1];

			# uniq the Paired reads, count only once for each fragment
			if ($sequencing =~ m/^P/) {
				my $key = '';
				if ($a[6] eq '=') # paired reads aligned to same reference
				{ 
					$a[6] = $a[2];
					if ($a[3] < $a[7]) { 
						$key = $a[0]."\t".$a[2]."\t".$a[3]."\t".$a[6]."\t".$a[7];
						$read_start = $a[3];
						$read_end   = $a[3] + $a[8] - 1;
					} else { 
						$key = $a[0]."\t".$a[6]."\t".$a[7]."\t".$a[2]."\t".$a[3];
						$read_start = $a[7];
						$read_end   = $a[7] + abs($a[8]) - 1;
					}
				} 
				else  # paired reads aligned to diff reference
				{
					$key = $a[0]."\t".$a[2]."\t".$a[3]."\t".$a[6]."\t".$a[7];
				}

				if ( defined $PE_uniq{$key} ) { next; }
				else { $PE_uniq{$key} = 1; }
			}

			my $bin_start_pos = int($read_start/$bin_size);
			my $bin_end_pos   = int($read_end/$bin_size);
			my %feature_pos = (
				$bin_start_pos-1 => 1, 
				$bin_start_pos =>1, 
				$bin_start_pos+1 => 1,
				$bin_end_pos-1 => 1,
				$bin_end_pos => 1,
				$bin_end_pos+1 => 1
			);

			my $feature = '';
			foreach my $fpos (sort {$a<=>$b} keys %feature_pos ) {
				$feature .= $feature_struc{$chr."\t".$fpos} if defined $feature_struc{$chr."\t".$fpos};
			}
			next unless $feature; # the read do not overlapped with feature

			my %uniq_feature;
			my $uniq_feature = '';
			chomp ($feature);
			my @fea = split(/\n/, $feature);
			foreach my $fea (@fea) {
				unless ( defined $uniq_feature{$fea} ) {
					$uniq_feature{$fea} = 1;
					$uniq_feature .= $fea."\n";
				}
			}

			# main: count reads number with strand-specific method
			# just count reads at this step
			my $fstat = "";
			if ($sequencing =~ m/^S/) {
				$fstat = "P" if ($flag eq '0'  || $flag eq '256');
				$fstat = "M" if ($flag eq '16' || $flag eq '272');
			} else {
				$fstat = "P" if defined $flag_plus{$flag};
				$fstat = "M" if defined $flag_minus{$flag};
			}
			next unless $fstat; # skip if do n

			chomp($uniq_feature); @fea = split(/\n/, $uniq_feature);
			foreach my $fea (@fea) {
				my @b = split(/\t/, $fea);
				if (($read_start >= $b[1] && $read_start <= $b[2]) || ($read_end >= $b[1] && $read_end <= $b[2]))
				{
					if ( defined $exp{$f}{$b[3]}{$fstat} ) { $exp{$f}{$b[3]}{$fstat}++; } 
					else { $exp{$f}{$b[3]}{$fstat}=1; }
				}
			}
		}
		$fh->close;
	}

	# output raw count
	my $sense_raw_table = "FeatureID";
	my $antisense_raw_table = "FeatureID";
	foreach my $f (@$files) {
		$sense_raw_table.="\t".$f;
		$antisense_raw_table.="\t".$f;
	}
	$sense_raw_table.="\n";
	$antisense_raw_table.="\n";

	my $fh = IO::File->new($feature_bed) || die $!;
	while(<$fh>)
	{
		chomp;
		next if $_ =~ m/^#/;
		my @a = split(/\t/, $_);
		my ($fid, $length, $strand) = ($a[3], $a[4], $a[5]);
		$sense_raw_table.=$fid;
		$antisense_raw_table.=$fid;

		foreach my $f (@$files) 
		{
			my ($ct_plus, $ct_minus) = (0,0);
			$ct_plus  = $exp{$f}{$fid}{'P'} if defined $exp{$f}{$fid}{'P'};
			$ct_minus = $exp{$f}{$fid}{'M'} if defined $exp{$f}{$fid}{'M'};
			my $ct = $ct_plus + $ct_minus;

			if ($sequencing =~ m/E$/) 
			{
				$sense_raw_table.= "\t".$ct;
			} 
			else # count raw count for strand specific read
			{
				if ( $library eq 'fr-firststrand' ) 
				{
					if ($strand eq "+") {
						$sense_raw_table.= "\t".$ct_minus;
						$antisense_raw_table.= "\t".$ct_plus;
					} else {
						$sense_raw_table.= "\t".$ct_plus;
						$antisense_raw_table.= "\t".$ct_minus;
					}
				} 
				else 
				{
					if ($strand eq "+") {
						$sense_raw_table.= "\t".$ct_plus;
						$antisense_raw_table.= "\t".$ct_minus;
					} else {
						$sense_raw_table.= "\t".$ct_minus;
						$antisense_raw_table.= "\t".$ct_plus;
					}
				}
			}
		}
		$sense_raw_table.="\n";
		$antisense_raw_table.="\n";
	}
	$fh->close;

	my $out1 = IO::File->new(">".$output_exp_raw) || die $!;
	print $out1 $sense_raw_table;
	$out1->close;

	if ($sequencing =~ m/S$/) {
		my $output_exp_raw_antisense = $output_prefix."_raw_count_anti.txt";
		my $out2 = IO::File->new(">".$output_exp_raw_antisense) || die $!;
		print $out2 $antisense_raw_table;
		$out2->close;
	}
}

=head2
 rnaseq_norm: normalization of expression file
=cut
sub rnaseq_norm
{
	my ($options, $file) = @_;
	my $usage = qq'
USAGE: $0 -t norm -f feature.bed -u report_tophat.txt project_raw_count.txt > project_rpkm.txt

';
	# check input file
	print $usage and exit unless defined $$file[0];
	my $exp_file = $$file[0];
	print "[ERR]no exp file $exp_file\n" and exit unless -s $exp_file;
	print "[ERR]no feature file\n$usage" and exit unless defined $$options{'f'};
	print "[ERR]no tophat report file\n$usage" and exit unless defined $$options{'u'};

	# check libsize and sample name
	my %sample_libsize;
	my $fh1 = IO::File->new($$options{'u'}) || die $!;
	<$fh1>;	# skip the title
	while(<$fh1>)
	{
		chomp;
		my @a = split(/\t/, $_);
		$sample_libsize{$a[0]} = $a[2];
	}
	$fh1->close;

	my $sname = `head -n 1 $exp_file`;
	chomp($sname);
	my @s = split(/\t/, $sname); shift @s;
	foreach my $s (@s) {
		print "[ERR]no libsize for sample: $s\n" and exit unless defined $sample_libsize{$s};
	}

	# load feature length from bed file
	my %feature_length;
	my $fh2 = IO::File->new($$options{'f'}) || die $!;
	while(<$fh2>)
	{
		chomp;
		next if $_ =~ m/^#/;
		# SL2.40ch00	16437	18189	Solyc00g005000.2.1	1693	+
		my @a = split(/\t/, $_);
		$feature_length{$a[3]} = $a[4];
	}
	$fh2->close;

	# normlization using RPKM/FPKM method
	my $fh3 = IO::File->new($exp_file) || die $!;
	my $title = <$fh3>; chomp($title);
	print $title."\n";
	my @title = split(/\t/, $title);
	while(<$fh3>)
	{
		chomp;
		next if $_ =~ m/^#/;
		my @a = split(/\t/, $_);
		print "[ERR]no feature length $a[0]\n" and exit unless defined $feature_length{$a[0]};
		my $length = $feature_length{$a[0]};
		my $output = $a[0];

		for(my $i=1; $i<@a; $i++) 
		{
			print "[ERR]no libsize for sample $title[$i]\n" and exit unless defined $sample_libsize{$title[$i]};
			my $lib_size = $sample_libsize{$title[$i]};
			my $rpkm = ($a[$i] * 1000 * 1000000) / ($length * $lib_size);
			$rpkm = sprintf("%.2f", $rpkm);
			$output.="\t".$rpkm;
		}

		print $output."\n";
	}
}

=head2
 rnaseq_corre: correlation analysis for rnaseq dataset
=cut
sub rnaseq_corre
{
	my ($options, $file) = @_;

	my $usage = qq'
USAGE: $0 -t corre [options] input_exp.txt > output_corre.txt

';

	print "[ERR]no exp file\n$usage" and exit unless defined $$file[0];
	my $exp_rpkm_file = $$file[0];
	print "[ERR]no exp file $exp_rpkm_file\n" and exit unless -s $exp_rpkm_file;

	my %sample_exp;
	my $fh = IO::File->new($exp_rpkm_file) || die $exp_rpkm_file;
	# parse first title line
	my $title = <$fh>; chomp($title);
	my @title = split(/\t/, $title);
	while(<$fh>)
	{
        	chomp;
		my @a = split(/\t/, $_);
		for(my $i=1; $i<@a; $i++)
		{
			my $tt = $title[$i];
			if (defined $sample_exp{$tt})   { $sample_exp{$tt}.="\t".$a[$i]; }
			else                            { $sample_exp{$tt} = $a[$i]; }
        	}
	}
	$fh->close;

	#################################################################
	# produce correlation for every line data                       #
	#################################################################

	print "Sample";
	shift @title;
	foreach my $tt ( @title ) { print "\t".$tt; }
	print "\n";

	foreach my $tt_i ( @title )
	{
        	print $tt_i;
		my @exp_i = split(/\t/, $sample_exp{$tt_i});

		foreach my $tt_j ( @title )
		{
                	my @exp_j = split(/\t/, $sample_exp{$tt_j});
			my $correlation = corr([@exp_i], [@exp_j]);
			print "\t".$correlation;
		}
		print "\n";
	}
}

=head2
 load_feature_bed : load feature informat with bed format
 data structure: 
 key: chr position(kb) (chr:1)
 value: structure of feature in the key position
=cut
sub load_feature_bed
{
	my ($feature_bed, $bin_size) = @_;

	# set bin size
	# change it base on genome size and density of feature
	# set it small will use more memory, set it big will take more time to locate read to feature
	# my $bin_size = 1000;

	# put feature to feature_struc hash
	my %feature_struc;
	my ($ref, $start, $end, $strand, $feature_id, $length, $kpos);
	my $fh = IO::File->new($feature_bed) || die $!;
	while(<$fh>)
	{
		chomp;
		next if $_ =~ m/^#/;
		my @a = split(/\t/, $_);
		($ref, $start, $end, $feature_id, $length, $strand) = ($a[0], $a[1]+1, $a[2], $a[3], $a[4], $a[5]);

		for(my $p = $start; $p < $end; $p = $p + $bin_size)
		{
			$kpos = int( $p / $bin_size );
			if ( defined $feature_struc{$ref."\t".$kpos} ) {
				$feature_struc{$ref."\t".$kpos}.= "$_\n";
			} else {
				$feature_struc{$ref."\t".$kpos} = $_."\n";
			}
		}
	}
	$fh->close;

	return %feature_struc;
}

#################################################################
# toolbox for rnaseq analysis					#
#################################################################
=head2
 rnaseq_unmap: get unmapped read for checking 
=cut
sub rnaseq_unmap
{
	my ($options, $file) = @_;
	
	my $usage = qq'
USAGE: $0 -t unmap [options] bam,read ... | bam,read_1,read2 ... 
	
	-u	[1-3] level for paired end mapped
		1. only paired aligned
		2. paired align + single align
		3. paired align + single align(paired)
';

	my $u = 2;
	$u = $$options{'u'} if (defined $$options{'u'} && $$options{'u'} > 0 && $$options{'u'} < 4);

	# check input files
	my @new_file;
	foreach my $set (@$file) {
		my @a = split(/,/, $set);
		print "[ERR]no bam file: $set\n" and exit unless ( $a[0] =~ m/\.bam$/ );
		print "[ERR]no bam file: $set\n" and exit unless -s $a[0];
		print "[ERR]no read1 file: $set\n" and exit unless -s $a[1];
		push(@new_file, "$a[0],$a[1]");
	}

	# main: get unmapped reads
	foreach my $set (@$file) 
	{
		my @a = split(/,/, $set);
		my $bam = $a[0];
		my $sam = $bam; $sam =~ s/\.bam$/\.sam/;
		run_cmd("samtools view $bam > $sam");

		# load reads flag to hash
		my %read_flag;
		my $fh1 = IO::File->new($sam) || die $!;
		while(<$fh1>)
		{
			chomp;
			next if $_ =~ m/^@/;
			my @a = split(/\t/, $_);
			if (defined $read_flag{$a[0]}) {
				my @f = split(/\t/, $read_flag{$a[0]});
				my %flg;
				foreach my $f (@f) { $flg{$f} = 1; }
				unless ( defined $flg{$a[1]} ) {
					$read_flag{$a[0]}.="\t".$a[1];
				}
			} else { 
				$read_flag{$a[0]} = $a[1]; 
			}
		}
		$fh1->close;

		if (scalar(@a) == 2)
		{
			my $unmap_file = $a[1].".unmap";
			my $out = IO::File->new(">".$unmap_file) || die $!;
			my $fh2 = IO::File->new($a[1]) || die $!;
			while(<$fh2>)
			{
				my $id = $_; chomp($id); my $format;
				if 	($id =~ m/^>/) { $format = 'fasta'; $id =~ s/^>//; }
				elsif	($id =~ m/^@/) { $format = 'fastq'; $id =~ s/^@//; }
				else	{ die "[ERR]seq format $id\n"; }
				my $tid = $id; $tid =~ s/ .*//;
				my $seq = <$fh2>; chomp($seq);
				if ($format eq 'fastq') { <$fh2>; <$fh2>; }
			
				if (defined $read_flag{$tid}) {
					my @f = split(/\t/, $read_flag{$tid});
					foreach my $flg (@f) 
					{
						my $flg_stat = parse_flag($flg);
						if ($flg_stat =~ m/unmapped/) { print $out ">$id\n$seq\n"; }
					}
				} else {
					print $out ">$id\n$seq\n";
				}
			}
			$fh2->close;
			$out->close;
		}
		elsif (scalar(@a) == 3)
		{
			# parse read flag according parameter U
			my %select_read;
			
			foreach my $r (sort keys %read_flag) 
			{
				my @fn = split(/\t/, $read_flag{$r});
				
				my $paired = 0;	# left, right, pair
				my $single = 0;
				foreach my $flag_num (@fn) 
				{
					my @flag_stat = parse_flag($flag_num);
					my $flag_stat = join("#", @flag_stat);
					if ($flag_stat =~ m/paired/) {
						if ($flag_stat =~ m/unmapped/) {
							$single = 'left' if $flag_stat =~ m/left/;
							$single = 'right' if $flag_stat =~ m/right/;
						} else {
							$paired = 1;
						}
					}
				}	
			}
		}
		else
		{
			print "[WARN]set has error: $set\n";
		}	
	}
}

=head2
 parse_flag: parse sam flag
=cut
sub parse_flag
{
	my $flag_num = shift;
	my @flag_stat = '';
	if ($flag_num & 0x1)    { push(@flag_stat, "Paired"); }
	if ($flag_num & 0x2)    { push(@flag_stat, "GoodAligned"); } 	
	if ($flag_num & 0x4)    { push(@flag_stat, "UnmappedSeg"); }
	if ($flag_num & 0x8)    { push(@flag_stat, "Unmapped 2ndSeg"); }
	if ($flag_num & 0x10)   { push(@flag_stat, "1stSEQ=RC"); }
	if ($flag_num & 0x20)   { push(@flag_stat, "2ndSEQ=RC"); }
	if ($flag_num & 0x40)   { push(@flag_stat, "Left"); }
	if ($flag_num & 0x80)   { push(@flag_stat, "Right"); }
	if ($flag_num & 0x100)  { push(@flag_stat, "Mhit"); }
	if ($flag_num & 0x200)  { push(@flag_stat, "lowQC"); }
	if ($flag_num & 0x400)  { push(@flag_stat, "Dup"); }
	if ($flag_num & 0x800)  { push(@flag_stat, "supplementary alignment"); } 
	return @flag_stat;
}

=head2
 parse_cigar: parse cigar to get mapped length 
=cut
sub parse_cigar
{
	my $cigar = shift;

	my $str_len = length($cigar);

	my $num = "";; my $mapped_length = 0;

	for(my $i=0; $i<$str_len; $i++)
	{
		my $str = substr($cigar, $i, 1);

		if ($str =~ m/\d+/)
		{
			$num = $num.$str;
		}
		elsif ($str eq "M" || $str eq "N" || $str eq "I")
		{
			$mapped_length = $mapped_length + $num;
			$num = "";
			
		}
		elsif ($str eq "D") 
		{
			$num = "";
		}
	}

	return $mapped_length;
}

=head2
 run_cmd: run other program using tophat
=cut
sub run_cmd
{
	my $cmd = shift;
	print $cmd."\n" and return(1) if $debug;
	system($cmd) && die "[ERR]cmd: $cmd\n";
}

=head2
 usage: print usage info
=cut
sub usage
{
	my $version = shift;
	print qq'
USAGE: $0 -t [tool] [options] input file

	align	align read to reference using tophat
	tport	summary tophat report info to table
	count	count aligned reads for each feature
	norm	normalization of expression file
	corre	generate correlation tables
	unmap	extract unampped reads (tophat could do it)
	denovo	denovo assembly using Trinity
	feature	generate feature file from contigs

'; 

	exit;
}

=head2
 pipeline: print pipeline
=cut
sub pipeline
{
	print qq'
>> RNASeq reference based pipeline:

1. align RNASeq read to reference (Tophat2)
   \$mRNAtool.pl -t align -s PS -d /home/database/tomato_genome -l fr-firststrand -p 24 -n 2 inputA_1.fastq,inputA_2.fastq ... inputN_1.fastq,inputN_2.fastq

2. count aligned reads for each feature
   \$mRNAtool.pl -t count -s PS -l fr-firststrand -f feature.bed -o project_exp inputA.bam ... inputN.bam

3. normalization
   \$mRNAtool.pl -t norm  -f feature.bed project_raw_count.txt > project_rpkm.txt

4. correlation
   \$mRNAtool.pl -t corre project_rpkm.txt > project_corr.txt

>> RNASeq denovo pipeline

1. assembly RNASeq read to contigs (Trinity)

2. align RNASeq read to reference (bowtie2/bwa)

3. count aligned reads for each contigs

4. 5. normalization and correlation like reference based 3. 4.

>> RNASeq debug utils

A. get unmapped reads


';
	exit;
}

