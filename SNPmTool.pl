#!/usr/bin/perl

=head1

 SNPmtool.pl -- call SNP for RNA-Seq 

 This pipeline are based on bwa alignment and pileup file generated by 
 samtools, C scripts are written by Linyong Mao for calling SNPs. Perl
 scripts are used for connecting all the steps to pipeline. 

 C script: Linyong Mao
 Filter SAM file: Honghe Sun
 perl script: Yi Zheng

 01/17/2014 filter SNP table 
 02/01/2014 parameter for remove multi-hit reads
 01/22/2014 fix bug for merge just one file, chrOrder 
 07/04/2013 init

=cut

use strict;
use warnings;
use FindBin;
use IO::File;
use Getopt::Std;

my $debug = 1;
my $version = 0.1;

my %options;
getopts('a:b:c:d:e:f:g:i:j:k:l:m:n:o:p:q:r:s:t:u:v:w:x:y:z:h', \%options);

unless (defined $options{'t'} ) { usage($version); }

if      ($options{'t'} eq 'identify')	{ snp_pipeline(\%options, \@ARGV); }    # 
elsif   ($options{'t'} eq 'filter1')	{ filter_RNASeq(@ARGV); }    		# 
elsif   ($options{'t'} eq 'filter2')	{ filter_indel(@ARGV); }    		#
elsif	($options{'t'} eq 'filterX')	{ snp_filter_table(\%options, @ARGV); }	# filter SNP combined table
elsif	($options{'t'} eq 'pipeline')	{ pipeline(); }
else	{ usage($version); }

#################################################################
# kentnf: subroutine						#
#################################################################

=head2
 snp_filter_table -- filter raw snp table
=cut
sub snp_filter_table
{
	my ($options, $files) = @_;
	
	my $input_file = $$files[0];
	die "[ERR]file not exist\n" unless -s $input_file;

	my $fh = IO::File->new($input_file) || die $!;
	my $title = <$fh>;
	print $title;
	while(<$fh>)
	{
		chomp;
		next if $_ =~ m/^#/;
		my @a = split(/\t/, $_);
		my $ref_id = shift @a;
		my $pos = shift @a;
		my $base = shift @a;
		foreach my $a (@a) {
			print $a."\n" and last if $a ne "N";
		}
	}
	$fh->close;
}

=head2
 snp_pipeline -- identify SNP from RNASeq dataset
=cut
sub snp_pipeline
{
	my ($options, $files) = @_;

	my $usage = qq'
USAGE: perl $0 -t identify -r reference [options]  input_RNASeq_list

	-p	thread
	-c	comparison file
	-e	run the script

* example of input RNASeq list
sampleName [tab] read_file1,read_file2 [tab] read_file3 ..... read_fileN

sampleName must diff with any read file name
read_file1,read_file2 are paired end reads
read_file3 are single end reads

* example of comparison_file
sampleNameA [tab] sampleNameB

';	
	print $usage and exit unless defined $$files[0];
	print "[ERR]no input RNAseq file\n" and exit unless -s $$files[0];
	my $input_list = $$files[0];

	print "[ERR]no reference\n" and exit unless defined $$options{'r'};
	my $genome = $$options{'r'};

	my $comparison_file = $$options{'c'} if defined $$options{'c'};
	
	#################################################################
	# init parameters and vars					#
	#################################################################
	my $ref_SNP_enable = 1;
	my $reSeqPrint = 1;

	my $debug = 1;
	my $cpu = 24;
	my $add_pileup = 1;

	#################################################################
	# check if the genome is indexed by bwa				#
	# check if the genome is indexed by samtools faidx		#
	# generate chrOrder file base one genome sequences		#
	#################################################################
	my $chrOrder_file = "chrOrder";
	check_genome($genome, $chrOrder_file);

	#################################################################
	# generate command to produce pileup files 			#
	#################################################################
	my (%cmd_pileup) = generate_pileup($input_list, $genome, $cpu, $debug);

	foreach my $cul (sort keys %cmd_pileup) { print $cmd_pileup{$cul}; }

	#################################################################
	# perform comparison analysis					#
	#################################################################
	my $script_bin = ${FindBin::RealBin}."/bin/SNPmao";

	if ( defined $comparison_file && -s $comparison_file)
	{
		#########################################################
		# load comparison file to hash				#
		# hash $cultivar					#
		# key: cultivar name; value: 1				#
		# hash $comparison					#
		# key: cultivarA \t cultivarB; value: 1			#
		# check cultivars exist in comparison file 		#
		#########################################################

		my ($cultivar, $comparison) = load_comparison($comparison_file);
		my $error = check_comparison($input_list, $cultivar); die if $error;

		foreach my $comparison (sort keys %$comparison)
		{
			my ($cultivarA, $cultivarB) = split(/\t/, $comparison);
			my ($pileupA, $pileupB) = ($cultivarA.".pileup", $cultivarB.".pileup");
			my $script = "$script_bin/combine2PileFiles";
			my $cmd_combine2PileFiles = "$script $pileupA $pileupB 0.9 0.8 $chrOrder_file 3";

			print $cmd_combine2PileFiles."\n";
			my $snp_file = $pileupA."_".$pileupB.".snp";
			my $snp_filter1 = $snp_file.".filter.txt";
			my $snp_filter2 = $snp_file.".filter.indel.txt";
			print "$0 -t filter1 $snp_file $snp_filter1\n";
			print "$0 -t filter2 $snp_filter1 $snp_filter2\n";

			unless ($debug) {	
				#system($cmd_combine2PileFiles) && die "Error in command: $cmd_combine2PileFiles\n";
				#(input_output)	
			}
		}
	}

	#################################################################
	# call SNPs between cultivar and reference			#
	#################################################################
	if ($ref_SNP_enable)
	{
		foreach my $cultivar (sort keys %cmd_pileup)
		{
			my $pileup = $cultivar.".pileup";
			my $script = "$script_bin/pileupFilter.AtoG";
			my $cmd_pileupFilter = "$script 0.9 0.8 3 $pileup";
			print $cmd_pileupFilter."\n";
			#system($cmd_pileupFilter) && die "Error in command: $cmd_pileupFilter\n";
		}
	}

	#################################################################
	# reSeqPrintSample virtual genome using SNP			#
	#################################################################
	if ($reSeqPrint)
	{
		foreach my $cultivar (sort keys %cmd_pileup)
		{
			my $pileup = $cultivar.".pileup";
			my $col = $cultivar.".1col";
			my $script = "$script_bin/reSeqPrintSample.indel.fast.strAssign.RNAseq.table";
			my $cmd_reSeqPrint = "$script $genome $col $pileup $cultivar 3 3 0.3";
			print $cmd_reSeqPrint."\n";
			#system($cmd_reSeqPrint) && die "Error in command: $cmd_reSeqPrint\n";
		}
	}
}

=head1 load_comparison
 load comparison and cultivar information to hash
=cut
sub load_comparison
{
	my $comparison_file = shift;
	my %comparison; my %cultivar;
	my $fh = IO::File->new($comparison_file) || die "Can not open comparison file $comparison_file $!\n";
	while(<$fh>)
	{
		chomp; 
		my @a = split(/\t/, $_);
		if ($_ =~ m/^#/) { next; }
		if (scalar(@a) != 2) { next; }
		$comparison{$_} = 1;
		$cultivar{$a[0]} = 1;
		$cultivar{$a[1]} = 1;
	}
	$fh->close;
	return (\%cultivar, \%comparison);
}

=head1 check_comparison
 check if the cultivars are consistent in comparison_file and input_list file
=cut
sub check_comparison
{
	my ($list_file, $cultivar) = @_;

	my %list_cultivar;
	my $fh = IO::File->new($list_file) || die "Can not open list file $list_file $!\n";
	while(<$fh>)
	{
		chomp;
		my @a = split(/\t/, $_);
		if ($_ =~ m/^#/) { next; }
		$list_cultivar{$a[0]} = 1;
	}
	$fh->close;

	my $error = 0;
	foreach my $cul (sort keys %$cultivar)
	{
		unless ( defined $list_cultivar{$cul} ) {   
			print "cultivar $cul in comparison file do not exist in list_file\n"; 
			$error = 1;
		}
	}
	return $error;
}

=head1 check_genome
 check if the genome is indexed by bwa
 check if the genome is indexed by samtools faidx
=cut
sub check_genome
{
	my ($genome, $chrOrder_file) = @_;
	my @files = ($genome.".amb", $genome.".ann", $genome.".bwt", $genome.".pac", $genome.".sa", $genome.".fai");
	foreach my $file (@files) 
	{
		unless (-s $file) { die "Error, please create genome index using bwa and samtools\n$file\n"; }
	}
	
	my $out = IO::File->new(">".$chrOrder_file) || die "Can not open chr Order file $chrOrder_file $!\n";
	my $fh = IO::File->new($genome) || die "Can not open genome file $genome $!\n";
	while(<$fh>)
	{
		chomp;
		if ($_ =~ m/^>/) 
		{
			my $chr = $_;
			$chr =~ s/^>//;
			$chr =~ s/ .*//ig;
			print $out $chr."\n";
		}
		else
		{
			next;
		}
	}
	$fh->close;
	$out->close;
}

=head1 generate_pileup
 generate pileup command
=cut
sub generate_pileup
{
	my ($list_file, $genome, $cpu, $debug) = @_;

	my $PE_filter = ${FindBin::RealBin}."/bin/SNPmao/filter_for_PEsnp.pl";
	my $SE_filter = ${FindBin::RealBin}."/bin/SNPmao/filter_for_SEsnp.pl";

	my %cmd_pileup;

	my $fh = IO::File->new($list_file) || die "Can not open input file: $list_file $!\n";
	while(<$fh>)
	{
		chomp;
		my @a = split(/\t/, $_);
		if ($_ =~ m/^#/) { next; }
		my $sample_name = $a[0];
		my $pileup_cmds = "";

		my @sort_bam;
		my @reads;

		# perform bwa alignment, generate sam, convert bam, sort, put it to hash
		my ($sai1, $sai2, $sai, $bam, $sort, $sort_bam);
		for(my $i=1; $i<@a; $i++)
		{
			@reads = split(/,/, $a[$i]);
			remove_dup($a[$i]);	# remove read duplication

			if ( scalar(@reads) == 2 )
			{

				my ($read1, $read2) = ($reads[0], $reads[1]);

				my $file_prefix1 = remove_file_suffix($read1);
				my $file_prefix2 = remove_file_suffix($read2);
				($sai1, $sai2, $bam, $sort, $sort_bam) = ($file_prefix1.".sai", $file_prefix2.".sai", $file_prefix1.".bam", $file_prefix1."_sort", $file_prefix1."_sort.bam");

				my $bwa_align_cmd1 = "bwa aln -t $cpu -n 0.02 -o 1 -e 2 -f $sai1 $genome $read1";
				my $bwa_align_cmd2 = "bwa aln -t $cpu -n 0.02 -o 1 -e 2 -f $sai2 $genome $read2";
				$pileup_cmds.=$bwa_align_cmd1."\n";
				$pileup_cmds.=$bwa_align_cmd2."\n";

				my $bwa_sam_cmd = "bwa sampe $genome $sai1 $sai2 $read1 $read2 | $PE_filter | samtools view -bS -o $bam -";
				$pileup_cmds.=$bwa_sam_cmd."\n";

				my $sort_cmd = "samtools sort $bam $sort";
				$pileup_cmds.=$sort_cmd."\n";
			
				push(@sort_bam, $sort_bam);
			}
			elsif ( scalar(@reads) == 1 )
			{
				my $read = $a[$i];
				my $file_prefix = remove_file_suffix($read);
				($sai, $bam, $sort, $sort_bam) = ($file_prefix.".sai", $file_prefix.".bam", $file_prefix."_sort", $file_prefix."_sort.bam");				

				my $bwa_align_cmd = "bwa aln -t $cpu -n 0.02 -o 1 -e 2 -f $sai $genome $read";
				$pileup_cmds.=$bwa_align_cmd."\n";
			
				my $bwa_sam_cmd = "bwa samse $genome $sai $read | $SE_filter | samtools view -bS -o $bam -";
				$pileup_cmds.=$bwa_sam_cmd."\n";
			
				my $sort_cmd = "samtools sort $bam $sort";
				$pileup_cmds.=$sort_cmd."\n";

				push(@sort_bam, $sort_bam);
			}
			else
			{
				print "Error in input sample files $!\n";
			}
		} 

		# merge all bam files
		my $all_bam = $sample_name.".bam";
		my $s_bam = join(" ", @sort_bam);
		my $sam_merge_cmd = "samtools merge -f $all_bam $s_bam";
		if (scalar(@sort_bam) == 1) { $sam_merge_cmd = "mv $s_bam $all_bam"; }
		$pileup_cmds.=$sam_merge_cmd."\n";

		# remove multi-hit reads

		# pileup all files
		my $all_pileup = $sample_name.".pileup";
		my $mpileup_cmd = "samtools mpileup -q 16 -Q 0 -d 10000 -f $genome $all_bam > $all_pileup";
		$pileup_cmds.=$mpileup_cmd."\n";

		$cmd_pileup{$sample_name} = $pileup_cmds;
	}
	$fh->close;

	return %cmd_pileup;
}

=head2
 remove_file_suffix -- remove fq fastq fa fasta gz
=cut
sub remove_file_suffix
{
	my $file_name = shift;
	$file_name =~ s/\.gz$//;
	$file_name =~ s/\.fastq//;
	$file_name =~ s/\.fq//;
	$file_name =~ s/\.fasta//;
	$file_name =~ s/\.fa//;
	return $file_name;
}

=head2
 filter_RNASeq : filter SNP result from RNASeq (input_file, output_file)
=cut 
sub filter_RNASeq
{
	my ($input_file, $output_file) = @_;

	my $usage = qq'
USAGE: $0 -t filter1 input_file output_file

* filter RNASeq SNP: 
  1) the start and end base is not count
  2) convert the pileup to simple format

';

	print $usage and exit unless defined $input_file;
	print $usage and exit unless defined $output_file;

	my $output_line = "type\tSNP\tchromosome\tposition\tRef base\ts1 coverage\ts1 base\ts2 coverage\ts2 base\n";

	my ($count1, $count2, $cov1, $cov2);
	my $fh = IO::File->new($input_file) || die $!;
	while(<$fh>) {
		chomp;
		my @a = split "\t";

		# remove SNP in read end
		$count1 = ($a[6] =~ tr/\^//);
		$count1 += ($a[6] =~ tr/\$//);
		$count2 = ($a[8] =~ tr/\^//);
 		$count2 += ($a[8] =~ tr/\$//);
		$cov1 = $a[5] - $count1;
		$cov2 = $a[7] - $count2;

		$a[1] =~ s/;//;
	
		$a[6] = repace_pileup($a[6]);
		$a[8] = repace_pileup($a[8]);
		$a[6] =~ s/\*/$a[4]/g;
        	$a[8] =~ s/\*/$a[4]/g;

        	if ($cov1 >= 4 && $cov2 >= 4) {
                	$output_line .= join("\t", @a)."\n";
        	}
	}

	if ( $input_file ne $output_file ) {
		my $out = IO::File->new(">$output_file") || die $!;
		print $out $output_line; 
		$out->close;
	} else {
		my $out = IO::File->new(">temp.RNASeq.snp.filter.txt") || die $!;
		print $out $output_line;
		$out->close;
		system("mv temp.RNASeq.snp.filter.txt $output_file");
	}
}
=head2
 repace_pileup: repace_pileup to ATCGN for mutation and -+ for indel
=cut
sub repace_pileup
{
	my $char = shift;
	$char =~ tr/agctn/AGCTN/;	# format the mismatch for forward and reverse strand
	$char =~ s/"//g;		# ?
	$char =~ s/\$//g;		# remove the end symbol (why do not remove corresponding base)
	$char =~ s/\^://g;		# ??? like remove start, but do not understand
	$char =~ s/\^F//g;		#
	$char =~ s/\^\d//g;		#
	$char =~ s/\^\)//g;		#
	$char =~ s/\.\+/\+/g;		# remove the match base before insertion and deletion 
	$char =~ s/\.\-/\-/g;		#
	$char =~ s/,\+/\+/g;		#
	$char =~ s/,\-/\-/g;		# 
	$char =~ s/\*/\-/g;		# ?????
	$char =~ s/,/\*/g;		# replace the comma (match reverse) to asterisk
	$char =~ s/\./\*/g;		# replace the dot (match forward) to asterisk
	$char =  " ".$char;		# for easy to import to excel
	#$char =~ s/\*/$a[4]/g;		# replace the asterisk to reference base
	return $char;
}

=head2
 filter_indel: filter indel result by indel Depth (input_file, $output_file)
=cut
sub filter_indel
{
	my ($input_file, $output_file) = @_;

	my $usage = qq'
USAGE: $0 -t filter2 input_file output_file

* filter indel lower than 90% depth

';	
	print $usage and exit unless defined $input_file;
	print $usage and exit unless defined $output_file;

	my %base_count;
        my $fh = IO::File->new($input_file) || die $!;
        my $output_line = <$fh>; # titile
        while(<$fh>)
        {
                chomp;
                my @a = split(/\t/, $_);
                $output_line.= $_."\n" and next if $a[0] eq 'M';
                my ($c1, $b1, $c2, $b2) = ($a[5], $a[6], $a[7], $a[8]);
                %base_count = ();
                $base_count{$b1} = $c1;
                $base_count{$b2} = $c2;
                my $label = 'removed';
                foreach my $base (sort keys %base_count) {
                        my $count = $base_count{$base};
                        my $count1 = $base =~ tr/+/+/;
                        my $count2 = $base =~ tr/-/-/;
                        if (($count1 / $count) >= 0.9 || ($count2 / $count) >= 0.9) {
                                $label = 'keep';
                        }
                        #print $_."\t$count\t$count1\t$count2\t$label\n";
                }
                $output_line.= $_."\t".$label."\n";
        }
        $fh->close;

	if ( $input_file ne $output_file ) {
                my $out = IO::File->new(">$output_file") || die $!;
                print $out $output_line;
                $out->close;
        } else {
                my $out = IO::File->new(">temp.RNASeq.snp.filter.txt") || die $!;
                print $out $output_line;
                $out->close;
                system("mv temp.RNASeq.snp.filter.txt $output_file");
        }
}

=head2
 remove duplication for read
=cut
sub remove_dup
{
	my $seq_files = shift;  # format, sample.fastq (single-end); sample_R1.fastq,sample_R2.fastq (paired-end)
	my $script = ${FindBin::RealBin}."/bin/SNPmao/removeRedundancy.pl";
	print "perl $script $seq_files\n";	
}

=head2
 usage: show to to use this pipeline
=cut
sub usage
{
	my $version = shift;
	my $usage = qq'
USAGE: $0 -t [tool] [options] input file

        identify  identify SNP from RNASeq
       	filter1   filter RNASeq SNP result
        filter2   filter RNASeq SNP indel

';
	print $usage;
	exit;
}

=head2
 pipeline: show how to use this pipeline
=cut
sub pipeline
{
	print qq'
===== A original SNP calling pipeline from Mao =====

1. remove redundancy reads using removeRedundancy.pl
   \$ perl removeRedundancy.pl input > output

2. align each sample to reference using bwa
   \$ bwa aln -t 24 -n 0.02 -o 1 -e 2 -f sample1.sai reference.fa sample1.fa
   \$ bwa samse reference.fa sample1.sai sample1.fa | filter_for_SEsnp.pl | samtools view -bS -o sample1.bam -
   \$ samtools sort sample1.bam sample1_sort

3. merge sorted bam files for each cultivar
   \$ samtools merge -f cultivar_A.bam sample1_sort.bam sample2_sort.bam ...... sampleN_sort.bam

4. generate pileup files for each cultivar
   \$ samtools mpileup -q 16 -Q 0 -d 10000 -f reference.fa cultivar_A.bam > cultivar_A.pileup

* choose one or more below step for next analysis.

5. generate virtual genome using SNP.
   \$ reSeqPrintSample.indel.fast.strAssign.RNAseq.table reference.fa cultivar_A.1col cultivar_A.pileup cultivar_A 3 3 0.3

6. call SNPs between cultivar and reference
   \$ pileupFilter.AtoG 0.9 0.8 3 cultivar_A.pileup

7. call SNPs between two cultivars
   \$ combine2PileFiles cultivar_A.pileup cultivar_B.pileup  0.9  0.8  ChrOrder  3

* the ChrOrder is the Chr ID, one ID per line.

===== B simple pipeline ===== 

* notice, this script only generate command to run

1. Prepare list file for pipeline

    * list file format
    cultivarA cultivarA_rep1.fa cultivarA_rep2.fa cultivarA_rep3.fa ... cultivarA_repN.fa
    cultivarB cultivarB_rep1.fa cultivarB_rep2.fa cultivarB_rep3.fa ... cultivarB_repN.fa

2. Prepare comparison file [option]

    * list file for comparison
    cultivarA cultivarB

3. Run Mao SNP pipeline

    \$ SNPmTool.pl list_file  comparison_file > run_cmd.sh

    Edit the run_cmd.sh file if required

    ./run_cmd.sh

';
	exit;
}


