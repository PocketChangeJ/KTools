#!/usr/bin/perl

=head1

 MaoSNP_pipeline.pl -- call SNP for RNA-Seq 

 This pipeline are based on bwa alignment and pileup file generated by 
 samtools, C scripts are written by Linyong Mao for calling SNPs. Perl
 scripts are used for connecting all the steps to pipeline. 

 C script: Linyong Mao
 Filter SAM file: Honghe Sun
 perl script: Yi Zheng

 
 02/01/2014 parameter for remove multi-hit reads
 01/22/2014 fix bug for merge just one file, chrOrder 
 07/04/2013 init

=cut

use strict;
use warnings;
use FindBin;
use IO::File;
use Getopt::Long;


sub snp_pipeline
{
	my ($options, $files) = @_;

	my $usage = qq'
USAGE: perl $0 -t SNP -r reference [options]  input_RNASeq_list

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
	# load comparison file to hash					#
	# hash $cultivar						#
	# key: cultivar name; value: 1					#
	# hash $comparison						#
	# key: cultivarA \t cultivarB; value: 1				#
	#								#
	# check cultivars exist in comparison file			#
	#################################################################
	my ($cultivar, $comparison) = load_comparison($comparison_file);
	my $error = check_comparison($input_list, $cultivar); die if $error;

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
	if (-s $comparison_file)
	{
		foreach my $comparison (sort keys %$comparison)
		{
			my ($cultivarA, $cultivarB) = split(/\t/, $comparison);
			my ($pileupA, $pileupB) = ($cultivarA.".pileup", $cultivarB.".pileup");
			my $script = ${FindBin::RealBin}."/bin/combine2PileFiles";
			my $cmd_combine2PileFiles = "$script $pileupA $pileupB 0.9 0.8 $chrOrder_file 3";
			print $cmd_combine2PileFiles."\n";
			#system($cmd_combine2PileFiles) && die "Error in command: $cmd_combine2PileFiles\n";
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
			my $script = ${FindBin::RealBin}."/bin/pileupFilter.AtoG";
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
			my $script = ${FindBin::RealBin}."/bin/reSeqPrintSample.indel.fast.strAssign.RNAseq.table";
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
		my ($bam, $sort, $sort_bam);
		for(my $i=1; $i<@a; $i++)
		{
			@reads = split(/,/, $a[$i]);

			if ( scalar(@reads) == 2 )
			{
				my ($read1, $read2) = ($reads[0], $reads[1]);
				my ($sai1, $sai2) = ($read1, $read2);

				if ($sai1 =~ m/\.gz$/) {  $sai1 =~ s/\.gz$//; }
				if ($sai2 =~ m/\.gz$/) {  $sai2 =~ s/\.gz$//; } 		
			
				$bam = $sai1; $sort = $sai1;

				$sai1 =~ s/\.fa$/\.sai/;
				$sai2 =~ s/\.fa$/\.sai/;
	                        $bam =~ s/\.fa$/\.bam/;
        	                $sort =~ s/\.fa$/_sort/;
                	        $sort_bam = $sort.".bam";

				my $bwa_align_cmd1 = "bwa aln -t $cpu -n 0.02 -o 1 -e 2 -f $sai1 $genome $read1";
				my $bwa_align_cmd2 = "bwa aln -t $cpu -n 0.02 -o 1 -e 2 -f $sai2 $genome $read2";
				$pileup_cmds.=$bwa_align_cmd1."\n";
				$pileup_cmds.=$bwa_align_cmd2."\n";

				my $bwa_sam_cmd = "bwa sampe $genome $sai1 $sai2 $read1 $read2 | filter_for_PEsnp.pl | samtools view -bS -o $bam -";
				$pileup_cmds.=$bwa_sam_cmd."\n";

				my $sort_cmd = "samtools sort $bam $sort";
				$pileup_cmds.=$sort_cmd."\n";
			
				push(@sort_bam, $sort_bam);
			}
			elsif ( scalar(@reads) == 1 )
			{
				my $read = $a[$i];
				my $sai = $read;
				if ($sai =~ m/\.gz$/) {  $sai =~ s/\.gz$//; }

				my ($bam, $sort) = ($sai, $sai);

				$sai =~ s/\.fa$/\.sai/;
				$bam =~ s/\.fa$/\.bam/;
				$sort =~ s/\.fa$/_sort/;
				$sort_bam = $sort.".bam";
				my $bwa_align_cmd = "bwa aln -t $cpu -n 0.02 -o 1 -e 2 -f $sai $genome $read";
				$pileup_cmds.=$bwa_align_cmd."\n";
			
				my $bwa_sam_cmd = "bwa samse $genome $sai $read | filter_for_SEsnp.pl | samtools view -bS -o $bam -";
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
 filter_RNASeq
=cut 
sub filter_RNASeq
{
	my ($input_file, $output_file) = @_;

	my $output_line = "type\tSNP\tchromosome\tposition\tRef base\ts1 coverage\ts1 base\ts2 coverage\ts2 base\n";

	my $fh = IO::File->new($input_file) || die $!;
while(<>) {
        chomp;
        @a = split "\t";
        $count1 = ($a[6] =~ tr/\^//);
        $count1 += ($a[6] =~ tr/\$//);
        $count2 = ($a[8] =~ tr/\^//);
        $count2 += ($a[8] =~ tr/\$//);
        $cov1 = $a[5] - $count1;
        $cov2 = $a[7] - $count2;
        $a[6] =~ tr/agctn/AGCTN/;
        $a[6] =~ s/"//g;
        $a[6] =~ s/\$//g;
        $a[6] =~ s/\^://g;
        $a[6] =~ s/\^F//g;
        $a[6] =~ s/\^\d//g;
        $a[6] =~ s/\^\)//g;
        $a[6] =~ s/\.\+/\+/g;
        $a[6] =~ s/\.\-/\-/g;
        $a[6] =~ s/,\+/\+/g;
        $a[6] =~ s/,\-/\-/g;
        $a[6] = " ".$a[6];
        $a[8] =~ tr/agctn/AGCTN/;
        $a[8] =~ s/"//g;
        $a[8] =~ s/\$//g;
        $a[8] =~ s/\^://g;
        $a[8] =~ s/\^F//g;
        $a[8] =~ s/\^\d//g;
        $a[8] =~ s/\^\)//g;
        $a[8] =~ s/\.\+/\+/g;
        $a[8] =~ s/\.\-/\-/g;
        $a[8] =~ s/,\+/\+/g;
        $a[8] =~ s/,\-/\-/g;
        $a[1] =~ s/;//;
        $a[6] =~ s/\*/\-/g;
        $a[6] =~ s/,/\*/g;
        $a[6] =~ s/\./\*/g;
        $a[6] =~ s/\*/$a[4]/g;
        $a[8] =~ s/\*/\-/g;
        $a[8] =~ s/,/\*/g;
        $a[8] =~ s/\./\*/g;
        $a[8] =~ " ".$a[8];
        $a[8] =~ s/\*/$a[4]/g;
        if ($cov1 >= 4 && $cov2 >= 4) {
                print join("\t", @a), "\n";
        }
}



}










sub pipeline
{
	print qq'
===== A original SNP calling pipeline from Mao =====

1. remove redundancy reads using removeRedundancy.pl
   $ perl removeRedundancy.pl input > output

2. align each sample to reference using bwa
   $ bwa aln -t 24 -n 0.02 -o 1 -e 2 -f sample1.sai reference.fa sample1.fa
   $ bwa samse reference.fa sample1.sai sample1.fa | filter_for_SEsnp.pl | samtools view -bS -o sample1.bam -
   $ samtools sort sample1.bam sample1_sort

3. merge sorted bam files for each cultivar
   $ samtools merge -f cultivar_A.bam sample1_sort.bam sample2_sort.bam ...... sampleN_sort.bam

4. generate pileup files for each cultivar
   $ samtools mpileup -q 16 -Q 0 -d 10000 -f reference.fa cultivar_A.bam > cultivar_A.pileup

* choose one or more below step for next analysis.

5. generate virtual genome using SNP.
   $ reSeqPrintSample.indel.fast.strAssign.RNAseq.table reference.fa cultivar_A.1col cultivar_A.pileup cultivar_A 3 3 0.3

6. call SNPs between cultivar and reference
   $ pileupFilter.AtoG 0.9 0.8 3 cultivar_A.pileup

7. call SNPs between two cultivars
   $ combine2PileFiles cultivar_A.pileup cultivar_B.pileup  0.9  0.8  ChrOrder  3

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

    $ SNPmTool.pl list_file  comparison_file > run_cmd.sh

    Edit the run_cmd.sh file if required

    ./run_cmd.sh

';

}


