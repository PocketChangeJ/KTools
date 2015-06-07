#!/usr/bin/perl

use strict;
use warnings;
use IO::File;
use Bio::SeqIO;
use Getopt::Std;

=head1 basic information
 Description of PB CCS
 1. the CCS reads will be named as "reads_of_insert"
 2. the CCS shorter than 300bp, and chimeric CSS will be discard, 
    the left reads were named as "isoseq_draft.fasta"
 3. the isoseq_draft.fasta will be classified as full length and 
    non-full length :
	isoseq_flnc.fasta
	isoseq_nfl.fasta

	5p end AAGCAGTGGTATCAACGCAGAGTACATGGG
	3p end GTACTCTGCGTTGATACCACTGCTT
	the 5p and 3p are same, and could be used for get strand
=cut

my %options;
getopts('a:b:c:d:e:f:g:i:j:k:l:m:n:o:p:q:r:s:t:u:v:w:x:y:z:h', \%options);
unless (defined $options{'t'} ) { usage(); } # set default tool

if	($options{'t'} eq 'formatCCS')	{ pb_format_ccs(\%options, \@ARGV); }
if	($options{'t'} eq 'align')	{ pb_align(); }
if      ($options{'t'} eq 'correctF')	{ pb_correctF(\%options, \@ARGV); }
else	{ usage(); }

#################################################################
# kentnf: subroutine						#
#################################################################

=head1 pb_correctF -- filter PB aligned reads
 1. filter out pb reads that does not need correction
	1.1 do not have any edit distance in best align
	1.2 with edit distance, but no presented in splicing site range
	1.3 also have error in splicing site range, but splicing site 
	    confirmed by illumima
 2. 
=cut
sub pb_correctF
{
	my ($options, $files) = @_;
	my $usage = qq'
USAGE: $0 -t correctF [options] pb_align.sam

* Before correction, below info need to be investigated
- min length of intron 
- max length of intron
- splice site frequency

Filter for ERR rate
	-e	max percentage of error		(default: 0.3)

Filter for junction filter
	-i	illumina reads alignment, SAM format
	-c	tophat aligned junctions, junction format converted using bed_to_junc
	-d	depth of Illumina junction (default: 1)
		[If set -d to 0, reads will be correct/filter if there none or
		 partially of junctions were supported by Illumina reads. But,
		 reads will not be correct/filter if 5bp continouse edit distance
		 present in the range of junction site.]


	-g	genome annotation, GTF format 	(required)
	-r	reference, fasta format		(required)
		[the reference must index by samtools]
	-a	min length of intron		(default: 40)
	-b	max length of intron		(default: 10000)
	-s	splice site 			(default: GT-AG,GC-AG) 
	
	-k	0. correct mode. (default) 

		1. check mode, only report mapping statistics info
		do no perform filter/correction

		2. only report soft clip reads

		report the statiscis of errors in junction site

';
	# parse required files
	my ($pb_sam, $gtf_file, $reference, $out_prefix);
	if (defined $$files[0] && defined $$options{'g'} && defined $$options{'r'} ) { 
		$pb_sam    = $$files[0]; 
		die "[ERR] input file suffix $pb_sam\n" unless $pb_sam =~ m/\.sam$/;
		$out_prefix = $pb_sam; 
		$out_prefix =~ s/\.sam$//;
		$gtf_file  = $$options{'g'};
		$reference = $$options{'r'};
		check_file($pb_sam, $gtf_file, $reference);
	}else { 
		print $usage; exit;
	}

	# set parameters
	my $ed_cutoff = 0.30;           # reads with more than 25% edit distance will be ignore
	$ed_cutoff = $$options{'e'} if (defined $$options{'e'} && $$options{'e'} >= 0);

	my ($min_intron, $max_intron) = (40, 10000);
	my %conserved_splice = (
		'GT-AG' => 1,
		'CT-AC' => 1,
		'GC-AG' => 1,
		'CT-GC' => 1
	);

	$min_intron = $$options{'e'} if (defined $$options{'a'} && $$options{'a'} >= 0);
	$max_intron = $$options{'b'} if (defined $$options{'b'} && $$options{'b'} >= $min_intron);
	if (defined $$options{'s'})
	{
		my @cs = split(/,/, $$options{'s'});
		foreach my $cs (@cs) {
			$conserved_splice{$cs} = 1;
			my $rcs = reverse $cs;
			$rcs = tr/ATCGNatcgn-/TAGCNtagcn-/;
			$conserved_splice{$rcs} = 1;
		}
	}

	my ($sr_sam, $sr_align_junc) = ('', '');
	$sr_sam 	= $$options{'i'} if defined $$options{'i'};
	$sr_align_junc 	= $$options{'c'} if defined $$options{'c'};
	my $sr_junc_depth = 1;          # junction depth for illumina reads
	$sr_junc_depth = $$options{'d'} if (defined $$options{'d'} && $$options{'d'} =~ m/^\d+$/);

	my $mode = 0;
	$mode = 1 if (defined $$options{'k'} && $$options{'k'} == 1);
	$mode = 2 if (defined $$options{'k'} && $$options{'k'} == 2);

	if ($mode == 1) {
		print "ID\tStrand\tRef\tPos\tLength\tSoftClipLen\tEditDist\tInsertion\tDeletion\tMismatch\n";
	}

	#########################################################
	# parameters pre defined				#
	#########################################################
	my $juction_flank_len = 2;      # search range of edit distance (ERR) around junction
	my $con_junc_ed = 5;            # cutoff for continous ERR 

	#########################################################
	# load junctions					#
	#########################################################
	# load illumina junctions by reads or by alignment software
	my %sr_junction;
	if (-s $sr_sam) {
		my $fh1 = IO::File->new($sr_sam) || die $!;
		while(<$fh1>)
		{
			chomp;
			next if $_ =~ m/^@/;
			my @a = split(/\t/, $_);
			next if $a[1] eq '4';
			my $apos = $a[3];
			my $cigar = $a[5];
			my ($junction, $insertion) = parse_cigar($cigar);
			next if scalar @$junction == 0;
			foreach my $j (@$junction) {
				my $jstart = $apos + $j->[1];
				my $jend   = $apos + $j->[2];
				if (defined $sr_junction{$a[2]."#".$jstart."#".$jend} ) {
					$sr_junction{$a[2]."#".$jstart."#".$jend}++;
				} else {
					$sr_junction{$a[2]."#".$jstart."#".$jend} = 1;
				}
			}
		}
		$fh1->close;
	} 

	if (-s $sr_align_junc) {


	}
	#print scalar(keys(%sr_junction)), " of junctions has been loaded from short read mapping\n"; #exit;

	# load junctions from gene annotation file
	my %gtf_junction;
	my %trans_info = parse_gtf($gtf_file);

		foreach my $tid (sort keys %trans_info) {
			my $chr = $trans_info{$tid}{'chr'};
			my @exon = split("\t", $trans_info{$tid}{'exon'});
			@exon = sort {$a<=>$b} @exon;
			next if @exon == 2;
			die "[ERR]exon num $tid\n" unless (scalar(@exon) % 2 == 0);
			shift @exon; pop @exon;
			for(my $i=0; $i<@exon; $i=$i+2) {
				my $jstart = $exon[$i]+1;
				my $jend = $exon[$i+1]-1;
				$gtf_junction{$chr."#".$jstart."#".$jend} = 1;
			}
		}
	#print scalar(keys(%gtf_junction)), " of junctions has been loaded from $$options{'g'}\n"; #exit;

	# store junction attribute information to hash
	# key1: read id
	# key2: RSEQ, raw read sequence
	#	CSEQ, correct read sequence
	#	SOFT, array of soft clip in end [left, right]
	#	HARD, array of hard clip in end [left, right]
	#	EDIT, edit distance
	#	JUNC, array of junctions
	# 	STR,  mapping strand 0 or 16
	#	CHR,  chr name 
	my %rid_attr;

	# store junction attribute information to hash
	# key1 : junction_key, 
	# key2 : PBS, PacBio Supported coverage
	#	 GTF, iTAG supported stat
	#	 SRS, illumina short reads converage
	#	 SPL, splicing base GT-AG
	#	 ERR, ERR count
	#	 CON, connected ERR count  
	#	 MAP: array of transcript ID:strand, separated by ","
	my %junc_attr;

	my ($count_hard_clip, $count_single_align) = (0, 0);

	# parse the pb sam file
	my ($ed, $md, $nm); # init var for edit distance (NM), MD
	my $fh0 = IO::File->new($pb_sam) || die $!;
	while(<$fh0>)
	{
		chomp;
		next if $_ =~ m/^@/;
		my @a = split(/\t/, $_);
		next if $a[1] eq '4';
		my $cigar = $a[5];
		
		# skip hard clip for it has two alignment in reads 
		$count_hard_clip++ and next if $cigar =~ m/H/;	
		$count_single_align++;
		
		# get ed md for compute edit distance, (mismatch and indel)
		if ($_ =~ m/NM:i:(\d+)/) { $ed = $1; } else { die "[ERR]in NM $_\n"; }
		if ($_ =~ m/MD:Z:(\S+)/) { $md = $1; } else { die "[ERR]in MD $_\n"; }

		# get reads align start point
		my $apos = $a[3];	

		# get soft clip
		my ($left_clip, $right_clip) = (0, 0);
		if ($cigar =~ m/^(\d+)S/) { $left_clip  = $1; }
		if ($cigar =~ m/(\d+)S$/) { $right_clip = $1; }

		# get junction, and insertion into array  
		# junction : array of junctions for each reads
		#	     member format: [ Left,Jstart,Jend ]
		#      | <- left junction base, 1 base
                # ++++++---------++++++
                #       |       | <- Jstart and Jend
		# insertion:
		my ($junction, $insertion) = parse_cigar($cigar);

		# correct sequence 
		# 1. remove soft clip
		# 2. remove insertion according result by parse_cigar
		# 3. add deletion and fix mismatch according to MD:Z 
		my $seq = $a[9];
		my $total_clip = $left_clip + $right_clip;
		$seq = substr($seq, $left_clip, (length($seq) - $total_clip)); # remove soft clip
		my ($correct_seq, $md_correct, $md_hash) = correct_pb_seq($seq, $insertion, $md);

		# get insertion, deletion, mismatch base num
		my $ins_num = 0;
		foreach my $ins (sort keys %$insertion) { $ins_num+=$$insertion{$ins}; }

		my ($del_num, $mis_num) = (0, 0);
		foreach my $k (sort keys %$md_hash) {
			if ($$md_hash{$k} =~ m/^\^/) { $del_num+=length($$md_hash{$k})-1; }
			else { $mis_num+=length($$md_hash{$k}); }
		}

		# report statistics info for infer the cutoff for quality control
		if ($mode == 1) {
			print "$a[0]\t$a[1]\t$a[2]\t$apos\t",length($seq),"\t$total_clip\t$ed\t$ins_num\t$del_num\t$mis_num\n";
			next;
		}

		# put the juction information to hash
		my @abs_junction;
		foreach my $j (@$junction) 
		{
			#      | <- left junction base, 1 base
			# ++++++---------++++++
			#       |       | <- Jstart and Jend
			my $left_jpos = $j->[0];
			my $jstart = $apos + $j->[1];
			my $jend   = $apos + $j->[2];
			push(@abs_junction, [$jstart, $jend]);
			#print "$left_jpos\t$apos\t$jstart\t$jend\n"; exit;
			
			# check if the pb junction support by illumina 
			my $junc_key = $a[2]."#".$jstart."#".$jend;

			#my $sr_support_stat = 0;
			#if ( $sr_junc_depth > 0 && defined $sr_junction{$junc_key} ) {
			#	$sr_support_stat = 1 if $sr_junction{$junc_key} >= $sr_junc_depth;
			#} 

			# check if the pb juction supported by annotation
			#my $gtf_support_stat = 0;
			#if (defined $gtf_junction{$junc_key}) {
			#	$gtf_support_stat = 1;
			#} 

			# generate attribute for junctions
			if (defined $junc_attr{$junc_key}{'PBS'})
			{
				$junc_attr{$junc_key}{'PBS'}++;
				$junc_attr{$junc_key}{'MAP'}.=",$a[0]:$a[1]";
			}
			else
			{
				# init ERR count
				$junc_attr{$junc_key}{'ERR'} = 0; # errors in window
				$junc_attr{$junc_key}{'CON'} = 0; # connect erros in window

				# get support information
				$junc_attr{$junc_key}{'PBS'} = 1; # coverage for junction supported by PB reads (PacBio Support)
				$junc_attr{$junc_key}{'GTF'} = 0;
				$junc_attr{$junc_key}{'SRS'} = 0;
				$junc_attr{$junc_key}{'GTF'} = 1 if defined $gtf_junction{$junc_key}; 
				$junc_attr{$junc_key}{'SRS'} = $sr_junction{$junc_key} if defined $sr_junction{$junc_key};
			
				# get basic information	
				$junc_attr{$junc_key}{'MAP'} = "$a[0]:$a[1]";
				my $intron_len = $jend-$jstart+1;
				$junc_attr{$junc_key}{'LEN'} = $intron_len;
			} 

			# check if the pb junction have mismatches
			my $connect_error = 0;
			my %connect_hash = (); # key, err_start_pos, value, error length
			my $p;
			my $err_pos = 0;
			my $err_stat = 0;
			for($p=$left_jpos-$juction_flank_len+1; $p<=$left_jpos+$juction_flank_len; $p++) 
			{
				$err_pos++;
				next if $p < 1;
				if (defined $$md_correct{$p}) {
					$err_stat = 1;
					$connect_error++;
					#print $a[0]."\t".$jstart."-".$jend."\t".$left_jpos."-".$p."-".$$md_correct{$p}."\n"; # mismatch error arround junction
				} else {
					$connect_hash{$err_pos-1} = $connect_error if $connect_error > 0;
					$connect_error = 0;
				}
			}
			$connect_hash{$err_pos-1} = $connect_error if $connect_error > 0;

			if ($err_stat == 1) {
				$junc_attr{$junc_key}{'ERR'}++;
			}

			# count the connected errors
			my $conn_filter_length_cutoff = 2;
			my $conn_filter_stat = 0;

			foreach my $p_end (sort keys %connect_hash) {
				my $len = $connect_hash{$p_end};
				my $p_start = $p_end - $len + 1;
				# check if the continous error cover the junction site
				if ($p_start < $left_jpos && $p_end > $left_jpos+1 && $len >= $conn_filter_length_cutoff) {
					$conn_filter_stat = 1;
				} 
			}

			if ($conn_filter_stat == 1) {
                                $junc_attr{$junc_key}{'CON'}++;
                        }
		}

		# rc correct_seq according to strand, it's necessary 
		$correct_seq = reverse_comp($correct_seq) if $a[1] == 16;
		my $raw_seq = $a[9]; $raw_seq = reverse_comp($a[9]) if $a[1] == 16;

		# get soft clip sequences
		my $left_soft_seq  = substr($a[9], 0, $left_clip);
                my $right_soft_seq = substr($a[9], -$right_clip);

		# save into to %rid_attr
		$rid_attr{$a[0]}{'CHR'}  = $a[2];
		$rid_attr{$a[0]}{'STR'}  = $a[1];
		$rid_attr{$a[0]}{'CSEQ'} = $correct_seq;
		$rid_attr{$a[0]}{'SOFT'} = [$left_clip, $left_soft_seq, $right_clip, $right_soft_seq];
		$rid_attr{$a[0]}{'EDIT'} = [$ed, $ins_num, $del_num, $mis_num];
		$rid_attr{$a[0]}{'JUNC'} = \@abs_junction;	

		# print to screen how many alignment has been processed (used for speed test)
		my $line_num = $.;
		if ($line_num % 1000 == 0) {
			my $t = $line_num / 1000;
			my $time = localtime();
			warn "[$time] $t x 1000 alignment finished\n";
		}
	}
	$fh0->close;
	#########################################################
	# end of parse PB reads alignment (SAM) file		#
	#########################################################

	$count_hard_clip = $count_hard_clip / 2;
	print "No. of single align: $count_single_align\n";
	print "No. of hard-clip: $count_hard_clip\n";
	print "No. of Read Attr: ", scalar(keys(%rid_attr)),"\n";
	print "No. of GTF Junction: ",scalar(keys(%gtf_junction)),"\n";
	print "No. of Short Read Junction: ",scalar(keys(%sr_junction)),"\n";
	print "No. of PB Read Junction: ",scalar(keys(%junc_attr)),"\n";

	# get splicing site for each junction
	retrive_splice_site($reference, \%junc_attr);

	# summary and report junction status
	if ($mode == 0) {
		my $report_junction_file = $out_prefix."_junc_report.txt";
		my $rjf = IO::File->new(">".$report_junction_file) || die $!;
		print $rjf "#Junction\tSPL\tLEN\tGTF\tPBS\tSRS\tERR\tCONER\tMAP\n";
		foreach my $jk (sort keys %junc_attr) {
			print $rjf $jk,"\t",$junc_attr{$jk}{'SPL'}
				,"\t",$junc_attr{$jk}{'LEN'},"\t",$junc_attr{$jk}{'GTF'},"\t",$junc_attr{$jk}{'PBS'}
				,"\t",$junc_attr{$jk}{'SRS'},"\t",$junc_attr{$jk}{'ERR'},"\t",$junc_attr{$jk}{'CON'},"\t",$junc_attr{$jk}{'MAP'},"\n";
		}
		$rjf->close;
	}

	# filter out corrected sequence file
	# PB reads filter	
	# Before filter, please combine the PB reads (stageB and MG), SR junctions(stageB and MG)	
	# A. remove all reads with hard clip (will perform fusion gene analysis in the further)	
	# B. remove all soft clip part, rename the reads ID for soft clip	
	#	example: B1k43249R101P-L200R300
	#	L: soft clip of 200bp from left
	#	R: soft clip of 300bp from right 
	#	output the clipped sequence for virus-mRNA analysis in the further
	# C. remove reads less than 300bp after soft clipping 	
	# D. For junctions	
	#	1. keep the junctions supported by GTF, illumina short Reads
	#	2. keep the junctions supported by more than 2 PB reads supported
	#	3. for single PB read supported junctions, keep the junction with conserved splicing site (GT-AG, GC-AG)
	#	4. for single PB reads supported unconsrved junctions, keep the juncsion without Errors in 4 base window
	#	from step2, remove juctions with longer or shorter intron length (junction not conserved with GTF)
	# E. filter reads with more Edit distance (>30% of aligned region)	
	# set output var
        my $pb_correct_seq = '';
	my $pb_softclip_seq = '';
	my $pb_filter_out_rpt1 = "RID\tLeftClip\tRightClip\tRealLength\tEdit\tINS\tDEL\tMIS\t%EDIT\n";
	my $pb_filter_out_rpt2 = "RID\tLeftClip\tRightClip\tRealLength\tEdit\tINS\tDEL\tMIS\t%EDIT\n";
	#my $pb_filter_out_rpt2 = "RID\tCHR\tSTR\tStart\tEnd\tLen\tPBS\tGTF\tSRS\tSplice\tERR\n";

	my ($count_short, $count_err_junc, $count_high_err, $count_low_err) = (0,0,0,0);
	foreach my $rid (sort keys %rid_attr) {
		my $correct_seq = $rid_attr{$rid}{'CSEQ'};
		my $correct_len = length($correct_seq);
		$count_short++ and next if length($correct_seq) < 300;
		
		my $soft = $rid_attr{$rid}{'SOFT'};
		my $edit = $rid_attr{$rid}{'EDIT'};
		my $junc = $rid_attr{$rid}{'JUNC'};	
		my $chr  = $rid_attr{$rid}{'CHR'};
		my $str  = $rid_attr{$rid}{'STR'};
		
		my $edit_four = join("\t", @$edit);
                my $edit_rate = sprintf("%.2f", (($edit->[0]/$correct_len)*100));

		my $with_err_junc = 0;
		my $junc_desc = '';
		foreach my $j (@$junc) {
			my $jstart = $j->[0];
			my $jend = $j->[1];
			my $jk = $chr."#".$jstart."#".$jend;
			die "[ERR]junction is not exist $jk\n" unless defined $junc_attr{$jk}{'PBS'};	
			my $length = $junc_attr{$jk}{'LEN'};
			my $pbs_ct = $junc_attr{$jk}{'PBS'};
			my $gtf_ct = $junc_attr{$jk}{'GTF'};
			my $srs_ct = $junc_attr{$jk}{'SRS'};
			my $splice = $junc_attr{$jk}{'SPL'};
			my $err_ct = $junc_attr{$jk}{'ERR'};
	
			$junc_desc.="$rid\t$chr\t$str\t$jstart\t$jend\t$length\t$pbs_ct\t$gtf_ct\t$srs_ct\t$splice\t$err_ct\n";

			if ($length >= $min_intron && $length <= $max_intron)
			{
				next if defined $conserved_splice{$splice};
				next if $gtf_ct == 1;
				next if $srs_ct >= $sr_junc_depth;
				next if $err_ct == 0;
				next if $pbs_ct > 1;
				next if ($pbs_ct == 1 && $edit_rate <=10);
			}
			else
			{
				next if $gtf_ct == 1;
				if ($pbs_ct > 1 && $srs_ct > 5 && defined $conserved_splice{$splice} && $err_ct == 0) {
					next;
				}
			}

			$with_err_junc = 1;
		}

		# reconstruct new transcript ID
		my $left_clip  = $soft->[0];
		my $left_seq   = $soft->[1];
		my $right_clip = $soft->[2];
		my $right_seq  = $soft->[3];
		
		$pb_softclip_seq.=">$rid-$str-L$left_clip\n$left_seq\n" if $left_clip > 30;
		$pb_softclip_seq.=">$rid-$str-R$right_clip\n$right_seq\n" if $right_clip > 30;
		
		my $new_id = $rid;
		if ($left_clip > 0 || $right_clip > 0) {
			$new_id.= "-".$str."L$left_clip"."R$right_clip";
		}

		if ($with_err_junc == 0) {
			if ($edit_rate <= 25)
			{
				$count_low_err++;
				$pb_correct_seq.=">$new_id\n$correct_seq\n";
				$pb_filter_out_rpt2.="$rid\t$left_clip\t$right_clip\t$correct_len\t$edit_four\t$edit_rate\n";
			}
			else
			{
				$count_high_err++;
				$pb_filter_out_rpt1.="$rid\t$left_clip\t$right_clip\t$correct_len\t$edit_four\t$edit_rate\n";
			}
		} else {
			$count_err_junc++;
			$pb_filter_out_rpt1.="$rid\t$left_clip\t$right_clip\t$correct_len\t$edit_four\t$edit_rate\n";
			#$pb_filter_out_rpt2.="$junc_desc";
		}
	}

	print "No. of Read short than 300: $count_short\n";
	print "       Read with low ERR (<=25%, kept): $count_low_err\n";
	print "       Read with high ERR (>25%, remove): $count_high_err\n";
	print "No. of Read with bad junc: $count_err_junc\n";

	# save output to file for next correction
	my $pb_correct_file = $out_prefix."_correct.fasta";
	my $pb_softclip_file = $out_prefix."_softclip.fasta";

	save_file($pb_correct_seq, $pb_correct_file);
	save_file($pb_softclip_seq, $pb_softclip_file);
	save_file($pb_filter_out_rpt1, $out_prefix."_RM_ReadRpt.txt");
	save_file($pb_filter_out_rpt2, $out_prefix."_KP_ReadRpt.txt");
	#save_file($pb_filter_out_rpt2, $out_prefix."_RM_JuncRpt.txt");
}

# child of pb_correctF
# parse cigar to get junction cite
#            |     | <--- will get this pos 
# +++++++++++-------+++++++
sub parse_cigar
{
	my $cigar = shift;

	my @junction = ();
	my %insertion = ();
	my $num = ""; 
	my $point0 = 0; # point for junction site in transcriptome
	my $point1 = 0;	# point for junction site in genome
	my $point2 = 0; # point for insertion in transcriptome

	for(my $i=0; $i<length($cigar); $i++)
	{
		my $str = substr($cigar, $i, 1);
		if ($str =~ m/\d+/) {
			$num = $num.$str;
		}
		elsif ($str eq "M")
		{
			$point0 = $point0 + $num;
			$point1 = $point1 + $num;
			$point2 = $point2 + $num;
			$num = "";
		}
		elsif ($str eq "N")
		{
			$point1 = $point1 + $num;
			push(@junction, [$point0, $point1-$num, $point1-1]);
			$num = "";
		}
		elsif ($str eq "D")
		{
			$point0 = $point0 + $num;
			$point1 = $point1 + $num;
			$num = "";
		}
		elsif ($str eq "I")
		{
			$insertion{$point2+1} = $num;
			$point2 = $point2 + $num;
			$num = "";
		}
		else 
		{
			$num = 0;
		}
	}
	return(\@junction, \%insertion);
}

# child of pb_correctF
# correct pb seq base on gmap alignment
sub correct_pb_seq
{
	my ($seq, $insert_hash, $md) = @_;

	# convert insert hash
	# befroe convert: key: insert_pos, value: insert length
	# after convert:  key: insert_base, value: 1
	#
	# Insert Pos:
	#     |  <- ins_pos
	# ++++II+++
	# 123456789
	#
	# before convert: key:5, value 2
	# after convert: key:5, value 1, key:6, value 1;
	my %insertion;
	foreach my $ins_pos (sort {$a<=>$b} keys %$insert_hash) {
		my $len_pos = $$insert_hash{$ins_pos};
		for(my $i=$ins_pos; $i<$ins_pos+$len_pos; $i++) {
			$insertion{$i} = 1;
		}
	}

	# remove insertion
	my $temp_seq1 = '';

	for(my $i=1; $i<=length($seq); $i++) {
		$temp_seq1 .= substr($seq, $i-1, 1) unless defined $insertion{$i};
	}

	# add deletion
	my @md_array;
	my $num = ''; my $pre_num = ''; my $md_str = '';
	for(my $i=0; $i<length($md); $i++) 
	{
		my $str = substr($md, $i, 1);
		if ($str =~ m/\d+/) {
			if ( $md_str ) {
				push(@md_array, [$pre_num, $md_str]);
				$md_str = '';
				$num = '';
			}
			$num = $num.$str;
		}else {
			$md_str = $md_str.$str;
			$pre_num = $num;
		}
	}
	
	my %md_hash;	# key, uncorrect seq pos; value, md value 
	my $md_pos = 0;
	my %md_correct;	# key, correct seq position; value, md value
	my $md_correct_pos = 0;
	foreach my $n (@md_array) {
		my $md_rel_pos = $n->[0]; 	# relative position
		my $md_content = $n->[1]; 	# correct base
		$md_pos = $md_pos+$md_rel_pos;	# abs position
		$md_correct_pos = $md_correct_pos + $md_rel_pos;

		if ($md_content =~ m/^\^/) {
			$md_hash{$md_pos}=$md_content;

			my @mm = split(//, $md_content); shift @mm;
			foreach my $m (@mm) {
				$md_correct_pos++;
				$md_correct{$md_correct_pos} = $m;
			}

		} else {
			my @mm = split(//, $md_content);
			for(my $i=0; $i<@mm; $i++) {
				$md_hash{$md_pos+1}=$mm[$i];
				#print "--",$md_pos+1,"\t",$mm[$i],"\n";
				$md_pos++;
			}

			foreach my $m (@mm){
				$md_correct_pos++;
				$md_correct{$md_correct_pos} = $m;
			}
		}
		# print "$md_correct_pos\t$md_content\n";
	}

	my $temp_seq2 = '';

	for(my $i=1; $i<=length($temp_seq1); $i++) 
	{
		my $base = substr($temp_seq1, $i-1, 1);	
		if (defined $md_hash{$i}) 
		{
			if ($md_hash{$i} =~ m/^\^/) {
				$temp_seq2.=$base;
				my $deletion = $md_hash{$i};
				$deletion =~ s/^\^//;
				$temp_seq2.=$deletion;
			}
		}

		if (defined $md_hash{$i})
		{
			unless ($md_hash{$i} =~ m/^\^/) 
			{
				$temp_seq2.=$md_hash{$i};
				#print "---".$i."\t$base->".$md_hash{$i}."\n";
			}
		} 
		else 
		{
			$temp_seq2.=$base;
		}
	}

	# print length($seq),"\n",$temp_seq1,"\n",length($temp_seq1),"\n";
	# print $temp_seq2,"\n",length($temp_seq2),"\n";
	return ($temp_seq2, \%md_correct, \%md_hash);
}

sub save_file
{
	my ($content, $file_name) = @_;
	my $fh = IO::File->new(">".$file_name) || die $!;
	print $fh $content;
	$fh->close;
}

sub reverse_comp
{
	my $seq = shift;
	my $revcomp = reverse($seq);
	$revcomp =~ tr/ACGTNacgtn/TGCANtgcan/;
	return $revcomp;
}

=head2 
 pb_format_ccs -- format ccs read to humand readable format
=cut
sub pb_format_ccs
{
	my ($options, $files) = @_;	
	my $usage = qq'
USAGE: $0 -t formatCCS [options] isoseq_draft.fasta 

	-l	library name (default: PB)

* output files PB.fasta

';
	print $usage and exit unless defined $$files[0];
	my $input_draft = $$files[0];
	die "[ERR]file not exist\n" unless -s $input_draft;

	my $lib = "PB";
	$lib = $$options{'l'} if defined $$options{'l'};

	my $out_seq = $lib.".fasta";
	die "[ERR]out file exist\n" if -s $out_seq;

	my $out = IO::File->new(">".$out_seq) || die $!;

	my ($id, $desc, $seq, $nid, $strand, $p5, $pA, $p3, $chimera, $full);
	my $in = Bio::SeqIO->new(-format=>'fasta', -file=>$input_draft);
	while(my $inseq = $in->next_seq) {
		$id = $inseq->id;
		$desc = $inseq->desc;
		$seq= $inseq->seq;		

		my @a = split(/\//, $id);
		die "[ERR]seqid $id\n" unless @a == 3;
		$nid = $lib.$a[1];

		if ($desc =~ m/strand=(\S+);fiveseen=(\d);polyAseen=(\d);threeseen=(\d);.*;primer=(\S+);chimera=(\S+)/) {
			($strand, $p5, $pA, $p3, $chimera) = ($1, $2, $3, $4, $6);
		} else {
			die "[ERR]seq desc $desc\n";
		}
		
		if 	($strand eq '+') { $strand = "F"; } 
		elsif 	($strand eq '-') { $strand = "R"; }
		else	{$strand = "N";}
		if ($chimera eq '0' ) {$chimera = "F"; } else { $chimera = "P"; }

		$nid.=$strand.$p5.$pA.$p3.$chimera;

		print $out ">$nid\n$seq\n";
	}
	$out->close;
}

=head2
 pb_align -- align pb mRNA reads to reference 
=cut
sub pb_align
{
	my $usage = qq'
# gmap location: /home/kentnf/software/gmap-2014-12-23/bin/gmap
# tomato refidx: /home/kentnf/06lncRNA/gmap/tomato_genome_v2_4

gmap -n 0 -D /home/kentnf/06lncRNA/gmap/tomato_genome_v2_4 -d gmap_db -f samse -t 31 --sam-use-0M input_fasta

# -n number of paths for each reads (0 means single alignment plus chimericalignments)
# -D  genome dir 
# -d  genomr database
# -f output format 
# -t CPU
# --sam-use-0M  Insert 0M in CIGAR between adjacent insertions and deletions Required by Picard, but can cause errors in other tools 
#
#
';
	print $usage; 
	exit;
}

=head2
 usage -- print usage information
=cut
sub usage
{
	my $usage = qq'
USAGE: $0 -t [tool]

	align		show align command used for pb reads
	formatCCS	format ccs for next analysis
	correctF	correct and filter the PB reads

';

	print $usage;
	exit;
}


=head2
 parse_gtf -- parse gtf file, return gtf information
=cut
sub parse_gtf
{
        my $input_file = shift;

        my %trans_info; # key: tid, chr, exon, gene, strand

        my $fh = IO::File->new($input_file) || die $!;
        while(<$fh>)
        {
                chomp;
                next if $_ =~ m/^#/;
                my @a = split(/\t/, $_);

                if ($a[2] eq 'exon')
                {
                        # analyzing the attribute 
                        my %attr_value;
                        my @b = split(/; /, $a[8]);
                        foreach my $b (@b) {
                                my @c = split(/ /, $b);
                                die "[ERR]attr $b in $a[8]\n" unless @c == 2;
                                $c[1] =~ s/"//ig;
                                $attr_value{$c[0]} = $c[1];
                        }

                        die "[ERR]No transcript id\n" unless defined $attr_value{'transcript_id'};
                        die "[ERR]No gene id\n" unless defined $attr_value{'gene_id'};
                        my ($tid, $gid) = ($attr_value{'transcript_id'}, $attr_value{'gene_id'});

                        if ( defined $trans_info{$tid}{'chr'} ) {
                                die "[ERR]inconsistency chr for $tid\n" if $trans_info{$tid}{'chr'} ne $a[0];
                        } else {
                                $trans_info{$tid}{'chr'} = $a[0];
                        }

                        if ( defined $trans_info{$tid}{'gid'} ) {
                                die "[ERR]inconsistency gid for $tid\n" if $trans_info{$tid}{'gid'} ne $gid;
                        } else {
                                $trans_info{$tid}{'gid'} = $gid;
                        }

                        if ( defined $trans_info{$tid}{'strand'} ) {
                                die "[ERR]inconsistency strand for $tid\n" if $trans_info{$tid}{'strand'} ne $a[6];
                        } else {
                                $trans_info{$tid}{'strand'} = $a[6];
                        }

                        if ( defined $trans_info{$tid}{'exon'}) {
                                $trans_info{$tid}{'exon'}.="\t".$a[3]."\t".$a[4];
                        } else {
                                $trans_info{$tid}{'exon'} = $a[3]."\t".$a[4];
                        }
                }
        }
        $fh->close;

        return %trans_info;
}

sub retrive_splice_site
{
	my ($reference, $junc_attr) = @_; 

	# convert junc_attr to junc_start hash
	# key: chr \t start
	# value: array of end; 
	my %junc_start;

	foreach my $jk (sort keys %$junc_attr) {
		my @a = split("#", $jk);
		my $js_key = $a[0]."\t".$a[1];
		push(@{$junc_start{$js_key}}, $a[2]);
	}

        my $in = Bio::SeqIO->new(-format=>'fasta', -file=>$reference);
        while(my $inseq = $in->next_seq)
        {
                my $chr = $inseq->id;
                my $seq = $inseq->seq;
                my $len = $inseq->length;

                # below start is the exon end
                # below end is the exon start
                for(1 .. $len)
                {
                        my $start = $_;

                        if ( defined $junc_start{$chr."\t".$start} )
                        {
                                my @end = @{$junc_start{$chr."\t".$start}};
                                foreach my $end (@end)
                                {
					my $jk = $chr."#".$start."#".$end;
                                        my $start_base = uc(substr($seq, $start-1, 2));
                                        my $end_base = uc(substr($seq, $end-2, 2));
                                        my $junction = $start_base."-".$end_base;
					$$junc_attr{$jk}{'SPL'} = $junction;
                                }
                        }
                }
        }
}

# check if file exist
sub check_file
{
	my @file = @_;
	foreach my $f (@file) {
		unless (-s $f) {
			print "[ERR]file not exist: $f\n";
			exit;
		}
	}
}


