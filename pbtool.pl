#!/usr/bin/perl

use strict;
use warnings;
use IO::File;
use Bio::SeqIO;
use Getopt::Std;

=head basic information
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

=head2
 pb_correct Filter -- filter PB aligned reads that do not need correct
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
USAGE: $0 -t correctF [options] pb_align.sam illumina_align.sam

	-e	max edit distance ratio (default: 0.2)
	-d	depth of Illumina junction (default: 1)
		If set -d to 0, reads will be correct/filter if there none or 
		partially of junctions were supported by Illumina reads. But, 
		reads will not be correct/filter if 5bp continouse edit distance 
		present in the range of junction site.

	-m	max mismatch ratio (default: 0.1)
	-c	max soft clip length (default: 100)
	-v	max soft clip coverage (default: 0.1)
	
	-k	0. correct mode. (default) 

		1. check mode, only report mapping statistics info
		do no perform filter/correction

';
	print $usage and exit unless defined $$files[0] && defined $$files[1];
	my $pb_sam = $$files[0];
	my $sr_sam = $$files[1];

	# set parameters
	my $sr_junc_depth = 1;		# junction depth for illumina reads	
	my $juction_flank_len = 5;	# search range of edit distance around junction
	my $con_junc_ed = 5;            # continous edit distance in junction range
	my $ed_cutoff = 0.20;		# reads with more than 25% edit distance will be ignore
	my $mismatch_cutoff = 0.10;	# reads with more than 20% mismatch will be ignore
	my $clip_cutoff = 100;		# reads with more than 100 soft clip will be ignore
	my $clip_cov    = 0.1;		# reads with more than 10% soft clip will be ignore

	$sr_junc_depth = $$options{'d'} if (defined $$options{'d'} && $$options{'d'} =~ m/^\d+$/);

	my $mode = 0;
	$mode = 1 if (defined $$options{'k'} && $$options{'k'} == 1);

	if ($mode == 1) {
		print "ID\tStrand\tRef\tPos\tLength\tSoftClipLen\tEditDist\tInsertion\tDeletion\tMismatch\n";
	}

	# load illumina junctions
	my %sr_junction;
	if ($sr_junc_depth > 0) {
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

	# set output
	my $pb_correct_seq = '';
	my $pb_raw_seq = '';
	my $pb_correct_sam = '';

	# parse the pb sam file
	my ($ed, $md, $nm);
	my $fh0 = IO::File->new($pb_sam) || die $!;
	while(<$fh0>)
	{
		chomp;
		if ($_ =~ m/^@/) { $pb_correct_sam.=$_."\n"; next; }
		my @a = split(/\t/, $_);
		next if $a[1] eq '4';
		my $cigar = $a[5];
		if ($_ =~ m/NM:i:(\d+)/) { $ed = $1; } else { die "[ERR]in NM $_\n"; }
		if ($_ =~ m/MD:Z:(\S+)/) { $md = $1; } else { die "[ERR]in MD $_\n"; }

		# get reads align start point
		my $apos = $a[3];	# start of align positon
		my ($left_clip, $right_clip) = (0, 0);
		if ($cigar =~ m/^(\d+)S/) { $left_clip  = $1; }
		if ($cigar =~ m/(\d+)S$/) { $right_clip = $1; }

		# get junction, and insertion into array  
		# junction : 
		# insertion:
		my ($junction, $insertion) = parse_cigar($cigar);

		# correct sequence 
		# 1. remove soft clip
		# 2. remove insertion according result by parse_cigar
		# 3. add deletion and fix mismatch according to MD:Z 
		my $seq = $a[9];
		my $total_clip = $left_clip + $right_clip;
		$seq = substr($seq, $left_clip, (length($seq) - $total_clip)); # remove soft clip
		my ($correct_seq, $md_correct) = correct_pb_seq($seq, $insertion, $md);

		# get insertion, deletion, mismatch base num
		my $ins_num = 0;
		foreach my $ins (sort keys %$insertion) { $ins_num+=$$insertion{$ins}; }

		my ($del_num, $mis_num) = (0, 0);
		foreach my $k (sort keys %$md_correct) {
			if ($$md_correct{$k} =~ m/^\^/) { $del_num+=length($$md_correct{$k})-1; }
			else { $mis_num+=length($$md_correct{$k}); }
		}

		# report statistics info
		if ($mode == 1) {
			print "$a[0]\t$a[1]\t$a[2]\t$apos\t",length($seq),"\t$total_clip\t$ed\t$ins_num\t$del_num\t$mis_num\n";
			next;
		}

		# find abnormal junction
		# 1. check if the pb junction support by illumina
		# 2. check edit distance in in junction

		my $unsupport_junc = 0;	# using this method, will report good single exon reads

		foreach my $j (@$junction) 
		{
			#      | <- left junction base, 1 base
			# ++++++---------++++++
			my $left_jpos = $j->[0];
			my $jstart = $apos + $j->[1];
			my $jend   = $apos + $j->[2];	
			# print "$left_jpos\t$jstart\t$jend\n";
			
			# check if the pb junction support by illumina 
			my $junc_key = $a[2]."#".$jstart."#".$jend;
			if ($sr_junc_depth > 0) 
			{
				if ( defined $sr_junction{$junc_key} ) {
					$unsupport_junc = 1 if $sr_junction{$junc_key} < $sr_junc_depth;
				} else {
					$unsupport_junc = 1;
				}
			}

			# check if the pb junction have mismatches
			my $connect_error = 0;
			my %connect_hash = ();
			my $p;
			for($p=$left_jpos-$juction_flank_len+1; $p<=$left_jpos+$juction_flank_len; $p++) 
			{
				if (defined $$md_correct{$p}) {
					$connect_error++;
					#print $a[0]."\t".$jstart."-".$jend."\t".$left_jpos."-".$p."-".$$md_correct{$p}."\n"; # mismatch error arround junction
				} else {
					$connect_hash{$p-1} = $connect_error if $connect_error > 0;
					$connect_error = 0;
				}
			}
			$connect_hash{$p-1} = $connect_error if $connect_error > 0;

			foreach my $p_end (sort keys %connect_hash) {
				my $len = $connect_hash{$p_end};
				my $p_start = $p_end - $len + 1;
				if ($p_start < $left_jpos && $p_end > $left_jpos+1 && $len >= $con_junc_ed) {
					$unsupport_junc = 1;
					#print "$a[0]\t$left_jpos --- $p_start -- $p_end\t$len\n";
				}
			}
		}

		# rc correct_seq according to strand
		$correct_seq = reverse_comp($correct_seq) if $a[1] == 16;
		my $raw_seq = $a[9]; $raw_seq = reverse_comp($a[9]) if $a[1] == 16;

		# === output result ===
		# Rule1:  
		# remove soft clip, hard clip read
		# remove reads according to edit distance, and number of mismatch
		# remove reads that not all junction supported by illumina
		if ($total_clip > $clip_cutoff || ($total_clip/length($correct_seq)) > $clip_cov || 
		    ($ed/length($correct_seq)) > $ed_cutoff || ($mis_num/length($correct_seq)) > $mismatch_cutoff || 
		    $cigar =~ m/H/ || $unsupport_junc == 1) 
		{
			$pb_raw_seq.= ">$a[0] $a[2]:$apos\n$raw_seq\n";
		} else {	
			$pb_correct_seq.= ">$a[0] $a[2]:$apos:$left_clip:$right_clip\n$correct_seq\n";
			$pb_correct_sam.=$_."\n";
		}
	}
	$fh0->close;

	exit if $mode == 1;

	# save output to file for next correction
	my $out_prefix = $pb_sam; $out_prefix =~ s/\.sam$//;
	my $pb_raw_file = $out_prefix."_raw.fasta";
	my $pb_correct_file = $out_prefix."_correct.fasta";
	my $pb_correct_samf = $out_prefix."_correct.sam";

	save_file($pb_raw_seq, $pb_raw_file);
	save_file($pb_correct_seq, $pb_correct_file);
	save_file($pb_correct_sam, $pb_correct_samf);
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
	return ($temp_seq2, \%md_correct);
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

	formatCCS	format ccs for next analysis
	correctF	correct and filter the PB reads


';

	print $usage;
	exit;
}



