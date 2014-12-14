#!/usr/bin/perl

=head1
 sRNAtool -- tools for sRNA data preparation
=cut

use strict;
use warnings;
use Bio::SeqIO;
use IO::File;
use Getopt::Std;
use FindBin;
use lib "$FindBin::RealBin/m";
use align;

my $version = 0.1;
my $debug = 0;

my %options;
getopts('a:b:c:d:e:g:i:j:k:l:m:n:o:p:q:r:s:t:v:w:x:y:z:fuh', \%options);
#getopts('t:i:o:p:s:c:fuh', \%options);
unless (defined $options{'t'}) { usage($version); }

if	($options{'t'} eq 'chkadp')	{ chkadp(\%options, \@ARGV);  }
elsif	($options{'t'} eq 'chkadp2')    { chkadp2(\%options, \@ARGV); }
elsif	($options{'t'} eq 'rmadp')	{ rmadp(\%options, \@ARGV); }
elsif	($options{'t'} eq 'convert')	{ convert(\%options); }
elsif	($options{'t'} eq 'lengthd')	{ lengthd(\%options); }
elsif	($options{'t'} eq 'unique' )	{ unique(\%options);  }
elsif   ($options{'t'} eq 'norm' )	{ norm(\%options);    }
elsif   ($options{'t'} eq 'normcut')	{ normcut(\%options); }
elsif   ($options{'t'} eq 'combine')	{ combine(\%options); }
elsif	($options{'t'} eq 'range')	{ range(\%options);   }
elsif	($options{'t'} eq 'pipeline')   { pipeline(\%options); }
else	{ usage($version); } 

#################################################################
# kentnf: subroutine						#
#################################################################
=head2
 rmadp -- remove adapter 
=cut
sub rmadp
{
	my ($options, $files) = @_;
	
	my $subUsage = qq'
USAGE: $0 -t rmadp [options] -s adapter_sequence  input_file1 ... input_fileN
	-s	3p adapter sequence [option if -p]
	-p	5p adapter sequence [option if -s]
	-l	3p adapter length (5-11) (default: 9)
	-d	distance between 3p adater and sRNA (default: 1)

* the 5p adapter will perform perfect match 

';

	# checking parameters
	print $subUsage and exit unless defined $$files[0];
	foreach my $f ( @$files ) { print "[ERR]no file $f\n" and exit unless -s $$files[0]; }

	my ($adp_3p, $adp_3p_len, $adp_3p_sub, $adp_5p, $adp_5p_len, $adp_5p_sub);
	$adp_3p_len = 9;
	$adp_5p_len = 5;
	if (defined $$options{'s'}|| defined $$options{'p'}) { } else { print $subUsage and exit; } 
	if (defined $$options{'s'}) {
		$adp_3p = $$options{'s'};
		print "[ERR]short adapter 3p\n" and exit if length($adp_3p) < 9;
		$adp_3p_len = $$options{'l'} if (defined $$options{'l'} && $$options{'l'} > 5);
		$adp_3p_sub = substr($adp_3p, 0, $adp_3p_len);
	} 
	if (defined $$options{'p'} ) {
		$adp_5p = $$options{'p'};
		print "[ERR]short adapter 5p\n" and exit if length($adp_5p) < 9;
		$adp_5p_sub = substr($adp_5p, -$adp_5p_len);
	} 
	
	my $distance = 1;
	$distance = $$options{'d'} if defined $$options{'d'} && $$options{'d'} >= 0;

	my $report_file = 'report_sRNA_trim.txt';
	print "[ERR]report file exist\n" if -s $report_file;
	my $report_info = "#sRNA\ttotal\t5Punmatch\t5Pnull\t5Pmatch\t3Punmatch\t3Pnull\t3Pmatch\tbaseN\tshort\tcleaned\tadp\tadp5\tadp3\n";

	foreach my $inFile ( @$files ) 
	{
		my $prefix = $inFile;
		$prefix =~ s/\.gz$//; $prefix =~ s/\.fastq$//; $prefix =~ s/\.fq$//; $prefix =~ s/\.fasta$//; $prefix =~ s/\.fa$//;

		my $outFile1 = $prefix.".trimmed".$distance;
		my $outFile2 = $prefix.".trimmed".$distance.".report.txt";
		#print "[ERR]out file exist\n" and exit if -s $outFile1;
		my $out1 = IO::File->new(">".$outFile1) || die $!;
		my $out2 = IO::File->new(">".$outFile2) || die $!;
		print $out2 "#ReadID\tlength\tstart\tend\t3p_edit_distacne\tlabel\n";

		my $fh;
		if ($inFile =~ m/\.gz$/) { 
			$fh = IO::File->new("gunzip -c $inFile | ") || die $!;
		} else { 
			$fh = IO::File->new($inFile) || die $!;
		}
	
		my ($format, $id1, $seq, $id2, $qul, $read_len);
		my ($total_num,$unmatch_5p_num,$null_5p_num,$match_5p_num,
			$unmatch_3p_num,$null_3p_num,$match_3p_num,
			$baseN_num,$short_num,
			$clean_num,$adp_clean,$adp5p_clean,$adp3p_clean)=(0,0,0,0,0,0,0,0,0,0,0,0,0);
		while(<$fh>)
		{
			$id1 = $_;	chomp($id1);
			$seq = <$fh>;	chomp($seq);	$seq = uc($seq);
			$read_len = length($seq);
	
			if      ($id1 =~ m/^>/) { $format = 'fasta'; $id1 =~ s/^>//; }
			elsif   ($id1 =~ m/^@/) { $format = 'fastq'; $id1 =~ s/^@//; }
			else    { die "[ERR]seq format $id1\n"; }

			# match 3' adapter to reads,
			# this method will find the best adapter
			my ($pos_3p, $match_ed);
			if (defined $adp_3p) {
				for (my $d=0; $d <=$distance; $d++) 
				{
					for (my $i=0; $i<(length($seq)-$adp_3p_len+1); $i++)
					{
						my $read_substr = substr($seq, $i, $adp_3p_len);
                        			my $edit_distance = hamming($read_substr,$adp_3p_sub);
                        			if ( $edit_distance <= $d ){
                                			$pos_3p = $i;
							last;
                        			}
					}
					$match_ed = $d;
					last if defined $pos_3p;
				}
			}
			$pos_3p = length($seq) unless defined $pos_3p;

			# match 5' adapter to reads
			my $pos_5p = 0;
			if (defined $adp_5p)
			{
				# find the position locate with 10 base seed
			

				# find the position locate with 5 base seed
				my ($match_start, $match_end, $match_len, $match_seq, $match_adp, $match_5p_ed, $pre_match_end);
				$match_start = 0;
				$match_end = 0;
				$pre_match_end = 0;
				while($seq =~ m/\Q$adp_5p_sub\E/g) {
					$match_end = pos($seq) - 1;
					#print "$seq, $adp_5p_sub\n$match_end\n";
					$match_len = $match_end - $match_start + 1;
					$match_seq = substr($seq, $match_start, $match_len);
					$match_adp = substr($adp_5p, -length($match_seq));
					#print "$match_seq, $match_adp\n";
					$match_5p_ed = hamming($match_seq, $match_adp);
					#print "$match_5p_ed\n";
					if ($match_5p_ed > 0) { last; }
					$match_start = $match_end + 1;
					$pre_match_end = $match_start;
				}

				# locate the best match 

				$pos_5p = $pre_match_end;	
			}

			# read id and qul from fastq file
			if ( $format eq 'fastq' ) {
				$id2 = <$fh>;   chomp($id2);	
				$qul = <$fh>;   chomp($qul);
			}

			# output result 
			my $label = '';
			if (defined $adp_5p) {
				if ($pos_5p == 0 ) {
					$label.=",5p_unmatch";
					$unmatch_5p_num++;
				} elsif ($pos_5p == length($seq)) {
					$label.=",5p_null";
					$null_5p_num++;
				} else {
					$label.=",5p_match"; 
					$match_5p_num++;
				}
			}

			if (defined $adp_3p) {
				if ( $pos_3p == 0 ) {
					$label.=",3p_null";
					$null_3p_num++;
				} elsif ($pos_3p == length($seq)) {
					$label.=",3p_unmatch";
					$match_ed = "NA";
					$unmatch_3p_num++;
				} else {
					$label.=",3p_match";
					$match_3p_num++;
				}
			}

			if ($pos_3p > $pos_5p)
			{
				my $trimmed_len = $pos_3p - $pos_5p;
				my $trimmed_seq = substr($seq, $pos_5p, $trimmed_len);
				my $trimmed_qul = substr($qul, $pos_5p, $trimmed_len) if $format eq 'fastq';
				my $baseN = $trimmed_seq =~ tr/N/N/;
				if ($baseN > 0) {
					$label.=",baseN";
					$baseN_num++;
				} elsif (length($trimmed_seq) < 15 ) {
					$label.=",short";
					$short_num++;
				} else {
					$clean_num++;
					if (defined $adp_5p && defined $adp_3p && $pos_5p > 0 && $pos_3p > 0) {
						$adp_clean++;
					} elsif (defined $adp_5p && $pos_5p > 0) {
						$adp5p_clean++;
					} elsif (defined $adp_3p && $pos_3p > 0) {
						$adp3p_clean++;
					}

					if ( $format eq 'fastq' ) {
						print $out1 "@".$id1."\n".$trimmed_seq."\n".$id2."\n".$trimmed_qul."\n";
					} else {
						print $out1 ">".$id1."\n".$trimmed_seq."\n";
					}
				}
			}

			$label =~ s/^,//;
			$total_num++;
			$match_ed = 'NA' unless defined $match_ed;
			print $out2 "$id1\t$read_len\t$pos_5p\t$pos_3p\t$match_ed\t$label\n";
 
		}
		$fh->close;
		$out1->close;
		$out2->close;
		$report_info.="$inFile\t$total_num\t".
			"$unmatch_5p_num\t$null_5p_num\t$match_5p_num\t".
			"$unmatch_3p_num\t$null_3p_num\t$match_3p_num\t".
			"$baseN_num\t$short_num\t$clean_num\t$adp_clean\t$adp5p_clean\t$adp3p_clean\n";
	}

	# report sRNA trim information
	my $outr = IO::File->new(">$report_file") || die $!;
	print $outr $report_info;
	$outr->close;
}

sub hamming($$) { length( $_[ 0 ] ) - ( ( $_[ 0 ] ^ $_[ 1 ] ) =~ tr[\0][\0] )  }

=head2
 chkadp2 -- check adapter using k-mer method
=cut
sub chkadp2
{
	my ($options, $files) = @_;

my $subUsage = qq'
USAGE: $0 chkadp2 [options] input_file
        -p      check 3p/5p adapter using kmer method [default:3]
        -y      number of reads processed (default:1e+5, min)

* the script will detect the 3p adapter by default, and it will detect 5p
  adapter by providing -p paramters, and the 3p detection will disable
* kmer=7 for 5p detection, and kmer=9 for 3p detection

';

	print $subUsage and exit unless $$files[0];
	my $inFile = $$files[0];
	die "[ERR]File not exist\n" unless -s $inFile;
	
	my $read_yeild = 1e+5;
	$read_yeild = int($$options{'y'}) if (defined $$options{'y'} && $$options{'y'} > 1e-5);

	my $strand = 3;
	$strand = 5 if (defined $$options{'p'} && $$options{'p'} ne '3');

	my ($min_len, $max_len) = (7, 10); # scan read length is from 5 to 10

	# load Illumina technical sequence form UniVec
        my %uv_seq;	# key: seq_ID, value: seq;
        my %uv_desc;	# key: seq_ID, value: desc;
	my %uv_match = (); # key: uv_id, value: number of match
	my %uv_match_stat = (); # number of adapter detected with different length
        my $univec = $FindBin::RealBin."/bin/adapters/UniVec";
        die "[ERR]cat not locate univec $univec\n" unless -s $univec;
        my $uv_in = Bio::SeqIO->new(-format=>'fasta', -file=>$univec);
        while(my $inseq = $uv_in->next_seq)
        {
                my $uv_id  = $inseq->id;
                my $uv_seq = $inseq->seq;
                my $uv_desc = $inseq->desc;
                next unless $uv_desc =~ m/Illumina/;
		$uv_seq{$uv_id} = $uv_seq;
                $uv_desc{$uv_id} = $uv_desc;
		$uv_match{$uv_id} = 0;
        }
	
	# compare UniVec with read;
	my $seq_ct = 0; my $format; 
	my $fh = IO::File->new($inFile) || die $!;
        while(<$fh>)
        {
                chomp;
                my $id = $_;
                if      ($id =~ m/^>/) { $format = 'fasta'; }
                elsif   ($id =~ m/^@/) { $format = 'fastq'; }
                else    { die "[ERR]seq format $id\n"; }
                my $seq = <$fh>; chomp($seq); $seq = uc($seq);

		# compare seq with univec, return location
		for(my $i=$max_len; $i>=$min_len; $i=$i-1) 
		{
			my $sub = substr($seq, 0, $i);	
			foreach my $uid (sort keys %uv_seq) 
			{
				my $useq = $uv_seq{$uid};
				while ($useq =~ m/\Q$sub\E/g) {
					if (pos($useq) eq length($useq)) {
						$uv_match{$uid}++;
						if (defined $uv_match_stat{$uid}{$i}) {
							$uv_match_stat{$uid}{$i}++;
						} else {
							$uv_match_stat{$uid}{$i} = 1;
						}
					}
				}
			}
		}

                if ($format eq 'fastq') { <$fh>; <$fh>; }
                $seq_ct++;
                last if $seq_ct == $read_yeild;
        }
        $fh->close;

	#output result
	print "=== adapter detection report for $inFile ===\n";
	foreach my $uid (sort keys %uv_seq) {
		if ($uv_match{$uid} > ($read_yeild/1000)) {
			my $match_stat = "";;
			for($min_len .. $max_len) {
				if ( defined $uv_match_stat{$uid}{$_} ) {
					$match_stat.=" $_:$uv_match_stat{$uid}{$_}";
				}
			}
			print "$uid\t$uv_match{$uid}\t$match_stat\t$uv_desc{$uid}\n$uv_seq{$uid}\n";
		}
	} 
}

=head2
 chkadp -- check adapter using k-mer method
=cut
sub chkadp
{
	my ($options, $files) = @_;
	
	my $subUsage = qq'
USAGE: $0 chkadp [options] input_file
	-p	check 3p/5p adapter using kmer method [default:3]
	-y	number of reads processed (default:1e+5, min)

* the script will detect the 3p adapter by default, and it will detect 5p
  adapter by providing -p paramters, and the 3p detection will disable
* kmer=7 for 5p detection, and kmer=9 for 3p detection

';
	print $subUsage and exit unless $$files[0];
	my $inFile = $$files[0];
	die "[ERR]File not exist\n" unless -s $inFile;

	my $read_yeild = 1e+5;
	$read_yeild = int($$options{'y'}) if (defined $$options{'y'} && $$options{'y'} > 1e-5);
	
	my $strand = 3;
	$strand = 5 if (defined $$options{'p'} && $$options{'p'} ne '3');

	my $min_len = 15; my $kmer_len = 9;
	$kmer_len = 7 if $strand == 5;

	my %kmer_ct;		# key: kmer; value: count;
	my $kmer_most;		# the most freq kmer
	my $kmer_most_freq = 0; # the count of most freq kmer

	my $seq_ct = 0; my $format;
	my $fh = IO::File->new($inFile) || die $!;
	while(<$fh>)
	{
		chomp;
		my $id = $_;
		if	($id =~ m/^>/) { $format = 'fasta'; }
		elsif	($id =~ m/^@/) { $format = 'fastq'; }
		else 	{ die "[ERR]seq format $id\n"; }
		my $seq = <$fh>; chomp($seq); $seq = uc($seq);

		my ($k_count, $k);
		if ($strand == 3) {
			$k_count = length($seq) - $min_len - $kmer_len + 1;
		} else { 
			$k_count = $min_len - $kmer_len + 1;
		}
			
		for(my $i=0; $i<$k_count; $i++) 
		{
			if ($strand == 3) {
				$k = substr($seq, $min_len + $i, $kmer_len);
			} else {
				$k = substr($seq, $i, $kmer_len);
			}
				
			if ( defined $kmer_ct{$k} ) 
			{
				$kmer_ct{$k}++;
				if ($kmer_ct{$k} > $kmer_most_freq) 
				{
					$kmer_most_freq = $kmer_ct{$k};
					$kmer_most = $k;
				}
			}
			else 
			{
				$kmer_ct{$k} = 1;
			}
		}

		if ($format eq 'fastq') { <$fh>; <$fh>; }
		$seq_ct++;
		last if $seq_ct == $read_yeild;
	}
	$fh->close;

	# load Illumina technical sequence form UniVec
	my %uv_kmer;
	my %uv_desc;
	my $univec = $FindBin::RealBin."/bin/adapters/UniVec";
	die "[ERR]cat not locate univec $univec\n" unless -s $univec;
	my $uv_in = Bio::SeqIO->new(-format=>'fasta', -file=>$univec);
	while(my $inseq = $uv_in->next_seq)
	{
		my $uv_id  = $inseq->id;
		my $uv_seq = $inseq->seq;
		my $uv_desc = $inseq->desc;
		next unless $uv_desc =~ m/Illumina/;

		$uv_desc{$uv_id} = $uv_desc;
		my $k_count = $inseq->length - $kmer_len + 1;		
		for(my $j=0; $j<$k_count; $j++) 
		{
			my $uvk = substr($uv_seq, $j, $kmer_len);
			if (defined $uv_kmer{$uvk}) {
				$uv_kmer{$uvk}.= ",$uv_id#$j";
			} else {
				$uv_kmer{$uvk} = "$uv_id#$j";
			}
		}
	}

	print "=== adapter detection report ===\n";
	print "No. of Illumina technical sequences in UniVec ". scalar(keys(%uv_desc))."\n";
	print "No. of $kmer_len-kmer in Illumina technical seq  ".scalar(keys(%uv_kmer))."\n";
	print "No. of read processed: $read_yeild\n";

	# output more than 1 result, and kmer to adapters 

	# process extend kmer most
	my @base = ("A", "T", "C", "G");

	print "=== potential adapters 1: $kmer_most, freq $kmer_most_freq ===\n";
	print "#left-kmer\tcount\tright-kmer\tcount\n";
	my ($pre_k_left, $pre_k_right) = ($kmer_most, $kmer_most);

	for(my $i=0; $i<$min_len; $i++)
	{
		my ($best_kmer_left, $best_kmer_count_left, $best_kmer_right, $best_kmer_count_right);
		$best_kmer_count_left = 0;
		$best_kmer_count_right = 0;

		my $subk_left  = substr($pre_k_left,  0, $kmer_len-1);
		my $subk_right = substr($pre_k_right, 1, $kmer_len-1);

		foreach my $b (@base)
		{
			my $k_left  = $b.$subk_left;
			if (defined $kmer_ct{$k_left} && $kmer_ct{$k_left} > $best_kmer_count_left) {
				$best_kmer_left = $k_left;
				$best_kmer_count_left = $kmer_ct{$k_left};
			}

			my $k_right = $subk_right.$b;
			if (defined $kmer_ct{$k_right} && $kmer_ct{$k_right} > $best_kmer_count_right) {
				$best_kmer_right = $k_right;
				$best_kmer_count_right = $kmer_ct{$k_right};
			}
			
		}
		$pre_k_left = $best_kmer_left;
		$pre_k_right = $best_kmer_right;

		my ($uv_id_left, $uv_id_right) = ('NA','NA');
		$uv_id_left = $uv_kmer{$best_kmer_left} if defined $uv_kmer{$best_kmer_left};
		$uv_id_right = $uv_kmer{$best_kmer_right} if defined $uv_kmer{$best_kmer_right};

		print "$best_kmer_left\t$best_kmer_count_left\t$best_kmer_right\t$best_kmer_count_right\t$uv_id_left\t$uv_id_right\n";
	}	
}

=head2
 convert -- convert table format (GEO database) to fasta format (default), or fasta format to table format
=cut
sub convert
{
	my $options = shift;

	my $subUsage = qq'
USAGE: $0 convert [options]
	-i	input file 
	-o	output prefix (defaul:sRNAseq)
	-p	prefix of out seqID (for table convert to fasta)
	-f	convert table to fasta (default:0) / fasta to table (1)

';

	print $subUsage and exit unless $$options{'i'}; 
	my $out_prefix = 'sRNAseq';
	$out_prefix = $$options{'o'} if $$options{'o'};


	my ($inFile, $outFile, $format, $prefix);
	$inFile = $$options{'i'};
	$prefix = $inFile; $prefix =~ s/\..*//;
	$prefix = $$options{'p'} if $$options{'p'};
	$format = 0;
	$format = 1 if $$options{'f'};

	my %read; my $num;
	if ($format) {
		$outFile = $out_prefix.".tab";
		die "[ERR]output file $outFile exist\n" if -s $outFile;
		my $out = IO::File->new(">".$outFile) || die $!;
		my $in = IO::File->new($inFile) || die $!;
		while(<$in>)
		{
			chomp;
			my $id = $_;
			my @a = split(/-/, $id);
			die "[ERR]sRNA ID: $id\n" unless @a == 2;
			die "[ERR]sRNA num: $a[1]\n" unless $a[1] > 0;
			my $rd = <$in>; chomp($rd);
			if ( defined $read{$a[0]} ) {
				$read{$rd} = $read{$rd} + $a[1];
				print "[WARN]Repeat Uniq sRNA Read $a[0]\n";
			} else {
				$read{$rd} = $a[1];
			}
		}
		$in->close;

		foreach my $r (sort keys %read)
		{
			print $out $r."\t".$read{$r}."\n";
		}
		$out->close;

	} else {
		$outFile = $out_prefix.".fasta";
		die "[ERR]output file $outFile exist\n" if -s $outFile;
		my $out = IO::File->new(">".$outFile) || die $!;
		my $in = IO::File->new($inFile) || die $!;
		while(<$in>)
		{
			chomp;
			my @a = split(/\t/, $_);
			if ( defined $read{$a[0]} ) {
				$read{$a[0]} = $read{$a[0]} + $a[1];
				print "[WARN]Repeat Uniq sRNA Read $a[0]\n"; 
			} else {
				$read{$a[0]} = $a[1];
			}
		}
		$in->close;

		foreach my $r (sort keys %read)
		{
        		$num++;
			my $count = $read{$r};
			print $out ">".$prefix."A".$num."-".$count."\n$r\n";
		}
		$out->close;
	}
	# foreach my $o (sort keys %$options) { print "$o\t$$options{$o}\n"; }
}

=head2
 unique -- convert the clean sRNA to unique (remove duplication)
=cut
sub unique
{
	my $options = shift;

	my $subUsage = qq'
USAGE: $0 unique [options]
        -i      input file 
        -u      input reads file is sRNA clean (default:0) / uniq (1) format

* convert the clean sRNA to unique
  convert the unique sRNA to clean if -u provide
* output example
  >id000001-10930   # 10930 is the number of sRNA
* support fastq file
* output read is sorted by number, could check the high expressed read

';

	print $subUsage and exit unless $$options{'i'};
	my $input_file = $$options{'i'};
	die "[ERR]cat not find input file $$options{'i'}\n" unless -s $$options{'i'};

	if ($$options{'u'}) { # convert uniq to clean(norm) format
		my $output_file = $input_file;
		$output_file =~ s/\..*$/_norm\.fasta/;
		die "[ERR]Output exist, $output_file\n" if -s $output_file;
		my $out = IO::File->new(">".$output_file) || die $!;

		my ($fh);
		if ($input_file =~ m/\.gz$/) {
			open($fh, '-|', "gzip -cd $input_file") || die $!;
		} else {
			open($fh, $input_file) || die $!;
		}

		while(<$fh>) 
		{
			chomp;
			my $id = $_;
			die "[ERR]seq format $id\n" unless $id =~ m/^>/;
			$id =~ s/^>//;
			my @a = split(/-/, $id);
			my $seq_num = pop @a;
			my $seq_id = join("-", @a);
        		my $seq = <$fh>; chomp($seq);

			for(my $i=1; $i<=$seq_num; $i++)
			{
				my $new_id = $seq_id."_$i";
				print $out ">".$new_id."\n".$seq."\n";
			}
		}

		close($fh);
		$out->close;

	} else {  # convert the clean to uniq 

		my $output_file = $input_file;
		$output_file =~ s/\..*$/_uniq\.fasta/;
		die "[ERR]Output exist, $output_file\n" if -s $output_file;
		
		my ($fh, $format, $total_seq_num, $length);
		if ($input_file =~ m/\.gz$/) {
			open($fh, '-|', "gzip -cd $input_file") || die $!;
		} else {
			open($fh, $input_file) || die $!;
		}

		# load seq/read count to hash
		my %seq_count;
		while(<$fh>)
		{
			chomp; 
			my $id = $_; $id =~ s/ .*//;
			$format = '';
			if ($id =~ m/^>/) { $format = 'fasta'; }
			elsif ($id =~ m/^@/) { $format = 'fastq'; }
			else { die "[ERR]seq fromat: $id\n"; }
		
			my $seq = <$fh>; chomp($seq);
			if ( defined $seq_count{$seq} ) {
				$seq_count{$seq}++;
			} else {
				$seq_count{$seq} = 1;
			}
			
			if ($format eq 'fastq') { <$fh>; <$fh>; }
			$total_seq_num++;
		}
		$fh->close;
		
		$length = length($total_seq_num);

		# sort by num for duplicate seq/read
		my %seq_count_sort;
		foreach my $sq (sort keys %seq_count) {
			my $count = $seq_count{$sq};
			if ($count > 1 ) {
				if (defined $seq_count_sort{$count}) {
					$seq_count_sort{$count}.= "\t".$sq;
				} else {
					$seq_count_sort{$count} = $sq;
				}
				delete $seq_count{$sq};
			} 
		}

		# output result
		my $seq_num = 0;
		my $out = IO::File->new(">".$output_file) || die $!;
		# --- output duplicate seq/read
		foreach my $ct (sort { $b<=>$a } keys %seq_count_sort) { 
			my @seq = split(/\t/, $seq_count_sort{$ct});
			foreach my $sq (@seq) {
				$seq_num++;
				my $zlen = $length - length($seq_num);
				my $z = "0"x$zlen;
				my $seq_id = "sRU".$z.$seq_num;
				print $out ">$seq_id-$ct\n$sq\n";
			}
		}	

		# --- output single seq/read
		foreach my $sq (sort keys %seq_count) {
			$seq_num++;
			my $zlen = $length - length($seq_num);
			my $z = "0"x$zlen;
			my $seq_id = "sRU".$z.$seq_num;
			print $out ">$seq_id-1\n$sq\n";
		}		

		$out->close;
	}
}

=head2
 lengthd -- get length distribution of sRNA 
=cut
sub lengthd
{
	my $options = shift;
	
	my $subUsage = qq'
USAGE: $0 lengthd [options]
        -i      input file 
        -u      input reads file is sRNA clean (default:0) / uniq (1) format

';

	print $subUsage and exit unless $$options{'i'};
	my $input_seq = $$options{'i'};
	my $key = $input_seq; $key =~ s/\..*$//;
	my $output_table = $key.".table";
	my $output_plots = $key.".pdf";
	my $output_image = $key.".png";

	my %length_dist;
	my $seq_num = 0;
	my ($seq_id_info, $seq_id, $seq_desc, $format, $sequence, $seq_length, $uniq_count);
	
	my $fh = IO::File->new($input_seq) || die $!;
	while(<$fh>)
	{
		chomp;
		$seq_id_info = $_;
		if      ($seq_id_info =~ m/^>/) { $format = 'fasta'; $seq_id_info =~ s/^>//; }
		elsif   ($seq_id_info =~ m/^@/) { $format = 'fastq'; $seq_id_info =~ s/^@//; }
		else    { die "[ERR]sRNA ID: $seq_id_info\n"; }
		($seq_id, $seq_desc) = split(/\s+/, $seq_id_info, 2);
		unless ($seq_desc) { $seq_desc = ""; }
		
		$sequence = <$fh>; chomp($sequence);
		$seq_length = length($sequence);

		if ($$options{'u'}) {
			my @nn = split(/-/, $seq_id);
			$uniq_count = $nn[scalar(@nn)-1];
			die "[ERR]sRNA count $seq_id_info, $seq_id, $uniq_count\n" if $uniq_count < 1;
			$seq_num = $seq_num + $uniq_count;

			if ( defined $length_dist{$seq_length} ) { $length_dist{$seq_length} = $length_dist{$seq_length} + $uniq_count; }
			else { $length_dist{$seq_length} = $uniq_count; }
		} else {
			$seq_num++;

			if ( defined $length_dist{$seq_length} ) { $length_dist{$seq_length}++; }
			else { $length_dist{$seq_length} = 1; }
		}

		if ($format eq 'fastq') { <$fh>; <$fh>; }
	}
	$fh->close;

	# output lengt distribution tables
	my $out = IO::File->new(">".$output_table) || die "Can not open output table file $output_table $!\n";
	foreach my $len (sort keys %length_dist) {
		my $freq = sprintf('%.4f', $length_dist{$len}/$seq_num);
		$freq = $freq * 100;
		print $out "$len\t$length_dist{$len}\t$freq\n";
	}
	$out->close;	

	# R code for length distribution
my $R_LD =<< "END";
a<-read.table("$output_table")
x<-a[,1]
y<-a[,2]
dat <- data.frame(fac = rep(x, y))
pdf("$output_plots",width=12,height=6)
barplot(table(dat)/sum(table(dat)), col="lightblue", xlab="Length(nt)", ylab="Frequency", main="Length distribution")
invisible(dev.off())
END

	open R,"|/usr/bin/R --vanilla --slave" or die $!;
	print R $R_LD;
	close R;	

	# convert pdf file to png
	my $cmd_convert = "convert $output_plots $output_image";
	system($cmd_convert) && die "[ERR]CMD: $cmd_convert\n";
}

=head2
 norm -- normalization of sRNA dataset
=cut
sub norm
{
	my $options = shift;
	
	my $subUsage = qq'
USAGE $0 norm [options] 
        -i      list of UNIQUE read
        -o      perfix of output files

* the input file in list MUST be UNIQUE format of sRNA
  the UNIQUE formart include read count in ID.

* the output files
[perfix]_sRNA_expr	raw expression
[perfix]_sRNA_libsize	library size
[perfix]_sRNA_expTPM	normalized exp
[perfix]_sRNA_seq	unique sRNA

';
	print $subUsage and exit unless ($$options{'i'} && $$options{'o'});
	
	my $list_uniq_read = $$options{'i'};
	my $prefix = $$options{'o'};
	my $output1 = $prefix."_sRNA_expr";
	my $output2 = $prefix."_sRNA_libsize";
	my $output3 = $prefix."_sRNA_expTPM";
	my $output4 = $prefix."_sRNA_seq";

	# put list of uniq small RNA reads to array
	my @list;
	my $fh = IO::File->new($list_uniq_read) || die "Can not open list file $list_uniq_read $!\n";
	while(<$fh>)
	{
		chomp;
		push(@list, $_);
		die "[ERR]cat not find uniq read file $_\n" unless -s $_;
	}
	$fh->close;

	# main expression and libsize value to hash
	my %uniq_read;
	my %libsize;	

	foreach my $file (@list)
	{
        	my $total_read = 0;
        	my $fu;
		if ($file =~ m/\.gz$/) {
			open ($fu,'-|', "gzip -cd $$file") || die $!;	# discard gzip IO, method from honghe
        	} else {
                	open($fu, $file) || die $!;
        	}

		while(<$fu>)
		{
                	chomp;
	                my $id = $_;
        	        my @a = split(/-/, $id);
	                my $exp = $a[scalar(@a)-1];
        	        my $seq = <$fu>; chomp($seq);
	                $uniq_read{$seq}{$file} = $exp;
        	        $total_read = $total_read + $exp;
	        }
	        close($fu);
	        $libsize{$file} = $total_read;
	}

	# output the libsize
	my $out2 = IO::File->new(">$output2") || die "Can not open output libsize $output2 $!\n";
	foreach my $k (sort keys %libsize) { print $out2 $k."\t".$libsize{$k}."\n"; }
	$out2->close;

	# output the expression 
	my $out1 = IO::File->new(">$output1") || die "Can not open output expression $output1 $!\n";
	my $out3 = IO::File->new(">$output3") || die "Can not open output expression TPM $output3 $!\n";
	my $out4 = IO::File->new(">$output4") || die "Can not open output small RNA $output4 $!\n";

	print $out1 "#ID\tUniqRead";
	print $out3 "#ID\tUniqRead";
	foreach my $file (@list) {
        		print $out1 "\t".$file;
			print $out3 "\t".$file;
	}
	print $out1 "\n";
	print $out3 "\n";

	my $seq_order = 0;
	my $length = length(scalar(keys(%uniq_read)));

	foreach my $seq (sort keys %uniq_read)
	{
		$seq_order++;
		my $zlen = $length - length($seq_order);
		my $z = "0"x$zlen;
		my $seq_id = "sR".$prefix.$z.$seq_order;

		print $out4 ">".$seq_id."\n".$seq."\n";

		my $line = $seq_id."\t".$seq;
		my $line_tpm = $seq_id."\t".$seq;

		foreach my $file (@list)
		{
			if ( defined $uniq_read{$seq}{$file} )
                	{
                        	$line.="\t$uniq_read{$seq}{$file}";

	                        my $tpm = $uniq_read{$seq}{$file} / $libsize{$file} * 1000000;
        	                $tpm = sprintf("%.2f", $tpm);
                	        $line_tpm.="\t$tpm";
                	}
                	else
                	{
                        	$line.="\t0";
                        	$line_tpm.="\t0";
                	}
        	}

	        print $out1 $line."\n";
	        print $out3 $line_tpm."\n";
	}
	$out1->close;
	$out3->close;
	$out4->close;

}

=head2
 normcut -- normalization of sRNA cutoff
=cut
sub normcut
{
	my $options = shift;
	
	my $subUsage = qq'
USAGE $0 normcut [options]
	-i	input file prefix
	-c	cutoff (default: 10 TPM)

* the input file should be:
[perfix]_sRNA_expr
[perfix]_sRNA_expTPM
[perfix]_sRNA_seq

';

	print $subUsage and exit unless $$options{'i'};
	my $file_prefix = $$options{'i'};
	my $cutoff = 10;
	$cutoff = $$options{'c'} if (defined $$options{'c'} && $$options{'c'} > 0);

	my $expr = $file_prefix."_sRNA_expr";
	my $norm = $file_prefix."_sRNA_expTPM";
	my $srna = $file_prefix."_sRNA_seq";

	my $out_expr = $file_prefix."_".$cutoff."TPM_sRNA_expr";
	my $out_norm = $file_prefix."_".$cutoff."TPM_sRNA_expTPM";
	my $out_sRNA = $file_prefix."_".$cutoff."TPM_sRNA_seq";

	my $out1 = IO::File->new(">".$out_expr) || die $!;
	my $out2 = IO::File->new(">".$out_norm) || die $!;
	my $out3 = IO::File->new(">".$out_sRNA) || die $!;

	my %id;
	my $fh2 = IO::File->new($norm) || die $!;
	my $title = <$fh2>; print $out2 $title;
	while(<$fh2>)
	{
        	chomp;
 		my @a = split(/\t/, $_);

		my $select = 0;
		for(my $i=2; $i<@a; $i++) {
               		if ( $a[$i] > $cutoff ) {
                        	$select = 1;
                	}
        	}

        	if ($select == 1) {
                	$id{$a[0]} = 1;
                	print $out2 $_."\n";
			print $out3 ">$a[0]\n$a[1]\n";
        	}
	}
	$fh2->close;
	$out2->close;
	$out3->close;

	my $fh1 = IO::File->new($expr) || die $!;
	my $t = <$fh1>; print $out1 $t;
	while(<$fh1>)
	{
        	chomp;
		my @a = split(/\t/, $_);
		if ( defined $id{$a[0]} ) { print $out1 $_."\n"; }
	}
	$fh1->close;
	$out1->close;
}

=head2
 combine -- combine sRNA replicate to one sample
=cut
sub combine
{
        my $options = shift;

        my $subUsage = qq'
USAGE $0 combine [options]
        -i      replicate1,replicate2,replicate3,...,repliateN
        -o      output_prefix (defaule: sRcombine)

* the input replicate should be uniq sRNA with fasta format
* the output file will be sRcombine.fasta, and sRNA ID start with sRcombine.

';

        print $subUsage and exit unless $$options{'i'};
	
	# check input files;
	my @inFiles = split(/,/, $$options{'i'});
	foreach my $f (@inFiles) {
		die "[ERR]File $f is not exist\n" unless -s $f;
	}

	# output file
	my $out_prefix = 'sRcombine';
	$out_prefix = $$options{'o'} if $$options{'o'};
	my $output_file = $out_prefix.".fasta";
	die "[ERR]output file $output_file exist\n" if -s $output_file;

	# put sRNA count to hash
	my %sRNA_count;

	foreach my $f (@inFiles) 
	{
		my $fh = IO::File->new($f) || die $!;
		while(<$fh>)
		{
			chomp;
			my $id = $_; 
			die "[ERR]seq id $id in file $f\n" unless $id =~ m/^>/;
			my @a = split(/-/, $id);
			my $ct = $a[scalar(@a)-1];
			die "[ERR]seq num $ct\n" if $ct < 1;
			my $seq = <$fh>; chomp($seq);

			if (defined $sRNA_count{$seq}) {
				$sRNA_count{$seq} = $sRNA_count{$seq} + $ct;
			} else {
				$sRNA_count{$seq} = $ct;
			}
		}
		$fh->close;
	}

	my $length = length(scalar(keys(%sRNA_count)));

	my $out = IO::File->new(">".$output_file) || die $!;
	my $s_num = 0;
	foreach my $sRNA (sort keys %sRNA_count)
	{
		my $ct = $sRNA_count{$sRNA};
		$s_num++;
		my $zlen = $length - length($s_num);
		my $z = "0"x$zlen;
		my $sid = $out_prefix.$z.$s_num."-".$ct;
		print $out ">$sid\n$sRNA\n";
	}
	$out->close;
}

=head2
 usage -- print usage information
=cut
sub usage
{
	print qq'
Program: sRNAtools (Tools for sRNA analysis)
Version: $version

USAGE: $0 <command> [options] 
Command: 
	chkadp		check adapter sequence (kmer method, better for unknown adapter)
	chkadp2		check adapter sequence (regexp method, better for known adapter)
	rmadp		remove sRNA adapter sequence
	convert		convert between table format and fastq/fasta format
	unique		convert between unique format and clean format
	norm	     	normalization (RPM)
	normcut		normalization cutoff	
	lengthd		length distribution of sRNA	
	combine		combine sRNA replicates
	
';
	exit;
}

=head2
 pipeline -- print pipeline for sRNA analysis
=cut
sub pipeline
{
	print qq'
A: sRNA clean, collesped, length distribution

B: microRNA identification

C: virus identification

D: Phasing RNA identification

E: mimic pairing

';

	exit;
}
