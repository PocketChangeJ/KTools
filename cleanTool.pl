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
elsif	($options{'t'} eq 'barcode')	{ clean_barcode(\%options, \@ARGV); }
elsif   ($options{'t'} eq 'pipeline')	{ pipeline(); }
else	{ usage($version); }
#################################################################
# kentnf: subroutine						#
#################################################################


=head2
 barcode -- split reads according to barcode, remove barcode
=cut
sub clean_barcode
{
	my ($options, $files) = @_;
	my $usage = qq'
USAGE $0 -t barcode -m mismatch(default:0) barcode_file input_reads

* format of barcode file
output_file_name1 [tab] bacode1
output_file_name2 [tab] bacode2
......

';
	# check input files
	print $usage and exit if (scalar(@$files) < 2);
	my ($barcode_file, $input_file) = ($$files[0], $$files[1]);
	print "[ERR]no file $barcode_file\n" and exit unless -s $barcode_file;
	print "[ERR]no file $input_file\n" and exit unless -s $input_file;

	my $cutoff = 0;
	$cutoff = $$options{'m'} if (defined $$options{'m'} && $$options{'m'} < 3);
	
	# define output file
	my $pre_name = $input_file;
	$pre_name =~ s/\.txt//g; $pre_name =~ s/\.gz//g;
	my $unmatch_file = $pre_name ."_unmatch.txt";
	my $ambiguous_file = $pre_name . "_ambiguous.txt";
	my $out1 = IO::File->new(">$unmatch_file") || die "Can't open unmatch output file\n";
	my $out2 = IO::File->new(">$ambiguous_file") || die "Can't open ambiguous output file\n";

        my $i = 0;
	my %tag_info;
	my @barcode;
	my $OUT;
	my $in1 = IO::File->new($barcode_file) || die "Can't open the barcode file\n";
	while(<$in1>)
	{
        	chomp;
		next if $_ =~ m/^#/;
		my @a = split(/\t/);
		$tag_info{$a[1]}{'file'} = $a[0];
		$tag_info{$a[1]}{'num'} = 0;
		push(@barcode, $a[1]);
		
		# Create file for each tag
		$i++;
		#my $out = $OUT.$i;
		my $out;
        	open($out, ">$a[0]") || die "can't create $a[0] $!\n";
		$tag_info{$a[1]}{'fh'} = $out;
	}
	$in1->close;

	# parse raw fastaq file
	my ($in2, $id1, $id2, $seq, $qul, $modified_seq, $modified_qul);
	my $tag_result = '';
	my ($unmatch_num, $ambiguous_num) = (0, 0);
	if ($input_file =~ m/\.gz$/) {
		open($in2, "gunzip -c $input_file |") || die "Can't open the fastq file\n";
	} else {
		open($in2, $input_file) || die "Can't open the fastq file\n";
	}

	while(<$in2>) 
	{
		$id1 = $_;      chomp($id1);
		$seq = <$in2>;   chomp($seq);
		$id2 = <$in2>;   chomp($id2);
		$qul = <$in2>;   chomp($qul);

		# compare seq with barcode
		foreach my $elements ( @barcode ) {
        		my $len = length($elements);
			next if length($seq) < $len;
			my $read_substr = substr($seq, 0, $len);

			# old method 
			if ($cutoff == 0) {
				my $m = 0;
                		for( my $i = 0; $i< $len; $i++ ) {
                        		if(substr($elements, $i, 1) eq substr($read_substr, $i, 1) ) {
						$m++;
                        		}
				}
				$tag_result.= "%".$elements if $m == $len;
			} 
			# new method, haming distance
			else 
			{
				my $mismatch = hamming($read_substr, $elements);
				$tag_result.= "%".$elements if $mismatch <= $cutoff;
			}
		}

		if($tag_result eq "") { 				# Print unmatch seq.
			print $out1 $id1."\n".$seq."\n".$id2."\n".$qul."\n";
			$unmatch_num++;
                } elsif($tag_result =~ /%[ACGTN]*%[ACGTN]*$/) {		# Print ambiguous seq.
                        print $out2 $id1."$tag_result\n".$seq."\n".$id2."\n".$qul."\n";
			$ambiguous_num++;
                } else {						# Print for other non-ambiguous tags
                        $tag_result =~ s/%//g;
                        my $fileH = $tag_info{$tag_result}{'fh'};
			$tag_info{$tag_result}{'num'}++;
                        $modified_seq = substr($seq,length($tag_result), );
                        $modified_qul = substr($qul,length($tag_result), );
                        print $fileH $id1."\n".$modified_seq."\n".$id2."\n".$modified_qul."\n";
                }
                $tag_result = "";
	}
	$out1->close;
	$out2->close;

	foreach my $elements (@barcode) {
        	my $fileH = $tag_info{$elements}{'fh'};
		close($fileH);
		print $tag_info{$elements}{'file'},"\t",$tag_info{$elements}{'num'},"\n";
	}

	print $ambiguous_file."\t".$ambiguous_num."\n";
	print $unmatch_file."\t".$unmatch_num."\n";
}

sub hamming($$) { length( $_[ 0 ] ) - ( ( $_[ 0 ] ^ $_[ 1 ] ) =~ tr[\0][\0] ) }

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

	trimo		trim adapter, low quality, and short reads using trimmomatic.	
	barcode		remove barcode and split file according to barcode

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
>> A. remove barcode
	\$ $0 -t barcode barcode_file input.fastq.gz

>> 1. clean mRNA
      	1.1 clean adapter, low quality, short reads
	\$ $0 -t trimo sample1_R1.fastq ... sampleN_R1.fastq | sample1_R1.fastq,sample1_R2.fastq ... sampleN_R1.fastq,sampleN_R2.fastq
	1.2 clean rRNA or other contanmination
	\$ $0 -t align -r reference -m 0 -v 0 -k 1 *.fastq
	
>> 2. clean sRNA/degradome
	2.1 clean adapter, short reads,
	2.2 clean rRNA or other contanmination

>> 3. clean DNA
	3.1 clean adapter, low quality, short reads

';
	exit;
}
