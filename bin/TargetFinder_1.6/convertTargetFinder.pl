#!/usr/bin/perl

=head1

 convertTargetFinder.pl -- convert TargetFinder.pl output format to our sRNA target output format

 TargetFinder.pl output format:
 query=sRNA3696, target=GSVIVT01013310001, score=0, range=797-817, strand=1

 target  5' GGAGAGCCUUACUGGAAAAAA 3'
            :::::::::::::::::::::
 query   3' CCUCUCGGAAUGACCUUUUUU 5'

 
 our sRNA target output format:
 AA000005-2      AT3G21295.1     1068    1090    3       UACGCAGGAAAAUGGUUUUAAAA AUCCGUUCUUUUAUGAAAAUUUU ||b|||o||||||ob|||||||| +

=cut
use strict;
use warnings;
use IO::File;

my $usage = qq'
perl convertTargetFinder.pl  output_folder_TargetFinder  cutoff_align_score  

';

my $input = shift || die $usage;
my $cutoff_align_score = shift || 4.5;

# check the input folder 
die "Error, can not find folder $input\n" unless -e $input;

# get output prefix 
my $output_prefix = $input; 
$output_prefix =~ s/_TargetAlign//;
$output_prefix.= "_targetfinder.txt";

# get alignment file list
my $align_list = `ls $input`; chomp $align_list;
my @align_file = split(/\n/, $align_list);

my $ofh = IO::File->new(">".$output_prefix) || die $!;

foreach my $file (@align_file)
{
	my $real_file = $input."/".$file;
	die "Error, can not find alignment file $real_file\n" unless -e $real_file;
	#print "$real_file, $output_prefix, $cutoff_align_score\n";

	my $line = convert_format($real_file, $output_prefix, $cutoff_align_score);

	if ($line)
	{
		print $ofh $line;
	}
}
$ofh->close;

#################################################################
# kentnf: subroutine						#
#################################################################

sub convert_format
{
	my ($align_file, $output_prefix, $cutoff_align_score) = @_;

	my $output_line = "";

	my $fh = IO::File->new($align_file) || die "Can not open alignment file $align_file $!\n";
	while(<$fh>)
	{
		my $info_line = $_; chomp($info_line);
		
		if ($info_line =~ m/No results for/) { last; }

		my @a = split(/, /, $info_line);

		my ($query, $target, $score, $range, $strand, $start, $end);
		if (scalar(@a) == 5)
		{
			($query, $target, $score, $range, $strand) = @a;

			$query =~ s/query=//;
			$target =~ s/target=//;
			$score =~ s/score=//;
			$range =~ s/range=//;
			$range =~ s/-/\t/;
			$strand =~ s/strand=//;
			if ($strand == 1) { $strand = "+"; } else { $strand = "-"; }
			
		}
		else
		{
			die "Error in info line $info_line\n";
		}

		<$fh>; # skip line
		my $target_seq = <$fh>; chomp($target_seq); my @mm = split(/\s+/, $target_seq);
		my $align_seq = <$fh>;  chomp($align_seq);
		my $query_seq = <$fh>;  chomp($query_seq);  my @nn = split(/\s+/, $query_seq);
		<$fh>; # skip line

		$target_seq = $mm[2]; 
		my $len = length($target_seq);
		$align_seq = substr($align_seq, 11, $len);
		$query_seq = $nn[2];
		$align_seq =~ s/ /b/ig;
		$align_seq =~ s/:/|/ig;
		$align_seq =~ s/\./o/ig;

		if ($score <= $cutoff_align_score)
		{
			$output_line.= "$query\t$target\t$range\t$score\t$target_seq\t$query_seq\t$align_seq\t$strand\n";
		}
	}
	$fh->close;

	return $output_line;
}
