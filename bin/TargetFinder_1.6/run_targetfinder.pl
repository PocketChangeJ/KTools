#!/usr/bin/perl

use strict;
use warnings;
use FindBin;

print "USAGE: $0 input_sRNA mRNA score\n" and exit if @ARGV < 2;

my $sRNA = $ARGV[0];
my $mRNA = $ARGV[1];
my $score = 4.5;
$score = $ARGV[2] if defined $ARGV[2];

target_finder($sRNA, $mRNA, $score);

sub target_finder
{
	my ($sRNA, $mRNA, $score) = @_;

	foreach my $f (($sRNA, $mRNA)) { 
		die "[ERR]file not exit, $f\n" unless -s $f;
	}

	my $targetfinder_bin = $FindBin::RealBin."/targetfinder.pl";
	my $fdr_targetfinder = $sRNA."_TargetAlign";
	my $tmp_fdr = $sRNA."_temp";

	if (-e $fdr_targetfinder ) { print "Folder $fdr_targetfinder exist! May contain old results !\n"; exit; }
	else { mkdir($fdr_targetfinder); }

	my $sq = IO::File->new($sRNA) || die "Can not open sRNA file $sRNA $!\n";
	while(<$sq>)
	{
		chomp;
		my $id = $_; $id =~ s/ .*//ig; $id =~ s/>//ig;
		my $seq = <$sq>; chomp($seq);
		my $align_file = $fdr_targetfinder."/".$id.".align";
		my $cmd_tf = "perl $targetfinder_bin $tmp_fdr -s $seq -d $mRNA -c $score -q $id > $align_file";
		print ($cmd_tf);
		system($cmd_tf) && die "[ERR]CMD $cmd_tf\n";
	}
	$sq->close;
}

