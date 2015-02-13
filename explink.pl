#!/usr/bin/perl

=head2
 explink.pl -- link expression tables
=cut

use strict;
use warnings;
use IO::File;
use Getopt::Std;

my $version = 0.1;
my $debug = 0;
my %options;
getopts('a:b:c:d:e:f:g:i:j:k:l:m:n:o:p:q:r:s:t:u:v:w:x:y:z:h', \%options);
unless (defined $options{'t'} ) { usage($version); }

if	($options{'t'} eq 'sRNA-mRNA')	{ elink_sRNA_mRNA(\%options, \@ARGV); } # parse multi dataset
else	{ usage($version); }

=head2
 elink_sRNA_mRNA -- link sRNA and mRNA result 
=cut
sub elink_sRNA_mRNA
{
	my ($options, $files) = @_;

	my $usage = qq'
USAGE: $0 -t sRNA-mRNA [options] sRNA_target_result sRNA_exp mRNA_exp

	-f   format of sRNA_target_result (default: fei)		

';
	print $usage and exit unless scalar @$files == 3;	
	my ($target_file, $sRNA_exp_file, $mRNA_exp_file) = @$files;

	foreach my $f (($target_file, $sRNA_exp_file, $mRNA_exp_file)) {
		die "[ERR]file not exist $f\n" unless -s $f;
	}

	# different sRNA target format will use different key to connect dataset
	my $format = 'fei'; my ($key_a, $key_b) = (7, 2);
	$format = $$options{'f'} if defined $$options{'f'};

	if ($format eq 'fei') {

	} elsif ($format eq 'sparta') {

	} elsif ($format eq 'targetfinder') {

	} else {
		die "[ERR]target format $format\n";
	}

	# using different column link tables for different format

	
	# load mRNA_exp to hash
	my %mRNA_exp = load_exp_table($mRNA_exp_file, 1);	

	# load sRNA_exp to hash
	my %sRNA_exp = load_exp_table($sRNA_exp_file, 1);

	# connect mRNA and sRNA exp by target result		
	my $fh = IO::File->new($target_file) || $!;	
	my $title = <$fh>; chomp($title);
	print $title."\t".$sRNA_exp{'title'}."\t".$mRNA_exp{'title'}."\n";
	while(<$fh>) 
	{
		chomp;
		my @a = split(/\t/, $_);
		die "[ERR]less column $_\n" unless (scalar @a >= max($key_a, $key_b)-1);
		my $sRNA_key = $a[$key_a-1];
		my $mRNA_key = $a[$key_b-1];
	
		# convert sRNA_key 
		if ($format eq 'fei') {
			$sRNA_key = reverse($sRNA_key);
			$sRNA_key =~ s/-//g;
			$sRNA_key =~ s/U/T/g;
		}

		# connect
		die "[ERR]undef $sRNA_key\n" unless defined $sRNA_exp{$sRNA_key};
		die "[ERR]undef $mRNA_key\n" unless defined $mRNA_exp{$mRNA_key};

		print "$_\t$sRNA_exp{$sRNA_key}\t$mRNA_exp{$mRNA_key}\n";
	}
	$fh->close;
}

# load expression table to hash
sub load_exp_table
{
	my ($file, $col) = @_;
	my %exp_table;
	my $fh = IO::File->new($file) || die $!;
	my $title = <$fh>; chomp($title);
	$exp_table{'title'} = $title;
	while(<$fh>) {
		chomp;
		my @a = split(/\t/, $_);
		$exp_table{$a[$col-1]} = $_;
	}
	$fh->close;
	return %exp_table;
}

# get max number 
sub max
{
	my @m = @_;
	my $max = 0;
	foreach my $m (@m) {
		$max = $m if $m > $max;
	}
	return $max;
}

=head2
  usage: print usage info
=cut
sub usage
{
	my $version = shift;
	print qq'
USAGE: $0 -t [tool] [options] input files

	sRNA-mRNA	link expression of sRNA and mRNA through target analysis

';
	exit;
}


