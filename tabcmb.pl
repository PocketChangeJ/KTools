#!/usr/bin/perl

use strict;
use warnings;
use IO::File;
use Getopt::Std;


my %options;
getopts('a:b:c:d:e:f:g:i:j:k:l:m:n:o:p:q:r:s:t:u:v:w:x:y:z:h', \%options);
unless (defined $options{'t'} ) { usage(); }


if 	( $options{'t'} eq 'reorder' )		{ tab_reorder(\%options, \@ARGV); }
if      ( $options{'t'} eq 'rotate' )          	{ tab_rotate(\%options, \@ARGV); }
else	{ usage(); }

sub tab_rotate
{
	my ($options, $files) = @_;
	my $usage = qq'
USAGE: $0 -t rotate input > output

';
	print $usage and exit unless defined $$files[0];
	my $input = $$files[0];
	die "[ERR]file not exist $input" unless -s $input;
	my @array;
	my $m = 0;
	my ($n, $max_n); $max_n = 0;
	my $fh = IO::File->new($input) || die $!;
	while(<$fh>)
	{
		chomp;
		my @a = split(/\t/, $_);
		$max_n = scalar(@a) if (scalar(@a) > $max_n);
		for(my $n=0; $n<@a; $n++) {
			$array[$n][$m] = $a[$n];
		}
		$m++;
	}
	$fh->close;

	for(my $i=0; $i<$max_n; $i++)
	{
		my $line = '';
		for(my $j=0; $j<$m; $j++)
		{
			$line.="\t".$array[$i][$j];
		}
		$line =~ s/^\t//;
		print $line."\n";
	}
}

sub tab_reorder
{
	my ($options, $files) = @_;
	my $usage = qq'
USAGE: $0 -t reorder  list tab

* list have key, tab will be order by the key of list 

';

	print $usage and exit unless @$files == 2;
	
	my $tab_col = 1;
	my $key_col = 1;

	my %tab;

	open(FH1, $$files[1]) || die $!;
	while(<FH1>)
	{
		chomp;
		next if $_ =~ m/^#/;
		my @a = split(/\t/, $_);
		$tab{$a[$tab_col-1]} = $_;
	}
	close(FH1);

	open(FH2, $$files[0]) || die $!;
	while(<FH2>)
	{
		chomp;
		next if $_ =~ m/^#/;
		my @a = split(/\t/, $_);
		my $k = $a[$key_col-1];
		
		if (defined $tab{$k}) {
			print $_."\t".$tab{$k}."\n";
		} else {
			print $_."\tNA\n";
		}
	}
	close(FH2);
}

sub usage
{
	my $usage = qq'
USAGE: $0 -t [options]
	reorder 

';
	print $usage and exit;
}

