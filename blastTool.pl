#!/usr/bin/perl

=head
 blastTool -- parse blast result
=cut

use strict;
use warnings;
use FindBin;
use IO::File;
use Bio::SearchIO;
use Getopt::Std;

my $version = 0.1;
my $debug = 0;
my %options;
getopts('a:b:c:d:e:f:g:i:j:k:l:m:n:o:p:q:r:s:t:u:v:w:x:y:z:h', \%options);
unless (defined $options{'t'} ) { usage($version); }

if	($options{'t'} eq 'b2table') 	{ blast_to_table(\%options, \@ARGV); } 	# parse multi dataset
elsif	($options{'t'} eq 'brh') 	{ blast_brh(\%options, \@ARGV); } 	# parse multi dataset
else	{ usage($version); }

#================================================================
# kentnf: subroutine						
#================================================================
sub blast_brh
{
	my ($options, $files) = @_;
	my $usage = qq'
USAGE $0 -t brh input_blast_table > output

* convert blast to table first

';
	print $usage and exit unless defined $$files[0];
	my $input_file = $$files[0];
	die "[ERR]input not exist\n" unless -s $input_file;

	my %id_best;
	my $in = IO::File->new($input_file) || die $!;
	while(<$in>)	
	{
		chomp;
		next if $_ =~ m/^#/;
		my @a = split(/\t/, $_);
		# query_name \\t hit_name \\t hit_description \\t score \\t significance
		if (defined $id_best{$a[0]}{'id'}) {
			if ($a[4] < $id_best{$a[0]}{'e'}) {
				$id_best{$a[0]}{'id'} = $a[1];
				$id_best{$a[0]}{'e'} = $a[4];
			}
		} else {
			$id_best{$a[0]}{'id'} = $a[1];
			$id_best{$a[0]}{'e'} = $a[4];
		}

		if (defined $id_best{$a[1]}{'id'}) {
			if ($a[4] < $id_best{$a[1]}{'e'}) {
				$id_best{$a[1]}{'id'} = $a[0];
				$id_best{$a[1]}{'e'} = $a[4];
			}
		} else {
			$id_best{$a[1]} = $a[0];
			$id_best{$a[1]}{'e'} = $a[4];
		}
	}
	$in->close;

	my %uniq_best;
	foreach my $id (sort keys %id_best) {
		my $id1 = $id_best{$id}{'id'};
		my $id2 = $id_best{$id1}{'id'};
		
		if ($id2 eq $id) 
		{
			if (defined $uniq_best{$id2."#".$id1} || defined $uniq_best{$id1."#".$id2})  {
				next;
			} else {
				print "$id2\t$id1\n";
				$uniq_best{$id2."#".$id1} = 1;
				$uniq_best{$id1."#".$id2} = 1;
			}
		}
	}
}
=head2
 blast_to_table : convert blast result to table
=cut
sub blast_to_table
{
	my ($options, $files) = @_;

	my $usage = qq'
USAGE: $0 -t b2table [options] blast_result > output_table
	
	-n max_hit_num[default:5]

* the output is the top 5 blast hits for each blast results
* format of output file
  query_name \\t hit_name \\t hit_description \\t score \\t significance 

';
	print $usage and exit unless defined $$files[0];
	my $input_file = $$files[0];
	die "[ERR]input not exist\n" unless -s $input_file;

	my $max_hit_num = 5;
	$max_hit_num = $$options{'n'} if (defined $$options{'n'} && $$options{'n'} > 0);

	my $evalue = 1e-5;
	$evalue = $$options{'e'} if (defined $$options{'e'} && $$options{'e'} > 0);

	my ($hit_num, $query_name, $query_length);
	my $in = Bio::SearchIO->new(-format=>'blast', -file=>$input_file);
	while(my $result = $in->next_result)
	{
		$hit_num = 0;
		$query_name = $result->query_name;
		$query_length = $result->query_length;

		while(my $hit = $result->next_hit)
		{
			if (($query_name ne $hit->name) && $hit_num < $max_hit_num && $hit->significance < $evalue)
			{
				$hit_num++;
				print 	$query_name."\t".
					$hit->accession."\t".
					$hit->description."\t".
					$hit->raw_score."\t".
					$hit->significance."\n";
			}
		}

		if ($hit_num == 0)
		{
			print $query_name."\tNo hit\n";
		}
	}
}

=head2
 usage: print usage information
=cut
sub usage
{
	my $version = shift;
	my $usage = qq'
USAGE: $0 -t [tools]

	b2table		convert blast result to table
	brh		best reciprocal hit

';
	print $usage; 
	exit;
}
