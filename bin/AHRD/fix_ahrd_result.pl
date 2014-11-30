#!/usr/bin/perl

=head2

 fix_ahrd_result.pl -- fix ahrd results with no hit 

=cut

use strict;
use warnings;
use Bio::SearchIO;

if (@ARGV < 2) {
	print qq'
USAGE: $0 input_ahrd.csv input_blast1 input_blast2 ... input_blastN

';
	exit;
}

my $evalue = 1e-5;

my $ahrd_result = $ARGV[0];

my @blast = @ARGV; shift @blast;

# key: protein/transcript ID
# value: hit

my %id_hit;
my $id_hit = \%id_hit;

foreach my $blast (@blast) 
{
	print "Parsing blast file: $blast ....\n";
	$id_hit = parse_blast($blast, $id_hit);
}

# parse ahrd
my $fh = IO::File->new($ahrd_result) || die $!;
for (1 .. 3) {
	my $line = <$fh>;
	print $line;
}

while(<$fh>)
{
	chomp;
	# Protein-Accession	Blast-Hit-Accession	AHRD-Quality-Code	Human-Readable-Description	Interpro-ID (Description)	Gene-Ontology-ID (Name)
	my @a = split(/\t/, $_);
	my $id = $a[0];
	if (defined $$id_hit{$id}) {
		
	} else {
		$a[3] = "No hit";
	}

	print join("\t", @a),"\n";
}
$fh->close;

print "done\n";

# kentnf: subroutine
sub parse_blast
{
	my ($blast, $id_hit) = @_;
	
	my ($hit_num, $query_name);

	my $in = Bio::SearchIO->new(-format=>'blast', -file=>$blast);
	while(my $result = $in->next_result)
	{
		$hit_num = 0;
		$query_name = $result->query_name;

		while(my $hit = $result->next_hit)
        	{
                	if (($query_name ne $hit->name) && $hit->significance < 1e-5  )
                	{
                        	$hit_num++;
			}
		}

		if ($hit_num > 0)
        	{
                	$$id_hit{$query_name} = 1;
        	}
	}

	return $id_hit;
}

