#!/usr/bin/perl

=head
 blastTool -- parse blast result
=cut

use strict;
use warnings;
use FindBin;
use IO::File;
use Bio::SeqIO;
use Bio::SearchIO;
use Getopt::Std;

my $version = 0.1;
my $debug = 0;
my %options;
getopts('a:b:c:d:e:f:g:i:j:k:l:m:n:o:p:q:r:s:t:u:v:w:x:y:z:h', \%options);
unless (defined $options{'t'} ) { usage($version); }

if	($options{'t'} eq 'b2table') 	{ blast_to_table(\%options, \@ARGV); }
elsif	($options{'t'} eq 'brh') 	{ blast_brh(\%options, \@ARGV); }
elsif	($options{'t'} eq 'unique')	{ blast_unique(\%options, \@ARGV); }
else	{ usage($version); }

#================================================================
# kentnf: subroutine						
#================================================================
=head
 blast_unique: remove redundancy sequence
=cut
sub blast_unique
{
	my ($options, $files) = @_;
	my $usage = qq'
USAGE $0 -t unique input_seq blast_m8 > output

* output file
	1. input_seq.unique
	2. input_seq.removed
	3. input_seq.rmDup_table.txt

* genrate unique seqs according self-blast m8 output
* only works for %100 identify and %100 alignment



';
	print $usage and exit unless (scalar @$files == 2);
	my ($input_seq, $blast_m8) = @$files;

	# put seq_info to hash
	my %seq_info;
	my $in = Bio::SeqIO->new(-format=>'fasta', -file=>$input_seq);
	while(my $inseq = $in->next_seq)
	{
		$seq_info{$inseq->id}{'seq'} = $inseq->seq;
		$seq_info{$inseq->id}{'len'} = $inseq->length;
	}

	# put removed seq to hash
	my %removed;
	my $fh = IO::File->new($blast_m8) || die $!;
	while(<$fh>)
	{
		chomp;
		next if $_ =~ m/^#/;
		my @a = split(/\t/, $_);
		# Qid, Sid, %iden, alignLen, mismatches, gap open, q.start, q.end, s.start, s.end, evalue, score
		die "[ERR]line:$_\n" unless @a == 12;
		my ($qid, $sid, $iden, $alen, $mismatch, $gap, $q_start, $q_end, $s_start, $s_end, $evalue, $score) = @a;
		my $qlen = $seq_info{$qid}{'len'};
		my $slen = $seq_info{$sid}{'len'};
		next if $qid eq $sid;
		next if ($iden < 100 || $mismatch > 0 || $gap > 0);
		next if ($alen < $qlen && $alen < $slen);

		if  ( $alen == $qlen && $q_start == 1 && $q_end == $qlen) {
			$removed{$qid} = $sid;
		}

		if  ( $alen == $slen && $s_start == 1 && $s_end == $slen ) {
                        $removed{$sid} = $qid;
                }

		if  ( $alen == $slen && $s_start == $slen && $s_end == 1 ) {
			$removed{$sid} = $qid;
		}
	}
	$fh->close;
	
	# output unique 
	my $out1 = IO::File->new(">".$input_seq.".removed") || die $!;
	my $out2 = IO::File->new(">".$input_seq.".unique") || die $!;
	my $out3 = IO::File->new(">".$input_seq.".rmDup_table.txt") || die $!;

	foreach my $id (sort keys %seq_info) {
		if (defined $removed{$id} ) {
			print $out1 ">".$id."\n".$seq_info{$id}{'seq'}."\n";
			print $out3 $id."\t".$removed{$id}."\n";
		} else {
			print $out2 ">".$id."\n".$seq_info{$id}{'seq'}."\n";
		}
	}

	$out1->close;
	$out2->close;
	$out3->close;
}

=head
 blast_brh: find best reciprocal hits in blast result
=cut
sub blast_brh
{
	my ($options, $files) = @_;
	my $usage = qq'
USAGE $0 -t brh input_blast_table > output

* convert blast to table first, then perform this analysis

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
				$id_best{$a[0]}{'id'}= $a[1];
				$id_best{$a[0]}{'e'} = $a[4];
				$id_best{$a[0]}{'s'} = $a[3];
			} elsif ($a[4] == $id_best{$a[0]}{'e'}) {
				if ($a[3] > $id_best{$a[0]}{'s'} ) {
					$id_best{$a[0]}{'id'}= $a[1];
					$id_best{$a[0]}{'e'} = $a[4];
					$id_best{$a[0]}{'s'} = $a[3];
				}
			}
		} else {
			$id_best{$a[0]}{'id'}= $a[1];
			$id_best{$a[0]}{'e'} = $a[4];
			$id_best{$a[0]}{'s'} = $a[3];
		}

		if (defined $id_best{$a[1]}{'id'}) {
			if ($a[4] < $id_best{$a[1]}{'e'}) {
				$id_best{$a[1]}{'id'}= $a[0];
				$id_best{$a[1]}{'e'} = $a[4];
				$id_best{$a[1]}{'s'} = $a[3];
			} elsif ($a[4] == $id_best{$a[1]}{'e'}) {
				if ($a[3] > $id_best{$a[1]}{'s'} ) {
					$id_best{$a[1]}{'id'}= $a[0];
					$id_best{$a[1]}{'e'} = $a[4];
					$id_best{$a[1]}{'s'} = $a[3];
				}
			}
		} else {
			$id_best{$a[1]}{'id'}= $a[0];
			$id_best{$a[1]}{'e'} = $a[4];
			$id_best{$a[1]}{'s'} = $a[3];
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
	unique		remove redundancy using self-blast m8 result

';
	print $usage; 
	exit;
}
