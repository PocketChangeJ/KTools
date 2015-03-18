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
else	{ usage(); }



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
 usage -- print usage information
=cut
sub usage
{
	my $usage = qq'
USAGE: $0 -t [tool]

	formatCCS	format ccs for next analysis

';

	print $usage;
	exit;
}



