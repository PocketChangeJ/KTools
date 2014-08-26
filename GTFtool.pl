#!/usr/bin/perl

=head1
 GTFtool -- parse GTF file
=cut
use strict;
use warnings;
use IO::File;
use Getopt::Std;

my $version = 0.1;
if (@ARGV < 1) { usage($version);}
my $tool = shift @ARGV;

my %options;
getopts('i:o:p:fuh', \%options);

if      ($tool eq 'stats')	{ stats(\%options); }
elsif   ($tool eq 'convert') { convert(\%options); }
else    { usage($version); }


#################################################################
# kentnf: subroutine						#
#################################################################
=head2
 stat -- generate statistics information
=cut
sub stats
{
	my $options = shift;
	my $subUsage = qq'
USAGE: $0 stats [options]
	-i      input file 

';
	print $subUsage and exit unless $$options{'i'};
	my ($inFile, $out_prefix);
	$inFile = $$options{'i'};

	my %trans_info = parse_gtf($$options{'i'});

	my %gene;
	my %trans;
	my %exon;
	my $chr;
	foreach my $tid (sort keys %trans_info)
	{
		$trans{$tid} = 1;
		$gene{$trans_info{$tid}{'gid'}} = 1;
		$chr = $trans_info{$tid}{'chr'};

		my @exon = split(/\t/, $trans_info{$tid}{'exon'});
		@exon = sort {$a<=>$b} @exon;

		for(my $i=0; $i<scalar(@exon)-1; $i=$i+2)
		{
			my $start = $exon[$i];
			my $end = $exon[$i+1];
			my $key = "$chr\t$start\t$end";
			$exon{$key} = 1;
		}
	}

	print $inFile,"\tGene:",scalar(keys(%gene)),"\tTrans:",scalar(keys(%trans)),"\tExon:",scalar(keys(%exon)),"\n";

}

=head2
 convert -- convert GTF to GFF,BED,TAB format
=cut
sub convert
{




}

=head2
 parse_gtf -- parse gtf file, return gtf information
=cut
sub parse_gtf
{
	my $input_file = shift;

	my %trans_info; # key: tid, chr, exon, gene, strand
	
	my $fh = IO::File->new($input_file) || die $!;
	while(<$fh>)
	{
		chomp;
		next if $_ =~ m/^#/;
		my @a = split(/\t/, $_);
		if ($a[2] eq 'exon') 
		{
			# analyzing the attribute 
			my %attr_value;
			my @b = split(/; /, $a[8]); 
			foreach my $b (@b) {
				my @c = split(/ /, $b);
				die "[ERR]attr $b in $a[8]\n" unless @c == 2;
				$c[1] =~ s/"//;
				$attr_value{$c[0]} = $c[1];
			}

			die "[ERR]No transcript id\n" unless defined $attr_value{'transcript_id'};
			die "[ERR]No gene id\n" unless defined $attr_value{'gene_id'};
			my ($tid, $gid) = ($attr_value{'transcript_id'}, $attr_value{'gene_id'});

			if ( defined $trans_info{$tid}{'chr'} ) {
				die "[ERR]inconsistency chr for $tid\n"	if $trans_info{$tid}{'chr'} ne $a[0];
			} else {
				$trans_info{$tid}{'chr'} = $a[0];
			}

			if ( defined $trans_info{$tid}{'gid'} ) {
				die "[ERR]inconsistency gid for $tid\n" if $trans_info{$tid}{'gid'} ne $gid;
			} else {
				$trans_info{$tid}{'gid'} = $gid;
			}

			if ( defined $trans_info{$tid}{'strand'} ) {
				die "[ERR]inconsistency strand for $tid\n" if $trans_info{$tid}{'strand'} ne $a[6];
			} else {
				$trans_info{$tid}{'strand'} = $a[6];
			}

			if ( defined $trans_info{$tid}{'exon'}) {
				$trans_info{$tid}{'exon'}.="\t".$a[3]."\t".$a[4];
			} else {
				$trans_info{$tid}{'exon'} = $a[3]."\t".$a[4];
			}
		}
	}
	$fh->close;

	return %trans_info;
}

=head2
 usage -- print usage information
=cut
sub usage
{
        print qq'

Program: GTFtools (Tools for GTF file)
Version: $version

USAGE: $0 <command> [options] 
Command:
	stats		statistics for GTF file 
        convert		convert GTF to BED/GFF/TAB format
	import		import GFF to GTF (gffread)
	extractList	extract GTF by list

* the gtf file must have exon feature, transcript_id, and gene_id attributes

';
        exit;
}

