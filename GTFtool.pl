#!/usr/bin/perl

=head1
 GTFtool -- parse GTF file
=cut
use strict;
use warnings;
use IO::File;
use Getopt::Std;

my $version = 0.1;

=head2
 convert 
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
	
	my $fh = IO::File->new($input_gtf) || die $!;
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

			


		}
	}
	$fh->close;

	return %gtf_info;
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
        convert		convert GTF to BED/GFF/TAB format
	import		import GFF to GTF (gffread)
	extractList	extract GTF by list
	stat		statistics for GTF file

* the gtf file must have exon feature, transcript_id, and gene_id attributes

';
        exit;
}

