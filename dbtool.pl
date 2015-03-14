#!/usr/bin/perl

=head
 dbtool.pl -- parse database for in house use 
=cut
use strict;
use warnings;
use IO::File;
use Getopt::Std;

my %options;
getopts('a:b:c:d:e:f:g:i:j:k:l:m:n:o:p:q:r:s:t:u:v:w:x:y:z:h', \%options);
unless (defined $options{'t'} ) { usage(); }

if	($options{'t'} eq 'uniprot2ec')		{ uniport2id(\@ARGV); }
elsif	($options{'t'} eq 'uniprot2go')		{ uniport2go(\@ARGV); }
else	{ usage(); }

#################################################################
# kentnf: subroutine						#
#################################################################

=head2
 uniprot2go -- convert uniport to GO
=cut
sub uniprot2go 
{
	my $input = shift;

}

=head2
 uniport2id -- convert uniport to fasta 
=cut
sub uniport2id
{
	my ($input) = @_;

	my $usage = qq'
USAGE: $0 -t uniprot2ec input_database(dat format)

';
	print $usage and exit unless defined $input;
	die "[ERR]file not exist\n" unless -s $input;


	my $db = 'sp';
	$db = 'tr' if $input =~ m/trembl/ig;

	my $output = $input.".id.txt";

	my ($id, $ac, @ec);
	open(OUT, ">".$output) || die $!;
	open(FH, $input) || die $!;
	while(<FH>) {
		chomp;
		if	(/^ID\s+(\S+)/)		{ $id = $1; }
		elsif	(/^AC\s+(\S+)/)		{ 
			unless ( $ac ) {
				$ac = $1; 
				$ac =~ s/;//ig; 
			}
		}
		elsif	(/^DE\s+EC=(\S+)/) { 
			my $ec = $1; $ec =~ s/;//ig; 
			push(@ec, $ec); 
		} 
		elsif	(/^\/\//) {
			die "[ERR]undef ID or ACC\n" unless ($id && $ac);
			my $ec_all = '';
			$ec_all = join(";", @ec); 
			if (scalar @ec > 0) {
				$ec_all = join(";", @ec);
				print OUT "$db|$ac|$id\t$ec_all\n";
			} else {
				print OUT "$db|$ac|$id\n";
			}
			$ac = ''; $id = ''; @ec = ();
		}
		else	{ next; }	
	}
	close(FH);
	close(OUT);
}

=head2
 usage -- print usage information
=cut
sub usage
{
	my $usage = qq'
USAGE: $0 -t tool [options]

	uniprot2ec	[]	prepare uniprot database for ec num table (For pathway analysis)
	uniprot2go	prepare uniprot database for go num (For GO analysis)

';

	print $usage;
	exit;
}








