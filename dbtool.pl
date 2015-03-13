#!/usr/bin/perl

=head
 dbtool.pl -- parse database for in house use 
=cut
use strict;
use warnings;






my $usage = qq'
USAGE: $0 input_database(dat format)

';

my $database = shift || die $usage;

my $db = 'sp';
$db = 'tr' if $database =~ m/trembl/ig;
uniport2id($database, $db);

#################################################################
# kentnf: subroutine						#
#################################################################

=head2
 convert uniport to fasta 
=cut
sub uniport2id
{
	my ($input, $db) = @_;
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
			$ec_all = join(";", @ec) if (scalar @ec > 0);
			print OUT "$db|$ac|$id\t$ec_all\n";	
			$ac = ''; $id = ''; @ec = ();
		}
		else	{ next; }	
	}
	close(FH);
	close(OUT);
}

