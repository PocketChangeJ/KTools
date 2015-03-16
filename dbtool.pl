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

if	($options{'t'} eq 'uniprot2ec')		{ uniprot2ec(\@ARGV); }
elsif	($options{'t'} eq 'uniprot2go')		{ uniprot2go(\@ARGV); }
else	{ usage(); }

#################################################################
# kentnf: subroutine						#
#################################################################

=head2
 uniprot2go -- convert uniport to GO
=cut
sub uniprot2go 
{
	my $files = shift;

	my $usage = qq'
USAGE: $0 -t uniprot2go idmapping_selected.tab > uniport_idmapping_GO_MMYYYY

* convert idmapping to GO mapping file
  then use the GO mapping file for GO and pathway annotation 

';
	print $usage and exit unless defined $$files[0];
	die "[ERR]file not exist\n" unless -s $$files[0];

	my $input = $$files[0];
	open(FH, $input) || die $!;
	while(<FH>)
	{
		chomp;
		my @a = split(/\t/, $_);
		next unless defined $a[6];
		next unless $a[6] =~ m/GO/;
		print "$a[0]|$a[1]\t$a[6]\n";
	}
	close(FH);
}

=head2
 uniport2ec -- convert uniport to EC num
=cut
sub uniprot2ec
{
	my $files = shift;

	my $usage = qq'
USAGE: $0 -t uniprot2ec input_database(dat format)

* convert to dat to ID and EC num
-t uniprot2ec uniprot_sp_plant.dat > uniprot_sp_plants.dat.id.txt
-t uniprot2ec uniprot_tr_plant.dat > uniprot_tr_plants.dat.id.txt

* the uniprot_sp/tr_plants.dat.id.txt could be used for extract info
  related plant species

* then the uniprot_sp/tr_plants.dat.id.txt only kept the record with
  EC ID for pathway analysis
  grep -P "\\t"  uniprot_sp/tr_plants.dat.id.txt | gzio > uniprot_sp/tr_plants.dat.id.txt.gz 

  
';
	print $usage and exit unless defined $$files[0];
	die "[ERR]file not exist\n" unless -s $$files[0];

	my $input = $$files[0];
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

	uniprot2ec	prepare uniprot database for ec num table (For pathway analysis)
			[input: uniprot_sp_plant.dat, uniprot_tr_plant.dat]
	uniprot2go	prepare uniprot database for go num (For GO analysis)
			[input: idmapping_selected.tab]

';

	print $usage;
	exit;
}








