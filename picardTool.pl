#!/usr/bin/perl

=head1
 picardTool -- learn how to huse picardTools
=cut
use strict;
use warnings;
use IO::File;
use FindBin;
use Getopt::Std;

# check soft ware setting
my $picard_path = find_software('picard');

my %options;
getopts('a:b:c:d:e:f:g:i:j:k:l:m:n:o:p:q:r:s:t:u:v:w:x:y:z:h', \%options);
unless ( defined $options{'t'} ) { usage(); }

if	( $options{'t'} eq 'CollectInsertSizeMetrics')	{ picard_CollectInsertSizeMetrics(\@ARGV); }
else	{ usage(); }

=head2
 picard_CollectInsertSizeMetrics
=cut
sub picard_CollectInsertSizeMetrics 
{
	my $files = shift;
	my $usage = qq'
perl $0 -t [picard tool name] input_file_list
';

	print $usage and exit unless defined $$files[0];

	foreach my $f (@$files) {
		my $out = $f.".insertsize.txt";
		my $hist = $f.".hist.pdf";
		print "java -jar $picard_path/picard.jar CollectInsertSizeMetrics I=$f O=$out H=$hist\n"; 
	}
}

=head2
 usage
=cut
sub usage
{
	my $usage = qq'
perl $0 -t [picard tool name] input_files

[ CollectInsertSizeMetrics: get insert size ] [ ... ]

';
	print $usage and exit;
}

=head2
 find_software
=cut 
sub find_software
{
	my $name = shift;
	my $config_file = $FindBin::RealBin."/bin/software.config";
	die "[ERR]file not exist $config_file\n" unless -s $config_file;

	my $path = '';
	my $fh = IO::File->new($config_file) || die $!;
	while(<$fh>) {
		chomp; next if $_ =~ m/^#/; 
		my @a = split(/=/, $_, 2);
		next if scalar @a < 2;
		$a[0] =~ s/^\s+//ig; $a[0] =~ s/\s+$//ig;
		$a[1] =~ s/^\s+//ig; $a[1] =~ s/\s+$//ig;
		$path = $a[1] if $a[0] eq $name;
	} 
	$fh->close;

	die "[ERR]please set software path in $config_file file\n" unless $path;
	return $path; 
}

