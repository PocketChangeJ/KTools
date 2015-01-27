#!/usr/bin/perl

package reads;

use strict;
use warnings;
use Moose;
use IO::File;

our @ISA = qw (Exporter);
our @EXPORT = qw (removeDup);

sub getPrefix
{
	my @files = @_;
}

=head2
my $usage = qq'
perl $0 read1 read2[option]
* support fasta and fastq format
* support single-end and paired-end reads
* For single-end reads, remove same reads
* For paired-end reads, remove same reads for both ends
';
=cut

sub removeDup
{
	my @files = @_;

	# check file exist
	my %link_uniq;
	my @suffix = ('.fa', '.fq', '.fasta', '.fastq');
	foreach my $f (@files) {
		die "[ERR]file not exist: $f\n" unless -s $f;
		my ($prefix, $suffix, $uniq_f);
		foreach my $s ( @suffix ) {
			$suffix = $s if $f =~ m/\Q$s\E$/;
		}
		die "[ERR]file suffix: $f\n" unless $suffix;

		$prefix = $f; $prefix =~ s/\Q$suffix\E$//;
		$uniq_f = $prefix."_uniq".$suffix;

		$link_uniq{$f} = $uniq_f;
	}

	# put reads to this hash
	my %uniq;

	if (scalar(@files) == 2)
	{
		my $out1 = $link_uniq{$files[0]};
		my $out2 = $link_uniq{$files[1]};
   
		my $fo1 = IO::File->new(">".$out1) || die "Can not open output file $out1 $!\n";
		my $fo2 = IO::File->new(">".$out2) || die "Can not open output file $out2 $!\n";

		my $fh1 = IO::File->new($files[0]) || die "Can not open input file $files[0] $!\n";
		my $fh2 = IO::File->new($files[1]) || die "Can not open input file $files[1] $!\n";

		my ($id, $seq, $uniq_read, $format);
		my %firstread;

		while(<$fh1>)
		{
			chomp;
			$id = $_; 
			$format = '';
			if      ( $id =~ m/^>/ ) { $format = 'fasta'; $id =~ s/^>//; }
			elsif   ( $id =~ m/^@/ ) { $format = 'fastq'; $id =~ s/^@//; }
			else    { die "Error at format for read $id.\n"; }
			$id =~ s/^>//; $id =~ s/\s+.*//ig;
			$seq = <$fh1>; chomp($seq);
			$firstread{$id} = $seq;
			if ( $format eq 'fastq' ) { <$fh1>; <$fh1>; }
		}

		while(<$fh2>)
		{
			chomp;
			$id = $_; 
			$format = '';
			if      ( $id =~ m/^>/ ) { $format = 'fasta'; $id =~ s/^>//; }
			elsif   ( $id =~ m/^@/ ) { $format = 'fastq'; $id =~ s/^@//; }
			else    { die "Error at format for read $id.\n"; }
			chomp($id); $id =~ s/^>//; $id =~ s/\s+.*//ig;

			$seq = <$fh2>; chomp($seq);
		
			unless (defined $firstread{$id}) { die "Error in $id\n"; }

			$uniq_read = $firstread{$id}."\t".$seq;

			# output uniq reads;
			unless ($uniq{$uniq_read})
			{
				$uniq{$uniq_read} = 1;
				print $fo1 ">$id\n$firstread{$id}\n";
				print $fo2 ">$id\n$seq\n";
			}

			if ( $format eq 'fastq' ) { <$fh2>; <$fh2>; }
		}
	
		$fh1->close;
		$fh2->close;
		$fo1->close;
		$fo2->close;

		return ($out1, $out2);
	}
	elsif (scalar(@files) == 1)
	{
		my $out = $link_uniq{$files[0]};
		my $fo = IO::File->new(">".$out) || die "Can not open output file $out $!\n";
		my $fh = IO::File->new($files[0]) || die "Can not open input file $files[0] $!\n";

		my ($id, $seq, $format);
		while(<$fh>)
		{
			chomp;
			$id = $_;
			$format = '';
			if 	( $id =~ m/^>/ ) { $format = 'fasta'; $id =~ s/^>//; }
			elsif 	( $id =~ m/^@/ ) { $format = 'fastq'; $id =~ s/^@//; }
			else	{ die "Error at format for read $id.\n"; }
	
			$seq = <$fh>; 
			chomp($seq);

			unless ( $uniq{$seq} )
			{
				$uniq{$seq} = 1;
				print $fo ">".$id."\n".$seq."\n";
			}

			if ( $format eq 'fastq' ) { <$fh>; <$fh>; }
		}
		$fh->close;
		$fo->close;

		return($out);
	}
	else 
	{
		die "[ERR]input file\n";
	}
}

1;
