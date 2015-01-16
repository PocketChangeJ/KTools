#!/usr/bin/perl

=head1

 remove Redundancy reads from single-end or paired-end reads
 support fasta and fastq format

 Yi Zheng

 07/13/2013

=cut

use strict;
use warnings;
use IO::File;

my $usage = qq'
perl $0 read1 read2[option]

* support fasta and fastq format
* support single-end and paired-end reads
* For single-end reads, remove same reads
* For paired-end reads, remove same reads for both ends

';

my $file1 = shift || die $usage;
my $file2 = shift;

my %uniq;

if (defined $file2 && -s $file2)
{
	my $out1 = $file1.".uniq";
        my $out2 = $file2.".uniq";
        my $fo1 = IO::File->new(">".$out1) || die "Can not open output file $out1 $!\n";
        my $fo2 = IO::File->new(">".$out2) || die "Can not open output file $out2 $!\n";

	my $fh1 = IO::File->new($file1) || die "Can not open input file $file1 $!\n";
	my $fh2 = IO::File->new($file2) || die "Can not open input file $file2 $!\n";

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
}
else
{
	my $out = $file1.".uniq";
	my $fo = IO::File->new(">".$out) || die "Can not open output file $out $!\n";
	my $fh = IO::File->new($file1) || die "Can not open input file $file1 $!\n";

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
}

