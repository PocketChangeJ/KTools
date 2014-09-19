#!/usr/bin/perl -w 
=head1 filter SAM file

 filter SAM file for next analysis
 Yi Zheng

=cut
use strict; 
use warnings;
use IO::File;
use Getopt::Long;

my $usage = qq'
perl $0 
	-i input_sam_file 
	-o output_sam_file 
	-h help info	

	[parameter for filter, only one can be set]
	-e add read by expression for small RNA uniq read
	-l filter by read length   ( >= 15 )
	-m filter by max hits ( <=6 )
	-n filter un-aligned reads
	-v filter by mismatch (1-3)
	-c count reads number for sRNA align
	-u filter multiple hits (see description)
	
* description of filter -u
Rules   1) if reads has uniq hits, keep it.
        2) if reads has multiply hits with same edit distance, remove all the hits for this reads.
        3) if reads has multiplu hist with diff edit distance, check how many hits for best edit distance,
           if only one, keep the best hits; if multiply, remove hits for this reads

';

my ($help, $sam_file_in, $sam_file_out, $add_byexp, $filter_byreadlength, $filter_multihits, $filter_unalignedread, $filter_editDistanceCutoff);
my $filter_uniq;
my $count_read;

GetOptions(
	"h"	=> \$help,
	"i=s"	=> \$sam_file_in,
	"o=s"	=> \$sam_file_out,
	
	"e"	=> \$add_byexp,
	"l=i"	=> \$filter_byreadlength,
	"m=i"	=> \$filter_multihits,
	"n"	=> \$filter_unalignedread,
	"v=i"	=> \$filter_editDistanceCutoff,
	"u"	=> \$filter_uniq,
	"c"	=> \$count_read,
);

die $usage if $help;
die $usage unless $sam_file_in;
die $usage unless $sam_file_in =~ m/\.sam$/;

if (defined $sam_file_out)
{
	die $usage unless $sam_file_out =~ m/\.sam$/;
}
else
{
	$sam_file_out ||= $sam_file_in;
	$sam_file_out =~ s/\.sam$/_out\.sam/;
}

#################################################################
# checking filter type						#
#################################################################
my $type_num = 0;
my @types = ($add_byexp, $filter_byreadlength, $filter_multihits, $filter_unalignedread, $filter_uniq, $count_read);
foreach my $type (@types) {
	if ($type) { $type_num++; }
}

if (defined $filter_editDistanceCutoff && $filter_editDistanceCutoff >= 0 && $filter_editDistanceCutoff <= 3) { $type_num++; }
if ($type_num != 1) { die "Error in number of filter type $usage\n"; }
if (defined $filter_byreadlength && $filter_byreadlength < 15 )  { die $usage; }
if (defined $filter_multihits && $filter_multihits < 1)  { die $usage; }

#################################################################
# main								#
#################################################################
my %outflag = make_flag_list(); 

if (defined $add_byexp) 		{ add_ByExp($sam_file_in, $sam_file_out); }
if (defined $filter_byreadlength)	{ filter_ByReadLength($sam_file_in, $sam_file_out, $filter_byreadlength); }
if (defined $filter_multihits)		{ filter_MultiHits($sam_file_in, $sam_file_out, $filter_multihits); }
if (defined $filter_unalignedread)	{ filter_UnAlignedRead($sam_file_in, $sam_file_out, \%outflag); }
if (defined $filter_editDistanceCutoff)	{ filter_EditDistanceCutoff($sam_file_in, $sam_file_out, $filter_editDistanceCutoff); }
if (defined $filter_uniq)		{ filter_uniq($sam_file_in, $sam_file_out); }
if (defined $count_read)		{ count_read($sam_file_in); }

#################################################################
# kentnf: subroutine						#
#################################################################

=head count read
  count number of aligned uniq sRNA reads
=cut
sub count_read
{
	my $sam_file_in = shift;


	my %uniq_read;
	my $seq_num = 0;
	
	my $in = IO::File->new($sam_file_in) || die $!;
	while(<$in>)
	{
		chomp;
		next if $_ =~ m/^@/;
		my @a = split(/\t/, $_);
		
		unless ( defined $uniq_read{$a[0]} )
		{
			my @b = split(/-/, $a[0]);
			my $num = pop @b;
			die "[ERR]seq num for $a[0]\n" if $num < 1;
			$seq_num = $seq_num + $num;
			$uniq_read{$a[0]} = 1;
		}
	}
	$in->close;

	print "File:$sam_file_in\tSeqNum:$seq_num\n";
}


=head1 add_ByExp

 add small RNA Read to sam by expression

=cut
sub add_ByExp
{
	my ($sam_file_in, $sam_file_out) = @_;

	my $out = IO::File->new(">".$sam_file_out) || die "Can not open output sam file $!\n";
	my $in = IO::File->new($sam_file_in) || die "Can not open input sam file $!\n";

	while(<$in>)
	{
		if (m/^@/) { print $out $_;}
		else
		{
			chomp;
			my @a = split(/\t/, $_, 2);
			my $read_id = $a[0];
			my @b = split(/-/, $read_id);
			if (scalar(@b) != 2) { die "Error in read ID for expression information, $_\n"; }
			
			if ($b[1] == 1)
			{
				print $out $_."\n";
			}
			else
			{
				for(my $i=1; $i<=$b[1]; $i++)
				{
					print $out $b[0]."_".$i."\t".$a[1]."\n";
				}
			}
		}
	}
}

=head1 filter_ByReadLength

 filter small RNA SAM by Read length (keep small RNA read alignment by length)

=cut
sub filter_ByReadLength
{
	my ($sam_file_in, $sam_file_out, $read_length) = @_;
	
	my $out = IO::File->new(">".$sam_file_out) || die "Can not open output sam file $!\n";
	my $in = IO::File->new($sam_file_in) || die "Can not open input sam file $!\n";

	while(<$in>)
	{
		if (m/^@/) { print $out $_;}
		else
		{
			chomp;
			my @a = split(/\t/, $_);
			if ( length($a[9]) eq $read_length ){ print $out $_."\n"; }
		}
	}

	$in->close;
	$out->close;
}

=head1 filter_uniq

 filter small RNA sam by MultiHits (remove MultiHits)

=cut
sub filter_uniq
{
        my ($sam_file_in, $sam_file_out, $hits) = @_;

        #########################################################
        # identify reads more than N hits                       #
        #########################################################
        my %hash;
	
        my $fin = IO::File->new($sam_file_in) || die "Can not open input file $sam_file_in \n";
        while(<$fin>)
        {
                chomp;
                unless ($_ =~ m/^@/)
                {
                        my @a = split(/\t/, $_);
			my $ed;
			if ($_ =~ m/NM:i:(\d+)/) { $ed = $1; }
			else { die "[ERR]Edit distance $_\n"; }

                        if (defined $hash{$a[0]}{$ed})	{ $hash{$a[0]}{$ed}++;   }
                        else				{ $hash{$a[0]}{$ed} = 1; }

			if (defined $hash{$a[0]}{'bestED'}) {
				$hash{$a[0]}{'bestED'} = $ed if $ed < $hash{$a[0]}{'bestED'};
			} else {
				$hash{$a[0]}{'bestED'} = $ed;
			}
                }
        }
        $fin->close;

        #########################################################
        # output result                                         #
        #########################################################
        my $out = IO::File->new(">".$sam_file_out) || die "Can not open output sam file $!\n";
        my $in = IO::File->new($sam_file_in) || die "Can not open input sam file $!\n";

        while(<$in>)
        {
                if (m/^@/) { print $out $_;}
                else
                {
                        chomp;
                        my @a = split(/\t/, $_);
			my $ed;
			if ($_ =~ m/NM:i:(\d+)/) { $ed = $1; }
			else { die "[ERR]Edit distance $_\n"; }

			if ( $hash{$a[0]}{'bestED'} == $ed && $hash{$a[0]}{$ed} == 1) 
			{
                                print $out $_."\n";
                        }
                }
        }

        $in->close;
        $out->close;
}

=head1 filter_MultiHits

 filter small RNA sam by MultiHits (remove MultiHits)

=cut
sub filter_MultiHits
{
	my ($sam_file_in, $sam_file_out, $hits) = @_;

	#########################################################
	# identify reads more than N hits			#
	#########################################################
	my %hash;
	my $fin = IO::File->new($sam_file_in) || die "Can not open input file $sam_file_in \n";
	while(<$fin>)
	{
		chomp;
		unless ($_ =~ m/^@/) 
		{
			my @a = split(/\t/, $_);
			if (defined $hash{$a[0]}) 	{ $hash{$a[0]}++;   }
			else 				{ $hash{$a[0]} = 1; }
		}
	}
	$fin->close;

	#########################################################
	# output result						#
	#########################################################
	my $out = IO::File->new(">".$sam_file_out) || die "Can not open output sam file $!\n";
	my $in = IO::File->new($sam_file_in) || die "Can not open input sam file $!\n";

	while(<$in>)
	{
		if (m/^@/) { print $out $_;}
		else
		{
			chomp;
			my @a = split(/\t/, $_);
			if ($hash{$a[0]} <= $hits )
			{
				print $out $_."\n";
			}
		}
	}

	$in->close;
	$out->close;
}

=head1 

=cut
sub filter_EditDistanceCutoff
{
	my ($sam_file_in, $sam_file_out, $filter_editDistanceCutoff) = @_;

	my $out = IO::File->new(">".$sam_file_out) || die "Can not open output sam file $!\n";
	my $in = IO::File->new($sam_file_in) || die "Can not open input sam file $!\n";
	while(<$in>)
	{
		if (m/^@/) {
			print $out $_;
		} else {
			chomp;
			my @a  = split(/\t/, $_);
			my $print = 0;
			for (my $col = 11; $col<scalar(@a); $col++) 
			{
				if ( $a[$col] =~ m/NM:i:(\d+)$/) 
				{
					my $edit_distance = $1;
					if ($edit_distance <= $filter_editDistanceCutoff) {
						$print = 1;
					}
				}
			}

			if ($print == 1) {
				print $out $_."\n";
			}
		}
	}
	$in->close;
	$out->close;
}

=head1 filter_UnAlignedRead

 filter small RNA SAM to remove unalinged reads

=cut
sub filter_UnAlignedRead
{
	my ($sam_file_in, $sam_file_out, $outflag) = @_;

	my $out = IO::File->new(">".$sam_file_out) || die "Can not open output sam file $!\n";
	my $in = IO::File->new($sam_file_in) || die "Can not open input sam file $!\n";
	while(<$in>)
	{
		if (m/^@/) {
                	print $out $_;
        	}else{
                	chomp;
                	my @ta = split(/\t/, $_);
                	$#ta > 10 or next;
                	$$outflag{$ta[1]} and next;
			# my $is_u = 0;
			# for my $tb (@ta[11..$#ta]) {
			#   $tb eq 'XT:A:U' and do { $is_u = 1; last; }; # Other alternative values are ":M", ":R", ":N". And I am not sure about ":M".
			# }
			# $is_u == 1 and print STDOUT "$_\n";
                	print $out "$_\n";
                	last;
        	}
	}

	while (<$in>) 
	{
        	chomp;
        	my @ta = split(/\t/, $_);
        	$#ta > 10 or next;
        	$$outflag{$ta[1]} and next;
		# my $is_u = 0;
		# for my $tb (@ta[11..$#ta]) {
		#	$tb eq 'XT:A:U' and do { $is_u = 1; last; };
		# }
		# $is_u == 1 and print STDOUT "$_\n";
        	print $out "$_\n";
	}

	$in->close;
	$out->close;
}

=head1 make_flag_list 

 Make FLAG list for SAM file

=cut
sub make_flag_list
{
	my %flag;
	for (0 .. 2047) {
		my $binum = unpack("B32", pack("N", $_)); 
		$flag{$_} = [ reverse( split(//, sprintf("%010d", $binum) ) ) ]; 
	}
	# BitPos	Flag	Chr	Description
	# 0	0x0001	p	the read is paired in sequencing
	# 1	0x0002	P	the read is mapped in a proper pair
	# 2	0x0004	u	the query sequence itself is unmapped
	# 3	0x0008	U	the mate is unmapped
	# 4	0x0010	r	strand of the query (1 for reverse)
	# 5	0x0020	R	strand of the mate
	# 6	0x0040	1	the read is the first read in a pair
	# 7	0x0080	2	the read is the second read in a pair
	# 8	0x0100	s	the alignment is not primary
	# 9	0x0200	f	the read fails platform/vendor quality checks
	# 10	0x0400	d	the read is either a PCR or an optical duplicate
	###### FLAG list prepared. 
	###### Set FLAGs in or out. 
	my (%inflag, %outflag); 
	for (sort { $a<=>$b } keys %flag) {
		my @ta = @{$flag{$_}}; 
		$ta[2] == 1 and $outflag{$_} = 1; 
	}
	###### Set FLAGs OK. 

	return %outflag;
}

