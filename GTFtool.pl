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
getopts('i:f:h', \%options);

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
	my $options = shift;
	my $subUsage = qq'
USAGE $0 convert [options]
	-i	input file
	-f	output format (bed/gff/tab)
	-t	feature type for bed: tid, gid, exon (default: tid)
	-a	add chromosome/scaffold seq feature to gtf

';
	print $subUsage and exit unless ( $$options{'i'} && $$options{'f'} ) ;
	my ($inFile, $out_format, $bed_type, $out_prefix);
	$inFile = $$options{'i'};
	$out_format = $$options{'f'};
	$out_prefix = $$options{'i'}; 
	$out_prefix =~ s/\.gtf$//ig;
	$bed_type = 'tid';
	$bed_type = $$options{'t'} if defined $$options{'t'};
	my $out_file = $out_prefix.".".$out_format;
	die "[ERR]output exist\n" if -s $out_file;

	my %trans_info = parse_gtf($$options{'i'});
	my %gene_info;
	# generate gene info
	my ($chr, $gid, $start, $end, $strand);
	foreach my $tid (sort keys %trans_info)
	{
		$chr = $trans_info{$tid}{'chr'};
		$strand = $trans_info{$tid}{'strand'};
		$gid = $trans_info{$tid}{'gid'};
		my @exon = split(/\t/, $trans_info{$tid}{'exon'});
		@exon = sort {$a<=>$b} @exon;
		$start = $exon[0];
		$end = $exon[scalar(@exon)-1];

		# generate gene info
		if (defined $gene_info{$gid}{'start'}) {
			$gene_info{$gid}{'start'} = $start if $start < $gene_info{$gid}{'start'};
		} else {
			$gene_info{$gid}{'start'} = $start;
		}
		
		if (defined $gene_info{$gid}{'end'}) {
                        $gene_info{$gid}{'end'} = $end if $end > $gene_info{$gid}{'end'};
                } else {
                        $gene_info{$gid}{'end'} = $end;
                }	

		if (defined $gene_info{$gid}{'chr'}) {
			die "[ERR]inconsistency chr for $gid $tid\n" if $chr ne $gene_info{$gid}{'chr'};
		} else {
			$gene_info{$gid}{'chr'} = $chr;
		}
		
		if (defined $gene_info{$gid}{'strand'}) {
			die "[ERR]inconsistency strand for $gid $tid\n" if $strand ne $gene_info{$gid}{'strand'};
		} else {
			$gene_info{$gid}{'strand'} = $strand;
		}

		if (defined $gene_info{$gid}{'tid'}) {
			$gene_info{$gid}{'tid'}.="\t".$tid;
		} else {
			$gene_info{$gid}{'tid'} = $tid;
		}
	}

	# convert to bed/tab/gff format
	my $out = IO::File->new(">".$out_file) || die $!;
	my %gid_uniq;
	foreach my $tid (sort keys %trans_info)
        {
                $chr = $trans_info{$tid}{'chr'};
                $strand = $trans_info{$tid}{'strand'};
                $gid = $trans_info{$tid}{'gid'};
                my @exon = split(/\t/, $trans_info{$tid}{'exon'});
                @exon = sort {$a<=>$b} @exon;
                $start = $exon[0];
                $end = $exon[scalar(@exon)-1];	

		if ($out_format eq 'bed') 
		{

		} 
		elsif ($out_format eq 'tab')
		{

		}
		elsif ($out_format eq 'gff')
		{
			unless (defined $gid_uniq{$gid}) {
				my $gstart = $gene_info{$gid}{'start'};
				my $gend = $gene_info{$gid}{'end'};
				$gid_uniq{$gid} = 1;
				print $out "$chr\tGTFtool\tgene\t$gstart\t$gend\t.\t$strand\t.\tID=$gid;Name=$gid;\n";
			}

			print $out "$chr\tGTFtool\ttranscript\t$start\t$end\t.\t$strand\t.\tID=$tid;Name=$tid;parent=$gid;\n";
			for(my $i=0; $i<scalar(@exon)-1; $i=$i+2)
			{
				my $e_start = $exon[$i];
				my $e_end = $exon[$i+1];
				print $out "$chr\tGTFtool\texon\t$e_start\t$e_end\t.\t$strand\t.\tID=$tid-exon;Name=$tid-exon;\n";
			}
		}
	}
	$out->close;
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
				$c[1] =~ s/"//ig;
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

%grades = (
  student1 => 90,
  student2 => 75,
  student3 => 96,
  student4 => 55,
  student5 => 76,
);

print "\nGRADES IN ASCENDING NUMERIC ORDER:\n";
foreach $key (sort hashValueAscendingNum (keys(%grades))) {
   print "\t$grades{$key} \t\t $key\n";
}

print "\nGRADES IN DESCENDING NUMERIC ORDER:\n";
foreach $key (sort hashValueDescendingNum (keys(%grades))) {
   print "\t$grades{$key} \t\t $key\n";
}

=cut
#sub hashValueAscendingNum {
#	$grades{$a} <=> $grades{$b};
#}

#sub hashValueDescendingNum {
#	$grades{$b} <=> $grades{$a};
#}

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

