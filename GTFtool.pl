#!/usr/bin/perl

=head1
 GTFtool -- parse GTF file
=cut
use strict;
use warnings;
use IO::File;
use Getopt::Std;

my $version = 0.1;
my %options;
getopts('a:b:c:d:e:f:g:i:j:k:l:m:n:o:p:q:r:s:t:u:v:w:x:y:z:h', \%options);
unless (defined $options{'t'} ) { usage($version); }

if      ($options{'t'} eq 'stats')	{ gtf_stats(\%options, \@ARGV); }
elsif   ($options{'t'} eq 'convert')	{ gtf_convert(\%options, \@ARGV); }
elsif	($options{'t'} eq 'extract')	{ gtf_extract(\%options, \@ARGV); }
elsif	($options{'t'} eq 'gffread')    { gtf_gffread(); }
elsif   ($options{'t'} eq 'togpd')    	{ gtf_to_gpd(\%options, \@ARGV); }
else    { usage($version); }

#################################################################
# kentnf: subroutine						#
#################################################################

=head2
  gtf_to_gpd -- convert gtf to gpd
=cut
sub gtf_to_gpd
{
	my ($options, $files) = @_;

	my $usage = qq'
USAGE: $0 -t togpd -i [gpd/gtf] input_file > ouput_file

example1 $0 -t togpd -i gpd gpd > gtf
example1 $0 -t togpd -i gtf gtf > gpd

';
	print $usage and exit unless defined $$files[0];
	my $input_file = $$files[0];
	die "[ERR]file not exist\n" unless -s $input_file;

	print $usage and exit unless defined $$options{'i'};
	
	if ($$options{'i'} eq 'gtf') 
	{
		my %trans_info = parse_gtf($input_file);
		foreach my $tid (sort keys %trans_info)
		{
			my $chr		= $trans_info{$tid}{'chr'};
			my $gid		= $trans_info{$tid}{'gid'};
			my $strand 	= $trans_info{$tid}{'strand'};
			my @exon	= split("\t",$trans_info{$tid}{'exon'});
			@exon = sort {$a<=>$b} @exon;
			die "[ERR]exon num\n" unless ((scalar @exon) % 2 == 0);
			my $exon_num = (scalar @exon) / 2;

			my @cds;
			if (defined $trans_info{$tid}{'cds'}) {
				@cds = split("\t",$trans_info{$tid}{'cds'});
			} else {
				@cds = @exon; 
			}

			my $start_e = $exon[0];
			my $end_e = $exon[scalar(@exon)-1];
			my $start_c = $cds[0];
			my $end_c = $cds[scalar(@cds)-1];
			my $exon_start = '';
			my $exon_end = '';

			for(my $i=0; $i<@exon; $i=$i+2)
			{
				$exon_start.=$exon[$i].",";
				$exon_end.=$exon[$i+1].",";
			}

			print "$gid\t$tid\t$chr\t$strand\t$start_e\t$end_e\t$start_c\t$end_c\t$exon_num\t$exon_start\t$exon_end\n";
		}
	} 
	elsif ($$options{'i'} eq 'gpd')
	{
		my $fh = IO::File->new($input_file) || die $!;
		while(<$fh>)
		{
			chomp;
			next if $_ =~ m/^#/;
			my @a = split(/\t/, $_);
			die "[ERR]col num: $_\n" if (scalar @a < 11);
			my ($gid, $tid, $chr, $strand, $e_start, $e_end) = ($a[0], $a[1], $a[2], $a[3], $a[9], $a[10]);
			my @m = split(/,/, $e_start);
			my @n = split(/,/, $e_end);
			die "[ERR]exon num\n" unless (scalar(@m) == scalar(@n));
			print "$chr\tGPD\ttranscript\t$a[4]\t$a[5]\t.\t$strand\t.\tgene_id \"$gid\"; transcript_id \"$tid\";\n";
			for(my $i=0; $i<@m; $i++)
			{
				print "$chr\tGPD\texon\t$m[$i]\t$n[$i]\t.\t$strand\t.\tgene_id \"$gid\"; transcript_id \"$tid\";\n";
			}	
		}
		$fh->close;
	}
	else
	{
		die "[ERR]parameter i $$options{'i'}\n";
	}
}

=head2
 gtf_extract -- extract gtf information by ID;
=cut
sub gtf_extract
{
	my ($options, $files) = @_;
	my $usage = qq'
USAGE: $0 -t extract [options] GTF list/BED

* listID could be gene ID, transcript ID;
* key for gene ID: gene_id "";
* key for transcript ID: transcript_id "";

';
	print $usage and exit unless scalar @$files == 2;
	my ($input_gtf, $list_file) = @$files;
	die "[ERR] file not exit $input_gtf" unless -s $input_gtf;
	die "[ERR] file not exit $list_file" unless -s $list_file;

	my $column = 1;
	$column = 4 if ($list_file =~ m/\.bed$/);

	# get listID from file
	my %list_id;
	my $fh1 = IO::File->new($list_file) || die $!;
	while(<$fh1>)
	{
	        chomp;
	        next if $_ =~ m/^#/;
	        my @a = split(/\t/, $_);
        	$list_id{$a[$column-1]} = 1;
	}
	$fh1->close;	

	my ($gid, $tid);

	my $fh2 = IO::File->new($input_gtf) || die $!;
	while(<$fh2>)
	{
		chomp;
		next if $_ =~ m/^#/;
		my @a = split(/\t/, $_);
		$gid = ''; $tid = '';
		if ($_ =~ m/gene_id "(\S+)"; /) { $gid = $1; }
		if ($_ =~ m/transcript_id "(\S+)"; /) { $tid = $1; }
		else { die "Error, can not find gid or tid $_\n"; }
		die "[ERR]Dup ID $_\n" if $gid eq $tid;
		print $_."\n" if $gid && defined $list_id{$gid};
		print $_."\n" if $tid && defined $list_id{$tid};
	}
	$fh2->close;
}

=head2
 stat -- generate statistics information
=cut
sub gtf_stats
{
	my ($options, $files) = @_;
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
 convert -- convert GTF to GFF,BED format
=cut
sub gtf_convert
{
	my ($options, $files) = @_;
	my $subUsage = qq'
USAGE $0 -t convert [options] input_gtf

	-o	output file (must be bed, or gff format, default: input.bed)
	-b	feature type for bed: transcript_id, gene_id, exon (default: transcript_id)
	-a	add chromosome/scaffold seq feature to gtf

';
	print $subUsage and exit unless defined $$files[0] ;
	my $input_gtf = $$files[0];
	die "[ERR]file not exist: $input_gtf\n" unless -s $input_gtf;
	die "[ERR]file format: $input_gtf\n" unless $input_gtf =~ m/\.gtf$/i;
	
	my $output_file = $input_gtf; $output_file =~ s/\.gtf$/\.bed/i;

	my $output_format = 'bed';
	$output_format = $$options{'o'} if defined $$options{'o'};

	die "[ERR]file format: $output_format\n" unless ($output_format eq "bed" || $output_format eq "gff");

	$output_file.=".".$output_format;
	die "[ERR]file exist: $output_file\n" if -s $output_file;

	my $out_format = 'bed';
	$out_format = 'gff' if $output_file =~ m/\.gff$/i;

	my $bed_type = 'transcript_id';
	$bed_type = $$options{'b'} if defined $$options{'b'};
	die "[ERR]bed feature type $bed_type\n" unless $bed_type =~ m/(transcript_id|gene_id|exon)/i;

	my %trans_info = parse_gtf($input_gtf);
	my %gene_info;
	# constract gene info according to trans info
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
	my $out = IO::File->new(">".$output_file) || die $!;
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
			if ($bed_type eq 'transcript_id') {
				print $out $chr,"\t",$start-1,"\t",$end,"\t",$tid,"\t.\t",$strand,"\n";	
			} elsif ( $bed_type eq 'gene_id') {
				unless (defined $gid_uniq{$gid}) {
					my $gstart = $gene_info{$gid}{'start'};
					my $gend = $gene_info{$gid}{'end'};
					$gid_uniq{$gid} = 1;
					print $out $chr,"\t",$gstart-1,"\t",$gend,"\t",$gid,"\t.\t",$strand,"\n";
				}			
			} else {
				my $exon_num = 0;
				for(my $i=0; $i<@exon; $i=$i+2) {
					$exon_num++;
					print $out $chr,"\t",$exon[$i]-1,"\t",$exon[$i+1],"\t",$tid.".exon".$exon_num,"\t.\t",$strand,"\n";
				}
			}
		} 
		elsif ($out_format eq 'gff')
		{
			unless (defined $gid_uniq{$gid}) {
				my $gstart = $gene_info{$gid}{'start'};
				my $gend = $gene_info{$gid}{'end'};
				$gid_uniq{$gid} = 1;
				print $out "$chr\tGTFtool\tgene\t$gstart\t$gend\t.\t$strand\t.\tID=$gid;Name=$gid;\n";
			}

			print $out "$chr\tGTFtool\tmRNA\t$start\t$end\t.\t$strand\t.\tID=$tid;Name=$tid;Parent=$gid;\n";
			for(my $i=0; $i<scalar(@exon)-1; $i=$i+2)
			{
				my $e_start = $exon[$i];
				my $e_end = $exon[$i+1];
				my $e_num = ($i/2)+1;
				print $out "$chr\tGTFtool\tCDS\t$e_start\t$e_end\t.\t$strand\t.\tID=$tid-exon$e_num;Name=$tid-exon$e_num;Parent=$tid;\n";
			}
		}
		else 
		{
			die "[ERR]output format\n";
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
 gffread -- print function for gffread
=cut
sub gtf_gffread
{
	print qq'
USAGE of gffread, which included in tophat package

1. convert GFF to GTF 
	gtffread -T -o output.gtf input.gff

2. extract transcript from GTF
	gffread -w transcript.fasta -g genome.fasta input.gtf

';

	exit;
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
        convert		convert GTF to BED/GFF format
	extract		extract GTF by list
	gffread		usage of gffread
	togpd		convert GTF to GPD format

* the gtf file must have exon feature, transcript_id, and gene_id attributes

';
        exit;
}

