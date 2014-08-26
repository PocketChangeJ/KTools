#!/usr/bin/perl

=head1
 sRNAtool -- tools for sRNA data preparation
=cut

use strict;
use warnings;
use IO::File;
use Getopt::Std;

my $version = 0.1;
if (@ARGV < 1) { usage($version);}
my $tool = shift @ARGV;

my %options;
getopts('i:o:p:fuh', \%options);

if	($tool eq 'convert') { convert(\%options); }
elsif	($tool eq 'lengthd') { lengthd(\%options); }
elsif	($tool eq 'unique' ) { unique(\%options);  }
elsif   ($tool eq 'norm' )   { norm(\%options);    }
elsif   ($tool eq 'normcut') { normcut(\%options); }
else	{ usage($version); } 

#################################################################
# kentnf: subroutine						#
#################################################################
=head2
 convert -- convert table format (GEO database) to fasta format (default), or fasta format to table format
=cut

sub convert
{
	my $options = shift;

	my $subUsage = qq'
USAGE: $0 convert [options]
	-i	input file 
	-o	output file
	-p	prefix of out seqID (for table convert to fasta)
	-f	convert table to fasta (default:0) / fasta to table (1)

';

	print $subUsage and exit unless $$options{'i'}; 

	my ($inFile, $outFile, $format, $prefix);
	$inFile = $$options{'i'};
	$prefix = $inFile; $prefix =~ s/\..*//;
	$prefix = $$options{'p'} if $$options{'p'};
	$format = 0;
	$format = 1 if $$options{'f'};
	my %read; my $num;
	if ($format) {
		my $in = IO::File->new($inFile) || die $!;
		while(<$in>)
		{
			chomp;
			my $id = $_;
			my @a = split(/-/, $id);
			die "[ERR]sRNA ID: $id\n" unless @a == 2;
			die "[ERR]sRNA num: $a[1]" unless $a[2] > 0;
			my $rd = <$in>; chomp($rd);
			if ( defined $read{$a[0]} ) {
				$read{$rd} = $read{$rd} + $a[1];
				print "[WARN]Repeat Uniq sRNA Read $a[0]\n";
			} else {
				$read{$rd} = $a[1];
			}
		}
		$in->close;

		foreach my $r (sort keys %read)
		{
			print $r."\t".$read{$r}."\n";
		}


	} else {
		my $in = IO::File->new($inFile) || die $!;
		while(<$in>)
		{
			chomp;
			my @a = split(/\t/, $_);
			if ( defined $read{$a[0]} ) {
				$read{$a[0]} = $read{$a[0]} + $a[1];
				print "[WARN]Repeat Uniq sRNA Read $a[0]\n"; 
			} else {
				$read{$a[0]} = $a[1];
			}
		}
		$in->close;

		foreach my $r (sort keys %read)
		{
        		$num++;
			my $count = $read{$r};
			print ">".$prefix."A".$num."-".$count."\n$r\n";
		}
	}
	# foreach my $o (sort keys %$options) { print "$o\t$$options{$o}\n"; }
}

=head2
 unique -- convert the clean sRNA to unique (remove duplication)
=cut
sub unique
{
	my $options = shift;

	my $subUsage = qq'
USAGE: $0 unique [options]
        -i      input file 
        -u      input reads file is sRNA clean (default:0) / uniq (1) format

* convert the clean sRNA to unique
  convert the unique sRNA to clean if -u provide
* output example
  >id000001-10930   # 10930 is the number of sRNA
* support fastq file
* output read is sorted by number, could check the high expressed read

';

	print $subUsage and exit unless $$options{'i'};
	my $input_file = $$options{'i'};
	die "[ERR]cat not find input file $$options{'i'}\n" unless -s $$options{'i'};

	if ($$options{'u'}) { # convert uniq to clean(norm) format
		my $output_file = $input_file;
		$output_file =~ s/\..*$/_norm\.fasta/;
		die "[ERR]Output exist, $output_file\n" if -s $output_file;
		my $out = IO::File->new(">".$output_file) || die $!;

		my ($fh);
		if ($input_file =~ m/\.gz$/) {
			open($fh, '-|', "gzip -cd $input_file") || die $!;
		} else {
			open($fh, $input_file) || die $!;
		}

		while(<$fh>) 
		{
			chomp;
			my $id = $_;
			die "[ERR]seq format $id\n" unless $id =~ m/^>/;
			$id =~ s/^>//;
			my @a = split(/-/, $id);
			my $seq_num = pop @a;
			my $seq_id = join("-", @a);
        		my $seq = <$fh>; chomp($seq);

			for(my $i=1; $i<=$seq_num; $i++)
			{
				my $new_id = $seq_id."_$i";
				print $out ">".$new_id."\n".$seq."\n";
			}
		}

		close($fh);
		$out->close;

	} else {  # convert the clean to uniq 

		my $output_file = $input_file;
		$output_file =~ s/\..*$/_uniq\.fasta/;
		die "[ERR]Output exist, $output_file\n" if -s $output_file;
		
		my ($fh, $format, $total_seq_num, $length);
		if ($input_file =~ m/\.gz$/) {
			open($fh, '-|', "gzip -cd $input_file") || die $!;
		} else {
			open($fh, $input_file) || die $!;
		}

		# load seq/read count to hash
		my %seq_count;
		while(<$fh>)
		{
			chomp; 
			my $id = $_; $id =~ s/ .*//;
			$format = '';
			if ($id =~ m/^>/) { $format = 'fasta'; }
			elsif ($id =~ m/^@/) { $format = 'fastq'; }
			else { die "[ERR]seq fromat: $id\n"; }
		
			my $seq = <$fh>; chomp($seq);
			if ( defined $seq_count{$seq} ) {
				$seq_count{$seq}++;
			} else {
				$seq_count{$seq} = 1;
			}
			
			if ($format eq 'fastq') { <$fh>; <$fh>; }
			$total_seq_num++;
		}
		$fh->close;
		
		$length = length($total_seq_num);

		# sort by num for duplicate seq/read
		my %seq_count_sort;
		foreach my $sq (sort keys %seq_count) {
			my $count = $seq_count{$sq};
			if ($count > 1 ) {
				if (defined $seq_count_sort{$count}) {
					$seq_count_sort{$count}.= "\t".$sq;
				} else {
					$seq_count_sort{$count} = $sq;
				}
				delete $seq_count{$sq};
			} 
		}

		# output result
		my $seq_num = 0;
		my $out = IO::File->new(">".$output_file) || die $!;
		# --- output duplicate seq/read
		foreach my $ct (sort { $b<=>$a } keys %seq_count_sort) { 
			my @seq = split(/\t/, $seq_count_sort{$ct});
			foreach my $sq (@seq) {
				$seq_num++;
				my $zlen = $length - length($seq_num);
				my $z = "0"x$zlen;
				my $seq_id = "sRU".$z.$seq_num;
				print $out ">$seq_id-$ct\n$sq\n";
			}
		}	

		# --- output single seq/read
		foreach my $sq (sort keys %seq_count) {
			$seq_num++;
			my $zlen = $length - length($seq_num);
			my $z = "0"x$zlen;
			my $seq_id = "sRU".$z.$seq_num;
			print $out ">$seq_id-1\n$sq\n";
		}		

		$out->close;
	}
}

=head2
 lengthd -- get length distribution of sRNA 
=cut
sub lengthd
{
	my $options = shift;
	
	my $subUsage = qq'
USAGE: $0 lengthd [options]
        -i      input file 
        -u      input reads file is sRNA clean (default:0) / uniq (1) format

';

	print $subUsage and exit unless $$options{'i'};
	my $input_seq = $$options{'i'};
	my $key = $input_seq; $key =~ s/\..*$//;
	my $output_table = $key.".table";
	my $output_plots = $key.".pdf";
	my $output_image = $key.".png";

	my %length_dist;
	my $seq_num = 0;
	my ($seq_id_info, $seq_id, $seq_desc, $format, $sequence, $seq_length, $uniq_count);
	
	my $fh = IO::File->new($input_seq) || die $!;
	while(<$fh>)
	{
		chomp;
		$seq_id_info = $_;
		if      ($seq_id_info =~ m/^>/) { $format = 'fasta'; $seq_id_info =~ s/^>//; }
		elsif   ($seq_id_info =~ m/^@/) { $format = 'fastq'; $seq_id_info =~ s/^@//; }
		else    { die "[ERR]sRNA ID: $seq_id_info\n"; }
		($seq_id, $seq_desc) = split(/\s+/, $seq_id_info, 2);
		unless ($seq_desc) { $seq_desc = ""; }
		
		$sequence = <$fh>; chomp($sequence);
		$seq_length = length($sequence);

		if ($$options{'u'}) {
			my @nn = split(/-/, $seq_id);
			$uniq_count = $nn[scalar(@nn)-1];
			die "[ERR]sRNA count $seq_id_info, $seq_id, $uniq_count\n" if $uniq_count < 1;
			$seq_num = $seq_num + $uniq_count;

			if ( defined $length_dist{$seq_length} ) { $length_dist{$seq_length} = $length_dist{$seq_length} + $uniq_count; }
			else { $length_dist{$seq_length} = $uniq_count; }
		} else {
			$seq_num++;

			if ( defined $length_dist{$seq_length} ) { $length_dist{$seq_length}++; }
			else { $length_dist{$seq_length} = 1; }
		}

		if ($format eq 'fastq') { <$fh>; <$fh>; }
	}
	$fh->close;

	# output lengt distribution tables
	my $out = IO::File->new(">".$output_table) || die "Can not open output table file $output_table $!\n";
	foreach my $len (sort keys %length_dist) {
		my $freq = sprintf('%.4f', $length_dist{$len}/$seq_num);
		$freq = $freq * 100;
		print $out "$len\t$length_dist{$len}\t$freq\n";
	}
	$out->close;	

	# R code for length distribution
my $R_LD =<< "END";
a<-read.table("$output_table")
x<-a[,1]
y<-a[,2]
dat <- data.frame(fac = rep(x, y))
pdf("$output_plots",width=12,height=6)
barplot(table(dat)/sum(table(dat)), col="lightblue", xlab="Length(nt)", ylab="Frequency", main="Length distribution")
invisible(dev.off())
END

	open R,"|/usr/bin/R --vanilla --slave" or die $!;
	print R $R_LD;
	close R;	

	# convert pdf file to png
	my $cmd_convert = "convert $output_plots $output_image";
	system($cmd_convert) && die "[ERR]CMD: $cmd_convert\n";
}

=head2
 norm -- normalization of sRNA dataset
=cut
sub norm
{
	my $options = shift;
	
	my $subUsage = qq'
USAGE $0 norm [options] 
        -i      list of UNIQUE read
        -o      perfix of output files

* the input file in list MUST be UNIQUE format of sRNA
  the UNIQUE formart include read count in ID.

* the output files
[perfix]_sRNA_expr	raw expression
[perfix]_sRNA_libsize	library size
[perfix]_sRNA_expTPM	normalized exp
[perfix]_sRNA_seq	unique sRNA

';
	print $subUsage and exit unless ($$options{'i'} && $$options{'o'});
	
	my $list_uniq_read = $$options{'i'};
	my $prefix = $$options{'o'};
	my $output1 = $prefix."_sRNA_expr";
	my $output2 = $prefix."_sRNA_libsize";
	my $output3 = $prefix."_sRNA_expTPM";
	my $output4 = $prefix."_sRNA_seq";

	# put list of uniq small RNA reads to array
	my @list;
	my $fh = IO::File->new($list_uniq_read) || die "Can not open list file $list_uniq_read $!\n";
	while(<$fh>)
	{
		chomp;
		push(@list, $_);
		die "[ERR]cat not find uniq read file $_\n" unless -s $_;
	}
	$fh->close;

	# main expression and libsize value to hash
	my %uniq_read;
	my %libsize;	

	foreach my $file (@list)
	{
        	my $total_read = 0;
        	my $fu;
		if ($file =~ m/\.gz$/) {
			open ($fu,'-|', "gzip -cd $$file") || die $!;	# discard gzip IO, method from honghe
        	} else {
                	open($fu, $file) || die $!;
        	}

		while(<$fu>)
		{
                	chomp;
	                my $id = $_;
        	        my @a = split(/-/, $id);
	                my $exp = $a[scalar(@a)-1];
        	        my $seq = <$fu>; chomp($seq);
	                $uniq_read{$seq}{$file} = $exp;
        	        $total_read = $total_read + $exp;
	        }
	        close($fu);
	        $libsize{$file} = $total_read;
	}

	# output the libsize
	my $out2 = IO::File->new(">$output2") || die "Can not open output libsize $output2 $!\n";
	foreach my $k (sort keys %libsize) { print $out2 $k."\t".$libsize{$k}."\n"; }
	$out2->close;

	# output the expression 
	my $out1 = IO::File->new(">$output1") || die "Can not open output expression $output1 $!\n";
	my $out3 = IO::File->new(">$output3") || die "Can not open output expression TPM $output3 $!\n";
	my $out4 = IO::File->new(">$output4") || die "Can not open output small RNA $output4 $!\n";

	print $out1 "#ID\tUniqRead";
	print $out3 "#ID\tUniqRead";
	foreach my $file (@list) {
        		print $out1 "\t".$file;
			print $out3 "\t".$file;
	}
	print $out1 "\n";
	print $out3 "\n";

	my $seq_order = 0;
	my $length = length(scalar(keys(%uniq_read)));

	foreach my $seq (sort keys %uniq_read)
	{
		$seq_order++;
		my $zlen = $length - length($seq_order);
		my $z = "0"x$zlen;
		my $seq_id = "sR".$prefix.$z.$seq_order;

		print $out4 ">".$seq_id."\n".$seq."\n";

		my $line = $seq_id."\t".$seq;
		my $line_tpm = $seq_id."\t".$seq;

		foreach my $file (@list)
		{
			if ( defined $uniq_read{$seq}{$file} )
                	{
                        	$line.="\t$uniq_read{$seq}{$file}";

	                        my $tpm = $uniq_read{$seq}{$file} / $libsize{$file} * 1000000;
        	                $tpm = sprintf("%.2f", $tpm);
                	        $line_tpm.="\t$tpm";
                	}
                	else
                	{
                        	$line.="\t0";
                        	$line_tpm.="\t0";
                	}
        	}

	        print $out1 $line."\n";
	        print $out3 $line_tpm."\n";
	}
	$out1->close;
	$out3->close;
	$out4->close;

}

=head2
 normcut -- normalization of sRNA cutoff
=cut
sub normcut
{
	my $options = shift;
	
	my $subUsage = qq'
USAGE $0 normcut [options]
	-i	input file prefix
	-c	cutoff (default: 10 TPM)

* the input file should be:
[perfix]_sRNA_expr
[perfix]_sRNA_expTPM
[perfix]_sRNA_seq

';

	print $subUsage and exit unless $$options{'i'};
	my $file_prefix = $$options{'i'};
	my $cutoff = 10;
	$cutoff = $$options{'c'} if (defined $$options{'c'} && $$options{'c'} > 0);

	my $expr = $file_prefix."_sRNA_expr";
	my $norm = $file_prefix."_sRNA_expTPM";
	my $srna = $file_prefix."_sRNA_seq";

	my $out_expr = $file_prefix."_".$cutoff."TPM_sRNA_expr";
	my $out_norm = $file_prefix."_".$cutoff."TPM_sRNA_expTPM";
	my $out_sRNA = $file_prefix."_".$cutoff."TPM_sRNA_seq";

	my $out1 = IO::File->new(">".$out_expr) || die $!;
	my $out2 = IO::File->new(">".$out_norm) || die $!;
	my $out3 = IO::File->new(">".$out_sRNA) || die $!;

	my %id;
	my $fh2 = IO::File->new($norm) || die $!;
	my $title = <$fh2>; print $out2 $title;
	while(<$fh2>)
	{
        	chomp;
 		my @a = split(/\t/, $_);

		my $select = 0;
		for(my $i=2; $i<@a; $i++) {
               		if ( $a[$i] > $cutoff ) {
                        	$select = 1;
                	}
        	}

        	if ($select == 1) {
                	$id{$a[0]} = 1;
                	print $out2 $_."\n";
			print $out3 ">$a[0]\n$a[1]\n";
        	}
	}
	$fh2->close;
	$out2->close;
	$out3->close;

	my $fh1 = IO::File->new($expr) || die $!;
	my $t = <$fh1>; print $out1 $t;
	while(<$fh1>)
	{
        	chomp;
		my @a = split(/\t/, $_);
		if ( defined $id{$a[0]} ) { print $out1 $_."\n"; }
	}
	$fh1->close;
	$out1->close;
}

=head2
 usage -- print usage information
=cut
sub usage
{
	print qq'

Program: sRNAtools (Tools for sRNA analysis)
Version: $version

USAGE: $0 <command> [options] 
Command: 
	rmadpter	remove adapter sequence
	convert		convert between table format and fastq/fasta format
	unique		convert between unique format and clean format
	norm	     	normalization (RPM)
	cutnorm		normalization cutoff	
	lengthd		length distribution of sRNA	

';
	exit;
}


