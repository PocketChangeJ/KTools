#!/usr/bin/perl

=head
 pathwayTool.pl -- tools for pathway analysis

 author: Yi Zheng
 2014-11-20 convert this script to pathwayTool.pl 
 2014-06-20 fix error in p value calculate
 2014-04-09 init
=cut

use strict;
use warnings;
use IO::File;
use Getopt::Std;

my %options;
getopts('a:b:c:d:e:f:g:i:j:k:l:m:n:o:p:q:r:s:t:u:v:w:x:y:z:h', \%options);
unless (defined $options{'t'} ) { usage(); }

if	($options{'t'} eq 'filter') { pwy_filter(@ARGV); }	# parse multi dataset
elsif	($options{'t'} eq 'unique') { pwy_uniq(@ARGV); }	# parse multi dataset
elsif	($options{'t'} eq 'enrich') { pwy_enrich(@ARGV); }
else	{ usage(); }

=head2
 pwy_filter: filter the pathway result for removing none-plant pathways
=cut
sub pwy_filter
{
	my $pathway_file = shift;

	my $pathway_kept   = $pathway_file.".kept";
	my $pathway_remove = $pathway_file.".remove";

	open(OUT1, ">".$pathway_kept) || die $!;
	open(OUT2, ">".$pathway_remove) || die $!;
	open(FH, $pathway_file) || die $!;
	while(<FH>)
	{
		chomp;
		next if $_ =~ m/^#/;
		my @a = split(/\t/, $_);
		if ($a[3] =~ m/\s+\((.+?)\)$/)
		{
			my $cc = $1;
			$cc =~ s/cytochrome c\) \(//;
			if($cc eq 'yeast' || $cc eq 'prokaryotic' || $cc eq 'obligate autotrophs' || $cc eq 'Gram-negative bacteria' ||
			   $cc eq 'Gram-positive bacteria' || $cc eq 'mammals' || $cc eq 'metazoan')
			{
				print OUT2 $_."\n";
			}
			else
			{
				print OUT1 $_."\n";
			}
		}
		else
		{
			print OUT1 $_."\n";
		}
	}
	close(FH);
	close(OUT1);
	close(OUT2);
}


=head2
 pwy_uniq: uniq the pathway report file
=cut
sub pwy_uniq
{
	my ($pathway_file) = shift;
	
	my $usage = qq'
USAGE: $0 -t unique pathway_file > uniq_pathway_file

';

	my %p;
	my $fh = IO::File->new($pathway_file) || die $!;
	while(<$fh>)
	{
        	chomp;
	        # MIN031488       Zeaxanthin epoxidase, chloroplastic     PWY-5945        zeaxanthin, antheraxanthin and violaxanthin interconversion
	        my @a = split(/\t/, $_);

	        if (defined $p{$a[2]}{'num'} ) {
	                $p{$a[2]}{'num'}++;
	        } else {
	                $p{$a[2]}{'num'} = 1;
	        }

	        $p{$a[2]}{'desc'} = $a[3];
	}
	$fh->close;

	foreach my $id (sort keys %p) {
        	print $id."\t".$p{$id}{'num'}."\t".$p{$id}{'desc'}."\n";
	}
}

=head2 
 pwy_enrich: pathway enrichment analysis
=cut
sub pwy_enrich
{
	my ($gene_list,$pathway_file) = @_;

	my $usage = qq'
USAGE: $0 -t enrich input_gene_list pathway_file

* the input gene list should be changed gene in DE analysis
* the pathway file should be output of pathways tools

';

	# load gene list (changed gene) to hash
	my %changed_gene;
	my $fh1 = IO::File->new($gene_list) || die $!;
	while(<$fh1>)
	{
		chomp;
		$changed_gene{$_} = 1;
	}
	$fh1->close;

	# load pathway to hash
	# key: pwy_id
	# value: pwy_name
	#
	# key: pwy_id
	# value: gene1 \t gene2 \t ... \t geneN
	# 
	# check pwy id and name uniq at same time
	my %pwy_name;
	my %pwy_gene;
	my %all_pwy_gene;
	my %all_pwy_changed_gene;

	my $fh2 = IO::File->new($pathway_file) || die $!;
	while(<$fh2>)
	{
		chomp;
		next if $_ =~ m/^#/;
		#gene_ID        gene_description        pathway_ID      pathway_name
		my @a = split(/\t/, $_);
		die "Error in line $_\n" unless scalar @a == 4;
		my ($gid, $g_desc, $pid, $p_name) = @a;

		# check pwy id and pwy name
		if (defined $pwy_name{$pid} )
		{
			die "Error in pwy $pid\n" if $p_name ne $pwy_name{$pid};
		}
		else
		{
			$pwy_name{$pid} = $p_name;
		}

		if (defined $pwy_gene{$pid})
		{
			$pwy_gene{$pid}.= "\t".$gid;
		}
		else
		{
			$pwy_gene{$pid} = $gid;
		}

		$all_pwy_gene{$a[0]} = 1;

		if ( defined $changed_gene{$a[0]} )
		{
			$all_pwy_changed_gene{$a[0]} = 1;
		}

	}
	$fh2->close;

	my $N = scalar(keys %all_pwy_gene);		# N: gene in all pathways
	my $n = scalar(keys %all_pwy_changed_gene);	# n: changed gene in pathways

	# uniq the gene in each pwy, then get pvalue of changed pathways
	my %uniq_gene;
	my @uniq_gene;
	my $temp_file = "temp_output";
	my $out1 = IO::File->new(">".$temp_file) || die $!;

	foreach my $pid (sort keys %pwy_gene)
	{
		%uniq_gene = ();
		my @gene = split(/\t/, $pwy_gene{$pid});

		foreach my $gid (@gene) { $uniq_gene{$gid} = 1; }
		@uniq_gene = sort keys %uniq_gene;

		my $M = scalar(@uniq_gene);		# M: gene in particular pathways
		my $x = 0;				# x: changed gene in particular pathways
		foreach my $gid (@uniq_gene) {
			if ( defined $changed_gene{$gid} ) {
				$x++;
			}
		}
	
		my $p_name = $pwy_name{$pid};

		# compute pvalue
		#################################################################################
		#  input format
		#  hypergeometric(N, n, M, x);
		#  N: gene in all pathways
		#  n: changed gene in pathways
		#  M: gene in particular pathways
		#  x: changed gene in particular pathways
		#
		#  check this link for confirmation:
		#  http://www.geneprof.org/GeneProf/tools/hypergeometric.jsp
		################################################################################

		if ($x > 0) {
			# the order should be N M n x
			# my $pvalue = hypergeometric($N, $n, $M, $x);
			my $pvalue = hypergeometric($N, $M, $n, $x);
			my $background = "$M out of $N genes";
			my $cluster = "$x out of $n genes";
			print $out1 "$pid\t$p_name\t$cluster\t$background\t$pvalue\n";
		}
	}
	$out1->close;

	# adjust p value to q value
	my $output_dat = $gene_list."_changed_pwy.table.txt";

	my $R_CODE =<< "END";
library(qvalue)
data<-read.delim(file="$temp_file", header=FALSE)
p<-data[,5]
qobj<-qvalue(p)
alldata<-cbind(data,qobj\$qvalue)
write.table(file="$output_dat", sep="\\t", alldata);
END

	#print $R_CODE; exit;
	open R,"|/usr/bin/R --vanilla --slave" or die $!;
	print R $R_CODE;
	close R;

	unlink($temp_file);


	# better to parse the output file format
}

#####################
###  Subroutines  ###
#####################

sub hypergeometric {
    my $n = $_[0]; # N Total number of genes in all the pathways
    my $np = $_[1];# M Total number of genes in a particular pathway
    my $k = $_[2]; # n Total number of changed genes (in the input list) from all the pathways
    my $r = $_[3]; # x total number of changed genes (in the input list) from the particular pathway
    my $nq;
    my $top;
    
    $nq = $n - $np;

    my $log_n_choose_k = lNchooseK( $n, $k );

    $top = $k;
    if ( $np < $k ) {
        $top = $np;
    }

    my $lfoo = lNchooseK($np, $top) + lNchooseK($nq, $k-$top);
    my $sum = 0;

    for (my $i = $top; $i >= $r; $i-- ) {
        $sum = $sum + exp($lfoo - $log_n_choose_k);

        if ( $i > $r) {
            $lfoo = $lfoo + log($i / ($np-$i+1)) +  log( ($nq - $k + $i) / ($k-$i+1)  )  ;
        }
    }
    return $sum;
}

sub lNchooseK {
    my $n = $_[0];
    my $k = $_[1];
    my $answer = 0;

    if( $k > ($n-$k) ){
        $k = ($n-$k);
    }

    for(my $i=$n; $i>($n-$k); $i-- ) {
        $answer = $answer + log($i);
    }

    $answer = $answer - lFactorial($k);
    return $answer;
}

sub lFactorial {
    my $returnValue = 0;
    my $number = $_[0];
    for(my $i = 2; $i <= $number; $i++) {
        $returnValue = $returnValue + log($i);
    }
    return $returnValue;
}


=head2
 usage: show usage information
=cut
sub usage
{
	my $usage = qq'
USAGE $0 -t [tool]

	filter	filter the pathway result to remove non-plant pathways
	unique	unique the pathway identified result
	enrich	enrichemnt analysis for gene list

Pipelines:
	\$ pathwayTool.pl -t filter input_pathway(the input file is generate by ptools)
	\$ pathwayTool.pl -t unique input_pathway.kept
	\$ pathwayTool.pl -t enrich gene_list input_pathway.kept

';

	print $usage;
	exit;
}

