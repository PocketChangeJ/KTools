#!/usr/bin/perl

use strict;
use warnings;
use diagnostics;

use IO::File;
use FindBin;

use GO::TermFinder;
use GO::AnnotationProvider::AnnotationParser;
use GO::OntologyProvider::OboParser;

use GO::TermFinderReport::Text;

use GO::Utils::File    qw (GenesFromFile);
use GO::Utils::General qw (CategorizeGenes);
use Getopt::Long;


my ($help, $annotationFile, $totalNum, $oboFile, $listFile, $pvalue_cutoff, $correction, $image_switch);
$correction = 'FDR';
$pvalue_cutoff = 0.05;
$oboFile = ${FindBin::RealBin}."/bin/gene_ontology_edit.obo";

GetOptions (
        "h|?|help"      => \$help,
        "i=s"   => \$listFile,
        "a=s"   => \$annotationFile,
        "n=s"   => \$totalNum,
        "o=s"   => \$oboFile,
        "p=s"   => \$pvalue_cutoff,
        "c=s"   => \$correction,
        "g"     => \$image_switch
);

=head2
 go_link_gene: link gene annotation to GO term according to sequence comparison
=cut
sub go_link_gene
{
my $usage =  qq'
perl $0 gene_blast parsed_go_mapping_file removed_go_term[option] > output

* 1. format of gene_blast file
* 2. one gene_ID will have more hit_ID, one result per line
* 3. including all gene_ID from blast results, if one gene 
* do not have blast result, it will print "GeneID \\t NA"
gene_ID \\t hit_ID

* format of parse_go_mapping_file;
* the hit ID may have no corresponding GO ID
hit_ID \\t GO_ID1 \\t GO_ID2 ... GO_IDN

* format of ouput files, if gene do not have GO ID, the only 
* gene ID will be prinit at the end
gene_ID \\t GO_ID1 # GO_ID2 # ... # GO_IDN
gene_ID \\n

This script will use a default list of removed GO term if user 
do not set the file.
Some GO terms are not related to plant species, and they should 
be removed from gene GO table;

';

my $input_gene_blast = shift || die $usage;
my $parsed_go_mapping = shift;
my $removed_go_term = shift;

# set default mapping file
unless ($parsed_go_mapping) {
	$parsed_go_mapping = "$FindBin::RealBin/mapping/uniport_idmapping_GO_07242013";
}
die "Can not find go mapping file: $parsed_go_mapping\n" unless -s $parsed_go_mapping;

# set default removed go term
unless ($removed_go_term) {
	$removed_go_term = "$FindBin::RealBin/bin/non_plant_GO";
	die "Can not find removed go term: $removed_go_term" unless -s $removed_go_term;
}

my %removed_go;
my $fhr = IO::File->new($removed_go_term) || die $!;
while(<$fhr>){
	chomp;
	next if $_ =~ m/^#/;
	my @a = split(/\t/, $_);
	$removed_go{$a[0]} = 1;
}
$fhr->close;
#print scalar(keys %removed_go), " of removed none-plant related GO\n";

# save hit id and GO to hash
# key: hit_id
# value: GO_ID
# * if one hit do not have GO ID, will not put it to hash
my $nn = 0;
my %hit_GO;
my $fh1;
if ($parsed_go_mapping =~ m/\.gz/)
{
        open($fh1, "<:gzip", $parsed_go_mapping) || die $!;
}
else
{
        open($fh1, $parsed_go_mapping) || die $!;
}

while(<$fh1>)
{
	chomp;
	my @a = split(/\t/, $_);

	my $hit_id = $a[0];

	my $go_id = "";
	for(my $i=1; $i<@a; $i++)
	{
		unless ( defined $removed_go{$a[$i]} ) {
			if ($a[$i] =~ m/GO:/) { $go_id.="\t".$a[$i]; }
		} else {
			$nn++;
		}
	}
	$go_id =~ s/^\t//;

	if ($go_id =~ m/GO:/)
	{
		if (defined $hit_GO{$hit_id})
		{
			$hit_GO{$hit_id}.= "\t".$go_id;
		}
		else
		{
			$hit_GO{$hit_id} = $go_id;	
		}
	}
}
$fh1->close;

#print "$nn of GO has been removed by non-plant_GO term\n";

# gene go id for each genes
# key: gene_ID
# value: GO_ID
#
# hash for all gene in blast result, if gene do not have blast hit ID, it will print NA behind geneID,
# key: gene_ID
# value: 1
my %gene_GO; my %gene_id;
my $fh2;
if ($input_gene_blast =~ m/\.gz/) 
{
	open($fh2, "<:gzip", $input_gene_blast) || die $!;
}
else 
{
	open($fh2, $input_gene_blast) || die $!;
}

while(<$fh2>)
{
	chomp;
	# gene_ID \t hit_ID
	my @a = split(/\t/, $_);

	if ($a[1] =~ m/tr\|/) {
		my @b = split(/\|/, $a[1]);
		$a[1] = $b[1];
	}

	#print "$_\t$a[1]\n";

	if (defined $hit_GO{$a[1]})
	{
		if ( defined $gene_GO{$a[0]} )
		{
			$gene_GO{$a[0]}.= "\t".$hit_GO{$a[1]};
		}
		else
		{
			$gene_GO{$a[0]} = $hit_GO{$a[1]};
		}
	}

	$gene_id{$a[0]} = 1;
}
$fh1->close;

# output result 
foreach my $gid (sort keys %gene_GO) 
{
	my $go_terms = $gene_GO{$gid};
	$go_terms =~ s/\t/#/ig;
	print $gid."\t".$go_terms."\n";
	delete $gene_id{$gid};
}

foreach my $gid ( sort keys %gene_id)
{
	print $gid."\n";
}	
}

=head2
 go_associate: generate go_associate file
=cut
sub go_associate
{
	my $usage  = qq'
USAGE : perl $0 
		-i input_file
		-o organism[default:Prometheus]
		-e gene_ontology_edit.obo[default: RealBin/bin/gene_ontology_edit.obo]

* input_file : all gene and corresponding GO information, format:
  geneID \\t GOID#GOID#..#GOID \\n

* the output file should named as:
  Prometheus_gene_GO
  Prometheus_associate_file

';

my ($input_gene_GO_file, $organism, $gene_ontology_file);
$gene_ontology_file 	= ${FindBin::RealBin}."/bin/gene_ontology_edit.obo";
$organism		= "Prometheus";

GetOptions(
	"i=s" => \$input_gene_GO_file,
	"o=s" => \$organism,
	"e=s" => \$gene_ontology_file
);

die $usage unless defined $input_gene_GO_file;
die "Error, in gene_ontology_edit.obo file\n" unless (-s $gene_ontology_file);
die "Error, the output file exist\n" if (-s $organism."_associate_file" || -s $organism."_gene_GO");

#################################################################
# main								#
#################################################################

# load gene and corresponding GO annotation to hash, then uniq them
# key: geneID
# value: GO_ID # GO_ID ... GO_ID # the GO_ID is uniq for each geneID
my %gene_GO = load_gene_GO($input_gene_GO_file);

# report
print scalar(keys(%gene_GO))." gene with annotation were loaded $!\n";

# check error
#foreach my $k (sort keys %gene_GO) { print $k."\t".$gene_GO{$k}."\n"; } die;
sub load_gene_GO
{
	my $input_gene_GO_file = shift;
	my %gene_GO;
	my $fh = IO::File->new($input_gene_GO_file) || die "Can not open input gene and GO file: $input_gene_GO_file $!\n";
	while(<$fh>)
	{
		chomp;
		my @a = split(/\t/, $_);
		##my @b = split(/#/, $a[1]);
		if (defined $gene_GO{$a[0]} ) {
			$gene_GO{$a[0]}.= "#".$a[1];
		} else {
			$gene_GO{$a[0]} = $a[1];
		}
	}
	$fh->close;

	my %uniq_gene_GO;
	foreach my $gid (sort keys %gene_GO)
	{
		if ($gene_GO{$gid})
		{
			my @a = split(/#/, $gene_GO{$gid});
			my %uniq_GO = (); foreach my $a (@a) { $uniq_GO{$a} = 1; }
			my $value = join("#", keys(%uniq_GO));
			$uniq_gene_GO{$gid} = $value;
		} 
		else
		{
			$uniq_gene_GO{$gid} = 1;
		}
	}

	return %uniq_gene_GO;
}

# load GO_cat info to hash from gene_ontology_edit.obo file
# key: GO_ID
# value: P or F or C
my %GO_cat = load_cat_GO($gene_ontology_file);

# report
print scalar(keys(%GO_cat))." number of GO and namespace were loaded into hash\n";

# check result
# foreach my $k (sort keys %GO_cat) { print $k."\t".$GO_cat{$k}."\n"; } die;

sub load_cat_GO
{
	my $gene_ontology_file = shift;

	my ($id, $namespace, $alt_id, $is_obsolete);

	my $tmp = IO::File->new(">TTTEEEMMMPPP") || die "Can not open temp file TTTEEEMMMPPP $!\n";
	my $fh = IO::File->new($gene_ontology_file) || die "Can not ontology file $gene_ontology_file $!\n";
	while(<$fh>)
	{
		chomp;
		print $tmp "\n$_" if /^id: GO/;
		print $tmp "\t$_" if /^namespace/;
		print $tmp "\t$_" if /^alt_id/;
		print $tmp "\t$_" if /^is_obsolete/;
	}
	$fh->close;
	$tmp->close;

	my %cat = ("biological_process" => "P", "molecular_function" => "F", "cellular_component" => "C");

	my %GO_cat;

	$tmp = IO::File->new("TTTEEEMMMPPP") || die "Can not open temp file TTTEEEMMMPPP $!\n";
	<$tmp>;
	while(<$tmp>)
	{
        	chomp;
        	my @a = split(/\t/, $_);

		die "Error in GO ID format: $a[0]\n" unless $a[0] =~ m/id: /;
		die "Error in GO namespace format: $a[1]\n" unless $a[1] =~ m/namespace: /;
        	$a[0] =~ s/id: //; 
		$a[1] =~ s/namespace: //;

        	if ($_ = m/is_obsolete/)
        	{
                	$GO_cat{$a[0]} = "obsolete";

                	for (my $i = 2; $i < scalar(@a); $i++)
                	{
                        	if ($a[$i] =~ /alt_id/)
                        	{
                                	$a[$i] =~ s/alt_id: //;
                                	$GO_cat{$a[$i]} = "obsolete";
                        	}
                	}
        	}
        	else
        	{
			die "Error, can not identify namespace for $a[0] : $a[1]\n" unless defined $cat{$a[1]};
                	$GO_cat{$a[0]} = "$a[0]//$cat{$a[1]}";
                	for (my $i = 2; $i < scalar(@a); $i++)
                	{
                        	if ($a[$i] =~ /alt_id/)
                        	{
                                	$a[$i] =~ s/alt_id: //;
                                	$GO_cat{$a[$i]} = "$a[0]//$cat{$a[1]}";
                        	}
                	}
        	}
	}
	$tmp->close;

	unlink("TTTEEEMMMPPP");

	return %GO_cat;
}


#########################################################################
# Produce gene_GO files and temp3 file
# temp3: File for produce Go_asso_file next step 
# gene_GO: unigene id and go id
#########################################################################

my $output_gene_GO = $organism."_gene_GO";
my $out1 = IO::File->new(">".$output_gene_GO) || die "Can not open gene and GO file $output_gene_GO $!\n";

my @temp_gene_with_GO;

foreach my $gid (sort keys %gene_GO)
{
	# skip gene without GO annotation
	if ($gene_GO{$gid} eq "1") { next; }

  	my @go_member = split(/#/, $gene_GO{$gid});

	# uniq the GO ID for each gene
	my %unique_go = ();
	foreach my $elem (@go_member) { $unique_go{$elem} = 1; }
	my @unique_go = sort keys %unique_go;

	# adjust GO member according to obsolete
	my %GO_unique = ();
	for(my $i=0; $i<@unique_go; $i++) 
	{
		my $go_id = $unique_go[$i];
		print "Warning, the GO ID ($go_id) does not defined in gene_ontology_edit.obo file\n" unless defined $GO_cat{$go_id}; 

 		if ($GO_cat{$go_id} eq "obsolete")
    		{
			#if the go id is obsolete, make it blank
      			$unique_go[$i] = "";
    		}
    		else
		{
			#die $go_member[$i];
			$unique_go[$i] = $GO_cat{$go_id};
			if ( $unique_go[$i] !~ /GO:0008150/ && $unique_go[$i] !~ /GO:0005575/ && $unique_go[$i] !~ /GO:0003674/ && 
			     $unique_go[$i] =~ /GO/ &&
			     !defined $GO_unique{$unique_go[$i]} )
      			{
        			my $gene_GO = $unique_go[$i];
        			$gene_GO =~ s/\/\/.*//;
        			print $out1 "$gid\t$gene_GO\t$organism\n";
      			}
      			else
      			{
        			$unique_go[$i] = "";
      			}
      			
			$GO_unique{$unique_go[$i]} = 1;
		}
	}

	# uniq the GO info in each gene and put it to temp hash
	%unique_go = ();
        foreach my $elem ( @unique_go ) { $unique_go{$elem} = 1 if $elem; }
        my @unique1 = sort keys %unique_go;
        if ( scalar(@unique1) > 0 ) {
                my $value = join("#", @unique1);
		$value =~ s/^#//ig;
		push(@temp_gene_with_GO, "$gid\t$value");
        }
}
$out1->close;

print "No. of gene has GO information: ", scalar(@temp_gene_with_GO), "\n";

#################################################################
# Produce GO_associate_file					#
# the file content is below.					#	
# Description:							#
# 	1. For every unigene, it must have 3 main cat		#
# 	2. Except 1 GOID, it has it's own GOID			#
#################################################################

# generate date
my ($day, $mon, $year) = (localtime)[3..5];
my $date = sprintf "%04d-%02d-%02d", 1900+$year, 1+$mon, $day;


# record how many gene has GO annotation;
# key: geneID
# value: 1
my %exist = ();

my $output_associate_file = $organism."_associate_file";
my $out2 = IO::File->new(">".$output_associate_file) || die "Can not open output associate file $output_associate_file $!\n";

foreach my $t (@temp_gene_with_GO)
{
	# make every unigene have 3 main cat and GO ID
	chomp($t);
	my @a = split(/\t/, $t);
	$exist{$a[0]} = 1;

	# add complementary category information to each gene,
	# for example: 
	# if geneA just has a GO that belong to F, we will add GO:0008150 and GO:0005575 on it in GO_associate_file
  	if ($a[1] !~ /F/)
	{ print $out2 "$organism\t$a[0]\t\t\tGO:0003674\tUniProt\tIEA\t\tF\t$a[0]\t\tgene\tTaxon:0000\t$date\tKNF\n";  }
	if ($a[1] !~ /P/)
  	{ print $out2 "$organism\t$a[0]\t\t\tGO:0008150\tUniProt\tIEA\t\tP\t$a[0]\t\tgene\tTaxon:0000\t$date\tKNF\n"; }
	if ($a[1] !~ /C/)
	{ print $out2 "$organism\t$a[0]\t\t\tGO:0005575\tUniProt\tIEA\t\tC\t$a[0]\t\tgene\tTaxon:0000\t$date\tKNF\n"; }

 	#print it's own GO ID
	my @go_member = split("#", $a[1]);
	for (my $i = 0; $i < @go_member; $i++)
	{
		if ($go_member[$i] =~ m/GO/)
 		{
			die "Error in format of GO info: $go_member[$i]\n" unless $go_member[$i] =~ m/\/\//;
      			$go_member[$i] =~ s/\/\//#/;
			my ($go_id, $go_cat) = split("#", $go_member[$i]);
			die "Error in format of GO info: $go_id, $go_cat\n" unless ($go_id && $go_cat);
      			print $out2 "$organism\t$a[0]\t\t\t$go_id\tUniProt\tIEA\t\t$go_cat\t$a[0]\t\tgene\tTaxon:0000\t$date\tKNF\n";
    		}
		else
		{
			die "Error in line $t\nIt has incorrect GO information";
		}
  	}
}
print scalar(keys(%exist))." gene has GO annotation\n";

#produce three main cat for each unigene:
#ID file have all unigene id;  hash %exist ids that have content; so !$exist is the left ids;
my $un_go_num = 0;

foreach my $gid (sort keys %gene_GO)
{
	if (!$exist{$gid})
  	{
    		print $out2 "$organism\t$gid\t\t\tGO:0003674\tUniProt\tIEA\t\tF\t$gid\t\tgene\tTaxon:0000\t$date\tKNF\n";
    		print $out2 "$organism\t$gid\t\t\tGO:0008150\tUniProt\tIEA\t\tP\t$gid\t\tgene\tTaxon:0000\t$date\tKNF\n";
    		print $out2 "$organism\t$gid\t\t\tGO:0005575\tUniProt\tIEA\t\tC\t$gid\t\tgene\tTaxon:0000\t$date\tKNF\n";
		$un_go_num++;
  	}

}
print "$un_go_num gene do not have GO annotation\n";

$out2->close;
}

=head2
 go_slim: perform GO slim analysis
=cut
sub go_slim
{
	my ($options, $usage) = @_;


my $usage = qq'
usage: perl $0 [options]
       perl $0 -a input_asso
       perl $0 -i listID -g background_asso

  -a input[GO_associate_file]
  -i list_of_IDs			[ChooseGene]
  -g GO_associate_file for list of IDs.	[Background]
  -o output_prefix [default: GSA]
  -s slim_file [defalut is goslim_plant.obo]

* -a is used for gene category base on plant goslim [required if no -i]
* -i and -g will generate enrichment p-value for each GO term [required if no -a]
  the list of ID must be a subset of background GO_associate_file.
* -s will use the default one [required if -i]

';

my ($help, $associate_in, $listID, $associate_list, $output_prefix, $plant_go_slim, $gene_ontology_obo);
GetOptions (
	"h|?|help"	=> \$help,	
	"a=s" 	=> \$associate_in,	 
	"i=s"	=> \$listID,	
	"g=s"	=> \$associate_list,
	"o=s"	=> \$output_prefix, 
	"s=s"	=> \$plant_go_slim,
	"n=s"	=> \$gene_ontology_obo
);

#################################################################
# init setting and checking input files				#
#################################################################

die $usage if $help;
$output_prefix ||= "GSA";

my ($out_slim, $out_slim_tab, $out_slim_all, $out_slim_all_tab);
if ($associate_in) { 
	$out_slim = $output_prefix."_goslim";
	$out_slim_tab = $output_prefix."_goslim_tab";
} elsif ($listID && $associate_list ) { 
	$out_slim = $output_prefix."_goslim";
	$out_slim_tab = $output_prefix."_goslim_tab";
	$out_slim_all = $output_prefix."_goslim_all";
	$out_slim_all_tab = $output_prefix."_goslim_all_tab";
} else { die $usage; }

$plant_go_slim ||=  ${FindBin::RealBin}."/bin/goslim_plant.obo";
$gene_ontology_obo ||= ${FindBin::RealBin}."/bin/gene_ontology_edit.obo";

unless (-s $gene_ontology_obo) { die "Error, can not find gene ontology obo file: $gene_ontology_obo\n"; }
unless (-s $plant_go_slim) { die "Error, can not find plant go slim file: $plant_go_slim\n"; }

my $map2slim_script = ${FindBin::RealBin}."/bin/map2slim";
unless(-s $map2slim_script) { die "Error in $!\n"; }

#################################################################
# main								#
#################################################################

if ($associate_in)
{	
	go_slim_analysis($associate_in, $plant_go_slim, $gene_ontology_obo, $out_slim, $out_slim_tab);
}
elsif ($listID && $associate_list)
{
	# create associate_in base on $listID and associate_list;
	$associate_in = "input_associate";
	
	# num_a: num of genes for list 
	# num_b: num of genes for all
	my ($num_a, $num_b) = generate_associate_byID($listID, $associate_list, $associate_in);
	go_slim_analysis($associate_in, $plant_go_slim, $gene_ontology_obo, $out_slim, $out_slim_tab);
	go_slim_analysis($associate_list, $plant_go_slim, $gene_ontology_obo, $out_slim_all, $out_slim_all_tab);
	compute_pvalue($num_a, $num_b, $out_slim_tab ,$out_slim_all_tab);

	unlink($associate_in);
}
else
{
	die "Error in main $usage\n";
}

#################################################################
# kentnf: subroutine						#
#################################################################
=head1 generate_associate_byID
 generate_associate_byID
=cut
sub generate_associate_byID
{
	my ($listID, $associate_list, $associate_in) = @_;

	my ($num_a, $num_b); my %all_ID;

	my %ID;
	my $fh = IO::File->new($listID) || die "Can not open list IDs file $listID $!\n";
	while(<$fh>)
	{
        	chomp;
        	$ID{$_} = 1;
	}
	$fh->close;

	my $out = IO::File->new(">".$associate_in) || die "Can not open file $associate_in\n";
	my $in  = IO::File->new($associate_list)   || die "Can not open file $associate_list\n";
	while(<$in>)
	{
		chomp;
		my @a = split(/\t/, $_);
		if (defined $ID{$a[1]})
		{
			print $out $_."\n";
		}

		unless (defined $all_ID{$a[1]}) { $num_b++; $all_ID{$a[1]} = 1; }

	}
	$in->close;
	$out->close;

	$num_a = scalar(keys(%ID));

	return ($num_a, $num_b);
}


=head1 go_slim_analysis 
 go_slim_analysis
=cut
sub go_slim_analysis
{
	my ($associate_in, $plant_go_slim, $gene_ontology_obo, $out_slim, $out_slim_tab) = @_;

	# perform GO slim analysis using map2slim script
	my $temp_file = "TTTMMMPPP";
	my $cmd_map2slim = "$map2slim_script $plant_go_slim $gene_ontology_obo $associate_in > $temp_file";
	print $cmd_map2slim."\n";
	system($cmd_map2slim) && die "Error in command: $cmd_map2slim\n";

	my %GO_gid;	# key: GO ID \t gene id; value: 1;
	%GO_gid = map2slim_result_to_hash($temp_file, $out_slim);
	unlink($temp_file);

	#########################################################
	# convert GO_gid hash to GO_num hash			#
	#########################################################
	my %GO_num;     # key: GO ID; value: number of genes
	foreach my $gg (sort keys %GO_gid) {
        	my ($gene_id, $go_id) = split(/\t/, $gg);
        	if (defined $GO_num{$go_id}) {
                	$GO_num{$go_id}++;
        	} else {
                	$GO_num{$go_id} = 1;
        	}
	}

	#########################################################
	# put go slim annotation to hash			#
	#########################################################
	my %GO_name;            # key: GO ID; value: GO name
	my %GO_namespace;       # key: GO ID; value: GO namespace
	my ($GO_name, $GO_namespace) = load_plant_go_slim($plant_go_slim);

	%GO_name = %$GO_name; %GO_namespace = %$GO_namespace;

	foreach my $a (sort keys %GO_name)
	{
		#print "$a\t$GO_name{$a}\n";
	}


	#########################################################
	# print the output result				#
	#########################################################
	my $out = IO::File->new(">$out_slim_tab") || die "Can not open file $out_slim_tab $!\n";
	foreach my $id (sort keys %GO_num)
	{
		if ($GO_num{$id} && $GO_namespace{$id} && $GO_name{$id})
		{
        		print $out $id."\t".$GO_num{$id}."\t".$GO_namespace{$id}."\t".$GO_name{$id}."\n";
		}
		else
		{
			print "Error\t$id\t$GO_num{$id}\n";
		}
	}
	$out->close;

}

=head1 
 map2slim_result_to_hash -- map2slim output to hash				
 output gene and corresponding GO info
=cut
sub map2slim_result_to_hash
{
	my ($map2slim_result, $output) = @_;;

	my %GO_gid = ();

	my ($go_id, $gene_id);

	my $out = IO::File->new(">$output") || die "Can not open output goslim: $output $!\n";
	my $fh = IO::File->new($map2slim_result) || die "Can not open map2slim result file: $map2slim_result\n";
	while(<$fh>)
	{
		chomp;
		my @a = split(/\t/, $_);
		($gene_id, $go_id) = ($a[1], $a[4]);

		# prase the output of map2slim
		if ($go_id ne "GO:0008150" && $go_id ne "GO:0005575" && $go_id ne "GO:0003674")
		{
			unless (defined $GO_gid{$gene_id."\t".$go_id})
			{
				$GO_gid{$gene_id."\t".$go_id} = 1;
				print $out $gene_id."\t".$go_id."\n";
			}
		}
	}
	$fh->close;
	$out->close;

	return %GO_gid;
}

=head1 
 load_plant_go_slim -- load plant go slim info to hash
 %go_name	# key: GO ID; value: GO name
 %go_namespace	# key: GO ID; value: GO namespace
=cut
sub load_plant_go_slim
{
	my $plant_go_slim = shift;

	my %GO_name; my %GO_namespace;

	my ($i, $j, $id, $name, $namespace);

	my $fh = IO::File->new($plant_go_slim) || die "Can not open plant GO slim file $plant_go_slim\n";
	while(<$fh>)
	{
		chomp;
		if ($_ =~ m/\[Term\]/ || $_ =~ m/\[Typedef\]/ )
		{
			if ($id && $name && $namespace)
			{
				$GO_name{$id} = $name;
				$GO_namespace{$id} = $namespace;
				$i++;
			}
		
			if ($_ =~ m/\[Term\]/) { $j++; }
			$id = ""; $name = ""; $namespace = "";
		}
		elsif ($_ =~ m/^id:\s(.*)/)
		{
			$id = $1;
		}
		elsif ($_ =~ m/^name:\s(.*)/)
		{
			$name = $1;
		}
		elsif ($_ =~ m/^namespace:\s(.*)/)
		{
			$namespace = $1;
		}
		else
		{
			next;
		}
	}
	$fh->close;

	# chech the number GO ID and GO name
	if ($i == $j) {
		#print "There are $i annotation in GO annotation $plant_go_slim\n";
	} else {
		print "There are $i annotation in GO annotation $plant_go_slim\n".
		      "But there are $j term in GO annotation $plant_go_slim\n";
	}

	return (\%GO_name, \%GO_namespace);
}

=head2 
 compute_pvalue: compute p-value for enrichment analysis
=cut
sub compute_pvalue
{
	my ($num_a, $num_b, $tab_a, $tab_b) = @_;

	# num in tab_b to hash
	my %hash_b;
	my $fh = IO::File->new($tab_b) || die "Can not open goslim table all $tab_b $!\n";
	while(<$fh>)
	{
		chomp;
		my @a = split(/\t/, $_);
		$hash_b{$a[0]} = $a[1];
	}
	$fh->close;

	my $tab_t = $tab_a."_temp";
	my $out = IO::File->new(">".$tab_t) || die "Can not open temp goslim table for list with p-value $tab_t $!\n";
	my $in  = IO::File->new($tab_a) || die "Can not open goslim table for list $tab_a $!\n";
	while(<$in>)
	{
		chomp;
		my ($go_id, $num_a_sub, $namespace, $name) = split(/\t/, $_);

		my ($population, $good, $bad, $sample, $select);

		$population = $num_b;
		$good = $hash_b{$go_id};
		$bad  = $num_b - $hash_b{$go_id};

		$sample = $num_a;
		$select = $num_a_sub;

		my $pvalue = hypergeom1($good, $bad, $sample, $select);

		print $out "$go_id\t$num_a_sub\t$num_a\t$hash_b{$go_id}\t$num_b\t$pvalue\t$namespace\t$name\n";
	}
	$in->close;
	$out->close;

	# adjust P value to Q value using R
	my $R_QV =<< "END";
a<-read.table("$tab_t", header = FALSE, sep = "\t")
p<-a[,6]
q<-p.adjust(p, method="BH", n=length(p))
write.table(q, file="qvalue.temp", sep = "\t", row.names = FALSE, col.names = FALSE)
END

	open R,"|/usr/bin/R --vanilla --slave" or die $!;
	print R $R_QV;
	close R;

	# combine the table into on
	my $tab_p = $tab_a."_pvalue";
	my $outp = IO::File->new(">".$tab_p) || die "Can not open goslim table for list with p-value $tab_p $!\n";
	my $in1 = IO::File->new($tab_t) || die "Can not open temp goslim table for list with p-value $tab_t $!\n";
	my $in2 = IO::File->new("qvalue.temp") || die "Can not open temp q-value $!\n";
	while(<$in1>)
	{
		chomp;
		my ($go_id, $num_a_sub, $a, $num_b_sub, $b, $pvalue, $namespace, $name) = split(/\t/, $_);
		my $qvalue = <$in2>;
		chomp($qvalue);
		unless($qvalue) { die "Error at qvalue file qvalue.temp $!\n"; }
		print $outp "$go_id\t$num_a_sub\t$a\t$num_b_sub\t$b\t$pvalue\t$qvalue\t$namespace\t$name\n";
	}
	$in1->close;
	$in2->close;
	$outp->close;
	unlink($tab_t);
	unlink("qvalue.temp");
}

=head2
 hypergeom1: hypergeometric distribution
=cut
sub hypergeom1 
{

 	# There are m "bad" and n "good" balls in an urn.
	# Pick N of them. The probability of i or more successful selections:
	# (m!n!N!(m+n-N)!)/(i!(n-i)!(m+i-N)!(N-i)!(m+n)!)
	# $m+n pop  $n successful
	# $N sample $i successful
	my ($n, $m, $N, $i) = @_;
    	my $loghyp1 = logfact1($m) +logfact1($n)+logfact1($N)+logfact1($m+$n-$N);
	my $loghyp2 = logfact1($i)+logfact1($n-$i)+logfact1($m+$i-$N)+logfact1($N-$i)+logfact1($m+$n);
	return exp($loghyp1 - $loghyp2);
}

sub logfact1 
{
	my $x = shift;
	my $ser = (   1.000000000190015
                + 76.18009172947146   / ($x + 2)
                - 86.50532032941677   / ($x + 3)
                + 24.01409824083091   / ($x + 4)
                - 1.231739572450155   / ($x + 5)
                + 0.12086509738661e-2 / ($x + 6)
                - 0.5395239384953e-5  / ($x + 7) );
	my $tmp = $x + 6.5;
	($x + 1.5) * log($tmp) - $tmp + log(2.5066282746310005 * $ser / ($x+1));
}

=head2
 go_enrich : enrichment analysis 
=cut

sub go_enrich
{
	my ($options, $files) = @_;

	my $usage = qq'
usage: perl $0 -t enrich [options] annotationFile listFile1 listFile2 ... listFileN

  -n	total number of genes ( automatically get this number from annotationFile)
  -o	gene_ontology_edit.obo ( will load it automatically )
  -c	correction method [bonferroni, none, simulation, FDR] ( default: FDR )
  -p	adjusted/raw p-value cutoff ( default: 0.05 )
  -g	disable/enable image output ( default: disable )
  -h|?|help help information

';
    exit;
}

	# checking input parameters


Usage if $help;
Usage unless $listFile;

	if ($correction eq "FDR" || $correction eq "none" || $correction eq "bonferroni" || $correction eq "simulation" ) {}
	else { Usage("Your correction method is not correct."); }

	# now check annotation file and get number of genes
	unless ($totalNum){ $totalNum = get_total_gene_num($annotationFile); }

	if ($oboFile !~ /\.obo$/){ Usage("Your obo file does not have a .obo extension."); }
	else { unless (-s $oboFile) { Usage("Can not find your obo file."); } }

	#################################################################
	# now set up the objects we need				#
	#################################################################
	my $process   = GO::OntologyProvider::OboParser->new(ontologyFile => $oboFile,
							     aspect       => 'P');
	my $component = GO::OntologyProvider::OboParser->new(ontologyFile => $oboFile,
							     aspect       => 'C');
	my $function  = GO::OntologyProvider::OboParser->new(ontologyFile => $oboFile,
							     aspect       => 'F');

	my $annotation = GO::AnnotationProvider::AnnotationParser->new(annotationFile=>$annotationFile);

	my $termFinderP = GO::TermFinder->new(annotationProvider=> $annotation,
					      ontologyProvider  => $process,
					      totalNumGenes     => $totalNum,
					      aspect            => 'P');

	my $termFinderC = GO::TermFinder->new(annotationProvider=> $annotation,
					      ontologyProvider  => $component,
					      totalNumGenes     => $totalNum,
					      aspect            => 'C');

	my $termFinderF = GO::TermFinder->new(annotationProvider=> $annotation,
					      ontologyProvider  => $function,
					      totalNumGenes     => $totalNum,
					      aspect            => 'F');

	my $report = GO::TermFinderReport::Text->new();

	# now go through each file
	my @listFiles = ($listFile);

	# foreach my $file (@listFiles)
	# {
    		print "Analyzing $file\n";
		my @genes = GenesFromFile($file); 

		my (@list, @notFound, @ambiguous);
		CategorizeGenes(annotation  => $annotation,
		    		genes       => \@genes,
		 		ambiguous   => \@ambiguous,
		    		unambiguous => \@list,
		    		notFound    => \@notFound);

		my $outfile = $file.".terms";
		my $outtable = $file.".table";

		my $fh = IO::File->new($outfile, q{>} ) || die "Cannot make $outfile : $!";

		print "Results being put in $outfile and $outtable\n";

		if (@list)
		{
			print $fh "The following gene(s) will be considered:\n\n";
			foreach my $gene (@list){ print $fh $gene, "\t", $annotation->standardNameByName($gene), "\n"; }
			print $fh "\n";
		}else{
			print $fh "None of the gene names were recognized\n";
			print $fh "They were:\n\n";
			print $fh join("\n", @notFound), "\n";
			$fh->close;
			next;
		}

		if (@ambiguous){
			# note, some of these ambiguous names would be perfectly fine
			# if put into GO::TermFinder if they are also standard names.
			# Currently the behavior of analyze.pl differs from the
			# default behavior of GO::TermFinder

		print $fh "The following gene(s) are ambiguously named, and so will not be used:\n";
		print $fh join("\n", @ambiguous), "\n\n";
	    	}

		if (@notFound){
			print $fh "The following gene(s) were not recognized, and will not be considered:\n\n";
			print $fh join("\n", @notFound), "\n\n";
		}

	
	    	foreach my $termFinder ($termFinderP, $termFinderC, $termFinderF)
		{
			# it's possible that the supplied number of genes on the
			# command line was less than indicated by the annotation
			# provider, and thus the TermFinder may have used a larger
			# number than was entered on the command line.

		my $totalNumGenesUsedInBackground = $termFinder->totalNumGenes;

		print $fh "Finding terms for ", $termFinder->aspect, "\n\n";

		my @pvalues;
		
		if ($correction eq "FDR") {
			@pvalues = $termFinder->findTerms(genes => \@list, calculateFDR => 1);
			print "FDR working\n";
		} else {
			@pvalues = $termFinder->findTerms(genes => \@list, correction => $correction, calculateFDR => 0);
			print "FDR disabled\n";
		}

		my $numHypotheses = $report->print(pvalues  => \@pvalues,
						   numGenes => scalar(@list),
						   totalNum => $totalNumGenesUsedInBackground,
						   cutoff   => $pvalue_cutoff,
						   fh       => $fh);

		# if they had no significant P-values

		if ($numHypotheses == 0)
		{
			print $fh "No terms were found for this aspect with a corrected P-value <= $pvalue_cutoff.\n";
		}

		print $fh "\n\n";
    	}
    	$fh->close;    

	parse_result($outfile, $outtable);
}


=head1 
 parse_result : parse go enrichment result 
=cut
sub parse_result
{
	my ($input, $output) = @_; 

	my $out = IO::File->new(">".$output) || die "Can not open go enrichment table file $output\n";
	my $in  = IO::File->new($input) || die "Can not open go enrichment term file $input\n";

	print $out "#GO ID\tGO term\tT\tCluster frequency\tBackground frequency\tRaw P-value\tCorrected P-value\tGenes annotated to the term\n";

	my ($aspect, $go_id, $term, $corrected_p, $raw_p, $m1, $m2, $m3, $m4, $freq1, $freq2, $genes, $cluster_freq, $background_freq);
	while(<$in>)
	{
		chomp;
		if ($_ =~ m/^Finding terms for (\D)/) 
		{
			$aspect = $1;
		}
		elsif ($_ =~ m/^GOID\s+(\S+)/)
		{
			$go_id = $1;
		}
		elsif ($_ =~ m/^TERM\s+(.*)/)
		{
			$term = $1;
		}
		elsif ($_ =~ m/^CORRECTED P-VALUE\s+(\S+)/)
		{
			$corrected_p = $1;
		}
		elsif ($_ =~ m/^UNCORRECTED P-VALUE\s+(\S+)/)
		{
			$raw_p = $1;
		}
		elsif ($_ =~ m/^NUM_ANNOTATIONS\s(\d+)\sof\s(\d+)\sin\sthe\slist\,\svs\s(\d+)\sof\s(\d+)/)
		{
			$m1 = $1; $m2 = $2; $m3 = $3; $m4 = $4;
			$freq1 = sprintf("%.2f", ($m1/$m2)); $freq1 = $freq1*100; $freq1.="%";
			$freq2 = sprintf("%.2f", ($m3/$m4)); $freq2 = $freq2*100; $freq2.="%";
		}
		elsif ($_ =~ m/^The genes annotated to this node are/)
		{
			$genes = <$in>; chomp($genes);
			$cluster_freq = "$m1 out of $m2 genes, $freq1";
			$background_freq = "$m3 out of $m4 genes, $freq2";
			print $out $go_id."\t".$term."\t".$aspect."\t".$cluster_freq."\t".$background_freq.
				   "\t".$raw_p."\t".$corrected_p."\t".$genes."\n";
		}
		else
		{
			next;
		}
	}
	$in->close;
	$out->close;
}

=head1 
 get_total_gene_num: get total gene number from input associate file
=cut
sub get_total_gene_num
{
	my $associate = shift;
	my %gid;
	my $fh = IO::File->new($associate) || die "Can not open input associate file: $associate \n";
	while(<$fh>)
	{
	    chomp;
	    unless ($_ =~ m/^!/)
	    {
		my @a = split(/\t/, $_);
		$gid{$a[1]} = 1;
	    }
	}
	$fh->close;	
	my $totalNum = scalar(keys(%gid));
	return $totalNum;
}

=head2
 usage -- print usage information
=cut
sub usage
{

}

=head2
 pipeline --
=cut
sub pipeline
{
	my $pipeline = qq'
update gene_ontology_edit.obo 02-24-2014

* There are two method to work with gene ontology obo file
        1. download Filtered ontology; the cross-products, inter-ontology links,
           and has_part relationships removed
        2. remove the links of parent with different namespace using remove_multi_parent.pl

GO Analysis pipelines

Place Plants GO associate file into this place

1. convert sequence aligment to GO annotation.

   due to various format of files, we use the simple method to perform this work with 3 steps

   1.1 sequence alignment

        using query sequences to search against databases such as UniPort, Pfam, InterProScan ...
        then parse the result to generate tab delimit file with two columns: seq_id \t hit_id1, hit_id2, ...

        parse blast alignment (UniProt)
       \$ seqTools/parse_blast.pl input_blast_output 5 > output

        parse HMMER3 alignment (Pfam Database)
       \$ seqTools/parse_hmmer.pl input_hmmer_output 5 > output

   1.2 generate gene and go annotation info

        \$ go_link_gene.pl parsed_gene_hit parsed_go_mapping_file

        parsed_gene_hit: gene_id \t hit_id \t ...
        parsed_go_mapping_file: hit_id \t GO id \t GO id .... \t GO id
        output: gene_id \t GO_ID # GO_ID # ...

        If you have parsed_gene_hit from different source, swissprot, trembl, pfam ....
        just combine them using cat

   1.3 convert gene and corresponding GO ID info to GO associate file using below method

        \$ go_generate_associate.pl input_gene_GO organism

        output: organism_associate_file, organism_GO

2. GOenrich -- GO enrichment analysis.

        \$ go_enrichment.pl -i list_gene_id -a associate_file

        usage: perl go_enrichment.pl [options]

          -i    listFile
          -a    annotationFile ( Usually associate file all genes )
          -n    total number of genes ( automatically get this number from annotationFile)
          -o    gene_ontology_edit.obo ( will load it automatically )
          -c    correction method [bonferroni, none, simulation, FDR] ( default: FDR )
          -p    adjusted/raw p-value cutoff ( default: 0.05 )
          -g    disable/enable image output ( default: disable )
          -h|?|help help information

        * if your associate file are generated according to gene_ontology_edit.obo version 1.2,
        * please use corresponding gene_ontology_edit.obo to perform GO slim and GO enrichment

3. GOslim -- category genes base on GO slim.

        \$ go_slim.pl -a input_associate_file

        output:


=== A Install GO::Parser ===

  The GOslim is base on the script map2slim, and the map2slim is
  base on GO::Parser, one of go-perl module. But the go-perl 0.14
  has error when installing by CPAN. So install 0.13 is the first
  step of GOslim.

  tar -xzvf go-perl-0.13.tar.gz
  cd go-perl-0.13
  perl Makefile.PL
  make test
  make install

  If you do not have root account. Please check the inlstall file
  on go-perl-0.13 folder. Install it on your home folder....

=== B Download GO annotation and mapping files ===

  Links for file download

  1. pfam
  http://www.geneontology.org/external2go/pfam2go

  2. interpro
  http://www.geneontology.org/external2go/interpro2go

  3. UniProt
  ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz

  Links for GO

  Full GO
  http://www.geneontology.org/ontology/obo_format_1_2/gene_ontology_ext.obo

  Filtered GO
  http://www.geneontology.org/ontology/obo_format_1_2/gene_ontology.1_2.obo

  Plant GO slims
  http://www.geneontology.org/GO_slims/goslim_plant.obo
';

	print $pipeline;
	exit;

}

