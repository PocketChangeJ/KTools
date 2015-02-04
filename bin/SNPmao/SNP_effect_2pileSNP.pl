#!/usr/bin/perl -w 
# Written by: Honghe Sun ( biosunhh@gmail.com ) 
# Method1:20120110 - each SNP is indenpendently from others, so only think about
#		     the effect by indenpendent SNP
#	Category of SNP effect
#	1. gene region : 1.1 up stream 500bp 1.2 down stream 500bp 1.3 gene region; 
#	2. intron region : 2.1 intron/exon border (2bp of exon border) 2.2 intron region;
#	3. coding region : 
#		3.1 start codon -> non-start codon
#		3.2 stop codon -> non-stop codon
#		3.3 Nonsynonymous substitution
#		3.4 non-stop codon -> stop codon
#		3.5 synonymous substitution 
#              	3.6 stop codon -> stop codon (synonymous substitution) 
#              	3.7 start codon -> start codon (synonymous substitution) 
#		3.8 undetermined with N in CDS
# edit20120118:           3.8 undetermined with N in CDS 
# edit20130712:           Edit to fit for the format of pileup2SNP format from linyong. 
# edit20150204:		  Edit to fit other plant species by Yi
# ClV6 genome intron length stat: SUM             MEAN                    MEDIAN  MIN     MAX     Count   NoNull
#                                 38918967        464.338157392383        199     11      9976    83816   83816

use strict; 
use warnings;
use IO::File;
use Bio::SeqIO;
use Getopt::Long;

my $usage = qq'
USAGE: $0 -cds input_CDS.fasta -gtf input.GTF -snp input_SNP_table.txt

* only accept GTF files

';

-t and !@ARGV and die $usage;

my ($cds_file, $gtf_file, $snp_file);
GetOptions (
	"cds=s" => \$cds_file,	# CDS sequence
	"gtf=s" => \$gtf_file,	# GTF file
	"snp=s" => \$snp_file	# SNP file
);

print $usage and exit unless (defined $cds_file && defined $gtf_file && defined $snp_file);

# Initial %codon;
my %codon; 
{
	my @aa = split(//, 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'); 
	my @ss = split(//, '---M---------------M---------------M----------------------------'); 
	my @l1 = split(//, 'TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG'); 
	my @l2 = split(//, 'TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG'); 
	my @l3 = split(//, 'TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG'); 
	for (my $i=0; $i<@aa; $i++) {
		my $setag = '-'; 
		$aa[$i] eq '*' and $setag = '*'; 
		$ss[$i] eq 'M' and $setag = 'M'; 
		my $bbb = join('', $l1[$i], $l2[$i], $l3[$i]); 
		$codon{$bbb} = [$setag, $aa[$i]]; 
	}
}

# load CDS sequence to hash
my %cds_seq;
my $in = Bio::SeqIO->new(-format=>'fasta', -file=>$cds_file);
while(my $inseq = $in->next_seq) {
	$cds_seq{$inseq->id} = $inseq->seq;
}

# set default feature, set it as parameters later
# key: transID
# subkey: chr, value: chr
# subkey: strand, value: +/-
# subkey: CDS/exon, value: region of CDS/exon

my $feature = 'CDS';
my %trans_info = parse_gtf($gtf_file, $feature);

# put trans info to hash
# key : chrID
# value: array; 
# 0 - transcript start, end (cds region)
# 1 - strand
# 2 - array of exon (cds region)
# 3 - transcript ID
# 4 - length of CDS seq
# 5 - sequence of CDS seq

my %anno; 
foreach my $tname (sort keys %trans_info) {
	die "[ERR]no cds seq $tname\n" unless defined $cds_seq{$tname};
	my $cds = $cds_seq{$tname};
	my $ref = $trans_info{$tname}{'chr'};
	my $str = $trans_info{$tname}{'strand'};
	my $c   = $trans_info{$tname}{$feature};
	my @cds_region = split(/\t/, $c);
	@cds_region = sort {$a<=>$b} @cds_region;
	my @se = ($cds_region[0], $cds_region[scalar(@cds_region)-1]);
	
	my @pack_region; # better to pack these region in sub parse_gtf
	if ($str eq '+') {
		for(my $i=0; $i<@cds_region; $i=$i+2) {
			push(@pack_region, [$cds_region[$i], $cds_region[$i+1]]);
		}
	} else {
		@cds_region = sort {$b<=>$a} @cds_region;
		for(my $i=0; $i<@cds_region; $i=$i+2) {
                        push(@pack_region, [$cds_region[$i+1], $cds_region[$i]]);
                }
	}
	push(@{$anno{$ref}}, [[@se], $str, [@pack_region], $tname, length($cds), $cds]);
}

# check anno
=head
foreach my $cid (sort keys %anno) {
	print $cid."\n";
	my $r = $anno{$cid};
	for my $r1 ( @$r ) {
		if ($r1->[1] eq '-') {
			print $r1->[3],"\t",$r1->[1],"\n";

			# print $r1->[2],"\n";
			foreach my $r2 ( @{$r1->[2]} ) {
				print $r2->[0]."\t".$r2->[1]."\n";
			}
			exit;
		}
	}
}
print "exit\n"; exit;
=cut

open(FH, $snp_file) || die $!;

SNP_LINE:
while (<FH>) {
	chomp; s/\s+$//; 
	/^\s*$/ and next;
	/^\s*#/ and next; 
#	/^chr(?:omosome|om)?\t/i and next; 
	my @ta = split(/\t/, $_); 
	my ($chr_id, $chr_pos) = @ta[2,3]; 
	my $refbase = undef(); 
	my $newbase = undef(); 
	my @td; # Storing the effections. 
	if ($ta[0] eq 'type') {
		print STDOUT join("\t", $_, 'Tag')."\n"; 
		next SNP_LINE; 
	} elsif ($ta[0] eq 'M') {
		$ta[1] =~ /^([ATGC])\->([ATGC]);?$/ or die "[Error] Format [$ta[1]] not known.\n$_\n"; 
		$refbase = $1; $newbase = $2; 
		$refbase = uc($refbase); $newbase = uc($newbase); 
		(defined $newbase and $newbase =~ /^[ATGC]$/) or die "Wrong Base:[$newbase]\n$_\n"; 
		my @tc = &parseSNP('', $anno{$chr_id}, $chr_pos, $newbase, $refbase); 
		for my $tcr (@tc) {
			push(@td, join(":", @$tcr)); 
		}
	} elsif ($ta[0] eq 'D') {
		push(@td, "D"); 
	} elsif ($ta[0] eq 'I') {
		push(@td, "I"); 
	} elsif ($ta[0] =~ m/^M([DI])$/) {
		my $td1 = ($1 eq "D") ? 'D' : 'I' ; 
		$ta[1] =~ /^([ATGC])\->([ATGC]);/ or die "[Error] Format of [$ta[1]] not known.\n"; 
		$refbase = $1; $newbase = $2; 
		$refbase = uc($refbase); $newbase = uc($newbase); 
		(defined $newbase and $newbase =~ /^[ATGC]$/) or die "Wrong Base:[$newbase]\n$_\n"; 
		my @tc = &parseSNP('', $anno{$chr_id}, $chr_pos, $newbase, $refbase); 
		for my $tcr (@tc) {
			push(@td, join(":", @$tcr)); 
		}
		$td[-1] .= ";$td1"; 
	} else {
		warn "Unknown type [$ta[0]]\n"; 
		push(@td, $ta[0]); 
	}
	print join("\t", $_, join(",", @td) ) . "\n"; 
}

# my ($cdsseqR, $gffR, $position, $newbase, $refbase) = @_; 
# gffR : [mRNA_S_E, Strand, CDS_SEs, GenID, cds_len, cds_seq]
# back : [effect_type1, effect_type2, ...]; 
# 显然, newbase是'+'链上的; 
sub parseSNP {
	my ($cdsseqR, $gffR, $position, $newbase, $refbase) = @_; 
	my $up_len = 500; 
	$up_len = 2000; # for 10kb upstream of genes; 
	my $down_len = 500; 
	$down_len = 2000; # for 10kb downstream of genes; 
#	defined $_[4] and $_[4]>=0 and $up_len = $_[4]; 
#	defined $_[5] and $_[5]>=0 and $down_len = $_[5]; 
	my $span = 2; 
	my @type; 
	GENE: 
	for my $r1 (@$gffR) {
		if ($r1->[1] eq '+') {
			$position < $r1->[0][0]-$up_len and next GENE; 
			$position > $r1->[0][1]+$down_len and next GENE; 
			if ($position < $r1->[0][0]) {
				push(@type, [$r1->[3], '1.1']); next GENE; 
			}elsif ($position > $r1->[0][1]) {
				push(@type, [$r1->[3], '1.2']); next GENE; 
			}else{
				# 在起止exon限定区间内; 
				my $cds_pos = 0; 
				for my $r2 (@{$r1->[2]}) {
					if ($position > $r2->[1]) {
						# intron2; 
						if ($position <= $r2->[1]+$span) {
							push(@type, [$r1->[3], '2.1']); next GENE; 
						}else{
							# else we should check next exon. 
							$cds_pos += ($r2->[1]-$r2->[0]+1); 
						}
					}elsif ($position < $r2->[0]) {
						# intron1; 
						if ($position >= $r2->[0]-$span) {
							push(@type, [$r1->[3], '2.1']); next GENE; 
						}else{ 
							# 此处需要配合 $position > $r2->[1] 后的 $position <= $r2->[1]+$span 共同使用才合适! 
							push(@type, [$r1->[3], '2.2']); next GENE; 
						}
					}else{
						# position fall in this exon now!
						$cds_pos += ($position-$r2->[0]+1); 
						# 判断cds_pos的类型; 
						my $tbase = $newbase; 
						my $ref_tbase = $refbase; 
						my @ttype = chk_cds_pos($cds_pos, $tbase, \$r1->[5], $r1->[4], $ref_tbase); 
						push(@type, [$r1->[3], @ttype]); next GENE; 
					}
				}# end for 
				scalar(@type) > 0 or die print_error($r1, $position, $newbase, $refbase); 
			}
		}elsif ($r1->[1] eq '-') {
			$position < $r1->[0][0]-$down_len and next GENE; 
			$position > $r1->[0][1]+$up_len and next GENE; 
			if ($position < $r1->[0][0]) {
				push(@type, [$r1->[3], '1.2']); next GENE; 
			}elsif ($position > $r1->[0][1]) {
				push(@type, [$r1->[3], '1.1']); next GENE; 
			}else{
				# 在起止exon限定区间内; 
				my $cds_pos = 0; 
				for my $r3 (@{$r1->[2]}) {
					if ($position < $r3->[0]) {
						# intron2; 
						if ($position >= $r3->[0]-$span) {
							push(@type, [$r1->[3], '2.1']); next GENE;
						}else{
							$cds_pos += ($r3->[1]-$r3->[0]+1); 
						}
					}elsif ($position > $r3->[1]) {
						# intron1; 
						if ($position <= $r3->[1]+$span) {
							push(@type, [$r1->[3], '2.1']); next GENE; 
						}else{
							push(@type, [$r1->[3], '2.2']); next GENE; 
						}
					}else{
						# position fall in this exon now!
						$cds_pos += ($r3->[1]-$position+1); 
						my $tbase = $newbase; 
						$tbase =~ tr/ATGC/TACG/; 
						my $ref_tbase = $refbase; 
						$ref_tbase =~ tr/ATGC/TACG/; 
						my @ttype = chk_cds_pos($cds_pos, $tbase, \$r1->[5], $r1->[4], $ref_tbase);
						push(@type, [$r1->[3], @ttype]); next GENE; 
					}
				}# end for 
				scalar(@type) > 0 or die print_error($r1, $position, $newbase, $refbase); 
			}
		}else{ # Failed to know '+'/'-'
			die "Failed to parse strand for $r1->[1]!\n"; 
		}
	}# end for gffR
	scalar(@type) == 0 and push(@type, ['NULL', '1.3']); 
	return (@type); 
}# end sub parseSNP 

# gffR : [mRNA_S_E, Strand, CDS_SEs, GenID, cds_len, cds_seq]
# back : [effect_type1, effect_type2, ...]; 

# mut应该是与seq同链的碱基类型; 

=head2
 input: cds_pos, tbase, $cds_seq, $cds_len, $ref_tbase
=cut
sub chk_cds_pos
{
	my ($p, $mut, $seqR, $len, $ref) = @_; 
	my ($newb3, $rawb3); 
	if ($p <= 3) {
		my $aa_pos = 1; 
		$rawb3 = $newb3 = substr($$seqR, 0, 3); 
		substr($newb3, $p-1, 1) = $mut; 
		substr($rawb3, $p-1, 1) = $ref; 
		if ($codon{$newb3}[0] eq '*') {
			return ('3.4', "$codon{$rawb3}[1]$aa_pos$codon{$newb3}[1]"); 
		}elsif ($codon{$newb3}[0] eq '-') {
			return ('3.1', "$codon{$rawb3}[1]$aa_pos$codon{$newb3}[1]"); 
		}elsif ($codon{$newb3}[1] ne $codon{$rawb3}[1]) {
			return ('3.7', "$codon{$rawb3}[1]$aa_pos$codon{$newb3}[1]"); 
		}else{
			return ('3.5', "$codon{$rawb3}[1]$aa_pos$codon{$newb3}[1]"); 
		}
	}elsif ($p > $len-3) {
		my $aa_pos = int($len/3); 
		$newb3 = $rawb3 = substr($$seqR, $len-3, 3); 
		$codon{$rawb3}[0] eq '*' or die "Tail BBB is $rawb3, not * but $codon{$rawb3}[1]\n"; 
		substr($newb3, 2-($len-$p), 1) = $mut; 
		substr($rawb3, 2-($len-$p), 1) = $ref; 
		if ($codon{$newb3}[0] ne '*') {
			return ('3.2', "*$aa_pos$codon{$newb3}[1]"); 
		}else{
			return ('3.6'); 
		}
	}else{
		my $aa_pos = int($p/3)+1; 
		my $frame = $p%3; 
		if ($frame == 0) {
			$frame = 3; 
			$aa_pos--; 
		}
		my $ssp = $p-$frame+1; 
		$rawb3 = $newb3 = substr($$seqR, $ssp-1, 3); 
		substr($newb3, $frame-1, 1) = $mut; 
		substr($rawb3, $frame-1, 1) = $ref; 
#defined  or die "$newb3\n$rawb3\n$mut\n$ref\n"; 
		if (!(defined $codon{$newb3}[0])) {
			if ($newb3 =~ /N/i) {
				return ('3.8', "X${aa_pos}X"); 
			}else{
				warn "$p, $mut, $rawb3, $newb3\n"; 
			}
		}
		if ($codon{$newb3}[0] eq '*') {
			return ('3.4', "$codon{$rawb3}[1]$aa_pos$codon{$newb3}[1]"); 
		}elsif ($codon{$newb3}[1] ne $codon{$rawb3}[1]) {
			return ('3.3', "$codon{$rawb3}[1]$aa_pos$codon{$newb3}[1]"); 
		}else{
			return ('3.5', "$codon{$rawb3}[1]$aa_pos$codon{$newb3}[1]"); 
		}
	}
}# end sub chk_cds_pos

=head2
 parse_gtf -- parse gtf file, return gtf information
=cut
sub parse_gtf
{
        my ($input_file, $feature) = @_;
        my %trans_info; # key: tid, chr, exon, gene, strand
        my $fh = IO::File->new($input_file) || die $!;
        while(<$fh>)
        {
                chomp;
                next if $_ =~ m/^#/;
                my @a = split(/\t/, $_);

                if ($a[2] eq $feature)
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
                                die "[ERR]inconsistency chr for $tid\n" if $trans_info{$tid}{'chr'} ne $a[0];
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

                        if ( defined $trans_info{$tid}{$feature}) {
                                $trans_info{$tid}{$feature}.="\t".$a[3]."\t".$a[4];
                        } else {
                                $trans_info{$tid}{$feature} = $a[3]."\t".$a[4];
                        }
                }
        }
        $fh->close;
        return %trans_info;
}

=head2
 print_error -- print line of array for error tracking
=cut
sub print_error
{
	my ($array, $position, $new_base, $ref_base) = @_;

	print "Err here\n";
	print "Start: ",$array->[0][0],"\n",
		"End: ", $array->[0][1],"\n",
		"Strand: ", $array->[1],"\n",
		"TransID: ", $array->[3],"\n";

	foreach my $r (@{$array->[2]}) {
		print $r->[0],"\t",$r->[1],"\n";
	}

	print "$position, $new_base, $ref_base\n";
}
