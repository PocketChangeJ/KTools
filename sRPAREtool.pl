#!/usr/bin/perl

=head
 combine the targetfinder, GSTAr, and fei target result into one for downstream analysis
=cut

use strict;
use warnings;
use IO::File;
use FindBin;
use Getopt::Std;

my $version = 0.1;

my %options;
getopts('a:b:c:d:e:f:g:i:j:k:l:m:n:o:p:q:r:s:t:u:v:w:x:y:z:h', \%options);
unless (defined $options{'t'} ) { usage($version); }

if	($options{'t'} eq 'GSTAr')		{ srna_GSTAr(\%options, \@ARGV); }		# target prediction using GSTAr
elsif	($options{'t'} eq 'targetfinder')	{ srna_targetfinder(\%options, \@ARGV); }	# target prediction using targetfinder
elsif	($options{'t'} eq 'feitarget')		{ srna_feitarget(\%options, \@ARGV); }		# target prediction using feitarget
elsif	($options{'t'} eq 'combine')		{ srna_target_combine(\%options, \@ARGV); } 	# combine the target pred results
elsif	($options{'t'} eq 'convertF2G')		{ srna_convertF2G(@ARGV); }			# convert fei to GSTAr format
elsif   ($options{'t'} eq 'filterGSTAr')	{ srna_filterGSTAr(@ARGV); }			# filter GSTAr result using score
elsif	($options{'t'} eq 'predcmb')		{ srna_pred_cmb(\%options, \@ARGV);}		# combine analysis 
else	{ usage($version); }

#################################################################
# kentnf: subroutine						#
#################################################################

=head2
 srna_filterGSTAr: filter GSTAr using score
=cut
sub srna_filterGSTAr
{
	my ($input, $score ,$output) = @_;
	my $usage = qq'
USAGE $0 -t filterGSTAr input score output
';
	print $usage and exit unless (defined $input && defined $score && defined $output);
	die "[ERR]file exist $input\n" unless -s $input;
	die "[ERR]file exist $input\n" if -s $output;	
	die "[ERR]score $score\n" unless $score > 0;

	my $fho = IO::File->new(">".$output) || die $!;
	my $fh = IO::File->new($input) || die $!;
	while(<$fh>)
	{
		chomp;
		my @a = split(/\t/, $_);
		if ($_ =~ m/^#/ || $_ =~ m/^Query/) {
			print $fho $_."\n";
			next;
		}

		if ($a[8] <= $score) {
			print $fho $_."\n";
		}
	}
	$fh->close;
	$fho->close;
}

=head2
 srna_target_pred : predict of sRNA target using combined method
=cut
sub srna_pred_cmb
{
	my ($option, $files) = @_;

	my $usage = qq'
USAGE: $0 -t predcmb [options] target_reference sRNA_uniq.fa

	-e	dRNA_reads (separate by "," for multiple input)
	-d	degradome density (separate by "," for multiple input)
        -p      number of threads (default 1)
        -c      score cutoff (default 4.5)
	-a 	score for filter GSTAr.txt (default 5.5)

* the target reference should be isoform of mRNA
* provide only one of -e or -d option

';
        print $usage and exit unless @$files == 2;
        my ($reference, $sRNA) = @$files;
        foreach my $f (($reference, $sRNA)) {
                die "[ERR]file not exit, $f\n" unless -s $f;
        }

        my $thread = 1;
        $thread = int($options{'p'}) if (defined $options{'p'} && $options{'p'} > 1);

        my $score = 4.5;
        $score = $options{'c'} if (defined $options{'c'} && $options{'c'} > 0);

	my @dRNA = (); my $mode;
	if (defined $options{'d'}) { 
		@dRNA = split(/,/, $options{'d'}); 
		$mode = "-d";
	} 
	elsif (defined $options{'e'}) {
		@dRNA = split(/,/, $options{'e'});
		$mode = "-e";
	}

	foreach my $d (@dRNA) { die "[ERR]file not exist, $d\n" unless -s $d; }

	# target prediction 
	my ($r1, $r2, $r3) = ($sRNA."_GSTAr.txt", $sRNA."_tf.txt", $sRNA."_fei.txt");
	my $cmd1 = "perl $0 -t GSTAr -p $thread $reference $sRNA";
	my $cmd2 = "perl $0 -t targetfinder -p $thread -c $score $reference $sRNA";
	my $cmd3 = "perl $0 -t feitarget -c $score $reference $sRNA";

	run_cmd($cmd1);
	run_cmd($cmd2);
	run_cmd($cmd3);

	# covert the fei format to GSTAr format
	my ($c1, $c2, $c3) = ($sRNA."_GSTAr.filter.txt", $sRNA."_tf.GSTAr.txt", $sRNA."_fei.GSTAr.txt");
	my $cmd4;
	if (defined $options{'a'}) {}
	else { $cmd4 = "cp $r1 $c1"; }
	my $cmd5 = "perl $0 -t convertF2G $r2 $c2";
	my $cmd6 = "perl $0 -t convertF2G $r3 $c3";

	run_cmd($cmd4);
	run_cmd($cmd5);
	run_cmd($cmd6);

	# combine the target result
	my $cmdA = "perl $0 -t combine $c2 $c3 $c1 > $sRNA\_cmbtarget.txt";
	run_cmd($cmdA);	
	
	# run cleaveland
	my $cl4 = $FindBin::RealBin."/bin/CleaveLand4-4.3/CleaveLand4.pl";
	if (scalar(@dRNA) > 0) {
		foreach my $d (@dRNA) {
			my $cmdB = "perl $cl4 $mode $d -g $sRNA\_cmbtarget.txt -n $reference -o $d\_CL4plot -p 0.05 > $d\_CL4.txt";
			run_cmd($cmdB);
		}
	}
}

=head2
 srna_target_combine: combine sRNA target table for cleaveland4 
=cut
sub srna_target_combine
{
	my ($option, $files) = @_;
	my $usage = qq'
USAGE: $0 -t combine targetfinder feitarget GSTAr > combine_target.txt

* the input file should be in order
';
	print $usage and exit unless @$files == 3;
	foreach my $f (@$files) {
		die "[ERR]file not exist: $f\n" unless $f;
	}

	# key queryID, target ID, start, end, value score;
	# if keys duplicate:
	#   1. using the best score (small is better)
	#   2. using the result of GSTAr.txt 
	my %target;

	# key: queryID, target ID, start, end, score;
	my %target_report;

	my $forder = 0;
	my $fname = '';
	my $head = '';
	foreach my $f (@$files)
	{
		$forder++;
        	$fname = "TargetFinder" if $forder == 1;
		$fname = "Fei" if $forder == 2;
 		$fname = "GSTAr" if $forder == 3;

		die "[ERR]no fname $fname\n" unless $fname;

		my $fh = IO::File->new($f) || die $!;
        	while(<$fh>)
        	{
                	chomp;
	                if ($_ =~ m/^#/ || $_ =~ m/^Query/) {
				$head.= $_."\n" if $fname eq "GSTAr";	
				next;
			}
	                my @a = split(/\t/, $_);
	                my ($qid, $tid, $start, $end, $score) = ($a[0], $a[1], $a[2], $a[3], $a[8]);
	                my $key1 = "$qid\t$tid\t$start\t$end";
	                my $key2 = "$qid\t$tid\t$start\t$end\t$score";

	                my $value2 = "$_";

	                if (defined $target{$key1})
	                {
	                        if ($score < $target{$key1})
	                        {
	                                my $pre_key2 = $key1."\t".$target{$key1};
	                                delete $target_report{$pre_key2};

	                                $target{$key1} = $score;
	                                $target_report{$key2} = $value2;
	                        }
	                        elsif ($score == $target{$key1})
	                        {
	                                $target{$key1} = $score;
	                                $target_report{$key2} = $value2; # using the last GSTAr result
	                        }
	                        else
        	                {
                	                #print "$f\t$key1\t$score\t$target{$key1}\n";
                        	}
                	} else {
                        	$target{$key1} = $score;
	                        $target_report{$key2} = $value2;
	                }
	        }
	        $fh->close;		
	}

	print $head;
	foreach my $key (sort keys %target_report)
	{
        	print $target_report{$key}."\n";
	}

	# sRsRV0007174    Solyc01g006710.2.1      3633    3656    3647    -34.50  -33.80  0.979710144927536       0.5     1-24,3656-3633  NA      ((((((((((((((((((((((((&))))))))))))))))))))))))       ACAUGUCCUAUAUGUCAUUUUUUU&AAAAAAAUGACAUGUAGGACAUGU

}

=head2
 sRNA target analysis software
 TargetFinder and Fei target could be used for plant miR target prediction without dRNA
 the GSTAr only works dRNA analysis
=cut
sub srna_targetfinder
{
	my ($options, $files) = @_;

	my $usage = qq'
USAGE: $0 -t targetfinder [options] target_reference sRNA_uniq.fa 

        -p      number of threads (default 1)
	-c	score cutoff (default 4.5)

* the target reference should be isoform of mRNA

';
        print $usage and exit unless @$files == 2;
        my ($reference, $sRNA) = @$files;
        foreach my $f (($reference, $sRNA)) {
                die "[ERR]file not exit, $f\n" unless -s $f;
        }

	my $thread = 1;
        $thread = int($options{'p'}) if (defined $options{'p'} && $options{'p'} > 1);

	my $score = 4.5;
	$score = $options{'c'} if (defined $options{'c'} && $options{'c'} > 0);

        # split files into sub files according to thread number
        my @sub_file = split_file($sRNA, $thread);

        # run targetfinder bin
        my $targetfinder = $FindBin::RealBin."/bin/TargetFinder_1.6/run_targetfinder.pl";
	my $convert2fei  = $FindBin::RealBin."/bin/TargetFinder_1.6/convertTargetFinder.pl";

        # constract script bash
        my $fh1 = IO::File->new(">tmp1_$sRNA.sh") || die $!;
	my $fh2 = IO::File->new(">tmp2_$sRNA.sh") || die $!;
        foreach my $subf (@sub_file) {
                print $fh1 "perl $targetfinder $subf $reference $score &\n";
		print $fh2 "perl $convert2fei $subf\_TargetAlign $score &\n";
        }
        $fh1->close;
	$fh2->close;

        # perfrom target analysis
        run_cmd("bash tmp1_$sRNA.sh");

        # get status of all the thread 
        # 0 -- finished
        # 1 -- still running 
        while(1) {
                my $s1 = process_ctrl("run_targetfinder.pl");
                last if $s1 == 0;
                sleep(20);
        }

	run_cmd("bash tmp2_$sRNA.sh");

	while(1) {
                my $s2 = process_ctrl("convertTargetFinder.pl");
                last if $s2 == 0;
                sleep(20);
        }

        # combine all results
	my $output_file = $sRNA."_tf.txt";
	my $fho = IO::File->new(">".$output_file) || die $!;

	foreach my $subf (@sub_file)
	{
		my $targetfinder_out = "$subf\_targetfinder.txt";
		die "[ERR] file not exist $targetfinder_out\n" unless -e $targetfinder_out;
	
		my $fhi = IO::File->new($targetfinder_out) || die $!;
		while(<$fhi>) {
			chomp;
			print $fho $_."\n";
		}
		$fhi->close;

		unlink($subf);      # remove sub files
		unlink($targetfinder_out); # remove sub files
	}
	$fho->close;

        unlink("tmp1_$sRNA.sh"); # remove temp command
	unlink("tmp2_$sRNA.sh"); # remove temp command
}

=head1
 srna_GSTAr : prediction of sRNA target using GSTAr
 Guide how to perfrom cleaveland analysis
        1. split sRNA sequence
                split -l input_sRNA_uniq.fa
        2. run GSTAr for each splited sRNA
                GSTAr.pl -q -t xa.fasta ref_transcripts.fasta > xa_GSTAr.txt
                ...
        3. then combine all GSTAr result into one file
                cat *_GSTAr.txt | grep -v "^#" | grep -v "^Query" > sRNA_GSTAr.txt
                * then add a correct GSTAr title for the combined sRNA_GSTAr.txt 
=cut
sub srna_GSTAr
{
	my ($options, $files) = @_;

	my $usage = qq'
USAGE: $0 -t GSTAr [options] target_reference sRNA_uniq.fa 

	-p	number of threads (default 1)

* the target reference should be isoform of mRNA
* the output file is: sRNA_uniq.fa_GSTAr.txt

';
	print $usage and exit unless @$files == 2;
	my ($reference, $sRNA) = @$files;
	die "[ERR] file not exit $reference\n" unless -s $reference;
	die "[ERR] file not exit $sRNA\n" unless -s $sRNA;

	my $thread = 1;
	$thread = int($options{'p'}) if (defined $options{'p'} && $options{'p'} > 1);
	
	# split files into sub files according to thread number
	my @sub_file = split_file($sRNA, $thread);
	
	# GSTAr bin
	my $gstar = $FindBin::RealBin."/bin/CleaveLand4-4.3/GSTAr_v1-0/GSTAr.pl";

	# constract script bash
	my $fht = IO::File->new(">tmp1_$sRNA.sh") || die $!;
	foreach my $subf (@sub_file) {
		print $fht "perl $gstar -q -t $subf $reference > $subf.gstar.txt &\n";
	}
	$fht->close;

	# perfrom target analysis
	run_cmd("bash tmp1_$sRNA.sh");

	# get status of all the thread 
	# 0 -- finished
	# 1 -- still running 
	while(1) {
		my $status = process_ctrl("GSTAr.pl");
		last if $status == 0;
		sleep(20);
	}

	# combine all results
	my $output_file = $sRNA."_GSTAr.txt";
	my $time = `date`; chomp($time);

	my $fho = IO::File->new(">".$output_file) || die $!;
	print $fho qq'# GSTAr version 1.0
# $time
# Queries: $sRNA
# Transcripts: $reference
# Minimum Free Energy Ratio cutoff (option -r): 0.65
# Sorted by: MFEratio
# Output Format: Tabular
Query\tTranscript\tTStart\tTStop\tTSlice\tMFEperfect\tMFEsite\tMFEratio\tAllenScore\tPaired\tUnpaired\tStructure\tSequence
';
	foreach my $subf (@sub_file) {
		my $gstar_out = "$subf.gstar.txt";
		die "[ERR] file not exist $gstar_out\n" unless -s $gstar_out;

		my $fhi = IO::File->new($gstar_out) || die $!;
		while(<$fhi>) {
			chomp;
			next if $_ =~ m/^#/;
			next if $_ =~ m/^Query/;
			print $fho $_."\n" ;
		}
		$fhi->close;
		
		unlink($subf);      # remove sub files	
		unlink($gstar_out); # remove sub files
	}
	$fho->close;

	unlink("tmp1_$sRNA.sh"); # remove temp command
}

=head2 
 process_ctrl: contronl the process 
=cut
sub process_ctrl
{
	my $task_name = shift;
	my $ps = `ps -f`; chomp($ps);
	my @p = split(/\n/, $ps);

	# get status of task 
	#  0 -- the process has been finished
	#  1 -- the process still running
	my $status = 0;
	shift @p;
	foreach my $p (@p) {
		$status = 1 if $p =~ m/\Q$task_name\E/; 
	}

	return $status;
}

=head2 
 split_file: split fasta file into different files
=cut

sub split_file
{
	my ($sRNA, $thread) = @_;

	# get number of total reads (count), and number of sub read (sub count)
	my $count = `wc -l $sRNA`; chomp $count; $count =~ s/ .*//;
	die "[ERR] sRNA seq num $count\n" unless ($count % 2 == 0);
	$count = $count / 2;
	my $sub_count = int($count/$thread);
	$sub_count = 3 if $sub_count <= 3;
	# print "$count\t$sub_count\n"; #exit;

	# put sub files to array
	my @sub_files;

	my $sub_order = 0;
	my $sub_file = "sub$sub_order"."_$sRNA";
	push(@sub_files, $sub_file);

	my $seq_num = 0;
	my $out = IO::File->new(">".$sub_file) || die $!;
	my $fh = IO::File->new($sRNA) || die $!;
	while(<$fh>)
	{
		my $id = $_; chomp($id);
		my $sq = <$fh>; chomp($sq);
		$seq_num++;

		if ($seq_num % $sub_count == 0) {
			print $out "$id\n$sq\n";
			$out->close;	# close prev files
		}
		elsif ($seq_num % $sub_count == 1)
		{
			if ($seq_num > $sub_count)
			{
				$sub_order++;	# start an new file
				$sub_file = "sub$sub_order"."_$sRNA";
				push(@sub_files, $sub_file);
				$out = IO::File->new(">".$sub_file) || die $!; 
			}
			print $out "$id\n$sq\n";
		}
		else
		{
			print $out "$id\n$sq\n";
		}
	}
	$fh->close;
	$out->close;
	return @sub_files;
}

=head1
 srna_feitarget : predict sRNA target using fei's pipeline 
=cut
sub srna_feitarget
{
	my ($options, $files) = @_;	

	my $usage = qq'
USAGE: $0 -t targetfinder [options] target_reference sRNA_uniq.fa

        -p      number of threads (default 1)
        -c      score cutoff (default 4.5)

* the target reference should be isoform of mRNA

';
        print $usage and exit unless @$files == 2;
        my ($reference, $sRNA) = @$files;
        foreach my $f (($reference, $sRNA)) {
                die "[ERR]file not exit, $f\n" unless -s $f;
        }

        my $thread = 1;
        $thread = int($options{'p'}) if (defined $options{'p'} && $options{'p'} > 1);

        my $score = 4.5;
        $score = $options{'c'} if (defined $options{'c'} && $options{'c'} > 0);
	
	my $script = $FindBin::RealBin."/bin/miRNA_target_pred_old.pl";

	my $output_file = $sRNA."_fei.txt";
	run_cmd("perl $script $sRNA $reference $reference $output_file");
}

#################################################################
# kentnf: convert						#
#################################################################

## convert target result from Fei to GSTAr
## Fei's format: 
## sRsRV0236704	Solyc03g116610.2.1	927	946	3	CGAAUUCCACAGUAACGUUU	ACAUAAAGUGUCAUUGCAAA	b|b|||b|||||||||||||	+
## GSTAr's format: 
## sRsRV0236704	Solyc03g116610.2.1	927	946	937	FNA	FNA	FNA	3	FNA	FNA	.(.(((.(((((((((((((&))))))))))))).))).).  CGAAUUCCACAGUAACGUUU&AAACGUUACUGUGAAAUACA
sub srna_convertF2G
{
	my @in = @_;
	my $usage = qq'
USAGE $0 input output

* the input is target pred result with fei format (targetfinder or fei)
* the output file is the format of GSTAr

';
	print $usage and exit unless @in == 2;
	my ($input, $output) = @in;
	print "[ERR]file not exit $input\n" unless -s $input;
	
	my $out = IO::File->new(">".$output) || die $!;
	my $in  = IO::File->new($input) || die $!;

	while(<$in>)
	{
        	chomp;
	        next if $_=~m/^#/;
	        my @a= split(/\t/, $_);
		die "[ERR]col num\n" unless @a == 9;
		my ($query, $target, $start, $stop, $score, $mrna, $srna, $align, $strand) = @a;

		$srna = reverse($srna);
		my $rev_align = reverse($align);
		my $cov_align = convert_stru($mrna, $align);
		my $cov_rev_align = convert_stru($srna, $rev_align);

		$cov_align =~ s/\|/\(/ig;
 		$cov_rev_align =~ s/\|/\)/ig;

		my $seq = $mrna."&".$srna;
		my $stru = $cov_align."&".$cov_rev_align;

        	my $splice = compute_splice_site($seq, $stop);

       		# output, the MEF info, paired, and unpaired were not present in Fei's result
		print $out "$query\t$target\t$start\t$stop\t$splice\tFNA\tFNA\tFNA\t$score\tFNA\tFNA\t$stru\t$seq\n";
	}
	$in->close;
	$out->close;
}

sub compute_splice_site 
{
        my ($seq, $tstop) = @_;

        my $q_seq;
        my $t_seq;

        if($seq =~ /^(\S+)\&(\S+)$/) {
                $t_seq = $1;
                $q_seq = $2;
        } else {
                die "ABORT: In sub-routine compute_slice_site, failed to parse alignment $seq\n";
        }

        my @q_char = split ('', $q_seq);
        my @t_char = split ('', $t_seq);

        my $real_t_pos = $tstop + 1;
        my $real_q_pos = 0;
        my $tch;

        foreach my $qch (@q_char)
        {
                $tch = pop @t_char;
                unless($tch eq "-") {
                        --$real_t_pos;
                }

                unless($qch eq "-") {
                        ++$real_q_pos;
                }

                if($real_q_pos == 10) {
                        return $real_t_pos;
                }
        }

        # Should not get here
       	die "ABORT: Failure in sub-routine compute_slice_site\n";
}

sub convert_stru
{
        my ($seq, $align) = @_;

        die "[ERR]seq not match align: $seq, $align\n" unless ( length($seq) == length($align) );
        my @s = split(//, $seq);
        my @a = split(//, $align);

        for(my $i=0; $i<@s; $i++) {
                $a[$i] = '-' if $s[$i] eq '-';
        }

        $align = join("", @a);
        $align =~ s/b/\./ig;
        $align =~ s/o/\./ig;
        return $align;
}

=head2 
=cut
sub usage
{
	my $usage = qq'
USAGE $0 -t tools [options]

	GSTAr		sRNA target prediction using GSTAr
	targetfinder	sRNA target prediction using targetfinder
	feitarget	sRNA target prediction using fei pipeline
	combine		combine the sRNA target prediction results
	
	convertF2G	convert target format from fei to GSTAr
	filterGSTAr	filter the GSTAr result using score

	predcmb		cleavage site prediction using combined target site

';
	print $usage;
	exit;
}

sub run_cmd
{
	my $cmd = shift;
	print $cmd."\n";
	system($cmd) && die "[ERR]CMD $cmd\n";
}
