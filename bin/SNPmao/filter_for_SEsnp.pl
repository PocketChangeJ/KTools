#!/usr/bin/perl -w 
use strict; 

-t and !@ARGV and die "perl $0 in_SE_aln.sam\n"; 

###### Make FLAG list 
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


while (<>) {
	if (m/^@/) {
		print; 
	}else{
		chomp; 
		my @ta = split(/\t/, $_); 
		$#ta > 10 or next; 
		$outflag{$ta[1]} and next; 
		print STDOUT "$_\n"; 
		last; 
	}
}
while (<>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	$#ta > 10 or next; 
	$outflag{$ta[1]} and next; 
	print STDOUT "$_\n"; 
}
