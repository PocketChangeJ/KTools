#!/usr/bin/perl

use strict;
use warnings;
use IO::File;
use FindBin;

#################################################################
# kentnf: subroutine						#
#################################################################


=head2
 trimo -- clean reads using trimmomatic
=cut
sub clean_trimo
{

}

=head2
 align -- clean reads through align reads
=cut
sub clean_align
{

}

=head2
 clean -- clean sRNA reads
=cut
sub clean_sRNA
{

}

=head2
 usage: print usage information
=cut
sub usage
{

	exit;
}

=head2
 usage: print pipeline or command 
=cut
sub pipeline
{
	print qq'
>> pipeline for $0
>> 1. clean mRNA
      	1.1 clean adapter, low quality, short reads
	1.2 clean rRNA or other contanmination
	
>> 2. clean sRNA/degradome
	2.1 clean adapter, short reads,
	2.2 clean rRNA or other contanmination

>> 3. clean DNA
	3.1 clean adapter, low quality, short reads
';
	exit;
}
