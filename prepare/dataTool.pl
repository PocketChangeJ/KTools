#!usr/bin/perl

=head
 save RNA-seq or DNA-seq raw dataset using cram format file and corresponding removed and trimmed read in clean process
 then extract the raw dataset from these file
=cut

use strict;
use warnings;
use IO::File;


=head2

=cut

=head2
 save data - will generated trimmed base, removed seq, aligned_cram, and unaligned_cream four files
	     and another reference info file will include reference file path, md5, and seq length
=cut
sub save_data
{
	
	my ($mapped_bam, $unmapped_bam, $removed_read, $trimed_base) = @_;
	
	# checking file
	print "[ERR]mapped file is not bam: $mapped_bam\n" and exit unless $mapped_bam =~ m/\.bam$/;
	



	# step1. file convert
	



	# step2. md5 check and write info file

# md5sum yeast.fasta	
# samtools faidx yeast.fasta
# samtools view -T yeast.fasta -C -o yeast.cram yeast.bam

}

=head2
 extract_data - will convert the four files to raw reads
                only the cram file match the reference info, data could be extracted
=cut
sub extract_data
{
	my ($mapped_cram, $unmapped_cram, $removed_read, $trimed_base) = @_;

	# step1. compare md5 between input file and compress info file
	
	# step2. file convert
	
	# step3. genrated the raw dataset
}

=head2
 run_cmd: run command
=cut
sub run_cmd
{

}

=head2
 pack_trimmomatic - pack the trimmomatic trim log and fastq file to tab format
 input file is raw fastq file and trim log info, the trimmed base were extracted 
 and combined with trim log file
=cut
sub pack_trimmomatic
{
	my ($fq_file, $trim_log) = shift;


	# load trim log file
	

}



