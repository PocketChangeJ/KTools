
PacBio RNASeq analyzing pipeline
================================

####1 generate CCS using iso-seq pipeline

#####1.1 Install smrtanalysis

check lib c version, must be greater than 2.5 

	$ whereis libc.so

set environment variables, and create fdr

	$ export SMRT_USER='kentnf'
	$ export SMRT_ROOT='/home/kentnf/software/smrtanalysis'
	$ export SMRT_GROUP='kentnf'
	$ mkdir $SMRT_ROOT
	$ chown -R $SMRT_ROOT:$SMRT_GROUP $SMRT_ROOT

download SMRT analysis from http://www.pacb.com/devnet, and install using

	$ bash smrtanalysis_2.3.0.140936.run -p smrtanalysis-patch_2.3.0.140936.p2.run --rootdir $SMRT_ROOT

>**update SMRT analysis**
> $SMRT_ROOT/admin/bin/smrtportald-initd stop 
> $SMRT_ROOT/admin/bin/smrtupdater smrtanalysis_2.3.0.140936.run
> $SMRT_ROOT/admin/bin/smrtportald-initd start
> $SMRT_ROOT/admin/bin/kodosd start

#####1.2 try to browse localhost:8080, create user 'administrator and set password

#####1.3 move files to /home/kentnf/software/smrtanalysis/userdata/inputs_dropbox

Heinz_B2_cDNA_1_2K_P6_021815_MB_115pM/F01_1.7085b7dd-7927-4a6f-8a76-4beb18c80880.tar.gz
Heinz_B2_cDNA_2_3K_P6_021815_MB_115pM/H01_1.fc8efdaa-8a75-43ca-a11c-980613282574.tar.gz
Heinz_B2_cDNA_3_6K_P6_021815_MB_115pM/F01_1.c5040f3a-1433-4ae0-985b-47a78ec86021.tar.gz
Heinz_MG_cDNA_1_2K_P6_021815_MB_115pM/G01_1.7fab5dfb-9c89-439b-8b3d-22dd4528ec91.tar.gz
Heinz_MG_cDNA_2_3K_P6_021815_MB_115pM/E01_1.bd0066ef-c821-42ea-af91-68dda201f9c1.tar.gz
Heinz_MG_cDNA_3_6K_P6_021815_MB_115pM/G01_1.149e3601-8a15-4b8f-96df-ce5f58cd4195.tar.gz
Heinz_MG_cDNA_5_50K_P6_022015_MB_125pM/C01_1.02ee5772-bd4c-47c4-8f78-809a919e814a.tar.gz

#####1.4 create reference 
	
	$ referenceUploader -c -p /path/to/repository -n GenomeName -f genome.fasta

#####1.5 perform iso-seq pipeline

	Design Job -> create new -> protecols [RS_isoseq.1] 

#####1.6 combine 4 cells

The result locate in smrtanalysis/userdata/jobs/016/016443/data folder.
We use the **isoseq_draft.fasta** for next analysis. The four cells (1-2k,
2-3k, 3-6k, 5-50k) were combined into one files. For easy to understand 
CCS reads ID, we change the format of CCS ID

> Example: B1k638F111F
> - B or M: the stage
> - 1k to 4K: 1k, 1-2k; 2k, 2-3k; 3k, 3-6k; 4k, 5-50k
> - 638: ZMW
> - R or F: CCS direction determined by pb adapter
> - 0 or 1: five seen
> - 0 or 1: polyA seen
> - 0 or 1: three seen
> - F or P: full length or partially

	$ pbtool.pl -t formatCCS -l B1k B2_cDNA_1-2k_isoseq_draft.fasta > B1k.fasta

Then combine 4 cells into one fasta file

>Statistics of PB CCS reads  
> - Heinz_B_SMRT_CCS_draft_newID.fasta -- SeqNum:282655; MaxLen: 46218; MinLen: 300; AvgLen: 3056.23; TotalBase: 863858829
> - Heinz_MG_SMRT_CCS_draft_newID.fasta -- SeqNum:287079; MaxLen: 43579; MinLen: 300; AvgLen: 3512.06; TotalBase: 1008239801

####2 correct PB reads

#####2.1 map combined CCS reads to reference 

	$ mkdir tomato_genome_v2_4
	$ gmap_build -d gmap_db -D tomato_genome_v2_4 /home/database/tomato_genome
	$ gmap -n 0 -D tomato_genome_v2_4 -d gmap_db -f samse -t 31 --sam-use-0M PB_B_map_raw_ID.fa

>Description of gamp parameters
> - -n  Maximum number of paths to show (default 5).  If you want a single alignment plus chimeric alignments, then set this to be 0.
> - -D  genome dir 
> - -d  genome database name
> - -f  output format 
> - -t  CPU
> - --sam-use-0M  Insert 0M in CIGAR between adjacent insertions and deletions Required by Picard, but can cause errors in other tools 

>Statistics of Mapping
> - B  stage: 147685 reads 149999 alignment, 81254 forward(0), 68745 reverse(16) 
> - MG stage: 113788 reads 115726 alignment, 60548 forward(0), 55178 reverse(16)

######2.2 filter mapped reads

> Filter Rules 
> - remove reads with hard clipping
> - keep reads which soft clipping length < 100bp and 10%
> - keep reads which edit distance < 20%  and  mismatch < 10%
> - keep reads all junction sites support by Illumina reads (min depth = 1)
> - keep reads which continuous edit distance should be less than 5 across junction site

> Problem
> The high quality reads without illumina support will be removed

	$ pbtool.pl -t correctF B_mapped.sam  B_SR.sam
	$ pbtool.pl -t correctF MG_mapped.sam MG_SR.sam

After filter, 80450 reads from B, 53108 reads from MG were filtered out because their mapping quality is good and do not need to correct.
So, the left reads need to correct is 67235, and 60680 for B and MG stages. This filter step will same more time on next correction step.
BTW, the correction will not always correct the reads 100%.

#####2.3 correct the reads

#####2.3.1 correct the reads using LSC

combine PAIRED short reads into one file.
	
	$ cat SR_R1.fq SR_R2.fq > SR_combine.fq

set below parameters in LSC config file (run.cfg), leave other parameters as default

	LR_pathfilename = data_path/Heinz_B_SMRT_CCS_draft.fasta
	LR_filetype = fa
	SR_pathfilename = data_path/SR_combine.fq
	SR_filetype = fq
	I_nonredundant = N	#Is this nonredundant SR data set
	SCD = 5			#Short-reads coverage depth (SCD)
	Nthread1 = 24		#for read alignment
	Nthread2 = 24		#for correction
	temp_foldername = temp
	output_foldername = output_all

Run LSC

	$ /path/to/runLSC.py run.cfg

>I use LSC to correct all the reads, then extract what I need by ID.
>For 67235 reads in stage B, 66849 of them were corrected by LSC. 386 reads were not in final correct result
>Manually check several reads in 386, mapping is not good for me. But I still add it to our file for next.
>The final corrected reads named as B_mapped_LSC.fasta

Map LSC corrected reads to reference

	$ gmap -n 0 -D tomato_genome_v2_4 -d gmap_db -f samse -t 31 --sam-use-0M B_mapped_LSC.fasta > B_mapped_LSC.sam
 
check the mapping statistics 

	$ pbtools.pl -t correctF -d 0 -k 1 B_mapped_LSC.sam B_SR.sam > B_mapped_LSC.report.txt

#####2.3.2 correct the reads using PBcR

We use all draft CCS as input to perform correction. From the report there are 146169 CCS reads were correct, and split into 204,539 subreads
The we check the correct reads overlapped with 67235 mapped but uncorrected. 60415 of them were corrected by PBcR, and split into 96266 subs.


#####2.3.3 correct the LSC reads using PBcR

	$ java -jar convertFastaAndQualToFastq.jar PB_B_map_raw_ID.fa > PB_B_map_raw_ID.fq 	# covnert fasta to fataq
	$ fastqToCA -insertsize 200.100 -libraryname HZ-4-B_paird -mates HZ-4-B_TCGAAG_paired_R1.fq,HZ-4-B_TCGAAG_paired_R2.fq	# fastq to CA
	$ /path/to/PBcR -length 50 -partitions 200 -l B_map_LSC -threads 40 -maxGap 10000 -s pacbio.spec -fastq PB_B_map_raw_ID.fq HZ-4-B_paird.frg > B_map_LSC_PBcR.log &

Input 66849 reads corrected by LSC, 65015 of the were further corrected by PBcR. But the corrected reads were split into 166,869 sub reads.
26361 reads were not split by PBcR, others were split from 2 to 44 subreads 

  Read  No. of sub
  26361 1
  14794 2
   9372 3
   5725 4
   3361 5
   1920 6
   1163 7
    726 8
    476 9
    344 10
    255 11
    183 12
    116 13
     98 14
     48 15
     27 16
     21 17
      9 18
      7 19
      4 20
      1 21
      1 22
      1 23
      1 30
      1 44

Next, map these splited reads to genome using gmap, check mapping statistics

	$ gmap -n 0 -D tomato_genome_v2_4 -d gmap_db -f samse -t 20 --sam-use-0M B_map_LSC_PBcR_corr.fasta
	$ pbtools.pl -t correctF -d 0 -k 1 B_map_LSC_PBcR_corr_gmap_mapped.sam B_SR.sam > B_mapped_LSC_PBcR.report.txt




