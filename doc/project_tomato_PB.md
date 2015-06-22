
PacBio RNASeq analyzing pipeline
================================

##1 generate CCS using iso-seq pipeline

###1.1 Install smrtanalysis

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

###1.2 try to browse localhost:8080, create user 'administrator' and set password

###1.3 move files to /home/kentnf/software/smrtanalysis/userdata/inputs_dropbox

Heinz_B2_cDNA_1_2K_P6_021815_MB_115pM/F01_1.7085b7dd-7927-4a6f-8a76-4beb18c80880.tar.gz
Heinz_B2_cDNA_2_3K_P6_021815_MB_115pM/H01_1.fc8efdaa-8a75-43ca-a11c-980613282574.tar.gz
Heinz_B2_cDNA_3_6K_P6_021815_MB_115pM/F01_1.c5040f3a-1433-4ae0-985b-47a78ec86021.tar.gz
??? file name of B 5-50k
Heinz_MG_cDNA_1_2K_P6_021815_MB_115pM/G01_1.7fab5dfb-9c89-439b-8b3d-22dd4528ec91.tar.gz
Heinz_MG_cDNA_2_3K_P6_021815_MB_115pM/E01_1.bd0066ef-c821-42ea-af91-68dda201f9c1.tar.gz
Heinz_MG_cDNA_3_6K_P6_021815_MB_115pM/G01_1.149e3601-8a15-4b8f-96df-ce5f58cd4195.tar.gz
Heinz_MG_cDNA_5_50K_P6_022015_MB_125pM/C01_1.02ee5772-bd4c-47c4-8f78-809a919e814a.tar.gz

###1.4 create reference 

```	
$ referenceUploader -c -p /path/to/repository -n GenomeName -f genome.fasta
```

###1.5 perform iso-seq pipeline

>Design Job -> create new -> protecols [RS_isoseq.1] 
>set minimum full passes to 0, and minimum predicted accuracy to 75
>set minimum sequence length to 300
>do not perfrom cluster or polish using ICE and Quiver

###1.6 combine 4 cells

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

```
$ pbtool.pl -t formatCCS -l B1k B2_cDNA_1-2k_isoseq_draft.fasta > B1k.fasta
```

Then combine 4 cells into one fasta file

>Statistics of PB CCS reads  
> - Heinz_B_SMRT_CCS_draft_newID.fasta -- SeqNum:282655; MaxLen: 46218; MinLen: 300; AvgLen: 3056.23; TotalBase: 863858829
> - Heinz_MG_SMRT_CCS_draft_newID.fasta -- SeqNum:287079; MaxLen: 43579; MinLen: 300; AvgLen: 3512.06; TotalBase: 1008239801

###1.7 map combined CCS reads to reference

```
$ mkdir tomato_genome_v2_4
$ gmap_build -d gmap_db -D tomato_genome_v2_4 /home/database/tomato_genome
$ gmap -n 0 -D tomato_genome_v2_4 -d gmap_db -f samse -t 31 --sam-use-0M PB_B_map_raw_ID.fa
```
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

###1.8 blast CCS reads to virus
We find that most of CCS are from virus using blast method. And they are mainly from 
Tomato mottle mosaic virus. So we use KF477193 (Tomato mottle mosaic virus isolate 
MX5, complete genome) to perform test.

```
$blastall -p blastn -i B_CCS_draft.fasta -d KF477193.fasta -m 8 -a 24 -F F -e 1e-5 -o B_KF477193.blast.m8.txt
```
In 282655 CCS from B, 116925 (41.36%) of them hit to KF477193. 

##2 correct PB reads
below include different correction test report for determine which 
method is good for PB ISO-seq correction

###2.1 correct the reads using LSC 

####correct the reads using LSC

combine 150 PAIRED short reads into one file.
```
$ cat SR_R1.fq SR_R2.fq > SR_combine.fq
```

set below parameters in LSC config file (run.cfg), leave other parameters as default

>LR_pathfilename = data_path/Heinz_B_SMRT_CCS_draft.fasta
>LR_filetype = fa
>SR_pathfilename = data_path/SR_combine.fq
>SR_filetype = fq
>I_nonredundant = N	# Is this nonredundant SR data set
>SCD = 5		# Short-reads coverage depth (SCD)
>thread1 = 24		# threads for read alignment
>Nthread2 = 10		# threads for correction
>temp_foldername = temp
>output_foldername = output_all

Run LSC

```
$ /path/to/runLSC.py run.cfg
```

Correct Test with different Short-reads coverage depth (SCD).
This test only for changed SCD to -1 at mode 2. In the mode 1,
the SCD is still 5

>SCD -1 for MG: SeqNum:271826; MaxLen: 43582; MinLen: 154; AvgLen: 3470.17; TotalBase: 943283641
>SCD 5 for MG:  SeqNum:271826; MaxLen: 43582; MinLen: 155; AvgLen: 3470.56; TotalBase: 943389079  

>In simple test, these two results are almost same, so the SCD will 
not work after mapping of mode 1

Corrected result
- LSC correct 273343 reads using short-reads
- 9312 reads were not corrected by LSC
- 428 of the uncorrected reads could be mapped to tomato genome

> need to discarded
>I use LSC to correct all the reads, then extract what I need by ID.
>For 67235 reads in stage B, 66849 of them were corrected by LSC. 386 reads were not in final correct result
>Manually check several reads in 386, mapping is not good for me. But I still add it to our file for next.
>The final corrected reads named as B_mapped_LSC.fasta

Map LSC corrected (include 9312 uncorrected) reads to reference

```
$ gmap -n 0 -D tomato_genome_v2_4 -d gmap_db -f samse -t 31 --sam-use-0M B_mapped_LSC.fasta > B_mapped_LSC.sam
```

Check the mapping statistics 

- There are 163238 alignment in B_mapped.sam
  - 157374 reads aligned one position of reference
  - 2932 reads aligned two position of reference, all of them have hard clip
compared 149999 alignment before correction, another 8.8% reads mapped after correction by LSC


```	
$ pbtools.pl -t correctF -d 0 -k 1 B_mapped_LSC.sam B_SR.sam > B_mapped_LSC.report.txt
```

###2.2 correct the reads using PBcR

####Install PBcR

Download wgs-assembler8.3rc1 from http://sourceforge.net/projects/wgs-assembler/files/
unzip the file wgs-8.3rc1-Linux_amd64.tar.bz2, then set the bin fdr to $PATH
__/home/kentnf/software/wgs-8.3rc1/Linux-amd64/bin_

####Prepare input files

convert PB fasta to fastq
```
$ java -jar ../convertFastaAndQualToFastq.jar B_LSC.fasta > B_LSC.fastq
```

convert illumina fastq to CA (must get insertsize)
```
$ bowtie -v 0 -k 1 --best -X 5000 -p 24 -S /home/database/01transcriptome/tomato_cdna \
  -1 HZ-3-MG_TCCCGA_paired_R1.fq -2 HZ-3-MG_TCCCGA_paired_R2.fq MG.sam
$ filter_SAM.pl -n -i MG.sam
$ java -Xmx20g -jar /home/kentnf/software/picard-tools-1.128/picard.jar \
   CollectInsertSizeMetrics I=B_out.sam O=B_out.insertSize H=B_out.insertSizeHistro.pdf

# the insert size for B : 218.725756    59.661029
# the insert sieze for MG : 226.77027   60.978336

$ fastqToCA -insertsize 218 60 -libraryname HeinzB -technology illumina -type sanger \
  -innie -mates HZ-4-B_TCGAAG_paired_R1.fq,HZ-4-B_TCGAAG_paired_R2.fq > HeinzB_illumina_PE.frg
```

####set pacbio.spec file
```
wget http://www.cbcb.umd.edu/software/PBcR/data/selfSampleData/pacbio.spec
```

># limit to 32GB. By default the pipeline will auto-detect memory and try to use maximum. This allow limiting it
>ovlMemory = 96
>ovlStoreMemory = 32000
>merylMemory = 32000
>assemble = 0

####check the input file size

> B_LSC.fasta SeqNum:282655; MaxLen: 46218; MinLen: 119; AvgLen: 3011.15; TotalBase: 851117306
> althought the minLen is 119, I would like to set the length to 300 for correction
> total base is 851,117,306, so If the transcript size is 40MB, the max corrected base is 22x

###4.5 Run PBcR
```
PBcR -length 300 -partitions 200 -l HB-cLP -threads 60 -s pacbio.spec -fastq B_LSC.fastq genomeSize=40000000 HeinzB_illumina_PE.frg > run.out 2>&1
```

###4.6 Connect the PBcR result (add uncorrected reads to final file)

####split the LOG(HB-cLP.log) file into 30
	$ perl PBcRconnect.pl HB-cLP.log
####run the PBcRconnect.pl for LOG simultaneously
	$ perl PBcRconnect.pl B_LSC.fasta HB-cLP.fasta LOG00 > LOG00.txt &

>After connect the PBcR result, add the uncorrect reads final reads for mapping

###4.7 align connected corrected reads to tomato genome
```
gmap -n 0 -D /home/kentnf/yb/software/tomato_genome_v2_4 -d gmap_db -f samse -t 40 \
   --sam-use-0M HB-cLP-conn-cmb.fasta > HB-cLP-conn-cmb.sam
```

###4.8 result description

We pefrom two sets of analysis using PBcR. One is correction of raw CCS reads (B_CCS_draft.fasta); 
and another is correction of LSC corrected reads. 

For raw CSS correction
- There are 158581 alignments
  - 152273 reads aligned one position of reference
  - 3154 reads aligned two position of reference, all of them have hard clip
compared 149999 alignment before correction, another 5.7% reads mapped after correction by PBcR, the 
correction is not better than LSC.

For LSC correction
- There are 167612 alignments
  - 156970 reads aligned one position of reference
  - 5321 reads aligned two position of reference, all of them have hard clip
compared 149999 alignment before correction, another 11.7% reads mapped after correction by LSC & PBcR

###2.3 correct using lordec

>__notice__ : I only use R1 for correction by input HZ-4-B_TCGAAG_paired_R1.fq,HZ-4-B_TCGAAG_paired_R2.fq
But it is ok for testing different parameters.

```
lordec-correct -2 HZ-4-B_TCGAAG_paired_R1.fq,HZ-4-B_TCGAAG_paired_R2.fq -k 19 -s 3 -i Heinz_B_SMRT_CCS_draft.fasta -o B_CCS_Lord_K19.fasta
```

Try to use different kmer for correction
Kmer	15	17	19	21	23
Single	145600	156208	159623	158092	156091
Double	2510	3445	4242	4377	4310

> the K=19 will generate best result

Try to use different solid threshold for correction
solid_threshold		2	3	4	5
Single			159221	159623	159736	159653
Double			4261	4242	4211	4261

> it seems the s=4 will genearate better result

Try to use different order for correction (Correct LSC and PBcR result)
(K=19, s=3)
method	Raw	LSC	PBcR	LSC+PBcR
Single	159623	161634	159419	161007
Double	4242	5389	4835	6179

> it seems the LSC+LoRDEC will get better result

__Note__
Please check the xls files for detail comparison of different strategy.

###2.4 Should I use All Illumina reads to correct PB reads

There are three replicates for sample B sequenced by Illumina. All of them
are strand-specific reads 
- Heinz_B_rep1_paired (B01), 4bp PE reads, 10859261 cleaned read pairs
- HZ-4-B_TCGAAG_paired (B02), 151bp PE reads, 30267228 cleaned read pairs
- Heinz_B_rep2 (B03), 48bp SE reads, 10777050 cleaned reads

####Using Kmer method

This method is quite simple. Generate 17mer for reads of each sample, then
compare the 17mer to figure the coverage among these samples.

- B01 has  66064418 17mer 
- B02 has 186305234 17mer 
- B03 has  32533278 17mer
- B02 could cover 81.12% (53589601) of B01 with 1 depth
- B02 could cover 82.45% (26826589) of B03 with 1 depth
- B01 could cover 66.45% (21621260) of B03 with 1 depth
- B02 could cover 98.36% (37027284/37643718) of B01 with 2 depth
- B02 could cover 98.19% (18100749/18434211) of B03 with 2 depth
- B01 could cover 90.90% (16756720/18434211) of B03 with 2 depth 
With solid = 2, B02 could cover B01 and B02 well. 

####Correlation analysis

Align B01 and B02 to genome using tophat with perfect match, and align B03
allowing 1 edit distance. Then count reads in gene region. 

>Sample  B02       B01     B03
>B02       1      0.98    0.97
>B01    0.98         1    0.99
>B03    0.97      0.99       1

####Junctions identified by tophat

convert junction bed to tab delimit format: bedt_to_juncs < junctions.bed > B01.junc,
then compare then using venn diagram

> B01, 79854 junction; B02, 118858 junction; B03 65923 junction
> They shared 57752 junctions, 
> B02 has 35886 specific junctions, 19803 shared with B01, 5637 shared with B03 
> B01 has 2214 specific junctions, accounting for 1.86% of 118858 in B02 
> B03 has 2449 specific junctions, accounting for 1.89% of 118858 in B02
> B01 and B03 has 285 (0.24%) specific junction compared to B02, 

Use only B02 will lose some result, but it is fine for just lose < 2% of junction
And this lose rate is same kmer analysis. 
**But I use all the SR reads for pbtool filter.

####SR data usage for correction

__corrected raw CCS using B02, B02+B01, B02+B01+B03__
As suggest by LoRDEC author, the min length set to 75 will get better result.
So I only use B02 dataset for correction.

####Correct reads using longer high Quality SR reads

As suggested by LoRDEC author, the quality of SR reads will affect the correction.
So I perform another cleaning for SR reads using Q30, and set min length to 75.
The Q30 were tested and compared with Q20 with better result

```
java -jar /data/home/kentnf/software/KTools/bin/trimmomatic0.32.jar PE -threads 24 -phred33 HZ-4-B_TCGAAG_paired_R1.fq HZ-4-B_TCGAAG_paired_R2.fq HZ-4-B_TCGAAG_paired_Q30_R1.fastq HZ-4-B_TCGAAG_single_Q30_R1.fastq HZ-4-B_TCGAAG_paired_Q30_R2.fastq HZ-4-B_TCGAAG_single_Q30_R2.fastq LEADING:3 TRAILING:3 SLIDINGWINDOW:4:30 MINLEN:75

java -jar /data/home/kentnf/software/KTools/bin/trimmomatic0.32.jar PE -threads 24 -phred33 HZ-3-MG_TCCCGA_paired_R1.fq  HZ-3-MG_TCCCGA_paired_R2.fq HZ-3-MG_TCCCGA_paired_Q30_R1.fastq HZ-3-MG_TCCCGA_single_Q30_R1.fastq HZ-3-MG_TCCCGA_paired_Q30_R2.fastq HZ-3-MG_TCCCGA_single_Q30_R2.fastq LEADING:3 TRAILING:3 SLIDINGWINDOW:4:30 MINLEN:75
```
>__clean report__
>TrimmomaticPE: Started with arguments: -threads 24 -phred33 HZ-3-MG_TCCCGA_paired_R1.fq HZ-3-MG_TCCCGA_paired_R2.fq HZ-3-MG_TCCCGA_paired_Q30_R1.fastq HZ-3-MG_TCCCGA_single_Q30_R1.fastq HZ-3-MG_TCCCGA_paired_Q30_R2.fastq HZ-3-MG_TCCCGA_single_Q30_R2.fastq LEADING:3 TRAILING:3 SLIDINGWINDOW:4:30 MINLEN:75
>Input Read Pairs: 21354750 Both Surviving: 11715187 (54.86%) Forward Only Surviving: 1936181 (9.07%) Reverse Only Surviving: 4365711 (20.44%) Dropped: 3337671 (15.63%)
TrimmomaticPE: Completed successfully

With Q30, more than half of paired, and ~30% of single high quality reads remained.
Here is the list of high quality reads used for correction

>HZ-4-B_TCGAAG_paired_Q30_R1.fastq
>HZ-4-B_TCGAAG_paired_Q30_R2.fastq
>HZ-4-B_TCGAAG_single_Q30_R1.fastq
>HZ-4-B_TCGAAG_single_Q30_R2.fastq

```
lordec-correct -2 B_list_Q30 -k 19 -s 3 -i B_LSC.fasta -o ZBQ30B.fasta
gmap -n 0 -D /home/kentnf/yb/tomato_genome_v2_4 -d gmap_db -f samse -t 64 --sam-use-0M ZBQ30B.fasta > ZBQ30B.sam
lordec-correct -2 MG_list_Q30 -k 19 -s 3 -i MG_LSC.fasta -o ZMGQ30B.fasta
gmap -n 0 -D /home/kentnf/yb/tomato_genome_v2_4 -d gmap_db -f samse -t 64 --sam-use-0M ZMGQ30B.fasta > ZMGQ30B.sam
```
method  MG	B(Q30)	B(Q20)	B(Q15) -- The Q15 only use Reverse Strand of PE reads
Single	126639	161815	161572	161634
Double	4810	5384	5484	5389

BQ20 count
ZBQ20B.fasta -- TotalRead:282655, Corrected:281925      TotalBase:838119417, correct: 640514046, un-correct: 197605371

__NOTE__
The Q30 generate a little better result than Q20. So we use Q30 for next analysis. 

###2.5 Filter the result using pbtool.

Before filter, please combine the PB reads (stageB and MG), SR junctions(stageB and MG, all SR samples)		
- A. remove all reads with hard clip (will perform fusion gene analysis in the further)		
- B. remove all soft clip part, rename the reads ID for soft clip		
 - example: B1k43249R101P-L200R300	
 - L: soft clip of 200bp from left	
 - R: soft clip of 300bp from right 	
 - output the clipped sequence for virus-mRNA analysis in the further	
- C. remove reads less than 300bp after soft clipping 		
- D. For junctions		
 - 1. keep the junctions supported by GTF, illumina short Reads	
 - If the intron length is from 40 to 10000, using below filters	
 - 2. keep the junctions supported by more than 2 PB reads supported	
 - 3. for single PB read supported junctions, keep the junction with conserved splicing site (GT-AG, GC-AG)	
 - 4. for single PB reads supported unconsrved junctions, keep the juncsion without Errors in 4 base window	
 - from step2, remove juctions with longer or shorter intron length (junction not conserved with GTF)	
 - 5. keep the longer and short junction 	
    - 5.1 supported by 1 pb reads, 5 SR reads at least
    - 5.2 with conserved splicing site
    - 5.3 without ERR
- E. filter reads with more Edit distance (>30% of aligned region)		

###2.6 summary

Why we does not use ICE and Quiver to generate polished cluster for CCS
- the numbler of output read is lower, may missing some AS events
  http://blog.pacificbiosciences.com/2014/10/data-release-whole-human-transcriptome.html

- there are still ERROR compared to genome 

The reads we generated

##3 build GTF

















