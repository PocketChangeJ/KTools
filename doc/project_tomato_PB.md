
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

##2 correct PB reads (pbtool)

###2.1 map combined CCS reads to reference 

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

###2.2 filter mapped reads

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

###3 correct the reads using LSC 

####3.1 correct the reads using LSC

combine PAIRED short reads into one file.
	
```
$ cat SR_R1.fq SR_R2.fq > SR_combine.fq
```

set below parameters in LSC config file (run.cfg), leave other parameters as default

>LR_pathfilename = data_path/Heinz_B_SMRT_CCS_draft.fasta
>LR_filetype = fa
>SR_pathfilename = data_path/SR_combine.fq
>SR_filetype = fq
>I_nonredundant = N	#Is this nonredundant SR data set
>SCD = 5		#Short-reads coverage depth (SCD)
>thread1 = 24		#for read alignment
>Nthread2 = 2		#for correction, will be slow if set two high
>temp_foldername = temp
>output_foldername = output_all

Run LSC

```
$ /path/to/runLSC.py run.cfg
```

Correct Test with different Short-reads coverage depth (SCD).
This test only for changed SCD to -1 at mode 2. In the mode 1,
the SCD is still 5

SCD -1 for MG: SeqNum:271826; MaxLen: 43582; MinLen: 154; AvgLen: 3470.17; TotalBase: 943283641
SCD 5 for MG:  SeqNum:271826; MaxLen: 43582; MinLen: 155; AvgLen: 3470.56; TotalBase: 943389079  

>These two results are almost same, so the SCD will not work after mapping of mode 1

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

##4 correct the reads using PBcR

###4.1Install PBcR

Download wgs-assembler8.3rc1 from http://sourceforge.net/projects/wgs-assembler/files/
unzip the file wgs-8.3rc1-Linux_amd64.tar.bz2, then set the bin fdr to $PATH
__/home/kentnf/software/wgs-8.3rc1/Linux-amd64/bin_

###4.2 Prepare input files

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

###4.3 set pacbio.spec file
```
wget http://www.cbcb.umd.edu/software/PBcR/data/selfSampleData/pacbio.spec
```

># limit to 32GB. By default the pipeline will auto-detect memory and try to use maximum. This allow limiting it
>ovlMemory = 96
>ovlStoreMemory = 32000
>merylMemory = 32000
>assemble = 0

###4.4 check the input file size

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

##5 correct using lordec

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

##5 Should I use All Illumina reads to correct PB reads

There are three replicates for sample B sequenced by Illumina. All of them
are strand-specific reads 
- Heinz_B_rep1_paired (B01), 4bp PE reads, 10859261 cleaned read pairs
- HZ-4-B_TCGAAG_paired (B02), 151bp PE reads, 30267228 cleaned read pairs
- Heinz_B_rep2 (B03), 48bp SE reads, 10777050 cleaned reads

###5.1 Using Kmer method

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

###5.2 correlation analysis

Align B01 and B02 to genome using tophat with perfect match, and align B03
allowing 1 edit distance. Then count reads in gene region. 

>Sample  B02       B01     B03
>B02       1      0.98    0.97
>B01    0.98         1    0.99
>B03    0.97      0.99       1

###5.3 junctions identified by tophat

convert junction bed to tab delimit format: bedt_to_juncs < junctions.bed > B01.junc,
then compare then using venn diagram

> B01, 79854 junction; B02, 118858 junction; B03 65923 junction
> They shared 57752 junctions, 
> B02 has 35886 specific junctions, 19803 shared with B01, 5637 shared with B03 
> B01 has 2214 specific junctions, accounting for 1.86% of 118858 in B02 
> B03 has 2449 specific junctions, accounting for 1.89% of 118858 in B02
> B01 and B03 has 285 (0.24%) specific junction compared to B02, 

Use only B02 will lose some result, but it is fine for just lose < 2% of junction
And this lose rate is same kmer analysis

###5.4 lordec method

corrected raw CCS using B02, B02+B01, B02+B01+B03



##6 filter by pbtool

Gene statistics information 

##7 isoform detection using IDP

##8 lncRNA analysis

##8 AS events

