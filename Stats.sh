########## Step 8: estimate the reads that will be filtered under the paramaters chosen in downstream analysis ##########

# Documentation:
# https://deeptools.readthedocs.io/en/develop/content/tools/estimateReadFiltering.html

##estimateReadFiltering: to estimate the reads that will be filtered under the paramaters chosen in downstream analysis.

## Usage example:
estimateReadFiltering -b paired_chr2L.bam \
--minMappingQuality 5 --samFlagInclude 16 \
--samFlagExclude 256 --ignoreDuplicates

## BamPEFragmentsize: This tool calculates the fragment sizes for read pairs given a BAM file from paired-end sequencing. Several regions are sampled depending on the size of the genome and number of processors to estimate the summary statistics on the fragment lengths. Properly paired reads are preferred for computation, i.e., it will only use discordant pairs if no concordant alignments overlap with a given region. The default setting simply prints the summary statistics to the screen.

## Usage example:
deepTools2.0/bin/bamPEFragmentSize \
-hist fragmentSize.png \
-T "Fragment size of PE RNA-seq data" \
--maxFragmentLength 1000 \
-b testFiles/RNAseq_sample1.bam testFiles/RNAseq_sample2.bam \
testFiles/RNAseq_sample3.bam testFiles/RNAseq_sample4.bam \
-samplesLabel sample1 sample2 sample3 sample4


#!/bin/bash # your header

module load samtools 

estimateReadFiltering -b file.bam file2.bam file3.bam --ignoreDuplicates --minMappingQuality 10 --binSize 10000 -o EstimateReadFilter_mapQ10.gh38.CR.txt 

bamPEFragmentSize -b file.bam file2.bam file3.bam --samplesLabel label1 label2 label3 --maxFragmentLength 1000 --histogram bamPEFragSize.CR.set1.gh38.png --table bamPEFragSize.CR.set1.gh38.txt


echo "start estimation of read filtering"
estimateReadFiltering -b ${SORT_PP_BAM} --ignoreDuplicates --minMappingQuality 10 --binSize 10000 -o EstimateReadFilter_${NAME}.txt --verbose -p 8
# Arguments:
#  --ignoreDuplicates	"reads that have the same orientation and start position will be considered only once"
# --minMappingQuality 10	mapping quality must be at least 10 to be included
# --binSize 10000	the bin sized used to estimate metrics
# --verbose	prints processing messages
# -p 8	number of processors to use to speed up the run
echo "finished estimate read filtering step"

date

echo "starting bamPEFragmentSize"
# Documentation:
# https://deeptools.readthedocs.io/en/develop/content/tools/bamPEFragmentSize.html
bamPEFragmentSize -b ${SORT_PP_BAM} --samplesLabel ${NAME} --maxFragmentLength 1000 --histogram bamPEFragSize.${NAME}.png --table bamPEFragSize.${NAME}.txt --verbose
echo "bamPEFragmentSize finished"

date

