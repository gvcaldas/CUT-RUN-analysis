############### Step 6: make bigwig file [deeptools} ##################
# Documentation:
# https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html
#Filter bam files/bamCoverage/normalize by sequencing depth only
# Here I am normalizing with RPGC : useful when visualizing libraries with different sequencing depths. 
# The RPGC normalization is equivalent to getting a 1x depth of coverage.

# What is RPGC?
# RPGC (1x sequencing depth ) : number of reads per bin /(total number of mapped reads * fragment length / effective genome size) 
# Denominator:  total number of mapped reads in the genome/ effective genome size = # reads per every basepair
		#reads per basepair x fragment length = # of reads expected for every chunk of Xbasepairs (fragment)
# RPGC = observed number of reads per bin / Denominator = ratio of actual reads in that bin normalized to expected number of reads per that bin size. 

## this is effective genome size for gh38: 2913022398

# Important Settings:
--ignoreDuplicates (some papers suggest to not do this for C&R or C&T data) : If set, reads that have the same orientation and start position will be considered only once. If reads are paired, the mate’s position also has to coincide to ignore a read. I believe that is OK to work with duplicates, as long as your library is not high on duplicates. For me, >40% duplicates is too much and I would prefer to ignore. But overall, I suggest to get a sense of your data by comparing removing duplicates or not.

--minMappingQuality :A mapping quality of 10 or less indicates that there is at least a 1 in 10 chance that the read truly originated elsewhere.

--centerReads : reads are centered with respect to the fragment length. For paired-end data, the read is centered at the fragment length defined by the two ends of the fragment.

--extendReads : Paired-end: Reads with mates are always extended to match the fragment size defined by the two read mates. Unmated reads, mate reads that map too far apart (>4x fragment length) or even map to different chromosomes are treated like single-end reads. The input of a fragment length value is optional. If no value is specified, it is estimated from the data (mean of the fragment size of all mate reads).

--smoothLength : The smooth length defines a window, larger than the binSize, to average the number of reads. For example, if the –binSize is set to 20 and the –smoothLength is set to 60, then, for each bin, the average of the bin and its left and right neighbors is considered.

## Additional parameters that I tested and seemed to work best: 

--minFragmentLength : The minimum fragment length needed for read/pair inclusion. (Default: 0)

--maxFragmentLength : The maximum fragment length needed for read/pair inclusion. (Default: 0)

## The user provided fragment length is only used as a fallback for singletons or mate reads that map too far apart (with a distance greater than four times the fragment length or are located on different chromosomes).

# Note: I have been processing with BinSize = 20, potentially, I could speed up with binsize = 50 , Plus --smooth length could be 100 or 150. However, Binsize 20 gives better resolution and I don't think smooth of 60 affects has an impact on the data.  

# in summary, these are the settings I think are the best. For data already filtered for discordant pairs.

bamCoverage -b file.st.bam --ignoreDuplicates --minMappingQuality 10 --centerReads --extendReads 200 --minFragmentLength 100 --maxFragmentLength 250 --binSize 20 --smoothLength 60 --normalizeUsing RPGC --effectiveGenomeSize 2913022398 -o file.bw



## example to batch: could be run as a loop.

#!/bin/bash # your header

bamCoverage -b file.bam --ignoreDuplicates --minMappingQuality 10 --centerReads --extendReads 200 --minFragmentLength 100 --maxFragmentLength 250 --binSize 20 --smoothLength 60 --normalizeUsing RPGC --effectiveGenomeSize 2913022398 -o file.bw

### bamCompare [deeptools]
##can be used to generate a bigWig or bedGraph file based on two BAM files that are compared to each other while being simultaneously normalized for sequencing depth.

bamCompare -b1 file.bam -b2 file2.bam --effectiveGenomeSize 2913022398 --ignoreDuplicates --minMappingQuality 10 --centerReads --minFragmentLength 100 --maxFragmentLength 250 --binSize 20 --smoothLength 60 --scaleFactorsMethod SES -o log2ratio.bw -o file.bw

## if want to run the whole pipeline by fastq
# Step 9: make bigwig file
echo "starting bigwig generation"

bamCoverage -b ${SORT_PP_BAM} --ignoreDuplicates --minMappingQuality 10 --centerReads --extendReads 200 --minFragmentLength 100 --maxFragmentLength 250 --binSize 20 --smoothLength 60 --normalizeUsing RPGC --effectiveGenomeSize 2652783500 --verbose -o  ${NAME}.bw
# Arguments
# --ignoreDuplicates	"reads that have the same orientation and start position will be considered only once"
# --minMappingQuality 10	"only reads that have a mapping quality score of at least this are considered"
# --centerReads	"For paired-end data, the read is centered at the fragment length defined by the two ends of the fragment."
# --extendReads 200	This option should be used with Chip-Seq data, where "fragments are known to map contiguously"
# --minFragmentLength 100 
# --maxFragmentLength 250
# --binSize 20	Size of bins used for bigwig/bedgraph output file
# --smoothLength 60	"if the –binSize is set to 20 and the –smoothLength is set to 60, then, for each bin, 
#	the average of the bin and its left and right neighbors is considered. Any value smaller than –binSize 
#	will be ignored and no smoothing will be applied"
# --normalizeUsing RPGC	"RPGC = reads per genomic content (1x normalization)""
# --effectiveGenomeSize 2652783500	effective genome size for the mouse mm10 reference genome from:
#	http://genomewiki.ucsc.edu/index.php/Mm10_Genome_size_statistics
echo "big wig file made"

