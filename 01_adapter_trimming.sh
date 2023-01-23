########## Step 1: Adapter trimming (using bbduk) ##########
# Documentation: https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/
# adapters.fa can be found at https://github.com/BioInfoTools/BBMap/blob/master/resources/adapters.fa

# if you want to batch this job as a loop (for all fastq files)
# create .sh file with any text editor in command line (or upload your file)

vim trimming.sh

##!/bin/bash - include here your bash header (cluster-specific)
# load environment or modules
for i in /data/your_path/FASTQ/*_1.fastq.gz ; do
name=$(basename ${i} _1.fastq.gz) ;
bbduk.sh -Xmx1g in1=${name}_1.fastq.gz in2=${name}_2.fastq.gz out1=/data/your_path/FASTQ/Trimmed/${name}_trim_1.fastq.gz out2=/data/your_path/FASTQ/Trimmed/${name}_trim_2.fastq.gz ref=/data/path_to_adapters/BBmap/bbmap/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo; 
done

## Alternative: if you want to run all pipeline per sample, this would be the part for trimming

echo "starting adapter trimming with bbduk"
bbduk.sh -Xmx1g in1=${R1} in2=${R2} out1=${TRIM_R1} out2=${TRIM_R2} ref=adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo
# Arguments:
# -Xmx1g forces bbduk to use 1 Gb of memory
# ktrim=r "once a reference kmer is matched in a read, that kmer and all the bases to the right will be trimmed, 
#	leaving only the bases to the left; this is the normal mode for adapter trimming"
# k=23 "the kmer size to use (must be at most the length of the adapters)"
# mink=11 "allows it to use shorter kmers at the ends of the read (for example, k=11 for the last 11 bases)"
# hdist = 1 is the number of allowed mismatches
echo "finished adapter trimming"


