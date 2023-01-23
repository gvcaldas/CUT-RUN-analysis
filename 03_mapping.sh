
########## Step 3: Align to reference genome ##########
# Note: check reference path variable at beginning of script!
# bowtie2 documentation: https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml

## Build the bowtie2 reference genome index if needed: bowtie2-build path/to/genome/genome.fa /path/to/bowtie2Index/index_name
# define Path to genome hg38 index and GTF.

# if you want to batch this job as a loop (for all fastq files)
# create .sh file with any text editor in command line (or upload your file)

vim 03_mapping.sh

##!/bin/bash - include here your bash header (cluster-specific)
# load environment or modules

module load bowtie2

for i in ./FASTQ/Trimmed/*1.fastq.gz; do 
name=$(basename ${i} _trim_1.fastq.gz);
bowtie2 --very-sensitive --no-mixed --no-discordant -I 10 -X 700 -x /path_to/Genome/Gh38 -1 ./FASTQ/Trimmed/${name}_trim_1.fastq.gz -2 ./FASTQ/Trimmed/${name}_trim_2.fastq.gz -S ./SAM/${name}.sam &> ./SAM/Summary/${name}_bt2.txt
done

# Arguments:
# -p 8	bowtie2 will use multiple processors/cores to speed up alignment step. 
# --very-sensitive	alignment preset option optimized for accuracy and sensitivity
# --no-mixed	bowtie2 won't try to find matches for individual mates if the pair can't be aligned
# --no-discordant	bowtie2 won't try to find discordant alignments if concordant alignments aren't found
# -I 10	"minimum fragment length for valid paired-end alignments"
# -X 700	"maximum fragment length for valid paired-end alignments"
# -x ${REF}	reference genome in bowtie2 format
#		how to make reference genome index: https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#the-bowtie2-build-indexer
# -1 and -2	forward and reverse read fastq files
# -S	name of output file (will be in .sam format)

## Alternative: if you want to run all pipeline per sample, this would be the part for mapping

# set up variables

# Path to reference genome
REF="/scratch/projects/path_to_reference_genome/Gh38"

# forward and reverse reads
R1=${NAME}_1.fastq
R2=${NAME}_2.fastq

# fastq files after adapter trimming step
TRIM_R1=${NAME}_trim_1.fastq.gz
TRIM_R2=${NAME}_trim_2.fastq.gz

# aligned files
SAM=${NAME}.sam

echo "starting reference alignment"
bowtie2 -p 8 --very-sensitive --no-mixed --no-discordant -I 10 -X 700 -x ${REF} -1 $TRIM_R1 -2 $TRIM_R2 -S ${SAM}
echo "finished aligning to reference"
