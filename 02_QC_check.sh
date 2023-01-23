

########## Step 2: QC check (using fastqc) ##########

# if you want to batch this job as a loop (for all fastq files)
# create .sh file with any text editor in command line (or upload your file)

vim QCcheck.sh

##!/bin/bash - include here your bash header (cluster-specific)
# load environment or modules

echo "running fastqc"
# Note: fastqc wants the adapters.fa in a different format than bbduk.
# To recreate "tab_adapters.fa", run:
# $ awk 'BEGIN{RS=">";OFS="\t"}NR>1{print $1,$2}' adapters.fa > tab_adapters.fa
# to convert fasta format to tab delimited format

for i in /data/your_path/FASTQ/trimmed/*_1.fastq.gz ; do
name=$(basename ${i} _1.fastq.gz) ;
fastqc -a tab_adapters.fa ${name}_1.fastq.gz ${name}_2.fastq.gz ;
done

## Alternative: if you want to run all pipeline per sample, this would be for QC check

# set up variables
# forward and reverse reads
R1=${NAME}_1.fastq
R2=${NAME}_2.fastq

# fastq files after adapter trimming step
TRIM_R1=${NAME}_trim_1.fastq.gz
TRIM_R2=${NAME}_trim_2.fastq.gz

fastqc -a tab_adapters.fa ${TRIM_R1} ${TRIM_R2}
echo "fastqc finished"
