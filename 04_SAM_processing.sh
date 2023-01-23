
########## Step 4: Convert sam to bam, filter properly paired reads and index ##########

# Run this in folder where .sam files are located
# if you want to batch this job as a loop (for all fastq files)
# create .sh file with any text editor in command line (or upload your file)

vim 04_sam2bam.sh

##!/bin/bash - include here your bash header (cluster-specific)
# load environment or modules

module load samtools
for file in *.sam; do
name="${file%.sam}"
samtools view -S -b "$name".sam > "$name".bam
done

# Filter to keep ONLY properly paired reads (and all concordant alignments, which means within expected insert size/distance and correct directions).

for file in *.bam; do
name="${file%.bam}"
samtools view -bf 0x2 "$name".bam > "$name".pp.bam
done

#The -f 0x2 option corresponds to the bitwise flags that specify that reads need to be properly paired. Proper pairing means reads are in Read1 forward, Read2 reverse orientation or Read1 reverse, Read2 forward orientation.
# $ samtools view -q 30 -f 0x2 -b -h in.bam > out.bam
# if mapping already specified filters to keep only properly paired, then this step is not needed. 

# Sort BAM files 
# to view a single sample: samtools view CRlibF.bam | head
# BAM file must be sorted such that the alignments occur in “genome order”. That is, ordered positionally based upon their alignment coordinates on each chromosome.

for file in *.pp.bam; do
name="${file%.bam}"
samtools sort "$name".bam > "$name".st.bam
done

# checking that sorted.SAM are truly sorted by name
# $ samtools view lib1.st.sam | head
# $ samtools view -H lib1.pp.st.bam
# @HD	VN:1.0	SO:coordinate #this is correct. The unsorted says unsorted

# Create .bam index file 

for file in *.st.bam; do
samtools index $file
done

# move all st.bam and bam index files to a separate folder for processing.


## Alternative: if you want to run all pipeline per sample, this would be the part for sam/bam processing

########## Step 4: Convert sam to bam ##########

# set up variables
# aligned files
SAM=${NAME}.sam
BAM=${NAME}.bam
PP_BAM=${NAME}.pp.bam # only properly paired reads
SORT_PP_BAM=${NAME}.pp.st.bam # sorted properly paired reads 

echo "converting sam to bam"
samtools view -S -b ${SAM} > ${BAM}
# -S	auto-detect input format
# -b 	output a bam file
echo "bam file made"

date

########## Step 5: Keep only properly paired reads ##########
echo "Keep only properly paired reads"
samtools view -bf 0x2 ${BAM} > ${PP_BAM}
# -bf 0x2	is a combination of two options:
# -b	 bam output file format
# -f 0x2	keeps only reads that are part of a pair that aligned in a paired-end fashion
echo "unpaired reads discarded"

date

########## Step 6: Sort properly paired bam file ##########
echo "sort properly paired bam file"
samtools sort ${PP_BAM} > ${SORT_PP_BAM}
echo "properly paired bam file sorted"

date

########## Step 7: Index sorted bam file ##########
echo "index sorted bam file"
samtools index ${SORT_PP_BAM}
echo "finished indexing bam"

date