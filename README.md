## CUT-RUN-analysis
Pipeline for analysis of CUT&RUN and CUT&TAG data - can also used for ChIPseq (SE, PE) after modifying few parameters.

## Summary of pipeline:

1. Trim all reads to remove adaptors and sequences with very low quality (bbduck or trimgalore). 
2. QC check (fastqc)
3. mapping (bowtie2) - (previous generation of genome indexes and getting latest genome and annotations)
4. optional: multiQC (to get report on general stats on all libraries)
5. filtering with samtools if needed (and indexing, sorting,etc)
6. statistics on libraries (deeptools) # fancier stats to look at quality of libraries. 
7. making bigwig files (deeptools) # normalization, excluding duplicates, others.
8. Plotting (deeptools)
9. visualization (higlass)

Intended for Biologists w/o little experience in Computational Biology - I provide all scripts in the right order for you to be able to run the entire analysis. I provide notes below on main parameters used.
I typically batch jobs on the cluster, so most scripts are intended for this. Running one by one is highly recommended if you are learning and this is your first time running an entire pipeline for analysis. It will get you familiarized with many things in comp.bio.
I have used this pipeline for C. elegans, mice and human libraries.

## Requirements:
Bowtie2 (Alignment and Filtering) - generates .BAM files and .SAM files
samtools (.sam .bam manipulation)
DeepTools ((Normalization, graphs) - generates .bw (bigwig) files
MACs (peak calling)
Homer (Normalization, peak identification, peak annonation) - generates .bw (bigwig) files
HiGlass (Visualization) 
trimgallore #optional
bbduck # I prefer to trim with bbduck. Much faster.
fastqc # optional
multiqc # optional. very useful if running too many samples.

## Installations:

conda install -c bioconda bowtie2 

conda install -c bioconda deeptools

Higlass installation: https://docs.higlass.io/

conda install -c maximinio macs3 (https://pypi.org/project/MACS3/)

## About Bowtie2 Alignment:
My prefered way:

bowtie2 --very-sensitive --no-mixed --no-discordant -I 10 -X 700 -x /path_to_genome/Gh38 -1 /path_to_FASTQ/Trimmed/*R1_val_1.fq.gz -2 /path_to_FASTQ/Trimmed/*_R2_val_2.fq.gz -S library_name.sam

Bowtie 2 runs a little faster in --no-mixed mode, but will only consider alignment status of pairs per se, not individual mates.
--no-discordant By default, bowtie2 looks for discordant alignments if it cannot find any concordant alignments. A discordant alignment is an alignment where both mates align uniquely, but that does not satisfy the paired-end constraints (--fr/--rf/--ff, -I, -X). This option disables that behavior.

Alternative: Map in Bowtie2 using default settings (sensitive) 

$ bowtie2 -X 700 --no-mixed --no-discordant ..... 

what are these settings?:
By default, Bowtie 2 performs end-to-end read alignment. That is, it searches for alignments involving all of the read characters (this is why my code does not specify --end-to-end).

From bowtie2 manual:

-D <int>
        Up to `<int>` consecutive seed extension attempts can "fail" before Bowtie 2
        moves on, using the alignments found so far.  A seed extension "fails" if it
        does not yield a new best or a new second-best alignment.  This limit is
        automatically adjusted up when -k or -a are specified.  Default: 15.

-R <int>
        `<int>` is the maximum number of times Bowtie 2 will "re-seed" reads with
        repetitive seeds. When "re-seeding," Bowtie 2 simply chooses a new set of reads
        (same length, same number of mismatches allowed) at different offsets and
        searches for more alignments.  A read is considered to have repetitive seeds if
        the total number of seed hits divided by the number of seeds that aligned at
        least once is greater than 300.  Default: 2.
-L <int>
        Sets the length of the seed substrings to align during multiseed alignment.
        Smaller values make alignment slower but more sensitive. Default: the
        `--sensitive` preset is used by default, which sets `-L` to 22 in
        `--end-to-end` mode and to 20 in `--local` mode.
-N <int>
        Sets the number of mismatches to allowed in a seed alignment during multiseed
        alignment. DEFAULT 0.

-i <func>
Sets a function governing the interval between seed substrings to use during multiseed alignment. 

[http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#presets-setting-many-settings-at-once]

## Notes:

A pair that aligns with the expected relative mate orientation and with the expected range of distances between mates is said to align “concordantly” [I DO NOT SPECIFY EXPECTED INSERT SIZE IN DEFAULT MAPPING MODE]. If both mates have unique alignments, but the alignments do not match paired-end expectations (i.e. the mates aren’t in the expected relative orientation, or aren’t within the expected distance range, or both), the pair is said to align “discordantly”. Discordant alignments may be of particular interest, for instance, when seeking structural variants.

By default, Bowtie 2 searches for both concordant and discordant alignments, though searching for discordant alignments can be disabled with the --no-discordant option. ## I have also mapped w/o using --no-discordant option and after .bam file generation, filtered out discordant pairs (I find easier to just remove from the beggining, but if for some reason those discordant reads may be needed, better to run it in this alternative way).
  
Bowtie2 defaults to a minimum and maximum insert size of 0 and 500. #which means this is my X = 500 using default settings. With X=700 I am allowing for mono, di and tri nucleosome sizes (150,300,450). 
# I have tested many parameters and basically using default mode is fine for mapping. When generating files for analysis I can add restrictions or not. I found that adding the restriction for fragment length 100-250bp does improve the signal, but there are no significant changes if the libraries are very good in # reads and complexity. I prefer to go with 700bp since for a lot of antibodies it is hard to obtain traces that are clean mononucleosomes. in other words, for several antibodies (histone marks) you almost always see a ladder. Therefore, I prefer to use X=700 and include that data.
I have seen in some papers people using X=1000. I almost never use this, but also have not tested properly to compare.

alternative mapping that could be done:

$ bowtie2 -X 1000 --no-mixed --no-discordant -x 
  
## About duplicates
  
from dovetail "default action in CUT&RUNTools is to retain duplicate reads, and users can choose to remove duplicates at their own discretion. We recommend users to be aware of low complexity of libraries with high duplication rates, as these may indicate a poor quality preparation" - I filter duplicates. I have tested and doing this does not change my data - but on average my libraries do not have high % duplicates (I tend to run 11-13 PCR cycles).
  my general suggestion is: if < 45% duplicates, then ignore them, otherwise do not (but think twice about this library- perhaps better to repeat with lower # PCR cycles, and/or more cells/tissue initial input).
  
## About Analysis using Deeptools
  
I feel Homer is limited in terms of statistical analysis and graphs when handling replicates. DeepTools has more options. Also, with Homer I tried visualizing raw data (which does not do any of the peak centering, extending, etc) and also normalized data (against input). The second one was my first set of bigwig files but I had the issue that in HiGlass the scale can't be force to show from zero above and my Homer results were in log2 scale (showing negatives). So I created the unnormalized version, raw scale. This looks better in Higlass and it looks like all the signals are similar (perhaps not all, but visually looks good). When looking at the input raw data looks messy though. However, the fact I am still getting all the peaks after normalization indicates normalization is not messing my data, and I like to normalize anyway. So if using Homer in the future I should prepare visualization files (Bedgraphs or Bigwig), normalized agains input but in raw scale (no Log2).
Now, Deeptools has the option of normalizing differently, and I liked that. Also, it has the option of generating the visualization files, only normalized by sequencing depth (using the RPGC method), not input, so it is possible to compare between experiments, plus you can center, extend reads, which visually is much better and makes a lot of sense to me since it makes interpretation easier. To create "unnormalized" (just by seq depth) files I used bamCoverage, and to generate normalized files against input I used bamCompare. 

