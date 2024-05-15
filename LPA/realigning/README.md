# Step LPA-2: Realigning

## External Software Requirement:

 - C++ part: same to normalize_mosdepth in main pipeline
 - SAMTOOLS 1.17 called from Linux environment variable "module load"


## Input:

- Reference of region of interest
   - fasta file, the output of extract_reference, the path is hard-coded in cpp file line 71

## External Data Requirement:

 - bam_list, used for loop over all the bam files (our bam saved under a child path)
    - provided in count_read

## Format:

 - cpp script: same to the CPP files in main pipeline,
 - SGE script / SLURM script (i.e. distributed computing scheduler script)
    - it's actually not distributed, but it took quite a long time to run.

## File:

 - realign_GRCh37.cpp
   - All the positions are hardcoded into the cpp file (see hardcoded_hg19_position.txt) because:
      1. the input stream is occupied in the task
      2. these don't need to be modified, it's just for LPA genes.
 - run_realign_GRCh37_SGE.sh/run_realign_GRCh37_SLURM.sh

## Running:

1. Please make sure the path to the output of extract_reference is correct, the path is hard-coded in cpp file line 71
2. compile the cpp script into binary first, and make sure the path at .sh file line 52 is the one you compiled.

For each .bam file:
3. samtools view 6:161032032-161068032 (LPA gene)
4. Only using flag 99,147,83,163 in samtools view.
  - It is different from the main pipeline (rest four 81/97/145/161 are not included).
5. sort the input by name and position
6. push the value into cpp binary file
7. The binary file run a read-depth based algorithm and compute a score for each repeat.
  -	If first position(161032593) is the smallest, it’s a 1B_KIV3
  -	If fifth position(index 4, 161054784) is the smallest it’s 1B_KIV2
  -	If first and fifth are equal and both smallest, it’s 1B_tied
  -	Otherwise if there’ a value smaller than 0 and 4, or all of them larger than 9999, it’s 1A

## Output:

A tab-split value txt file with 4 columns

## Corresponding original file:

code/genotyping_opt/LPA/realign.cpp
code/genotyping_opt/LPA/run_realign.sh
