# Step LPA-1: Extract Reference

## External Software Requirement:

 - Bedtools/2.30

## External Data Requirement:

 - Human Genomic reference file, too large to included in the repository:
    - source file using GRCh38_full_analysis_set_plus_decoy_hla.no_chr.fa for hg38 and we tried both 1000 genome and hs37d5 reference.
    - It is confirmed that the hg19 and hg38 reference produce identical reference on the region.
    - It is confirmed that g1k_v37 and hs37d5 produce identical output when working on 6:161032031-161068032
 - bed file: indicating the necessary position need
    - ref_hg38.bed is from source file
    - ref_hg37.bed is converted from hg38, and change "chr6" to '6'


## Format:

Linux Bash script (it's running relatively fast and it's unnecessary to put it into scheduler)

## File:

 - extract_reference.sh

## Running:

 - simply extract a part of reference by *bedtools getfasta*

## Output:

 - fasta file, reference for region of interest.

## Corresponding original file:

code/genotyping_opt/LPA/ref_hg38.bed
code/genotyping_opt/LPA/extract_KIV2_sequence.sh
