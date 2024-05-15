# Step LPA-3: Computing the diploid copy number

## External Software Requirement:

None

## Input:

 - Aligning result: txt file, the output of *realigning*
 - find_neighbors result: txt.gz file, the output of main pipeline *find_neighbors*

## Format:

Linux Bash script (it's running relatively fast and it's unnecessary to put it into scheduler)

## File:

 - run_compute_dipCN_GRCh37.sh

## Running:

This step is equivalent to main pipeline *neighbors_normalization* but slightly different:

For each type of KIV2 repeat ("1B_KIV3","1B_notKIV3","1B","1A"):
   - find top 200 neighbors(rather than 300 in the main pipeline *neighbors_normalization*)
   - normalize the count from realigning by **both individual scale across repeats in this person** and **repeat-wise mean scale across neighbors**.

## Output:

4 individual files from OUTPUT_PREFIX:

 - OUTPUT_PREFIX.exon1A.dipCN.txt
 - OUTPUT_PREFIX.exon1B.dipCN.txt
 - OUTPUT_PREFIX.exon1B_KIV3.dipCN.txt
 - OUTPUT_PREFIX.exon1B_notKIV3.dipCN.txt

## Corresponding original file:

code/genotyping_opt/LPA/run_compute_dipCN.sh
