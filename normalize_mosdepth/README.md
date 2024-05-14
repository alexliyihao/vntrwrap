# Normalize Mosdepth

## External Software Requirement (see external_source/cpp_external.txt):

 - g++
 - [boost 1.58](https://www.boost.org/users/history/version_1_58_0.html)
 - [Eagle 2.4.1](https://github.com/poruloh/Eagle)(only the source code, not running)

## External Data Requirement:

 - Repeat Mask:
   - Mukamel provided a Hg38 version and it can be translated into Hg37
      - repeat_mask_list**.ucsc_bed files
   - We are currently using from https://genome.ucsc.edu/cgi-bin/hgTables, GRCh37 chr6:1-171,115,067 group == repeats, which is way larger than the hg38 version
      - GRCh37_repeat_mask_chr_6_cleaned.ucsc_bed

## Format:

C++ script, the g++ compile and running is at the top of the script

## File:

 - normalize_mosdepth_inflow_rewritten.cpp fixed to chr 6, check line 123-130

## Running:

 - A very weird normalization process...

## Corresponding original file:

code/WES_read_depth/normalize_mosdepth.cpp

## Note:

 - fixed a bug that the first read has 100x as much as others weights during the mean_depths computing
 - *pyextern* is an attempt of python argparse, not tested yet
