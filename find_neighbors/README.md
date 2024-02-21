# Find Neighbors

## External Software Requirement:

 - C++ part: same to normalize_mosdepth

## External Data Requirement:

None

## Format:

C++ script, the g++ compile and running is at the top of the script

## File:
 - find_neighbor.cpp

## Running:
 - A naive pairwise distance computing process - see the notes sent before.

## Corresponding original file:

code/WES_read_depth/find_neighbors

## Besides:

 - fixed a bug that the first read has 100 as much as others weights during the mean_depths computing
 - *pyextern* is an attempt of python argparse, not tested yet
