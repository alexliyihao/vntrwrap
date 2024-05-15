# Step 5: Neighbor-wise normalization

## External Software Requirement:

None

## External Data Requirement:

 - IBD2R file provided in the source file (actually the IBD2R is not necessary, we only used first 3 columns)

## Format:

Shell script (no need for distributed computing)

## File:

  - GRCh37_run_neighbor_norm_VNTR_reads_20240103_old_mask.sh used in our pipeline
  - GRCh37_run_neighbor_norm_VNTR_reads_cleaned.sh should be better but not tested

## Running:

1. Read the IBD2r file, for each GENE get row, get Chr, pos start and pos end,
2. create a name VNTR=prefix_BP_START_BP_END
3. For all the first 300 neighbors, compute mean(scales/dist)
4. normalized_reads = (count from VNTR_counts)/(personal scale from neighbor file)/mean(scales/dist)

## Corresponding original file:

code/WES_read_depth/run_neighbor_norm_VNTR_reads.sh
