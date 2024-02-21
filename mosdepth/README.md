# Mosdepth

## External Software Requirement:

 - [mosdepth_v0.2.5](https://github.com/brentp/mosdepth/releases?page=2)

##External Data Requirement:

 - Human Genomic reference file, too large to included in the repository:
    - source file using GRCh38_full_analysis_set_plus_decoy_hla.no_chr.fa for hg38
    - I'm currently using human_g1k_v37.fasta
 - bam_list, used for loop over all the bam files (our bam saved under a child path)
    - provided in count_read

## Format:

SGE script (i.e. distributed computing)

## File:

 - run_mosdepth.sh used in out pipeline, fixed the chr to 6, check line 77,78

## Running:

1.	Split the BAM list into small batches(using batch_size = 25)
2.  For the file in each batch run mosdepth

## Corresponding original file:

code/WES_read_depth/run_mosdepth.sh

## Notes:

Some mosdepth results are all zero and caused un-preventable NaN in the later pipeline. Please consider drop them
(I mannually dropped one and ran a batch of 9 bams, looking good)
