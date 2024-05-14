# Mosdepth

## External Software Requirement:

 - [mosdepth_v0.2.5](https://github.com/brentp/mosdepth/releases?page=2)

## External Data Requirement:

 - Human Genomic reference file, too large to included in the repository:
    - source file using GRCh38_full_analysis_set_plus_decoy_hla.no_chr.fa for hg38
    - The SGE version using human_g1k_v37.fasta (line 32-34)
    - The SLURM version using hs37d5 reference (line 32-33)
    - It is confirmed that g1k_v37 and hs37d5 produce identical output when working on Chr6
 - bam_list, used for loop over all the bam files (our bam saved under a child path)
    - provided in count_read

## Format:

SGE script / SLURM script (i.e. distributed computing scheduler script)

## File:

 - run_mosdepth_SGE.sh used in out pipeline on SGE scheduler,
 - run_mosdepth_SLURM.sh used in out pipeline on SLURM scheduler,
 - fixed the chr to 6, check line 82-86 for either script

## Running:

1.	Split the BAM list into small batches(using batch_size = 25)
2.  For the file in each batch run mosdepth

## Corresponding original file:

code/WES_read_depth/run_mosdepth.sh

## Notes:

Some mosdepth results are all zero and caused un-preventable NaN in the later pipeline. Please consider dropping them
(ran a "unfilled" batch won't affect the output)
