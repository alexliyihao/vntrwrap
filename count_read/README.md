# Count the reads

## External Software Requirement:

 - SAMTOOLS 1.17 called from Linux environment variable "module load"

## External Data Requirement:

 - IBD2R file provided in the source file (actually the IBD2R is not necessary, we only used first 3 columns)
 - bam_list, used for loop over all the bam files (our bam saved under a child path)

## Format:

SGE script / SLURM script (i.e. distributed computing scheduler script)

## File:

 - count_read_chr_6_only_SGE.sh used in SGE
 - count_read_chr_6_only_SLURM.sh used in SLURM

## Running:

Calling samtools view with following flags(Details can be seen [here](https://broadinstitute.github.io/picard/explain-flags.html)): 81, 83, 97, 99, 145(mate of 97), 147(mate of 99), 161(mate of 81), 163(mate of 83)

## Corresponding original file:

code/WES_read_depth/run_download_count_VNTR_reads.sh
