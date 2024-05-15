REFERENCE_PATH="/mnt/vast/hpc/bvardarajan_lab/LPA_analysis/VNTR_pipeline/dataset/hs37d5/hs37d5.fa"
BED_PATH="/mnt/vast/hpc/bvardarajan_lab/LPA_analysis/VNTR_pipeline/pipeline_columbia_local/sources/genotyping_opt/ref_hg37.bed"
OUTPUT_PATH="ref_hg37_hs37d5.fasta"

module load Bedtools/2.30
bedtools getfasta \
    -fi $REFERENCE_PATH \
    -bed $BED_PATH \
    > $OUTPUT_PATH
