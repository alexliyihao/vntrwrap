#!/bin/bash -l

## Section 1: Header ----------------------------------------------------

# Specify name to be used to identify this run
#SBATCH --job-name=realign_WHICAP

# Memory requirement
#SBATCH --mem=2G

# Change directory to the current
#SBATCH --chdir=./

# Specify that bash shell should be used to process this script
#SBATCH #!/bin/bash

# Specify the outerr file
#SBATCH --output=/mnt/vast/hpc/bvardarajan_lab/LPA_analysis/VNTR_pipeline/SLURM_logs/VNTR_realign/%x/%j.out
#SBATCH --error=/mnt/vast/hpc/bvardarajan_lab/LPA_analysis/VNTR_pipeline/SLURM_logs/VNTR_realign/%x/%j.err

# Partition, use yours
#SBATCH --partition=GEN

## Section 2: Path Settings & module loading ----------------------------------------------------

module load SAMTOOLS/1.17

SOURCE_BAM_PATH="/mnt/mfs/hgrcgrid/data/whicap/WHICAP_WES/BAM/washeiDragenBamsList/washeiBamsUpdate2/BQSR/bqsrRGbam/"
SOURCE_LIST="/mnt/mfs/hgrcgrid/shared/LPA_analysis/coassin_pipeline/data_inflow/bam_list_VNTR_N=3915.txt"
OUTPUT="/mnt/mfs/hgrcgrid/shared/LPA_analysis/VNTR_pipeline/analysis_results/estimate_KIV2_length/GRCh37_counts.1B_KIV3_KIV2_tied.1A.txt"

## Section 3: Running -----------------------------------------

start_time=$(date +%s)

rm -f $OUTPUT

# read the whole list
readarray -t INPUTFILES < $SOURCE_LIST
echo "Setting: ${#INPUTFILES[@]} files detected"

for FILE in "${INPUTFILES[@]}"
do
    # get bam path
    BAM_PATH=$SOURCE_BAM_PATH/$bam_file
    # get the corresponding bai path
    BAI_PATH=${BAM_PATH%?}i
    # This ID is along with "hispanic/" prefix!
    ID=$FILE
    if [ -f $BAM_PATH -a -f $BAI_PATH ]; then
        echo -n $ID >> $OUTPUT
        samtools view $BAM_PATH 6:161032032-161068032 --verbosity 0 \
	    | awk '($2==99||$2==147||$2==83||$2==163) {print $1,$4,$6,$10,$11}' \
	    | sort -k1,1 -k2,2n \
      | ./realign_GRCh37 >> $OUTPUT
    fi
done
echo "Successfully completed VNTR alignment, time taken: $(($(date +%s)-$start_time))"
