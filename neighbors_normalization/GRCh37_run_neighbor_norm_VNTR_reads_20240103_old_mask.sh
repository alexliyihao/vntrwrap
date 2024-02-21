COL=1

while read LINE
do
    COL=$[COL+1]
    ARR=( $LINE )

    CHR=${ARR[0]}
    BP_START_HG37=${ARR[1]}
    BP_END_HG37=${ARR[2]}
    # modify the name here, it's just a name
    VNTR=hg37_20240103_${CHR}_${BP_START_HG37}_${BP_END_HG37}

    echo $COL $VNTR

    N_NBR=300

    awk -v Nnbr=$N_NBR -v col=$COL -v chr=$CHR '
ARGIND==1 {reads[$1]=$col}
ARGIND==2 {
  sum=0; num=0;
  for (i=1; i<=Nnbr; i++) {
      sum+=reads[$(3*i)]/$(3*i+1);
      num++;
  }
  norm_reads = reads[$1]/$2/(sum/num);
  print $1,norm_reads
}' \
    /mnt/mfs/hgrcgrid/shared/LPA_analysis/VNTR_pipeline/analysis_results/VNTR_counts/VNTR_counts_chr6_GRCh37_second_run.txt \
    <(zcat /mnt/mfs/hgrcgrid/shared/LPA_analysis/VNTR_pipeline/analysis_results/find_neighbor/GRCh37_find_neighbor_20240102_old_mask.zMax2.txt.gz) \
    > /mnt/mfs/hgrcgrid/shared/LPA_analysis/VNTR_pipeline/analysis_results/neighbor_norm_VNTR_READ/$VNTR.dipCN.txt

done < <(tail -n +2 /mnt/mfs/hgrcgrid/shared/LPA_analysis/VNTR_pipeline/pipeline_columbia_local/sources/WES_read_depth/GRCh37_chr6_possible_coding_vntr_regions.IBD2R_gt_0.25.txt)
