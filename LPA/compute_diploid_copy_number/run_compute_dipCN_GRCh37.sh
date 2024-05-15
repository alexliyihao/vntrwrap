COUNT_RESULT="/mnt/mfs/hgrcgrid/shared/LPA_analysis/VNTR_pipeline/analysis_results/estimate_KIV2_length/GRCh37_counts.1B_KIV3_KIV2_tied.1A.txt"
FIND_NEIGHBOR_RESULT="/mnt/mfs/hgrcgrid/shared/LPA_analysis/VNTR_pipeline/analysis_results/find_neighbor/GRCh37_find_neighbor_20240102_new_mask.zMax2.txt.gz"
OUTPUT_PREFIX="/mnt/mfs/hgrcgrid/shared/LPA_analysis/VNTR_pipeline/analysis_results/estimate_KIV2_length/GRCh37_covid/GRCh37_covid_LPA"
N_NBR=200

for EXON in 1B_KIV3 1B_notKIV3 1B 1A
do
    awk -v Nnbr=$N_NBR -v exon=$EXON '
BEGIN {print "ID","dipCN"}
ARGIND==1 {
  if (exon=="1B_KIV3") reads[$1]=$2;
  if (exon=="1B_notKIV3") reads[$1]=$3+$4;
  if (exon=="1B") reads[$1]=$2+$3+$4;
  if (exon=="1A") reads[$1]=$5;
}
ARGIND==2 && reads[$1]!="" {
  sum=0; num=0;
  for (i=1; i<=Nnbr; i++)
    if (reads[$(3*i)]!="") {
      sum+=reads[$(3*i)]/$(3*i+1);
      num++;
    }
  print $1,reads[$1]/$2/(sum/num); #*38.5-2;
}' \
    $COUNT_RESULT \
    <(zcat $FIND_NEIGHBOR_RESULT) \
    > ${OUTPUT_PREFIX}.exon$EXON.dipCN.txt
done
