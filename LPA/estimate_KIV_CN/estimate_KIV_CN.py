import pandas as pd

def extract_VNTR_estimate(
    prefix,
    source_path = "/mnt/vast/hpc/bvardarajan_lab/LPA_analysis/VNTR_pipeline/analysis_results/",
    test = True):

    #exon1A1 = pd.read_csv(
    #    source_path + f"estimate_KIV2_length_PSV/{prefix}.exon1A1.dipCN.txt",
    #    index_col = 0, sep = " ")
    #exon1A2 = pd.read_csv(
    #    source_path + f"estimate_KIV2_length_PSV/{prefix}.exon1A2.dipCN.txt",
    #    index_col = 0, sep = " ")
    exon1A = pd.read_csv(
        source_path + f"estimate_KIV2_length/{prefix}/{prefix}_LPA.exon1A.dipCN.txt",
        index_col=0,
        sep=" ")
    exon1B = pd.read_csv(
        source_path + f"estimate_KIV2_length/{prefix}/{prefix}_LPA.exon1B.dipCN.txt",
        index_col=0,
        sep=" ")
    complete = pd.concat(
        [exon1A.rename(columns = {"dipCN": "exon1A"}),
         #exon1A1.rename(columns = {"dipCN": "exon1A1"}),
         #exon1A2.rename(columns = {"dipCN": "exon1A2"}),
         exon1B.rename(columns = {"dipCN": "exon1B"})], join = "inner", axis = 1)
    # in our pipeline we don't have fine phasing as Mukamel's, so this test
    complete["dip_estimate"] = 34.9 * complete["exon1A"] + 5.2 * complete["exon1B"] - 1
    complete["estimate"] = complete["dip_estimate"]/2
    return complete
