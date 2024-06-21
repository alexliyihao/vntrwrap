# vntrwrap

This repository is an implementation(wrapper) for [*Protein-coding repeat polymorphisms strongly shape diverse human phenotypes paper code*](https://www.science.org/doi/10.1126/science.abg8289)'s code repo for easier re-use.

This pipeline is designed for the collaborated study project on WHICAP cohort [paper: TOBE ADDED]() of Reyes-Soffer Lab, Department of Preventive Medicine, Columbia University Irving Medical Center(CUIMC) and Badri Vardarajan Lab, The Gertrude H. Sergievsky Center, Department of Neurology CUIMC. We are working on **Hg19** rather than the **Hg38** references and originally dealing with LPA KIV2 repeats on Chr 6. So it's not exactly working the same as the source pipeline, a part of them has been modified, and the LPA-specific pipeline can be found at *LPA* path.

Current version is the raw code intended to be run on SGE or SLURM cluster with a minimal documentation for running.

## Steps:

1. count_read
2. mosdepth
3. normalize_mosdepth
4. find_neighbors
5. neighbors_normalization

## Updates:

2024/05/14: added SLURM version due to the scheduler transition of CUIMC neurology cluster from SGE to SLURM.

## TODO:

1. Adding makefile
2. Adding an argparse interface (either from cpp or python).
3. Adding phasing and inputation (We have no spare valid access to UK biobank running right now)

## Contact

If need any help or explanation, want to collaborate, or can help with any parts, please email Yihao Li (yl4326@cumc.columbia.edu), and CC both Dr. Gissette Reyes-Soffer (gr2104@cumc.columbia.edu) and Dr. Badri Vardarajan (bnv2103@cumc.columbia.edu).
