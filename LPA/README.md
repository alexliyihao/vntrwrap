# LPA sub-pipeline

This part is working on LPA genes' KIV-2 repeats in paper [TODO: add link](https:?) with Dr.Reyes-Soffer Gissette and Dr.Badri Vardarajan.


To work on LPA, this sub-pipeline should be run after the main pipeline step 4 (find_neighbors), and main pipeline step 5 (neighbors_normalization) can be omitted

## Steps:

- LPA-1: extract_reference
- LPA-2: realigning
- LPA-3: compute_diploid_copy_number
- LPA-4: **(deprecated and not included)** the PSV phasing imputation
    - it's not working quite well on our data on either the first sub-step or both sub-steps and we decided to deprecate the step for now.
- LPA-5: estimate_KIV_CN
