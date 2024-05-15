# Step LPA-5: estimate KIV-CN

## External Software Requirement:

  - Python
  - Pandas

## Format:

Python script with a function
 - You don't have to use python here, we use python simply because our analysis is majorly operated on python
 - See *Running* section for the formula and run it in whatever way you need.

## File:

 - estimate_KIV_CN.py

## Running:

 - extract the 1A output from {prefix}_LPA.exon1A.dipCN.txt and the 1B output from {prefix}_LPA.exon1B.dipCN.txt
 - estimate = 34.9 * exon1A + 5.2 * exon1B - 1

## Output:

 - pd.DataFrame

## Corresponding original file:

None: we are not using the PSV correction step (whose should have been the step LPA-4) thus simply compute the weighted estimate is file.
