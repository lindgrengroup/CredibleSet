# CredibleSet
R script for credible sets

This script works with 2 input files:

1 - Regions for credible sets determined by chromosomes and base pair positions e.g. 1-Mb regions. One row for each locus.

2 - Genome-wide summary statistics.

You can also determine the alpha in the R command line. e.g. 0.99 for 99% credible sets. Please see the examples given below the script.


Output is a file with full sets and all the statistics and another file with a summary table of the credible sets showing the number of variants and the distance for comparison.
