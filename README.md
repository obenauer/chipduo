# ChipDuo
ChipDuo detects peak differences in ChIP-seq and ATAC-seq experiments.  The input files should 
be aligned BAM files representing replicates of 2 groups: some kind of treatment (t) 
and a set of untreated controls (c).  "Treatment" is used loosely here to mean any kind 
of experimental samples: they could be drug-treated, or mutants, or knockout mice, or 
have different ages, for example.  ChipDuo divides the BAM files into millions of windows 
of the human or mouse genome (or any other species with chromosomes defined in the BAM files) 
and compares the read counts and read depths in each window.  To reduce false positives, 
every replicate in the treatment group must be higher (or lower) than every replicate in 
the control group to be considered significant.  The BAMs are normalized by their total 
mapped reads before being compared, and significant results are written out with their log2 
fold changes, raw p-values, adjusted p-values, and genes nearest the peaks.
