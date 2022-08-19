# illumuna-error

Next generation sequencing methods usually have noticeable error rate. Pair-end sequencing can generate reverse complementary read pair with overlap. Researchers could identify and correct the system errors of sequencing by analyse with the overlaps. Besides, investigating the position with what characteristic could lead to high incidence of system error is benefit to gain accurate and reliable data.

This is a python script uses aligned BAM file as input to merge up paired end read by overlap while correct the system errors. At the same time, it could record and generate statistical analysis (include visualization) of system error.

Usage (arguments within squared brackets are optional):
  
  python |illumina_error_pair_read.py| --input_bamfile <input.bam> --ref <reference.fasta> [--index <index.bai>] [--output <output filename>]
  
If running the script on local, these are required packages or tools:
1. Samtools 
2. Pysam (version 0.18.0) (https://github.com/pysam-developers/pysam)
3. WebLogo (version 3.7.12) (http://weblogo.berkeley.edu/)
4. R script
5. ggplot2 (version 3.3.5)
6. RColorBrewer (version 1.1-2)
7. Reshape2 (1.4.4)
