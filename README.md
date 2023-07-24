# cutnrun_repeats
snakemake pipeline to align cutnrun data to repetitive regions of the genome:
1. trim cutnrun reads w/o deduping
2. pair reads using bbmerge
3. align to t2tv2 genome with bowtie and remove multimappers using samtools
4. align to ecoli genome to calculate normalization factor
5. feed scaled data to macs2
6. call peaks, calculate FE
7. intersect peaks with unique 100-mers to get high confidence dataset

   IMP. NOTE: Scale factor calculation is slightly different from cutntag with SEACR for peak calling. Calculation of normalization factor in this pipeline utilizes total mapped reads instead of a "constant factor", since MACS2 subsequently scales for total reads
