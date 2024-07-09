## Paternal obesity-induced changes in sperm chromatin accessibility and have mild effect on offspring metabolic health.
DOI: https://doi.org/10.1016/j.heliyon.2024.e34043 


R code used to perform peak differential accessibility analysis and prepare inputs for HOMER TF motif analysis.

#### Generate `consensus_counts.tsv` and `consensus_overlaps.tsv` for ATAC-seq analysis:

```
#create consensus matrix from all peak callings
cat *peaks.narrowPeak | sort -k1,1 -k2,2n | bedtools merge > consensus_pe.bed

#create counts per peak of consensus matrix per sample
consensus_coverage.R -t 20 consensus_pe.bed <BAM files> > consensus_counts.tsv

#create presense of peak (0/1) of consensus matrix per sample
consensus_coverage.R -t 20 -m 1 consensus_pe.bed *peaks.narrowPeak > consensus_overlaps.tsv
```

For downstream analyses check the R scripts included.
