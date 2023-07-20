## Chromosomal positioning and epigenetic architecture influence DNA methylation patterns triggered by galactic cosmic radiation.

## R Scripts:
1. DNA.methylation.R
2. Genes.for.hic.astronauts.R
3. Microarray.R
4. RNAseq.R

R Scripts contain all necessary code to perform the analysis and generate all figures included in the article entilted: "Chromosomal positioning and epigenetic architecture influence DNA methylation patterns triggered by galactic cosmic radiation" published in XXX.

## Files necessary to run code:
1. GSE108187_processed_matrix.csv (you can download it here: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE108187)
2. GSE108187_colnames.xlsx
3. DMPs.annotations.all.probes.2000.rdata (you can generate it via script DNA.methylation.R - it takes several hours and 200MB of disk space)
4. DMRs.newest.GRanges.rdata
5. DMRs.newest.annotations.rdata
6. Den.newest.GRanges.rdata

The rest of files needed for certain step are avaiable on publicly available respositories and are decribied in detail in the article. 

## Table of contents:

### DNA.methylation.R (Running time when all needed files are provided ~ 20 minutes)
1. Install and load packages
2. Import and prepare primary files
3. Differential methylation analysis
4. Pick significant DMPs
5. DMRs Building
6. Simulate DMRs building with rows permutation
7. Load and prepare lists of DMRs
8. Acute methylation change versus control between DMPs
9. % of DMPs in DMRs
10. DMPs which changed at least of 0.58 or -0.58 M Value versus control
11. Detailed methylation patterns
12. Select DMRs, DMPs in or out of DMRs and all DMPs for further analysis
13. Annotate data.raw and prepare for further normalization
14. Annotate Genic Locations
15. Annotate to CpG Islands
16. Chronic methylation change
17. Chromosomes with the highest DMPs frequency
18. Assign Hi-C compartments to the data
19. Data preparation for Hi-C analysis
20. Order and plot data.raw chromosome position
21. Chromosomal frequencies in each nucleus layer
22. DMPs frequency and DNA methylation change within nucleus layers
23. Load and prepare histone modifications datasets for methylation analysis
24. Assign histone modifications to Hi-C layers
25. Plot data.raw probes into histone modifications picks and Hi-C layers
26. Data preparation for regulatory regions (promoter and enhancers) analysis
27. Add Hi-C to regulatory regions
28. DNA Methylation change in regulatory regions
29. Baseline methylation level in regulatory regions
30. Add Hi-C to regulatory regions (without DMPs overlapping)
