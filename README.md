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
3. DMPs.annotations.all.probes.2000.rdata (you can generate it via in second step of DNA.methylation.R script - it takes several hours and 200MB of disk space)
4. DMRs.newest.GRanges.rdata
5. DMRs.newest.annotations.rdata
6. Bed.newest.GRanges.rdata

The rest of files needed for certain step are avaiable on publicly available respositories and are decribied in detail in the article. 

## Table of contents:

### I. DNA.methylation.R (Running time when all needed files are provided ~ 20 minutes)
1. Install and load packages
2. Import and prepare primary files
3. Differential methylation analysis
4. Pick significant DMPs
5. DMRs Building
6. Simulate DMRs building with rows permutation (Optional step - can take several hours to run based on number of CPUs and repetitions)
7. Load and prepare lists of DMRs
8. Acute methylation change versus control between DMPs **(Figure_1D; Figure_S1B; Figure_S1C; Figure_1D)**
9. % of DMPs in DMRs (Figure_S1H)
10. DMPs which changed at least of 0.58 or -0.58 M Value versus control (Figure_1E; Figure_1D; Figure_S1D)
11. Detailed methylation patterns **(Figure_S1E; Figure_S1F)**
12. Select DMRs, DMPs in or out of DMRs and all DMPs for further analysis
13. Annotate data.raw and prepare for further normalization
14. Annotate Genic Locations **(Figure_1G)**
15. Annotate to CpG Islands **(Figure_1H; Figure_1I; Figure_S1J; Figure_S1K; Figure_S1L)**
16. Chronic methylation change **(Figure_4A; Figure_4B; Figure_S4A)**
17. Chromosomes with the highest DMPs frequency **(Figure_S2B; Figure_S2C; Figure_S2D; Figure_S2E; Figure_S2F; Figure_S2G)**
18. Assign Hi-C compartments to the data
19. Data preparation for Hi-C analysis
20. Order and plot data.raw chromosome position **(Figure_2B; Figure_S2A)**
21. Chromosomal frequencies in each nucleus layer **(Figure_2D; Figure_S2I; Figure_S2J; Figure_S2K; Figure_S2L; Figure_S2M)**
22. DMPs frequency and DNA methylation change within nucleus layers **(Figure_2C; Figure_2E; Figure_S2H; Figure_S2N)**
23. Load and prepare histone modifications datasets for methylation analysis
24. Assign histone modifications to Hi-C layers **(Figure_3B; Figure_3C; Figure_3D; Figure_3E; Figure_3F; Figure_S3B; Figure_S3C; Figure_S3D)**
25. Plot data.raw probes into histone modifications picks and Hi-C layers **(Figure_3A; Figure_S4A)**
26. Data preparation for regulatory regions (promoter and enhancers) analysis
27. Add Hi-C to regulatory regions **(Figure_S3F)**
28. DNA Methylation change in regulatory regions **(Figure_3G; Figure_S3I; Figure_S3J; Figure_S3K; Figure_S3L; Figure_S3G)**
29. Baseline methylation level in regulatory regions
30. Add Hi-C to regulatory regions (without DMPs overlapping) **(Figure_S3E)**

We suggest running it chronologically as some steps may have consequences in further ones.

### II. Genes.for.hic.astronatus.R (Running time ~ 3 minutes)
This step creates space.genes.rdata file necessary for RNAseq.R script

### III. Microarray.R (Running time ~ 2 minutes)
1. Load packages
2. Microarray for Fe - Heart **(Figure_4F; Figure_4G)**
3. Microarray for Si - Breast **(Figure_S3A; Figure_S3B)**
4. Microarray for Fe - Heart _ Protons **(Figure_S3C)**
5. Microarray for Si - Breast _ Gamma Rays **(Figure_S5E; Figure_S5F)**
6. Microarray for Cs137 - Blood **(Figure_S5G; Figure_S5H)**

### IV. RNAseq.R (Running time ~ minutes)
1. Load packages
2. RNA-seq - Fe liver **(Figure_4D; Figure_4E)**
3. RNA-seq - Co57 retina **(Figure_S5I; Figure_S5J)**
4. RNA-seq - Astronauts - blood **(Figure_4G; Figure_4H)**
5. Astronauts - Hi-C **(Figure_4F; Figure_S4M)**
