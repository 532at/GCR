# 1. Load packages ####
library(data.table)
library(reshape2)
library(ggplot2)
library(dplyr)
library(tibble)
library(readxl)
library(RColorBrewer)
library(rlist)
library(datasets)
library(tidyr)
library(robustbase)
library(Kendall)
library(org.Hs.eg.db)

# 2. RNA-seq - Fe liver ####

getwd()
wpath <- "/mnt/work1/adrian/Nasa/"
setwd(wpath)
options(width=200)
set.seed(123)

counts.raw <- data.frame(data.table::fread("GLDS-294_rna-seq_GSE136165_Expression_Data.csv.xls", header=TRUE, sep=",", check.names = F, stringsAsFactors = F))

rownames(counts.raw) <- counts.raw[,1]

counts.raw <- counts.raw[,-1]

samples_tab <- data.frame("sample"=colnames(counts.raw))

## create a table with sample information / phenotypic information / any other information on your samples which might be relevant for samples grouping,
## variability of the expression data etc. Most importantly this table should contain the biological coeficient you want to test, e.g. CTRL vs. PT.
samples_tab <- data.frame("sample"=colnames(counts.raw)) 
samples <- colnames(counts.raw)

samples_tab$type <- c(rep("con",15), rep("rad",15))
samples_tab$time <- rep(c(rep("1",3), rep("2",3), rep("4",3), rep("9",3), rep("12",3)),2)
samples_tab[,4] <- unite(samples_tab[,2:3], ID, sep = "_", remove = T)

##### FILTER DATA
library(edgeR)
## How many genes (rows) sum up to 0 in terms of their raw expression? These are presentend in a form of TRUE/FALSE values
table(rowSums(counts.raw)==0)

## Here we check how stringent can we be when filtering the data. We will use the miniumum number of samples and minimum expression (count per million)
## We iterate through numbers from 0 to 10 to check how this influences the final number of genes in our data:
stats.minexpr <- sapply(c(0:10), function(x) sum(rowSums(cpm(counts.raw[rowSums(counts.raw)!=0,])>=1)>=x))
names(stats.minexpr) <- paste0(">=",c(0:10))
print(stats.minexpr)

## These are our final thresholds - at least 5 samples where expression is at least 1 CPM
minsamples <- 5; mincpm <- 1
## we store this information in a vector of TRUE/FALSE values for every gene in the raw expression table
keep <- rowSums(cpm(counts.raw)>=mincpm)>=minsamples

counts.filt <- counts.raw[keep,] # Reduce the counts to only the once fullfilling our ctiteria (so genes with TRUE label). We don not filter by column - sample

##### VISUALIZE DATA
table(samples_tab$ID)
## define a variable that will hold colors we want to use for plots. Here we can define colors as numbers, or define them with literal color names
my_cols <- c(1:10)
## alternatively:
#my_cols <- c("green","cyan","blue","red","purple")
names(my_cols) <- c("con_1","con_2","con_4","con_9", "con_12",
                    "rad_1","rad_2","rad_4","rad_9", "rad_12") # assign names to your colors. Here names are sample types
print(my_cols)

## for plotting, it is best to use log values. Because when transforming to log we cannot have values equal to zero, hence we need to add a so-called pseudo count
## Here wa use a pseudo count equall to 1
counts.filt.pseudo <- log2(counts.filt+1)

## plot distribution of counts per sample
par(mfrow=c(1,1))
boxplot(counts.filt.pseudo, las=2, cex=.5, cex.axis=.5, main="Raw counts, filtered, log") # we make pseudocounts of the raw data

## Visualize our data as MDS plot - multidimensional scaling plot of distances between gene expression profiles. Remember, MDS is not the same as PCA plot
plotMDS(counts.filt.pseudo, las=2, cex=1, col=my_cols[samples_tab$ID])


##### EDGER OBJECT 

## Here we set up an EdgeR object that contains expression data and sample description. We group samples per sample type (CTRL, PT, etc.)
dgeFull <- DGEList(counts.filt, group=samples_tab$ID)
dgeFull$samples$Individual <- samples_tab$sample
dgeFull$samples

## Normalize values with TMM - weighted trimmed mean of the log expression ratios (trimmed mean of M values)
dgeFull <- calcNormFactors(dgeFull, method="TMM")
dgeFull$samples # check how the column called 'norm.factors' has changed as compared to the previous execution of the same line
## Here we want to define the order of our sample groups - otherwise they are ordered alphabetically
dgeFull$samples$group <- factor(dgeFull$samples$group, levels=c("con_1", "con_2", "con_4", "con_9", "con_12","rad_1","rad_2","rad_4","rad_9", "rad_12"))


##### NORMALIZED COUNTS

## get a table of normalized counts and a table of log-transformed normalized counts
counts.norm <- cpm(dgeFull)
counts.norm.pseudo <- log2(counts.norm + 1)


## Plot, side by side, boxplots of raw and normalized per sample expression counts
par(mfrow=c(1,2), mar=c(4,2,2,1))
boxplot(counts.filt.pseudo, las=2, cex=.5, cex.axis=.5, main="Raw counts, filtered, log")
boxplot(counts.norm.pseudo, las=2, cex=.5, cex.axis=.5, main="Normalized counts, filtered, log")

## Plot, side by side, MDS plots of raw and normalized data
par(mfrow=c(1,2), mar=c(4,3,2,1), cex.axis=.8, cex.lab=.8)
plotMDS(counts.norm.pseudo, cex=1, col=my_cols[samples_tab$ID], main="Raw counts, filtered, log")
plotMDS(counts.filt.pseudo, cex=1, col=my_cols[samples_tab$ID], main="Normalized counts, filtered, log")


## PCA
library(FactoMineR)
## Here we run PCA analysis to yet again visualize our data, but also to check how it groups based on the biological type. We can also identify here any potential problems, outliers etc.
## Remember - for PCA we need a transposed table of normalized counts - genes as columns and samples as rows
norm_t <- t(counts.norm)
#norm_t[1:5,1:5]
my.pca <- PCA(log(norm_t+1), graph=FALSE)

library(factoextra)
# PCA plot using first two axes. We group/color samples based on the sample type. We use a previously defined set of colors.
pca.plot1 <- fviz_pca_ind(my.pca, habillage=dgeFull$samples$group, invisible="quali", repel=2, labelsize=3, pointsize=3, title="", palette=my_cols) +
  theme(legend.position="bottom")
# PCA plot using first axis 1 and axis 3
pca.plot2 <- fviz_pca_ind(my.pca, habillage=dgeFull$samples$group, invisible="quali", repel=2, labelsize=3, pointsize=3, title="", palette=my_cols, axes = c(1,3)) +
  theme(legend.position="bottom")

print(pca.plot1)
print(pca.plot2)


##### DIFFERENTIAL TEST
## To do the comparisons (DE tests) we need to define a design matrix for our data, as well as a set of contrasts (comparisons) to be tested
designMat <- model.matrix(~0 + dgeFull$samples$group) # here, we can add any other important information to our model, such as gender, age, sequencing batch, etc.
print(designMat)
## modify column names and row names of the design matrix so it is easier to read
colnames(designMat) <- c("con_1", "con_2", "con_4", "con_9", "con_12","rad_1","rad_2","rad_4","rad_9", "rad_12")
rownames(designMat) <- rownames(dgeFull$samples)
print(designMat)


contMatrix <- makeContrasts(con_1-rad_1, con_1-rad_2, con_1-rad_4, con_1-rad_9, con_1-rad_12,
                            rad_1-rad_2, rad_1-rad_4, rad_1-rad_9, rad_1-rad_12,
                            levels=designMat)
print(contMatrix)


## Now we need to estimate the dispersion parameter for our negative binomial model
dgeFull <- estimateDisp(dgeFull, designMat)
dgeFull <- estimateGLMCommonDisp(dgeFull, designMat)
dgeFull <- estimateGLMTrendedDisp(dgeFull, designMat)
dgeFull <- estimateGLMTagwiseDisp(dgeFull, designMat)
## Fit a negative binomial generalized log-linear model to the read counts for each gene.
## Conduct genewise statistical tests for a given coefficient or coefficient contrast.
fit <- glmFit(dgeFull, designMat)

## Fit a quasi-likelihood negative binomial generalized log-linear model to count data.
## Conduct genewise statistical tests for a given coefficient or contrast.
fitQLF <- glmQLFit(dgeFull, designMat)

## How can we make this for all contrasts in one go? Here we use a custom function and we loop throug all contrasts with the list-apply (lapply) function:

## Loop contrasts and apply LRT
LRTtests <- lapply(colnames(contMatrix), function(x) {
  glmLRT(fit, contrast=contMatrix[,x]) %>% topTags(n=nrow(dgeFull))
})

## Loop contrasts and apply QLF
QLFtests <- lapply(colnames(contMatrix), function(x) {
  glmQLFTest(fitQLF, contrast=contMatrix[,x]) %>% topTags(n=nrow(dgeFull))
})

## assign names to the list with results of both tests. The names are simply names of contrasts (column names of contMatrix table)
names(LRTtests) <- names(QLFtests) <- colnames(contMatrix)

LRTtests.sel <- lapply(LRTtests, function(x) subset(x$table, FDR<0.05 & abs(logFC) >= 0))

QLFtests.sel <- lapply(QLFtests, function(x) subset(x$table, FDR<0.05 & abs(logFC) >= 0))
#& abs(logFC) >= 1
## let's additionally sort our DEGs based on log fold change of expression:
LRTtests.sel <- lapply(LRTtests.sel, function(x) x[order(x$logFC),])
QLFtests.sel <- lapply(QLFtests.sel, function(x) x[order(x$logFC),])

## Here we print a small summary table of how many DEGs were found per contrast and per type of test
rbind("LRT"=sapply(LRTtests.sel, nrow), "QLF"=sapply(QLFtests.sel, nrow))

########### Further analysis on selected genes (FDR, logFC)

##### CONVERT 
library(org.Mm.eg.db)

LRTtests.sel.genes <- lapply(LRTtests.sel, function(x) {x$Gene <- toupper(mapIds(
  org.Mm.eg.db,
  keys = rownames(x),
  column = 'SYMBOL',
  keytype = 'ENSEMBL'));x})

##### CONVERT without FDR threshold

LRTtests.noFDR <- lapply(lapply(lapply(LRTtests, function(x) x$table), function(x) x[order(x$logFC),]), function(x) {x$Gene <- toupper(mapIds(
  org.Mm.eg.db,
  keys = rownames(x),
  column = 'SYMBOL',
  keytype = 'ENSEMBL'));x})


Fe.liver <- filter(LRTtests.sel.genes$`con_1 - rad_1`, !Gene == "NA")[,c(1,6)]
colnames(Fe.liver)[1] <- "logFC_1"

Fe.liver.2 <- left_join(Fe.liver, LRTtests.noFDR$`con_1 - rad_2`[,c(1,6)], by="Gene")[,c(1,3,2)]
colnames(Fe.liver.2)[2] <- "logFC_2"

Fe.liver.4 <- left_join(Fe.liver.2, LRTtests.noFDR$`con_1 - rad_4`[,c(1,6)], by="Gene")[,c(1,2,4,3)]
colnames(Fe.liver.4)[3] <- "logFC_4"

Fe.liver.9 <- left_join(Fe.liver.4, LRTtests.noFDR$`con_1 - rad_9`[,c(1,6)], by="Gene")[,c(1,2,3,5,4)]
colnames(Fe.liver.9)[4] <- "logFC_9"

Fe.liver.12 <- left_join(Fe.liver.9, LRTtests.noFDR$`con_1 - rad_12`[,c(1,6)], by="Gene")[,c(1,2,3,4,6,5)]
colnames(Fe.liver.12)[5] <- "logFC_12"

rownames(Fe.liver.12) <- Fe.liver.12$Gene

Fe.heatmap <- Fe.liver.12[,1:5]

colnames(Fe.heatmap) <- c("1", "2", "4", "9", "12")

#Figure_4D
ggplot(melt(filter(Fe.heatmap, Fe.heatmap[,1] > 0)), aes(factor(variable), value)) +
  geom_boxplot(fill = "lightblue", alpha = 0.75, width=0.5) +
  geom_hline(yintercept=0, linetype="dashed", color = "blue", size=1, alpha = 0.5) +
  theme_classic() +
  ylim(-6,6) +
  labs(title="Up-regulated genes in liver\ntissue (Fe 0.2 Gy)", x="Time after exposure (months)", y="logFC", color="Sample") +
  theme(
    axis.title.x = element_text(size = 8, face="bold"),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    axis.title.y = element_text(size = 8, face="bold"),
    title = element_text(size = 8, face = "bold"))

#Figure_4E
ggplot(melt(filter(Fe.heatmap, Fe.heatmap[,1] < 0)), aes(factor(variable), value)) +
  geom_boxplot(fill = "lightblue", alpha = 0.75, width=0.5) +
  geom_hline(yintercept=0, linetype="dashed", color = "blue", size=1, alpha = 0.5) +
  theme_classic() +
  ylim(-6,6) +
  labs(title="Down-regulated genes in liver\ntissue (Fe 0.2 Gy)", x="Time after exposure (months)", y="logFC", color="Sample") +
  theme(
    axis.title.x = element_text(size = 8, face="bold"),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    axis.title.y = element_text(size = 8, face="bold"),
    title = element_text(size = 8, face = "bold"))

#Trend
anova(lm(colMedians(as.matrix(filter(Fe.heatmap, Fe.heatmap[,1] < 0)))~as.vector(1:5)))[[5]][1]
anova(lm(colMedians(as.matrix(filter(Fe.heatmap, Fe.heatmap[,1] > 0)))~as.vector(1:5)))[[5]][1]


# 3. RNA-seq - Co57 retina ####

getwd()
wpath <- "/mnt/work1/adrian/Nasa/"
setwd(wpath)
options(width=200)
set.seed(123)

counts.raw.info <- data.frame(data.table::fread("GLDS-203-samples.csv", header=TRUE, sep=",", check.names = F, stringsAsFactors = F))

counts.raw.info <- counts.raw.info[,2:3]
counts.raw.info$Sample.Name <- gsub("C57-6J", "C57.6J", counts.raw.info$Sample.Name)

counts.raw <- data.frame(data.table::fread("GLDS-203_rna_seq_Normalized_Counts.csv", header=TRUE, sep=",", check.names = F, stringsAsFactors = F))

col.names <- left_join(data.frame("Sample.Name"=colnames(counts.raw[-1])), counts.raw.info, by="Sample.Name")

rownames(counts.raw) <- counts.raw[,1]

counts.raw <- counts.raw[,-1]

colnames(counts.raw) <- col.names$Comment..Original.Submitted.Sample.Name

colnames(counts.raw)

counts.raw <- counts.raw[,c(1:16,26:31)]

samples_tab <- data.frame("sample"=colnames(counts.raw))

## create a table with sample information / phenotypic information / any other information on your samples which might be relevant for samples grouping,
## variability of the expression data etc. Most importantly this table should contain the biological coeficient you want to test, e.g. CTRL vs. PT.
samples_tab <- data.frame("sample"=colnames(counts.raw)) 
samples <- colnames(counts.raw)

samples_tab$type <- c(rep("rad",16), rep("con",6))
samples_tab$time <- c(rep("m1",6), rep("m4",5), rep("d7",11))
samples_tab[,4] <- unite(samples_tab[,2:3], ID, sep = "_", remove = T)

##### FILTER DATA

## How many genes (rows) sum up to 0 in terms of their raw expression? These are presentend in a form of TRUE/FALSE values
table(rowSums(counts.raw)==0)

## Here we check how stringent can we be when filtering the data. We will use the miniumum number of samples and minimum expression (count per million)
## We iterate through numbers from 0 to 10 to check how this influences the final number of genes in our data:
stats.minexpr <- sapply(c(0:10), function(x) sum(rowSums(cpm(counts.raw[rowSums(counts.raw)!=0,])>=1)>=x))
names(stats.minexpr) <- paste0(">=",c(0:10))
print(stats.minexpr)

## These are our final thresholds - at least 5 samples where expression is at least 1 CPM
minsamples <- 5; mincpm <- 1
## we store this information in a vector of TRUE/FALSE values for every gene in the raw expression table
keep <- rowSums(cpm(counts.raw)>=mincpm)>=minsamples

counts.filt <- counts.raw[keep,] # Reduce the counts to only the once fullfilling our ctiteria (so genes with TRUE label). We don not filter by column - sample
dim(counts.filt)


##### VISUALIZE DATA
table(samples_tab$ID)
## define a variable that will hold colors we want to use for plots. Here we can define colors as numbers, or define them with literal color names
my_cols <- c(1:4)
## alternatively:
#my_cols <- c("green","cyan","blue","red","purple")
names(my_cols) <- c("con_d7","rad_d7","rad_m1","rad_m4") # assign names to your colors. Here names are sample types
print(my_cols)

## for plotting, it is best to use log values. Because when transforming to log we cannot have values equal to zero, hence we need to add a so-called pseudo count
## Here wa use a pseudo count equall to 1
counts.filt.pseudo <- log2(counts.filt+1)

## plot distribution of counts per sample
par(mfrow=c(1,1))
boxplot(counts.filt.pseudo, las=2, cex=.5, cex.axis=.5, main="Raw counts, filtered, log") # we make pseudocounts of the raw data

## Visualize our data as MDS plot - multidimensional scaling plot of distances between gene expression profiles. Remember, MDS is not the same as PCA plot
plotMDS(counts.filt.pseudo, las=2, cex=1, col=my_cols[samples_tab$ID])


##### EDGER OBJECT 

## Here we set up an EdgeR object that contains expression data and sample description. We group samples per sample type (CTRL, PT, etc.)
dgeFull <- DGEList(counts.filt, group=samples_tab$ID)
dgeFull$samples$Individual <- samples_tab$sample
dgeFull$samples

## Normalize values with TMM - weighted trimmed mean of the log expression ratios (trimmed mean of M values)
dgeFull <- calcNormFactors(dgeFull, method="TMM")
dgeFull$samples # check how the column called 'norm.factors' has changed as compared to the previous execution of the same line
## Here we want to define the order of our sample groups - otherwise they are ordered alphabetically
dgeFull$samples$group <- factor(dgeFull$samples$group, levels=c("con_d7","rad_d7","rad_m1","rad_m4"))


##### NORMALIZED COUNTS

## get a table of normalized counts and a table of log-transformed normalized counts
counts.norm <- cpm(dgeFull)
counts.norm.pseudo <- log2(counts.norm + 1)


## Plot, side by side, boxplots of raw and normalized per sample expression counts
par(mfrow=c(1,2), mar=c(4,2,2,1))
boxplot(counts.filt.pseudo, las=2, cex=.5, cex.axis=.5, main="Raw counts, filtered, log")
boxplot(counts.norm.pseudo, las=2, cex=.5, cex.axis=.5, main="Normalized counts, filtered, log")

## Plot, side by side, MDS plots of raw and normalized data
par(mfrow=c(1,2), mar=c(4,3,2,1), cex.axis=.8, cex.lab=.8)
plotMDS(counts.norm.pseudo, cex=1, col=my_cols[samples_tab$ID], main="Raw counts, filtered, log")
plotMDS(counts.filt.pseudo, cex=1, col=my_cols[samples_tab$ID], main="Normalized counts, filtered, log")


## PCA

## Here we run PCA analysis to yet again visualize our data, but also to check how it groups based on the biological type. We can also identify here any potential problems, outliers etc.
## Remember - for PCA we need a transposed table of normalized counts - genes as columns and samples as rows
norm_t <- t(counts.norm)
#norm_t[1:5,1:5]
my.pca <- PCA(log(norm_t+1), graph=FALSE)

# PCA plot using first two axes. We group/color samples based on the sample type. We use a previously defined set of colors.
pca.plot1 <- fviz_pca_ind(my.pca, habillage=dgeFull$samples$group, invisible="quali", repel=2, labelsize=3, pointsize=3, title="", palette=my_cols) +
  theme(legend.position="bottom")
# PCA plot using first axis 1 and axis 3
pca.plot2 <- fviz_pca_ind(my.pca, habillage=dgeFull$samples$group, invisible="quali", repel=2, labelsize=3, pointsize=3, title="", palette=my_cols, axes = c(1,3)) +
  theme(legend.position="bottom")

print(pca.plot1)
print(pca.plot2)


##### DIFFERENTIAL TEST
## To do the comparisons (DE tests) we need to define a design matrix for our data, as well as a set of contrasts (comparisons) to be tested
designMat <- model.matrix(~0 + dgeFull$samples$group) # here, we can add any other important information to our model, such as gender, age, sequencing batch, etc.
print(designMat)
## modify column names and row names of the design matrix so it is easier to read
colnames(designMat) <- c("con_d7","rad_d7","rad_m1","rad_m4")
rownames(designMat) <- rownames(dgeFull$samples)
print(designMat)


contMatrix <- makeContrasts(con_d7-rad_d7, con_d7-rad_m1, con_d7-rad_m4,
                            rad_d7-rad_m1, rad_d7-rad_m4,
                            levels=designMat)
print(contMatrix)


## Now we need to estimate the dispersion parameter for our negative binomial model
dgeFull <- estimateDisp(dgeFull, designMat)
dgeFull <- estimateGLMCommonDisp(dgeFull, designMat)
dgeFull <- estimateGLMTrendedDisp(dgeFull, designMat)
dgeFull <- estimateGLMTagwiseDisp(dgeFull, designMat)
## Fit a negative binomial generalized log-linear model to the read counts for each gene.
## Conduct genewise statistical tests for a given coefficient or coefficient contrast.
fit <- glmFit(dgeFull, designMat)

## Fit a quasi-likelihood negative binomial generalized log-linear model to count data.
## Conduct genewise statistical tests for a given coefficient or contrast.
fitQLF <- glmQLFit(dgeFull, designMat)

## How can we make this for all contrasts in one go? Here we use a custom function and we loop throug all contrasts with the list-apply (lapply) function:

## Loop contrasts and apply LRT
LRTtests <- lapply(colnames(contMatrix), function(x) {
  glmLRT(fit, contrast=contMatrix[,x]) %>% topTags(n=nrow(dgeFull))
})

## Loop contrasts and apply QLF
QLFtests <- lapply(colnames(contMatrix), function(x) {
  glmQLFTest(fitQLF, contrast=contMatrix[,x]) %>% topTags(n=nrow(dgeFull))
})

## assign names to the list with results of both tests. The names are simply names of contrasts (column names of contMatrix table)
names(LRTtests) <- names(QLFtests) <- colnames(contMatrix)

LRTtests.sel <- lapply(LRTtests, function(x) subset(x$table, FDR<0.05 & abs(logFC) >= 0))

QLFtests.sel <- lapply(QLFtests, function(x) subset(x$table, FDR<0.05 & abs(logFC) >= 0))
#& abs(logFC) >= 1
## let's additionally sort our DEGs based on log fold change of expression:
LRTtests.sel <- lapply(LRTtests.sel, function(x) x[order(x$logFC),])
QLFtests.sel <- lapply(QLFtests.sel, function(x) x[order(x$logFC),])

## Here we print a small summary table of how many DEGs were found per contrast and per type of test
rbind("LRT"=sapply(LRTtests.sel, nrow), "QLF"=sapply(QLFtests.sel, nrow))

########### Further analysis on selected genes (FDR, logFC)

##### CONVERT 
library(org.Mm.eg.db)

LRTtests.sel.genes <- lapply(LRTtests.sel, function(x) {x$Gene <- toupper(mapIds(
  org.Mm.eg.db,
  keys = rownames(x),
  column = 'SYMBOL',
  keytype = 'ENSEMBL'));x})

##### CONVERT without FDR threshold

LRTtests.noFDR <- lapply(lapply(lapply(LRTtests, function(x) x$table), function(x) x[order(x$logFC),]), function(x) {x$Gene <- toupper(mapIds(
  org.Mm.eg.db,
  keys = rownames(x),
  column = 'SYMBOL',
  keytype = 'ENSEMBL'));x})


Co.retina <- filter(LRTtests.sel.genes$`con_d7 - rad_d7`, !Gene == "NA")[,c(1,6)]
colnames(Co.retina)[1] <- "logFC_d7"

Co.retina.2 <- left_join(Co.retina, LRTtests.noFDR$`con_d7 - rad_m1`[,c(1,6)], by="Gene")[,c(1,3,2)]
colnames(Co.retina.2)[2] <- "logFC_m1"

Co.retina.4 <- left_join(Co.retina.2, LRTtests.noFDR$`con_d7 - rad_m4`[,c(1,6)], by="Gene")[,c(1,2,4,3)]
colnames(Co.retina.4)[3] <- "logFC_m4"

rownames(Co.retina.4) <- Co.retina.4$Gene

Co.heatmap <- Co.retina.4[,-4]


colnames(Co.heatmap) <- c("7 days", "1 month", "4 months")

#Figure_S5I
ggplot(melt(filter(Co.heatmap, Co.heatmap[,1] > 0)), aes(factor(variable), value)) +
  geom_boxplot(fill = "lightblue", alpha = 0.75, width=0.5) +
  geom_hline(yintercept=0, linetype="dashed", color = "blue", size=1, alpha = 0.5) +
  theme_classic() +
  labs(title="Up-regulated genes in retina\ntissue (Co 0.04 Gy)", x="Time after exposure", y="logFC", color="Sample") +
  theme(
    axis.title.x = element_text(size = 8, face="bold"),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    axis.title.y = element_text(size = 8, face="bold"),
    title = element_text(size = 8, face = "bold"))

#Figure_S5J
ggplot(melt(filter(Co.heatmap, Co.heatmap[,1] < 0)), aes(factor(variable), value)) +
  geom_boxplot(fill = "lightblue", alpha = 0.75, width=0.5) +
  geom_hline(yintercept=0, linetype="dashed", color = "blue", size=1, alpha = 0.5) +
  theme_classic() +
  labs(title="Down-regulated genes in retina\ntissue (Co 0.04 Gy)", x="Time after exposure", y="logFC", color="Sample") +
  theme(
    axis.title.x = element_text(size = 8, face="bold"),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    axis.title.y = element_text(size = 8, face="bold"),
    title = element_text(size = 8, face = "bold"))

#Trend
anova(lm(colMedians(as.matrix(filter(Co.heatmap, Co.heatmap[,1] < 0)))~as.vector(1:3)))[[5]][1]
anova(lm(colMedians(as.matrix(filter(Co.heatmap, Co.heatmap[,1] > 0)))~as.vector(1:3)))[[5]][1]



# 4. RNA-seq - Astronauts - blood ####

getwd()
wpath <- "/mnt/work1/adrian/Nasa/"
setwd(wpath)
options(width=200)
set.seed(123)

counts.raw <- data.frame(read_xlsx("GLDS-530_rna-seq_TGB_050_1_2_64samples_11group_totalcount_all0removed_scalingnormalized_SEM.xlsx"))

rownames(counts.raw) <- counts.raw[,1]

counts.raw <- counts.raw[,-1]

counts.raw <- counts.raw[,c(3,5,7,9,11,13,15,17,19,21)]
colnames(counts.raw) <- gsub("*...Normalized.means", "", colnames(counts.raw))

## create a table with sample information / phenotypic information / any other information on your samples which might be relevant for samples grouping,
## variability of the expression data etc. Most importantly this table should contain the biological coeficient you want to test, e.g. CTRL vs. PT.
samples_tab <- data.frame("sample"=colnames(counts.raw)) 
samples <- colnames(counts.raw)

samples_tab$type <- c(rep("con",2), rep("rad",8))
samples_tab$time <- c(1,1,2,2,3,3,4,4,5,5)
samples_tab[,4] <- unite(samples_tab[,2:3], ID, sep = "_", remove = T)

##### FILTER DATA
library(edgeR)
## How many genes (rows) sum up to 0 in terms of their raw expression? These are presentend in a form of TRUE/FALSE values
table(rowSums(counts.raw)==0)

## Here we check how stringent can we be when filtering the data. We will use the miniumum number of samples and minimum expression (count per million)
## We iterate through numbers from 0 to 10 to check how this influences the final number of genes in our data:
stats.minexpr <- sapply(c(0:10), function(x) sum(rowSums(cpm(counts.raw[rowSums(counts.raw)!=0,])>=1)>=x))
names(stats.minexpr) <- paste0(">=",c(0:10))
print(stats.minexpr)

## These are our final thresholds - at least 5 samples where expression is at least 1 CPM
minsamples <- 5; mincpm <- 1
## we store this information in a vector of TRUE/FALSE values for every gene in the raw expression table
keep <- rowSums(cpm(counts.raw)>=mincpm)>=minsamples

counts.filt <- counts.raw[keep,] # Reduce the counts to only the once fullfilling our ctiteria (so genes with TRUE label). We don not filter by column - sample
dim(counts.filt)


##### VISUALIZE DATA
table(samples_tab$ID)
## define a variable that will hold colors we want to use for plots. Here we can define colors as numbers, or define them with literal color names
my_cols <- c(1:5)
## alternatively:
#my_cols <- c("green","cyan","blue","red","purple")
names(my_cols) <- c("con_1","rad_2","rad_3","rad_4", "rad_5") # assign names to your colors. Here names are sample types
print(my_cols)

## for plotting, it is best to use log values. Because when transforming to log we cannot have values equal to zero, hence we need to add a so-called pseudo count
## Here wa use a pseudo count equall to 1
counts.filt.pseudo <- log2(counts.filt+1)

## plot distribution of counts per sample
par(mfrow=c(1,1))
boxplot(counts.filt.pseudo, las=2, cex=.5, cex.axis=.5, main="Raw counts, filtered, log") # we make pseudocounts of the raw data

## Visualize our data as MDS plot - multidimensional scaling plot of distances between gene expression profiles. Remember, MDS is not the same as PCA plot
plotMDS(counts.filt, las=2, cex=1, col=my_cols[samples_tab$ID])


##### EDGER OBJECT 

## Here we set up an EdgeR object that contains expression data and sample description. We group samples per sample type (CTRL, PT, etc.)
dgeFull <- DGEList(counts.filt, group=samples_tab$ID)
dgeFull$samples$Individual <- samples_tab$sample
dgeFull$samples

## Normalize values with TMM - weighted trimmed mean of the log expression ratios (trimmed mean of M values)
dgeFull <- calcNormFactors(dgeFull, method="TMM")
dgeFull$samples # check how the column called 'norm.factors' has changed as compared to the previous execution of the same line
## Here we want to define the order of our sample groups - otherwise they are ordered alphabetically
dgeFull$samples$group <- factor(dgeFull$samples$group, levels=c("con_1","rad_2","rad_3","rad_4", "rad_5"))


##### NORMALIZED COUNTS

## get a table of normalized counts and a table of log-transformed normalized counts
counts.norm <- cpm(dgeFull)
counts.norm.pseudo <- log2(counts.norm + 1)


## Plot, side by side, boxplots of raw and normalized per sample expression counts
par(mfrow=c(1,2), mar=c(4,2,2,1))
boxplot(counts.filt.pseudo, las=2, cex=.5, cex.axis=.5, main="Raw counts, filtered, log")
boxplot(counts.norm.pseudo, las=2, cex=.5, cex.axis=.5, main="Normalized counts, filtered, log")

## Plot, side by side, MDS plots of raw and normalized data
par(mfrow=c(1,2), mar=c(4,3,2,1), cex.axis=.8, cex.lab=.8)
plotMDS(counts.norm.pseudo, cex=1, col=my_cols[samples_tab$ID], main="Raw counts, filtered, log")
plotMDS(counts.filt.pseudo, cex=1, col=my_cols[samples_tab$ID], main="Normalized counts, filtered, log")


## PCA
library(FactoMineR)
## Here we run PCA analysis to yet again visualize our data, but also to check how it groups based on the biological type. We can also identify here any potential problems, outliers etc.
## Remember - for PCA we need a transposed table of normalized counts - genes as columns and samples as rows
norm_t <- t(counts.norm)
#norm_t[1:5,1:5]
my.pca <- PCA(log(norm_t+1), graph=FALSE)

library(factoextra)
# PCA plot using first two axes. We group/color samples based on the sample type. We use a previously defined set of colors.
pca.plot1 <- fviz_pca_ind(my.pca, habillage=dgeFull$samples$group, invisible="quali", repel=2, labelsize=3, pointsize=3, title="", palette=my_cols) +
  theme(legend.position="bottom")
# PCA plot using first axis 1 and axis 3
pca.plot2 <- fviz_pca_ind(my.pca, habillage=dgeFull$samples$group, invisible="quali", repel=2, labelsize=3, pointsize=3, title="", palette=my_cols, axes = c(1,3)) +
  theme(legend.position="bottom")

print(pca.plot1)
print(pca.plot2)


##### DIFFERENTIAL TEST
## To do the comparisons (DE tests) we need to define a design matrix for our data, as well as a set of contrasts (comparisons) to be tested
designMat <- model.matrix(~0 + dgeFull$samples$group) # here, we can add any other important information to our model, such as gender, age, sequencing batch, etc.
print(designMat)
## modify column names and row names of the design matrix so it is easier to read
colnames(designMat) <- c("con_1","rad_2","rad_3","rad_4", "rad_5")
rownames(designMat) <- rownames(dgeFull$samples)
print(designMat)


contMatrix <- makeContrasts(con_1-rad_2, con_1-rad_3, con_1-rad_4, con_1-rad_5,
                            levels=designMat)

print(contMatrix)


## Now we need to estimate the dispersion parameter for our negative binomial model
dgeFull <- estimateDisp(dgeFull, designMat)
dgeFull <- estimateGLMCommonDisp(dgeFull, designMat)
dgeFull <- estimateGLMTrendedDisp(dgeFull, designMat)
dgeFull <- estimateGLMTagwiseDisp(dgeFull, designMat)
## Fit a negative binomial generalized log-linear model to the read counts for each gene.
## Conduct genewise statistical tests for a given coefficient or coefficient contrast.
fit <- glmFit(dgeFull, designMat)

## Fit a quasi-likelihood negative binomial generalized log-linear model to count data.
## Conduct genewise statistical tests for a given coefficient or contrast.
fitQLF <- glmQLFit(dgeFull, designMat)

## How can we make this for all contrasts in one go? Here we use a custom function and we loop throug all contrasts with the list-apply (lapply) function:

## Loop contrasts and apply LRT
LRTtests <- lapply(colnames(contMatrix), function(x) {
  glmLRT(fit, contrast=contMatrix[,x]) %>% topTags(n=nrow(dgeFull))
})

## Loop contrasts and apply QLF
QLFtests <- lapply(colnames(contMatrix), function(x) {
  glmQLFTest(fitQLF, contrast=contMatrix[,x]) %>% topTags(n=nrow(dgeFull))
})

## assign names to the list with results of both tests. The names are simply names of contrasts (column names of contMatrix table)
names(LRTtests) <- names(QLFtests) <- colnames(contMatrix)

LRTtests.sel <- lapply(LRTtests, function(x) subset(x$table, FDR<0.05 & abs(logFC) >= 0))

QLFtests.sel <- lapply(QLFtests, function(x) subset(x$table, FDR<0.05 & abs(logFC) >= 0))
#& abs(logFC) >= 1
## let's additionally sort our DEGs based on log fold change of expression:
LRTtests.sel <- lapply(LRTtests.sel, function(x) x[order(x$logFC),])
QLFtests.sel <- lapply(QLFtests.sel, function(x) x[order(x$logFC),])

## Here we print a small summary table of how many DEGs were found per contrast and per type of test
rbind("LRT"=sapply(LRTtests.sel, nrow), "QLF"=sapply(QLFtests.sel, nrow))

########### Further analysis on selected genes (FDR, logFC)

##### CONVERT 
LRTtests.sel <- lapply(lapply(LRTtests.sel, function(x) {x$Gene <- rownames(x);x}), function(x) x[,c(6,1:5)])
LRTtests.sel <- lapply(LRTtests.sel, function(x) {rownames(x) <- NULL;x})

LRTtests.noFDR <- lapply(lapply(lapply(lapply(LRTtests, function(x) data.frame(x)), function(x) {x$Gene <- rownames(x);x}), function(x) x[,c(6,1:5)]),
                       function(x) {rownames(x) <- NULL;x})

Astro <- dplyr::filter(LRTtests.sel$`con_1 - rad_2`, !Gene == "NA")[,c(1,2)]
colnames(Astro)[2] <- "logFC_1"

Astro.2 <- left_join(Astro, LRTtests.noFDR$`con_1 - rad_3`[,c(1,2)], by="Gene")
colnames(Astro.2)[3] <- "logFC_2"

Astro.3 <- left_join(Astro.2, LRTtests.noFDR$`con_1 - rad_4`[,c(1,2)], by="Gene")
colnames(Astro.3)[4] <- "logFC_3"

Astro.4 <- left_join(Astro.3, LRTtests.noFDR$`con_1 - rad_5`[,c(1,2)], by="Gene")
colnames(Astro.4)[5] <- "logFC_4"

rownames(Astro.4) <- Astro.4$Gene

Astro.4 <- Astro.4[,2:5]

colnames(Astro.4) <- c("5-30days\ninFlight", "60-120days\ninFlight", "3-30days\npostFlight", "60-120days\npostFlight")

#Figure_4G
ggplot(melt(dplyr::filter(Astro.4, Astro.4[,1] > 0)), aes(factor(variable), value)) +
  geom_boxplot(fill = "lightblue", alpha = 0.75, width=0.5) +
  geom_hline(yintercept=0, linetype="dashed", color = "blue", size=1, alpha = 0.5) +
  theme_classic() +
  labs(title="Up-regulated genes in human blood", x="Time after exposure (days)", y="logFC (versus 112-56days pre-flight)") +
  theme(
    axis.title.x = element_text(size = 8, face="bold"),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    axis.title.y = element_text(size = 8, face="bold"),
    title = element_text(size = 8, face = "bold"))

#Figure_4H
ggplot(melt(dplyr::filter(Astro.4, Astro.4[,1] < 0)), aes(factor(variable), value)) +
  geom_boxplot(fill = "lightblue", alpha = 0.75, width=0.5) +
  geom_hline(yintercept=0, linetype="dashed", color = "blue", size=1, alpha = 0.5) +
  theme_classic() +
  labs(title="Down-regulated genes in human blood", x="Time after exposure (days)", y="logFC (versus 112-56days pre-flight)") +
  theme(
    axis.title.x = element_text(size = 8, face="bold"),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    axis.title.y = element_text(size = 8, face="bold"),
    title = element_text(size = 8, face = "bold"))

#Trend
anova(lm(colMedians(as.matrix(dplyr::filter(Astro.4, Astro.4[,1] < 0)))~as.vector(1:4)))[[5]][1]
anova(lm(colMedians(as.matrix(dplyr::filter(Astro.4, Astro.4[,1] > 0)))~as.vector(1:4)))[[5]][1]
anova(lm(colMedians(as.matrix(Astro.4))~as.vector(1:4)))[[5]][1]


# 5. Astronauts - Hi-C ####
space.genes <- list.load("space.genes.rdata")

genes.hic <- lapply(space.genes, function(x) ("ENSEMBL"=as.character(x)))

genes.hic <- lapply(genes.hic, function(x) {AnnotationDbi::select(
  org.Hs.eg.db,
  keys = x,
  column = 'SYMBOL', 
  keytype = 'ENSEMBL')})

genes.hic <- lapply(genes.hic, function(x) data.frame("Gene"=x[,2]))

LRTtests.sel$`con_1 - rad_2` #Genes DE
LRTtests.noFDR #All genes

genes.hic.DE <- lapply(lapply(genes.hic, function(x) left_join(x, LRTtests.sel$`con_1 - rad_2`, by="Gene")), function(x) na.omit(x))
genes.hic.all <- lapply(lapply(genes.hic, function(x) left_join(x, LRTtests.noFDR$`con_1 - rad_2`, by="Gene")), function(x) na.omit(x))


#Frequency of DE genes in layers
normal.pen <- data.frame("ALLgenes"=t(data.frame(bind_rows(lapply(genes.hic.all[-1], function(x) nrow(x))))),
                         "DEgenes"=t(data.frame(bind_rows(lapply(genes.hic.DE[-1], function(x) nrow(x))))))

normal.pen.pct <- data.frame(apply(normal.pen,2,function(x){x/sum(x)}))

normal.pen.pct$Enrichment <- as.numeric(normal.pen.pct[,2]/normal.pen.pct[,1])
normal.pen.pct$Layers <- rownames(normal.pen.pct)

#Gene Expression Change
mean.logFC.DE <- data.frame(bind_rows(lapply(genes.hic.DE[-1], function(x) data.frame("mean_logFC"=mean(x[,2])))))
mean.logFC.all <- data.frame(bind_rows(lapply(genes.hic.all[-1], function(x) data.frame("mean_logFC"=mean(x[,2])))))

fig <- data.frame("mean.logFC"=rbind(mean.logFC.DE, mean.logFC.all), "Layer"=c("L1","L2","L3","L4","L5"), "Genes"=c(rep("DE",5), rep("all",5)))

#Figure_4F
ggplot(fig[1:5,], aes(x=Layer, y=mean_logFC)) +
  geom_point(size=4) +
  geom_line(aes(group=interaction(Genes)), size=2, alpha=0.5) +
  theme_classic() +
  labs(title="Pre-flight versus\nDuring-flight (5-30days from lunch)", x="Nucleus layer", y="Absolute mean logFC change versus control") +
  theme(
    axis.title.x = element_text(size = 10, face="bold"),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 10, face="bold"),
    title = element_text(size = 10, face = "bold"))

#Figure_S4M
ggplot(fig[6:10,], aes(x=Layer, y=mean_logFC)) +
  geom_point(size=4) +
  geom_line(aes(group=interaction(Genes)), size=2, alpha=0.5) +
  theme_classic() +
  labs(title="Pre-flight versus\nDuring-flight (5-30days from lunch)", x="Nucleus layer", y="Absolute mean logFC change versus control") +
  theme(
    axis.title.x = element_text(size = 10, face="bold"),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 10, face="bold"),
    title = element_text(size = 10, face = "bold"))


