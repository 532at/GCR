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
library(edgeR)
library(FactoMineR)
library(mogene10sttranscriptcluster.db)
library(factoextra)

# 2. Microarray for Fe - Heart ####

getwd()
wpath <- "/mnt/work1/adrian/Nasa/Microarray2/"
setwd(wpath)
options(width=200)
set.seed(123)

annotTable <- AnnotationDbi::select(
  mogene10sttranscriptcluster.db,
  keys = keys(mogene10sttranscriptcluster.db),
  column = c('PROBEID', 'SYMBOL', 'ENTREZID', 'ENSEMBL'),
  keytype = 'PROBEID')

annotTable <- annotTable[!is.na(annotTable$SYMBOL),]

array2.info <- read.csv("GLDS-109-samples.csv", header=TRUE, sep=",", check.names = F, stringsAsFactors = F)

# GSM1684884; GSM1684875 - removed after pca
array2.processed <- do.call(cbind, list(data.frame(data.table::fread("GSM1684876_sample_table.txt", header=TRUE, sep="\t", check.names = F, stringsAsFactors = F)),
                                        "GSM1684877"=data.frame(data.table::fread("GSM1684877_sample_table.txt", header=TRUE, sep="\t", check.names = F, stringsAsFactors = F))[,2],
                                        "GSM1684878"=data.frame(data.table::fread("GSM1684878_sample_table.txt", header=TRUE, sep="\t", check.names = F, stringsAsFactors = F))[,2],
                                        "GSM1684879"=data.frame(data.table::fread("GSM1684879_sample_table.txt", header=TRUE, sep="\t", check.names = F, stringsAsFactors = F))[,2],
                                        "GSM1684880"=data.frame(data.table::fread("GSM1684880_sample_table.txt", header=TRUE, sep="\t", check.names = F, stringsAsFactors = F))[,2],
                                        "GSM1684881"=data.frame(data.table::fread("GSM1684881_sample_table.txt", header=TRUE, sep="\t", check.names = F, stringsAsFactors = F))[,2],
                                        "GSM1684882"=data.frame(data.table::fread("GSM1684882_sample_table.txt", header=TRUE, sep="\t", check.names = F, stringsAsFactors = F))[,2],
                                        "GSM1684883"=data.frame(data.table::fread("GSM1684883_sample_table.txt", header=TRUE, sep="\t", check.names = F, stringsAsFactors = F))[,2],
                                        "GSM1684885"=data.frame(data.table::fread("GSM1684885_sample_table.txt", header=TRUE, sep="\t", check.names = F, stringsAsFactors = F))[,2],
                                        "GSM1684886"=data.frame(data.table::fread("GSM1684886_sample_table.txt", header=TRUE, sep="\t", check.names = F, stringsAsFactors = F))[,2],
                                        "GSM1684887"=data.frame(data.table::fread("GSM1684887_sample_table.txt", header=TRUE, sep="\t", check.names = F, stringsAsFactors = F))[,2],
                                        "GSM1684888"=data.frame(data.table::fread("GSM1684888_sample_table.txt", header=TRUE, sep="\t", check.names = F, stringsAsFactors = F))[,2]))

colnames(array2.processed)[1:2] <- c("PROBEID","GSM1684876")
array2.processed$PROBEID <- as.character(array2.processed$PROBEID)
rownames(array2.processed) <- array2.processed[,1]
array2.processed <- array2.processed[,-1]

## create a table with sample information / phenotypic information / any other information on your samples which might be relevant for samples grouping,
## variability of the expression data etc. Most importantly this table should contain the biological coeficient you want to test, e.g. CTRL vs. PT.
samples_tab <- data.frame("sample"=colnames(array2.processed)) 
samples <- colnames(array2.processed)

samples_tab$type <- c(rep("rad",10), rep("con",2))
samples_tab$time <- c(rep("1",2), rep("3",3), rep("7",2), rep("14",1), rep("28",2), rep("1",2))
samples_tab[,4] <- unite(samples_tab[,2:3], ID, sep = "_", remove = T)
samples_tab$ID <- factor(samples_tab$ID, levels=c("con_1","rad_1","rad_3","rad_7","rad_14", "rad_28"), labels=c("con_1","rad_1","rad_3","rad_7","rad_14", "rad_28"))

##### FILTER DATA

## Check the rowMeans
array2.medians <- matrixStats::rowMedians(as.matrix(array2.processed))

par(mfrow=c(1,1))
hist(array2.medians)

## Filter those with median less than 4
array2.filtered <- filter(array2.processed, matrixStats::rowMedians(as.matrix(array2.processed)) > 3)

## Removing multiple mappings

anno_grouped <- group_by(array2.filtered, "PROBEID"=rownames(array2.filtered))
anno_summarized <- 
  dplyr::summarize(anno_grouped, no_of_matches = n_distinct(PROBEID))

sum(anno_summarized$no_of_matches > 1)

##### VISUALIZE DATA
table(samples_tab$ID)
## define a variable that will hold colors we want to use for plots. Here we can define colors as numbers, or define them with literal color names
my_cols <- c(1:6)
## alternatively:
#my_cols <- c("green","cyan","blue","red","purple")
names(my_cols) <- c("con_1","rad_1","rad_3","rad_7","rad_14", "rad_28") # assign names to your colors. Here names are sample types
print(my_cols)

## plot distribution of counts per sample
par(mfrow=c(1,1))
boxplot(array2.filtered, las=2, cex=.5, cex.axis=.5, main="Raw counts, filtered, log") # we make pseudocounts of the raw data

## Visualize our data as MDS plot - multidimensional scaling plot of distances between gene expression profiles. Remember, MDS is not the same as PCA plot
plotMDS(array2.filtered, las=2, cex=1, col=my_cols[samples_tab$ID])

##### EDGER OBJECT 
## Here we set up an EdgeR object that contains expression data and sample description. We group samples per sample type (CTRL, PT, etc.)
dgeFull <- DGEList(array2.filtered, group=samples_tab$ID)
dgeFull$samples$Individual <- samples_tab$sample
dgeFull$samples

## Here we want to define the order of our sample groups - otherwise they are ordered alphabetically
dgeFull$samples$group <- factor(dgeFull$samples$group, levels=c("con_1","rad_1","rad_3","rad_7","rad_14", "rad_28"))

## PCA
## Here we run PCA analysis to yet again visualize our data, but also to check how it groups based on the biological type. We can also identify here any potential problems, outliers etc.
## Remember - for PCA we need a transposed table of normalized counts - genes as columns and samples as rows
norm_t <- t(array2.filtered)
#norm_t[1:5,1:5]
my.pca <- PCA(norm_t, graph=FALSE)

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
colnames(designMat) <- c("con_1","rad_1","rad_3","rad_7","rad_14", "rad_28")
rownames(designMat) <- rownames(dgeFull$samples)
print(designMat)

contMatrix <- makeContrasts(con_1-rad_1, con_1-rad_3, con_1-rad_7, con_1-rad_14, con_1-rad_28,
                            rad_1-rad_3, rad_1-rad_7, rad_1-rad_14, rad_1-rad_28,
                            levels=designMat)
print(contMatrix)
print(designMat)

fit <- lmFit(array2.filtered, designMat)
fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)

topTable(fit2, coef=1, adjust="BH")
results <- decideTests(fit2)
vennDiagram(results)

topTable(fit2, number=10, adjust="BH", coef=c(1:5))

array.results <- list("con_1-rad_1"=topTable(fit2, coef=1, adjust="BH", number=24396),
                      "con_1-rad_3"=topTable(fit2, coef=2, adjust="BH", number=24396),
                      "con_1-rad_7"=topTable(fit2, coef=3, adjust="BH", number=24396),
                      "con_1-rad_14"=topTable(fit2, coef=4, adjust="BH", number=24396),
                      "con_1-rad_28"=topTable(fit2, coef=5, adjust="BH", number=24396),
                      "rad_1-rad_3"=topTable(fit2, coef=6, adjust="BH", number=24396),
                      "rad_1-rad_7"=topTable(fit2, coef=7, adjust="BH", number=24396),
                      "rad_1-rad_14"=topTable(fit2, coef=8, adjust="BH", number=24396),
                      "rad_1-rad_28"=topTable(fit2, coef=9, adjust="BH", number=24396))

##### CONVERT 
array.results <- lapply(lapply(lapply(lapply(lapply(array.results, function(x) {x$PROBEID <- rownames(x);x}), function(x) left_join(x, annotTable[,1:2], by="PROBEID")), function(x) distinct(x)),
                        function(x) filter(x, !SYMBOL == "NA")), function(x) {x[,9] <- unite(x[,c(7:8)], ID, sep="_", remove=T);x})

array.sel.genes <- filter(filter(array.results$`con_1-rad_3`, adj.P.Val <=0.05), !SYMBOL == "NA")

# Remove duplicates by PROBEID and SYMBOL leaving those with the lowest P Val
array.sel.genes.ordered <- array.sel.genes[order(array.sel.genes$adj.P.Val, decreasing = FALSE),]
array.sel.genes.filtered <- array.sel.genes.ordered[!duplicated(array.sel.genes.ordered$PROBEID),]
array.sel.genes.filtered <- array.sel.genes.filtered[!duplicated(array.sel.genes.filtered$SYMBOL),]


array.heart <- bind_cols(lapply(lapply(array.results, function(x) left_join(array.sel.genes.filtered[,c(1,9)], x, by="ID")), function(x) x[,c(3,7,10)]))[,c(4,5,7,10,13,15)]
array.heart.ordered <- array.heart[order(array.heart$adj.P.Val...5, decreasing = FALSE),]
array.heart.filtered <- array.heart.ordered[!duplicated(array.heart.ordered$SYMBOL...15),]
array.heart.filtered <- array.heart.filtered[complete.cases(array.heart.filtered),][,-2]

colnames(array.heart.filtered) <- c("logFC_3", "logFC_7", "logFC_14", "logFC_28", "Gene")
rownames(array.heart.filtered) <- array.heart.filtered$Gene
array.heart.filtered <- array.heart.filtered[,-5]


array.heart.filtered$gene <- rownames(array.heart.filtered)

par(mfrow=c(1,1), mar=c(4,4,4,4))
colnames(array.heart.filtered) <- c("3", "7", "14", "28", "Gene")

#Figure_4F
ggplot(melt(filter(array.heart.filtered, array.heart.filtered[,1] > 0)), aes(factor(variable), value)) +
  geom_boxplot(fill = "lightblue", alpha = 0.75, width=0.5) +
  geom_hline(yintercept=0, linetype="dashed", color = "blue", size=1, alpha = 0.5) +
  theme_classic() +
  ylim(-2.5,2.5) +
  labs(title="Up-regulated genes in heart\ntissue (Fe 0.15 Gy)", x="Time after exposure (days)", y="logFC", color="Sample") +
  theme(
    axis.title.x = element_text(size = 8, face="bold"),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    axis.title.y = element_text(size = 8, face="bold"),
    title = element_text(size = 8, face = "bold"))

#Figure_4G
ggplot(melt(filter(array.heart.filtered, array.heart.filtered[,1] < 0)), aes(factor(variable), value)) +
  geom_boxplot(fill = "lightblue", alpha = 0.75, width=0.5) +
  geom_hline(yintercept=0, linetype="dashed", color = "blue", size=1, alpha = 0.5) +
  theme_classic() +
  ylim(-2.5,2.5) +
  labs(title="Down-regulated genes in heart\ntissue (Fe 0.15 Gy)", x="Time after exposure (days)", y="logFC", color="Sample") +
  theme(
    axis.title.x = element_text(size = 8, face="bold"),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    axis.title.y = element_text(size = 8, face="bold"),
    title = element_text(size = 8, face = "bold"))

#Trend
anova(lm(colMedians(apply(as.matrix(filter(array.heart.filtered, array.heart.filtered[,1] < 0))[,-5], 2, as.numeric))~as.vector(1:4)))[[5]][1]
anova(lm(colMedians(apply(as.matrix(filter(array.heart.filtered, array.heart.filtered[,1] > 0))[,-5], 2, as.numeric))~as.vector(1:4)))[[5]][1]

# 3. Microarray for Si - Breast ####

getwd()
wpath <- "/mnt/work1/adrian/Nasa/Microarray_Breast//"
setwd(wpath)
options(width=200)
set.seed(123)

annotTable <- AnnotationDbi::select(
  mogene10sttranscriptcluster.db,
  keys = keys(mogene10sttranscriptcluster.db),
  column = c('PROBEID', 'SYMBOL', 'ENTREZID', 'ENSEMBL'),
  keytype = 'PROBEID')[!is.na(annotTable$SYMBOL),]

array3.info <- read.csv("GLDS-80-samples.csv", header=TRUE, sep=",", check.names = F, stringsAsFactors = F)[,-3]
array3.info <- filter(array3.info, !array3.info$`Factor Value: Ionizing Radiation` == "Gamma Radiation")

# GSM1176961; GSM1176965 - removed after pca
array3.processed <- do.call(cbind, list(data.frame(data.table::fread("GSM1176957_sample_table.txt", header=TRUE, sep="\t", check.names = F, stringsAsFactors = F)),
                                        "GSM1176958"=data.frame(data.table::fread("GSM1176958_sample_table.txt", header=TRUE, sep="\t", check.names = F, stringsAsFactors = F))[,2],
                                        "GSM1176959"=data.frame(data.table::fread("GSM1176959_sample_table.txt", header=TRUE, sep="\t", check.names = F, stringsAsFactors = F))[,2],
                                        "GSM1176960"=data.frame(data.table::fread("GSM1176960_sample_table.txt", header=TRUE, sep="\t", check.names = F, stringsAsFactors = F))[,2],
                                        "GSM1176962"=data.frame(data.table::fread("GSM1176962_sample_table.txt", header=TRUE, sep="\t", check.names = F, stringsAsFactors = F))[,2],
                                        "GSM1176963"=data.frame(data.table::fread("GSM1176963_sample_table.txt", header=TRUE, sep="\t", check.names = F, stringsAsFactors = F))[,2],
                                        "GSM1176964"=data.frame(data.table::fread("GSM1176964_sample_table.txt", header=TRUE, sep="\t", check.names = F, stringsAsFactors = F))[,2],
                                        "GSM1176969"=data.frame(data.table::fread("GSM1176969_sample_table.txt", header=TRUE, sep="\t", check.names = F, stringsAsFactors = F))[,2],
                                        "GSM1176970"=data.frame(data.table::fread("GSM1176970_sample_table.txt", header=TRUE, sep="\t", check.names = F, stringsAsFactors = F))[,2],
                                        "GSM1176971"=data.frame(data.table::fread("GSM1176971_sample_table.txt", header=TRUE, sep="\t", check.names = F, stringsAsFactors = F))[,2],
                                        "GSM1176972"=data.frame(data.table::fread("GSM1176972_sample_table.txt", header=TRUE, sep="\t", check.names = F, stringsAsFactors = F))[,2],
                                        "GSM1176973"=data.frame(data.table::fread("GSM1176973_sample_table.txt", header=TRUE, sep="\t", check.names = F, stringsAsFactors = F))[,2],
                                        "GSM1176974"=data.frame(data.table::fread("GSM1176974_sample_table.txt", header=TRUE, sep="\t", check.names = F, stringsAsFactors = F))[,2],
                                        "GSM1176982"=data.frame(data.table::fread("GSM1176982_sample_table.txt", header=TRUE, sep="\t", check.names = F, stringsAsFactors = F))[,2],
                                        "GSM1176983"=data.frame(data.table::fread("GSM1176983_sample_table.txt", header=TRUE, sep="\t", check.names = F, stringsAsFactors = F))[,2],
                                        "GSM1176984"=data.frame(data.table::fread("GSM1176984_sample_table.txt", header=TRUE, sep="\t", check.names = F, stringsAsFactors = F))[,2],
                                        "GSM1176985"=data.frame(data.table::fread("GSM1176985_sample_table.txt", header=TRUE, sep="\t", check.names = F, stringsAsFactors = F))[,2],
                                        "GSM1176986"=data.frame(data.table::fread("GSM1176986_sample_table.txt", header=TRUE, sep="\t", check.names = F, stringsAsFactors = F))[,2],
                                        "GSM1176987"=data.frame(data.table::fread("GSM1176987_sample_table.txt", header=TRUE, sep="\t", check.names = F, stringsAsFactors = F))[,2]))

colnames(array3.processed)[1:2] <- c("PROBEID","GSM1176957")
array3.processed$PROBEID <- as.character(array3.processed$PROBEID)
rownames(array3.processed) <- array3.processed[,1]
array3.processed <- array3.processed[,-1]

## create a table with sample information / phenotypic information / any other information on your samples which might be relevant for samples grouping,
## variability of the expression data etc. Most importantly this table should contain the biological coeficient you want to test, e.g. CTRL vs. PT.
samples.breast_tab <- data.frame("sample"=colnames(array3.processed)) 
samples.breast <- colnames(array3.processed)

samples.breast_tab$type <- c(rep("rad",5), rep("con",2), rep("rad",6), rep("rad",6))
samples.breast_tab$dose <- c(rep("30CG", 3), rep("10CG", 2), rep("0CG", 2), rep("30CG", 3), rep("10CG", 3), rep("30CG", 3), rep("10CG", 3))
samples.breast_tab$time <- c(rep("1",7), rep("4",6), rep("12",6))
samples.breast_tab[,5] <- unite(samples.breast_tab[,2:4], ID, sep = "_", remove = T)
samples.breast_tab$ID <- factor(samples.breast_tab$ID, levels=c("con_0CG_1","rad_10CG_1","rad_30CG_1","rad_10CG_4","rad_30CG_4", "rad_10CG_12", "rad_30CG_12"),
                                labels=c("con_0CG_1","rad_10CG_1","rad_30CG_1","rad_10CG_4","rad_30CG_4", "rad_10CG_12", "rad_30CG_12"))

##### FILTER DATA

## Check the rowMeans
library(matrixStats)

array3.medians <- matrixStats::rowMedians(as.matrix(array3.processed))

par(mfrow=c(1,1))
hist(array3.medians)

## Filter those with median less than 4
array3.filtered <- filter(array3.processed, matrixStats::rowMedians(as.matrix(array3.processed)) > 3)

## Removing multiple mappings

anno_grouped <- group_by(array3.filtered, "PROBEID"=rownames(array3.filtered))
anno_summarized <- 
  dplyr::summarize(anno_grouped, no_of_matches = n_distinct(PROBEID))

sum(anno_summarized$no_of_matches > 1)

##### VISUALIZE DATA
table(samples.breast_tab$ID)
## define a variable that will hold colors we want to use for plots. Here we can define colors as numbers, or define them with literal color names
my_cols <- c(1:7)
## alternatively:
#my_cols <- c("green","cyan","blue","red","purple")
names(my_cols) <- c("con_0CG_1","rad_10CG_1","rad_30CG_1","rad_10CG_4","rad_30CG_4", "rad_10CG_12", "rad_30CG_12") # assign names to your colors. Here names are sample types
print(my_cols)

## plot distribution of counts per sample
par(mfrow=c(1,1))
boxplot(array3.filtered, las=2, cex=.5, cex.axis=.5, main="Raw counts, filtered, log") # we make pseudocounts of the raw data

## Visualize our data as MDS plot - multidimensional scaling plot of distances between gene expression profiles. Remember, MDS is not the same as PCA plot
plotMDS(array3.filtered, las=2, cex=1, col=my_cols[samples.breast_tab$ID])


##### EDGER OBJECT 

## Here we set up an EdgeR object that contains expression data and sample description. We group samples per sample type (CTRL, PT, etc.)
dgeFull <- DGEList(array3.filtered, group=samples.breast_tab$ID)
dgeFull$samples$Individual <- samples.breast_tab$sample
dgeFull$samples

## Here we want to define the order of our sample groups - otherwise they are ordered alphabetically
dgeFull$samples$group <- factor(dgeFull$samples$group, levels=c("con_0CG_1","rad_10CG_1","rad_30CG_1","rad_10CG_4","rad_30CG_4", "rad_10CG_12", "rad_30CG_12"))

## PCA

## Here we run PCA analysis to yet again visualize our data, but also to check how it groups based on the biological type. We can also identify here any potential problems, outliers etc.
## Remember - for PCA we need a transposed table of normalized counts - genes as columns and samples as rows
norm_t <- t(array3.filtered)
#norm_t[1:5,1:5]
my.pca <- PCA(norm_t, graph=FALSE)

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
colnames(designMat) <- c("con_0CG_1","rad_10CG_1","rad_30CG_1","rad_10CG_4","rad_30CG_4", "rad_10CG_12", "rad_30CG_12")
rownames(designMat) <- rownames(dgeFull$samples)
print(designMat)


contMatrix <- makeContrasts(con_0CG_1-rad_10CG_1, con_0CG_1-rad_30CG_1, 
                            con_0CG_1-rad_10CG_4, con_0CG_1-rad_30CG_4,
                            con_0CG_1-rad_10CG_12, con_0CG_1-rad_30CG_12,
                            rad_10CG_1-rad_10CG_4, rad_30CG_1-rad_30CG_4,
                            rad_10CG_1-rad_10CG_12, rad_30CG_1-rad_30CG_12,
                            levels=designMat)
print(contMatrix)
print(designMat)


fit <- lmFit(array3.filtered, designMat)
fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)

topTable(fit2, coef=1, adjust="BH")
results <- decideTests(fit2)

topTable(fit2, number=10, adjust="BH", coef=c(1:5))


array.results <- list("con_0CG_1-rad_10CG_1"=topTable(fit2, coef=1, adjust="BH", number=29377),
                      "con_0CG_1-rad_30CG_1"=topTable(fit2, coef=2, adjust="BH", number=29377),
                      "con_0CG_1-rad_10CG_4"=topTable(fit2, coef=3, adjust="BH", number=29377),
                      "con_0CG_1-rad_30CG_4"=topTable(fit2, coef=4, adjust="BH", number=29377),
                      "con_0CG_1-rad_10CG_12"=topTable(fit2, coef=5, adjust="BH", number=29377),
                      "con_0CG_1-rad_30CG_12"=topTable(fit2, coef=6, adjust="BH", number=29377),
                      "rad_10CG_1-rad_10CG_4"=topTable(fit2, coef=7, adjust="BH", number=29377),
                      "rad_30CG_1-rad_30CG_4"=topTable(fit2, coef=8, adjust="BH", number=29377),
                      "rad_10CG_1-rad_10CG_12"=topTable(fit2, coef=9, adjust="BH", number=29377),
                      "rad_30CG_1-rad_30CG_12"=topTable(fit2, coef=10, adjust="BH", number=29377))

##### CONVERT 
array.results.10Gy <- lapply(lapply(lapply(lapply(lapply(array.results, function(x) {x$PROBEID <- rownames(x);x}), function(x) left_join(x, annotTable[,1:2], by="PROBEID")), function(x) distinct(x)),
                               function(x) filter(x, !SYMBOL == "NA")), function(x) {x[,9] <- unite(x[,c(7:8)], ID, sep="_", remove=T);x})[c(1,3,5)]

array.results.30Gy <- lapply(lapply(lapply(lapply(lapply(array.results, function(x) {x$PROBEID <- rownames(x);x}), function(x) left_join(x, annotTable[,1:2], by="PROBEID")), function(x) distinct(x)),
                               function(x) filter(x, !SYMBOL == "NA")), function(x) {x[,9] <- unite(x[,c(7:8)], ID, sep="_", remove=T);x})[c(2,4,6)]
names(array.results.10Gy)

array.sel.genes.10Gy <- filter(filter(array.results.10Gy$`con_0CG_1-rad_10CG_1`, adj.P.Val <=0.05), !SYMBOL == "NA")
array.sel.genes.30Gy <- filter(filter(array.results.30Gy$`con_0CG_1-rad_30CG_1`, adj.P.Val <=0.05), !SYMBOL == "NA")

# 30Gy
# Remove duplicates by PROBEID and SYMBOL leaving those with the lowest P Val
array.sel.genes.ordered.30Gy <- array.sel.genes.30Gy[order(array.sel.genes.30Gy$adj.P.Val, decreasing = FALSE),]
array.sel.genes.filtered.30Gy <- array.sel.genes.ordered.30Gy[!duplicated(array.sel.genes.ordered.30Gy$PROBEID),]
array.sel.genes.filtered.30Gy <- array.sel.genes.filtered.30Gy[!duplicated(array.sel.genes.filtered.30Gy$SYMBOL),]

array.heart.30Gy <- bind_cols(lapply(lapply(array.results.30Gy, function(x) left_join(array.sel.genes.filtered.30Gy[,c(1,9)], x, by="ID")), function(x) x[,c(3,7,10)]))[,c(1,2,4,7,9)]

array.heart.ordered.30Gy <- array.heart.30Gy[order(array.heart.30Gy$adj.P.Val...2, decreasing = FALSE),]
array.heart.filtered.30Gy <- array.heart.ordered.30Gy[!duplicated(array.heart.ordered.30Gy$SYMBOL...9),]
array.heart.filtered.30Gy <- array.heart.filtered.30Gy[complete.cases(array.heart.filtered.30Gy),][,-2]

colnames(array.heart.filtered.30Gy) <- c("logFC_1", "logFC_4", "logFC_12", "Gene")
rownames(array.heart.filtered.30Gy) <- array.heart.filtered.30Gy$Gene
array.heart.filtered.30Gy <- array.heart.filtered.30Gy[,-4]

colnames(array.heart.filtered.30Gy) <- c("1", "4", "12")

#Figure_S3A
ggplot(melt(filter(array.heart.filtered.30Gy, array.heart.filtered.30Gy[,1] > 0)), aes(factor(variable), value)) +
  geom_boxplot(fill = "lightblue", alpha = 0.75, width=0.5) +
  geom_hline(yintercept=0, linetype="dashed", color = "blue", size=1, alpha = 0.5) +
  theme_classic() +
  ylim(-1.5,1.5) +
  labs(title="Up-regulated genes in breast\ntissue (Si 0.3 Gy)", x="Time after exposure (weeks)", y="logFC", color="Sample") +
  theme(
    axis.title.x = element_text(size = 8, face="bold"),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    axis.title.y = element_text(size = 8, face="bold"),
    title = element_text(size = 8, face = "bold"))

#Figure_S3B
ggplot(melt(filter(array.heart.filtered.30Gy, array.heart.filtered.30Gy[,1] < 0)), aes(factor(variable), value)) +
  geom_boxplot(fill = "lightblue", alpha = 0.75, width=0.5) +
  geom_hline(yintercept=0, linetype="dashed", color = "blue", size=1, alpha = 0.5) +
  theme_classic() +
  ylim(-1.5,1.5) +
  labs(title="Down-regulated genes in breast\ntissue (Si 0.3 Gy)", x="Time after exposure (weeks)", y="logFC", color="Sample") +
  theme(
    axis.title.x = element_text(size = 8, face="bold"),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    axis.title.y = element_text(size = 8, face="bold"),
    title = element_text(size = 8, face = "bold"))

#Trend
anova(lm(colMedians(as.matrix(filter(array.heart.filtered.30Gy, array.heart.filtered.30Gy[,1] < 0)))~as.vector(1:3)))[[5]][1]
anova(lm(colMedians(as.matrix(filter(array.heart.filtered.30Gy, array.heart.filtered.30Gy[,1] > 0)))~as.vector(1:3)))[[5]][1]

# 4. Microarray for Fe - Heart _ Protons ####

getwd()
wpath <- "/mnt/work1/adrian/Nasa/Microarray_Heart_Proton/"
setwd(wpath)
options(width=200)
set.seed(123)

annotTable <- AnnotationDbi::select(
  mogene10sttranscriptcluster.db,
  keys = keys(mogene10sttranscriptcluster.db),
  column = c('PROBEID', 'SYMBOL', 'ENTREZID', 'ENSEMBL'),
  keytype = 'PROBEID')

annotTable <- annotTable[!is.na(annotTable$SYMBOL),]

array4.info <- read.csv("GLDS-117-samples.csv", header=TRUE, sep=",", check.names = F, stringsAsFactors = F)
array4.info[,c(1,2,8,9)]

# GSM1684884; GSM1684875 - removed after pca
array4.processed <- do.call(cbind, list(data.frame(data.table::fread("GSM1684890_sample_table.txt", header=TRUE, sep="\t", check.names = F, stringsAsFactors = F)),
                                        "GSM1684891"=data.frame(data.table::fread("GSM1684891_sample_table.txt", header=TRUE, sep="\t", check.names = F, stringsAsFactors = F))[,2],
                                        "GSM1684892"=data.frame(data.table::fread("GSM1684892_sample_table.txt", header=TRUE, sep="\t", check.names = F, stringsAsFactors = F))[,2],
                                        "GSM1684893"=data.frame(data.table::fread("GSM1684893_sample_table.txt", header=TRUE, sep="\t", check.names = F, stringsAsFactors = F))[,2],
                                        "GSM1684894"=data.frame(data.table::fread("GSM1684894_sample_table.txt", header=TRUE, sep="\t", check.names = F, stringsAsFactors = F))[,2],
                                        "GSM1684895"=data.frame(data.table::fread("GSM1684895_sample_table.txt", header=TRUE, sep="\t", check.names = F, stringsAsFactors = F))[,2],
                                        "GSM1684896"=data.frame(data.table::fread("GSM1684896_sample_table.txt", header=TRUE, sep="\t", check.names = F, stringsAsFactors = F))[,2],
                                        "GSM1684897"=data.frame(data.table::fread("GSM1684897_sample_table.txt", header=TRUE, sep="\t", check.names = F, stringsAsFactors = F))[,2],
                                        "GSM1684898"=data.frame(data.table::fread("GSM1684898_sample_table.txt", header=TRUE, sep="\t", check.names = F, stringsAsFactors = F))[,2],
                                        "GSM1684899"=data.frame(data.table::fread("GSM1684899_sample_table.txt", header=TRUE, sep="\t", check.names = F, stringsAsFactors = F))[,2],
                                        "GSM1684900"=data.frame(data.table::fread("GSM1684900_sample_table.txt", header=TRUE, sep="\t", check.names = F, stringsAsFactors = F))[,2],
                                        "GSM1684901"=data.frame(data.table::fread("GSM1684901_sample_table.txt", header=TRUE, sep="\t", check.names = F, stringsAsFactors = F))[,2],
                                        "GSM1684902"=data.frame(data.table::fread("GSM1684902_sample_table.txt", header=TRUE, sep="\t", check.names = F, stringsAsFactors = F))[,2],
                                        "GSM1684903"=data.frame(data.table::fread("GSM1684903_sample_table.txt", header=TRUE, sep="\t", check.names = F, stringsAsFactors = F))[,2]))

colnames(array4.processed)[1:2] <- c("PROBEID","GSM1684890")
array4.processed$PROBEID <- as.character(array4.processed$PROBEID)
rownames(array4.processed) <- array4.processed[,1]
array4.processed <- array4.processed[,-1]

## create a table with sample information / phenotypic information / any other information on your samples which might be relevant for samples grouping,
## variability of the expression data etc. Most importantly this table should contain the biological coeficient you want to test, e.g. CTRL vs. PT.
samples_tab <- data.frame("sample"=colnames(array4.processed)) 
samples <- colnames(array4.processed)

samples_tab$type <- c(rep("rad",12), rep("con",2))
samples_tab$time <- c(rep("1",2), rep("3",3), rep("5",3), rep("12",2), rep("26",2), rep("1",2))
samples_tab[,4] <- unite(samples_tab[,2:3], ID, sep = "_", remove = T)
samples_tab$ID <- factor(samples_tab$ID, levels=c("con_1","rad_1","rad_3","rad_5","rad_12", "rad_26"), labels=c("con_1","rad_1","rad_3","rad_5","rad_12", "rad_26"))
samples_tab
##### FILTER DATA

## Check the rowMeans
array4.medians <- matrixStats::rowMedians(as.matrix(array4.processed))

par(mfrow=c(1,1))
hist(array4.medians)

## Filter those with median less than 4
array4.filtered <- filter(array4.processed, matrixStats::rowMedians(as.matrix(array4.processed)) > 3)

## Removing multiple mappings

anno_grouped <- group_by(array4.filtered, "PROBEID"=rownames(array4.filtered))
anno_summarized <- 
  dplyr::summarize(anno_grouped, no_of_matches = n_distinct(PROBEID))

sum(anno_summarized$no_of_matches > 1)

dim(array2.filtered)

##### VISUALIZE DATA
table(samples_tab$ID)
## define a variable that will hold colors we want to use for plots. Here we can define colors as numbers, or define them with literal color names
my_cols <- c(1:6)
## alternatively:
#my_cols <- c("green","cyan","blue","red","purple")
names(my_cols) <- c("con_1","rad_1","rad_3","rad_5","rad_12", "rad_26") # assign names to your colors. Here names are sample types
print(my_cols)

## plot distribution of counts per sample
par(mfrow=c(1,1))
boxplot(array4.filtered, las=2, cex=.5, cex.axis=.5, main="Raw counts, filtered, log") # we make pseudocounts of the raw data

## Visualize our data as MDS plot - multidimensional scaling plot of distances between gene expression profiles. Remember, MDS is not the same as PCA plot
plotMDS(array4.filtered, las=2, cex=1, col=my_cols[samples_tab$ID])


##### EDGER OBJECT 

## Here we set up an EdgeR object that contains expression data and sample description. We group samples per sample type (CTRL, PT, etc.)
dgeFull <- DGEList(array4.filtered, group=samples_tab$ID)
dgeFull$samples$Individual <- samples_tab$sample
dgeFull$samples

## Here we want to define the order of our sample groups - otherwise they are ordered alphabetically
dgeFull$samples$group <- factor(dgeFull$samples$group, levels=c("con_1","rad_1","rad_3","rad_5","rad_12", "rad_26"))

## PCA

## Here we run PCA analysis to yet again visualize our data, but also to check how it groups based on the biological type. We can also identify here any potential problems, outliers etc.
## Remember - for PCA we need a transposed table of normalized counts - genes as columns and samples as rows
norm_t <- t(array4.filtered)
#norm_t[1:5,1:5]
my.pca <- PCA(norm_t, graph=FALSE)

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
colnames(designMat) <- c("con_1","rad_1","rad_3","rad_5","rad_12", "rad_26")
rownames(designMat) <- rownames(dgeFull$samples)
print(designMat)


contMatrix <- makeContrasts(con_1-rad_1, con_1-rad_3, con_1-rad_5, con_1-rad_12, con_1-rad_26,
                            rad_1-rad_3, rad_1-rad_5, rad_1-rad_12, rad_1-rad_26,
                            levels=designMat)
print(contMatrix)
print(designMat)


fit <- lmFit(array4.filtered, designMat)
fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)

topTable(fit2, coef=1, adjust="BH")
results <- decideTests(fit2)

topTable(fit2, number=10, adjust="BH", coef=c(1:5))

array.results <- list("con_1-rad_1"=topTable(fit2, coef=1, adjust="BH", number=24361),
                      "con_1-rad_3"=topTable(fit2, coef=2, adjust="BH", number=24361),
                      "con_1-rad_5"=topTable(fit2, coef=3, adjust="BH", number=24361),
                      "con_1-rad_12"=topTable(fit2, coef=4, adjust="BH", number=24361),
                      "con_1-rad_26"=topTable(fit2, coef=5, adjust="BH", number=24361),
                      "rad_1-rad_3"=topTable(fit2, coef=6, adjust="BH", number=24361),
                      "rad_1-rad_5"=topTable(fit2, coef=7, adjust="BH", number=24361),
                      "rad_1-rad_12"=topTable(fit2, coef=8, adjust="BH", number=24361),
                      "rad_1-rad_26"=topTable(fit2, coef=9, adjust="BH", number=24361))


hist(array.results$`con_1-rad_3`$P.Value)
sum(array.results$`con_1-rad_1`$adj.P.Val <= 0.05)

##### CONVERT 
array.results <- lapply(lapply(lapply(lapply(lapply(array.results, function(x) {x$PROBEID <- rownames(x);x}), function(x) left_join(x, annotTable[,1:2], by="PROBEID")), function(x) distinct(x)),
                               function(x) filter(x, !SYMBOL == "NA")), function(x) {x[,9] <- unite(x[,c(7:8)], ID, sep="_", remove=T);x})

array.sel.genes <- filter(filter(array.results$`con_1-rad_3`, adj.P.Val <=0.05), !SYMBOL == "NA")

# Remove duplicates by PROBEID and SYMBOL leaving those with the lowest P Val
array.sel.genes.ordered <- array.sel.genes[order(array.sel.genes$adj.P.Val, decreasing = FALSE),]
array.sel.genes.filtered <- array.sel.genes.ordered[!duplicated(array.sel.genes.ordered$PROBEID),]
array.sel.genes.filtered <- array.sel.genes.filtered[!duplicated(array.sel.genes.filtered$SYMBOL),]


array.heart <- bind_cols(lapply(lapply(array.results, function(x) left_join(array.sel.genes.filtered[,c(1,9)], x, by="ID")), function(x) x[,c(3,7,10)]))[,c(4,5,7,10,13,15)]
head(array.heart)
array.heart.ordered <- array.heart[order(array.heart$adj.P.Val...5, decreasing = FALSE),]
array.heart.filtered <- array.heart.ordered[!duplicated(array.heart.ordered$SYMBOL...15),]
array.heart.filtered <- array.heart.filtered[complete.cases(array.heart.filtered),][,-2]

colnames(array.heart.filtered) <- c("logFC_3", "logFC_5", "logFC_12", "logFC_26", "Gene")
rownames(array.heart.filtered) <- array.heart.filtered$Gene
array.heart.filtered <- array.heart.filtered[,-5]


colnames(array.heart.filtered) <- c("3", "5", "12", "26")

#Figure_S3C
ggplot(melt(filter(array.heart.filtered, array.heart.filtered[,1] < 0)), aes(factor(variable), value)) +
  geom_boxplot(fill = "lightblue", alpha = 0.75, width=0.5) +
  geom_hline(yintercept=0, linetype="dashed", color = "blue", size=1, alpha = 0.5) +
  theme_classic() +
  labs(title="Down-regulated genes in retina\ntissue (H 1 Gy)", x="Time after exposure (days)", y="logFC", color="Sample") +
  theme(
    axis.title.x = element_text(size = 8, face="bold"),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    axis.title.y = element_text(size = 8, face="bold"),
    title = element_text(size = 8, face = "bold"))

#Trend
anova(lm(colMedians(as.matrix(filter(array.heart.filtered, array.heart.filtered[,1] < 0)))~as.vector(1:4)))[[5]][1]
anova(lm(colMedians(as.matrix(filter(array.heart.filtered, array.heart.filtered[,1] > 0)))~as.vector(1:4)))[[5]][1]

# 5. Microarray for Si - Breast _ Gamma Rays ####

getwd()
wpath <- "/mnt/work1/adrian/Nasa/Microarray_Breast//"
setwd(wpath)
options(width=200)
set.seed(123)

annotTable <- AnnotationDbi::select(
  mogene10sttranscriptcluster.db,
  keys = keys(mogene10sttranscriptcluster.db),
  column = c('PROBEID', 'SYMBOL', 'ENTREZID', 'ENSEMBL'),
  keytype = 'PROBEID')[!is.na(annotTable$SYMBOL),]

array5.info <- read.csv("GLDS-80-samples.csv", header=TRUE, sep=",", check.names = F, stringsAsFactors = F)[,-3]
array5.info <- filter(array3.info, !array3.info$`Factor Value: Ionizing Radiation` == "High-LET Si Particle Radiation (350 MeV/amu)")

# GSM1176961; GSM1176965 - removed after pce
array5.processed <- do.call(cbind, list(data.frame(data.table::fread("GSM1176954_sample_table.txt", header=TRUE, sep="\t", check.names = F, stringsAsFactors = F)),
                                        "GSM1176955"=data.frame(data.table::fread("GSM1176958_sample_table.txt", header=TRUE, sep="\t", check.names = F, stringsAsFactors = F))[,2],
                                        "GSM1176956"=data.frame(data.table::fread("GSM1176959_sample_table.txt", header=TRUE, sep="\t", check.names = F, stringsAsFactors = F))[,2],
                                        "GSM1176963"=data.frame(data.table::fread("GSM1176963_sample_table.txt", header=TRUE, sep="\t", check.names = F, stringsAsFactors = F))[,2],
                                        "GSM1176964"=data.frame(data.table::fread("GSM1176964_sample_table.txt", header=TRUE, sep="\t", check.names = F, stringsAsFactors = F))[,2],
                                        "GSM1176966"=data.frame(data.table::fread("GSM1176969_sample_table.txt", header=TRUE, sep="\t", check.names = F, stringsAsFactors = F))[,2],
                                        "GSM1176967"=data.frame(data.table::fread("GSM1176970_sample_table.txt", header=TRUE, sep="\t", check.names = F, stringsAsFactors = F))[,2],
                                        "GSM1176968"=data.frame(data.table::fread("GSM1176971_sample_table.txt", header=TRUE, sep="\t", check.names = F, stringsAsFactors = F))[,2],
                                        "GSM1176979"=data.frame(data.table::fread("GSM1176972_sample_table.txt", header=TRUE, sep="\t", check.names = F, stringsAsFactors = F))[,2],
                                        "GSM1176980"=data.frame(data.table::fread("GSM1176973_sample_table.txt", header=TRUE, sep="\t", check.names = F, stringsAsFactors = F))[,2],
                                        "GSM1176981"=data.frame(data.table::fread("GSM1176974_sample_table.txt", header=TRUE, sep="\t", check.names = F, stringsAsFactors = F))[,2]))

colnames(array5.processed)[1:2] <- c("PROBEID","GSM1176954")
array5.processed$PROBEID <- as.character(array5.processed$PROBEID)
rownames(array5.processed) <- array5.processed[,1]
array5.processed <- array5.processed[,-1]

## create a table with sample information / phenotypic information / any other information on your samples which might be relevant for samples grouping,
## variability of the expression data etc. Most importantly this table should contain the biological coeficient you want to test, e.g. CTRL vs. PT.
samples.breast_tab <- data.frame("sample"=colnames(array5.processed)) 
samples.breast <- colnames(array3.processed)

samples.breast_tab$type <- c(rep("rad",3), rep("con",2), rep("rad",3), rep("rad",3))
samples.breast_tab$dose <- c(rep("100CG", 3), rep("0CG", 2), rep("100CG", 3), rep("100CG", 3))
samples.breast_tab$time <- c(rep("1",5), rep("4",3), rep("12",3))
samples.breast_tab[,5] <- unite(samples.breast_tab[,2:4], ID, sep = "_", remove = T)
samples.breast_tab$ID <- factor(samples.breast_tab$ID, levels=c("con_0CG_1","rad_100CG_1","rad_100CG_4","rad_100CG_12"),
                                labels=c("con_0CG_1","rad_100CG_1","rad_100CG_4","rad_100CG_12"))

##### FILTER DATA

## Check the rowMeans
array5.medians <- matrixStats::rowMedians(as.matrix(array5.processed))

par(mfrow=c(1,1))
hist(array5.medians)

## Filter those with median less than 4
array5.filtered <- filter(array5.processed, matrixStats::rowMedians(as.matrix(array5.processed)) > 3)

## Removing multiple mappings

anno_grouped <- group_by(array5.filtered, "PROBEID"=rownames(array5.filtered))
anno_summarized <- 
  dplyr::summarize(anno_grouped, no_of_matches = n_distinct(PROBEID))

sum(anno_summarized$no_of_matches > 1)


##### VISUALIZE DATA
table(samples.breast_tab$ID)
## define a variable that will hold colors we want to use for plots. Here we can define colors as numbers, or define them with literal color names
my_cols <- c(1:4)
## alternatively:
#my_cols <- c("green","cyan","blue","red","purple")
names(my_cols) <- c("con_0CG_1","rad_100CG_1","rad_100CG_4","rad_100CG_12") # assign names to your colors. Here names are sample types
print(my_cols)

## plot distribution of counts per sample
par(mfrow=c(1,1))
boxplot(array5.filtered, las=2, cex=.5, cex.axis=.5, main="Raw counts, filtered, log") # we make pseudocounts of the raw data

## Visualize our data as MDS plot - multidimensional scaling plot of distances between gene expression profiles. Remember, MDS is not the same as PCA plot
plotMDS(array5.filtered, las=2, cex=1, col=my_cols[samples.breast_tab$ID])


##### EDGER OBJECT 

## Here we set up an EdgeR object that contains expression data and sample description. We group samples per sample type (CTRL, PT, etc.)
dgeFull <- DGEList(array5.filtered, group=samples.breast_tab$ID)
dgeFull$samples$Individual <- samples.breast_tab$sample
dgeFull$samples

## Here we want to define the order of our sample groups - otherwise they are ordered alphabetically
dgeFull$samples$group <- factor(dgeFull$samples$group, levels=c("con_0CG_1","rad_100CG_1","rad_100CG_4","rad_100CG_12"))

## PCA

## Here we run PCA analysis to yet again visualize our data, but also to check how it groups based on the biological type. We can also identify here any potential problems, outliers etc.
## Remember - for PCA we need a transposed table of normalized counts - genes as columns and samples as rows
norm_t <- t(array5.filtered)
#norm_t[1:5,1:5]
my.pca <- PCA(norm_t, graph=FALSE)

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
colnames(designMat) <- c("con_0CG_1","rad_100CG_1","rad_100CG_4","rad_100CG_12")
rownames(designMat) <- rownames(dgeFull$samples)
print(designMat)


contMatrix <- makeContrasts(con_0CG_1-rad_100CG_1, con_0CG_1-rad_100CG_4, con_0CG_1-rad_100CG_12,
                            rad_100CG_1-rad_100CG_4, rad_100CG_1-rad_100CG_12,
                            levels=designMat)
print(contMatrix)
print(designMat)


fit <- lmFit(array5.filtered, designMat)
fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)

topTable(fit2, coef=1, adjust="BH")
results <- decideTests(fit2)

topTable(fit2, number=10, adjust="BH", coef=c(1:5))


array.results <- list("con_0CG_1-rad_100CG_1"=topTable(fit2, coef=1, adjust="BH", number=23217),
                      "con_0CG_1-rad_100CG_4"=topTable(fit2, coef=2, adjust="BH", number=23217),
                      "con_0CG_1-rad_100CG_12"=topTable(fit2, coef=3, adjust="BH", number=23217),
                      "rad_100CG_1-rad_100CG_4"=topTable(fit2, coef=4, adjust="BH", number=23217),
                      "rad_100CG_1-rad_100CG_12"=topTable(fit2, coef=5, adjust="BH", number=23217))


hist(array.results$`con_0CG_1-rad_100CG_4`$adj.P.Val)
sum(array.results$`con_0CG_1-rad_100CG_4`$adj.P.Val <= 0.05)

##### CONVERT 
array.results.100Gy <- lapply(lapply(lapply(lapply(lapply(array.results, function(x) {x$PROBEID <- rownames(x);x}), function(x) left_join(x, annotTable[,1:2], by="PROBEID")), function(x) distinct(x)),
                                    function(x) filter(x, !SYMBOL == "NA")), function(x) {x[,9] <- unite(x[,c(7:8)], ID, sep="_", remove=T);x})

array.sel.genes.100Gy <- filter(filter(array.results.100Gy$`con_0CG_1-rad_100CG_4`, adj.P.Val <=0.05), !SYMBOL == "NA")

# 100Gy
# Remove duplicates by PROBEID and SYMBOL leaving those with the lowest P Val
array.sel.genes.ordered.100Gy <- array.sel.genes.100Gy[order(array.sel.genes.100Gy$adj.P.Val, decreasing = FALSE),]
array.sel.genes.filtered.100Gy <- array.sel.genes.ordered.100Gy[!duplicated(array.sel.genes.ordered.100Gy$PROBEID),]
array.sel.genes.filtered.100Gy <- array.sel.genes.filtered.100Gy[!duplicated(array.sel.genes.filtered.100Gy$SYMBOL),]

array.heart.100Gy <- bind_cols(lapply(lapply(array.results.100Gy, function(x) left_join(array.sel.genes.filtered.100Gy[,c(1,9)], x, by="ID")), function(x) x[,c(3,7,10)]))[,c(1,2,4,7,9)]

array.heart.ordered.100Gy <- array.heart.100Gy[order(array.heart.100Gy$adj.P.Val...2, decreasing = FALSE),]
array.heart.filtered.100Gy <- array.heart.ordered.100Gy[!duplicated(array.heart.ordered.100Gy$SYMBOL...9),]
array.heart.filtered.100Gy <- array.heart.filtered.100Gy[complete.cases(array.heart.filtered.100Gy),][,-2]

colnames(array.heart.filtered.100Gy) <- c("logFC_1", "logFC_4", "logFC_12", "Gene")
rownames(array.heart.filtered.100Gy) <- array.heart.filtered.100Gy$Gene
array.heart.filtered.100Gy <- array.heart.filtered.100Gy[,-4]


colnames(array.heart.filtered.100Gy) <- c("1", "4", "12")

#Figure_S5E
ggplot(melt(filter(array.heart.filtered.100Gy, array.heart.filtered.100Gy[,1] > 0)), aes(factor(variable), value)) +
  geom_boxplot(fill = "lightblue", alpha = 0.75, width=0.5) +
  geom_hline(yintercept=0, linetype="dashed", color = "blue", size=1, alpha = 0.5) +
  theme_classic() +
  ylim(-2,2) +
  labs(title="Up-regulated genes in breast\ntissue (Gamma 1 Gy)", x="Time after exposure (weeks)", y="logFC", color="Sample") +
  theme(
    axis.title.x = element_text(size = 8, face="bold"),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    axis.title.y = element_text(size = 8, face="bold"),
    title = element_text(size = 8, face = "bold"))

#Figure_S5F
ggplot(melt(filter(array.heart.filtered.100Gy, array.heart.filtered.100Gy[,1] < 0)), aes(factor(variable), value)) +
  geom_boxplot(fill = "lightblue", alpha = 0.75, width=0.5) +
  geom_hline(yintercept=0, linetype="dashed", color = "blue", size=1, alpha = 0.5) +
  theme_classic() +
  ylim(-2,2) +
  labs(title="Down-regulated genes in breast \ntissue (Gamma 1 Gy)", x="Time after exposure (weeks)", y="logFC", color="Sample") +
  theme(
    axis.title.x = element_text(size = 8, face="bold"),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    axis.title.y = element_text(size = 8, face="bold"),
    title = element_text(size = 8, face = "bold"))

#Trend
anova(lm(colMedians(as.matrix(filter(array.heart.filtered.100Gy, array.heart.filtered.100Gy[,1] < 0)))~as.vector(1:3)))[[5]][1]
anova(lm(colMedians(as.matrix(filter(array.heart.filtered.100Gy, array.heart.filtered.100Gy[,1] > 0)))~as.vector(1:3)))[[5]][1]

# 6. Microarray for Cs137 - Blood ####

getwd()
wpath <- "/mnt/work1/adrian/Nasa/Microarray_Blood//"
setwd(wpath)
options(width=200)
set.seed(123)

annotTable <- AnnotationDbi::select(
  mogene10sttranscriptcluster.db,
  keys = keys(mogene10sttranscriptcluster.db),
  column = c('PROBEID', 'SYMBOL', 'ENTREZID', 'ENSEMBL'),
  keytype = 'PROBEID')

annotTable <- annotTable[!is.na(annotTable$SYMBOL),]

array6.info <- read.csv("GLDS-159-samples.csv", header=TRUE, sep=",", check.names = F, stringsAsFactors = F)
array6.info$ID <- gsub("WholeBlood_", "", array6.info$`Comment: Sample_title`)
array6.info <- setorder(array6.info, `Comment: Sample_title`)[c(7:12,25:39,43:48),]
array6.info$ID
array6.info[,c(1,12)]

array6.processed <- do.call(cbind, list(data.frame(data.table::fread("GSM1274150_sample_table.txt", header=TRUE, sep="\t", check.names = F, stringsAsFactors = F)),
                                        "GSM1274151"=data.frame(data.table::fread("GSM1274151_sample_table.txt", header=TRUE, sep="\t", check.names = F, stringsAsFactors = F))[,2],
                                        "GSM1274152"=data.frame(data.table::fread("GSM1274152_sample_table.txt", header=TRUE, sep="\t", check.names = F, stringsAsFactors = F))[,2],
                                        "GSM1274153"=data.frame(data.table::fread("GSM1274153_sample_table.txt", header=TRUE, sep="\t", check.names = F, stringsAsFactors = F))[,2],
                                        "GSM1274154"=data.frame(data.table::fread("GSM1274154_sample_table.txt", header=TRUE, sep="\t", check.names = F, stringsAsFactors = F))[,2],
                                        "GSM1274155"=data.frame(data.table::fread("GSM1274155_sample_table.txt", header=TRUE, sep="\t", check.names = F, stringsAsFactors = F))[,2],
                                        "GSM1274186"=data.frame(data.table::fread("GSM1274186_sample_table.txt", header=TRUE, sep="\t", check.names = F, stringsAsFactors = F))[,2],
                                        "GSM1274187"=data.frame(data.table::fread("GSM1274187_sample_table.txt", header=TRUE, sep="\t", check.names = F, stringsAsFactors = F))[,2],
                                        "GSM1274188"=data.frame(data.table::fread("GSM1274188_sample_table.txt", header=TRUE, sep="\t", check.names = F, stringsAsFactors = F))[,2],
                                        "GSM1274189"=data.frame(data.table::fread("GSM1274189_sample_table.txt", header=TRUE, sep="\t", check.names = F, stringsAsFactors = F))[,2],
                                        "GSM1274190"=data.frame(data.table::fread("GSM1274190_sample_table.txt", header=TRUE, sep="\t", check.names = F, stringsAsFactors = F))[,2],
                                        "GSM1274191"=data.frame(data.table::fread("GSM1274191_sample_table.txt", header=TRUE, sep="\t", check.names = F, stringsAsFactors = F))[,2],
                                        "GSM1274162"=data.frame(data.table::fread("GSM1274162_sample_table.txt", header=TRUE, sep="\t", check.names = F, stringsAsFactors = F))[,2],
                                        "GSM1274163"=data.frame(data.table::fread("GSM1274163_sample_table.txt", header=TRUE, sep="\t", check.names = F, stringsAsFactors = F))[,2],
                                        "GSM1274164"=data.frame(data.table::fread("GSM1274164_sample_table.txt", header=TRUE, sep="\t", check.names = F, stringsAsFactors = F))[,2],
                                        "GSM1274165"=data.frame(data.table::fread("GSM1274165_sample_table.txt", header=TRUE, sep="\t", check.names = F, stringsAsFactors = F))[,2],
                                        "GSM1274166"=data.frame(data.table::fread("GSM1274166_sample_table.txt", header=TRUE, sep="\t", check.names = F, stringsAsFactors = F))[,2],
                                        "GSM1274167"=data.frame(data.table::fread("GSM1274167_sample_table.txt", header=TRUE, sep="\t", check.names = F, stringsAsFactors = F))[,2],
                                        "GSM1274168"=data.frame(data.table::fread("GSM1274168_sample_table.txt", header=TRUE, sep="\t", check.names = F, stringsAsFactors = F))[,2],
                                        "GSM1274169"=data.frame(data.table::fread("GSM1274169_sample_table.txt", header=TRUE, sep="\t", check.names = F, stringsAsFactors = F))[,2],
                                        "GSM1274170"=data.frame(data.table::fread("GSM1274170_sample_table.txt", header=TRUE, sep="\t", check.names = F, stringsAsFactors = F))[,2],
                                        "GSM1274174"=data.frame(data.table::fread("GSM1274174_sample_table.txt", header=TRUE, sep="\t", check.names = F, stringsAsFactors = F))[,2],
                                        "GSM1274175"=data.frame(data.table::fread("GSM1274175_sample_table.txt", header=TRUE, sep="\t", check.names = F, stringsAsFactors = F))[,2],
                                        "GSM1274176"=data.frame(data.table::fread("GSM1274176_sample_table.txt", header=TRUE, sep="\t", check.names = F, stringsAsFactors = F))[,2],
                                        "GSM1274177"=data.frame(data.table::fread("GSM1274177_sample_table.txt", header=TRUE, sep="\t", check.names = F, stringsAsFactors = F))[,2],
                                        "GSM1274178"=data.frame(data.table::fread("GSM1274178_sample_table.txt", header=TRUE, sep="\t", check.names = F, stringsAsFactors = F))[,2],
                                        "GSM1274179"=data.frame(data.table::fread("GSM1274179_sample_table.txt", header=TRUE, sep="\t", check.names = F, stringsAsFactors = F))[,2]))

colnames(array6.processed)[1:2] <- c("PROBEID","GSM1274150")
array6.processed$PROBEID <- as.character(array6.processed$PROBEID)
rownames(array6.processed) <- array6.processed[,1]
array6.processed <- array6.processed[,-1]
array6.processed <- array6.processed[complete.cases(array6.processed),] #30k probes has NAs

## create a table with sample information / phenotypic information / any other information on your samples which might be relevant for samples grouping,
## variability of the expression data etc. Most importantly this table should contain the biological coeficient you want to test, e.g. CTRL vs. PT.
samples_tab <- data.frame("sample"=colnames(array6.processed)) 
samples <- colnames(array6.processed)

samples_tab$type <- c(rep("rad",12), rep("con",3), rep("rad",12))
samples_tab$time <- c(rep("20",6), rep("30",6), rep("3",3), rep("3",6), rep("5",6))
samples_tab$dose <- c(rep("9.5Gy",6), rep("9.59Gy",6), rep("0Gy",3), rep("2.7Gy",6), rep("4.1Gy",6))
samples_tab[,5] <- unite(samples_tab[,2:4], ID, sep = "_", remove = T)
samples_tab$ID <- factor(samples_tab$ID, levels=c("con_3_0Gy","rad_3_2.7Gy","rad_5_4.1Gy","rad_20_9.5Gy","rad_30_9.59Gy"), labels=c("con_3_0Gy","rad_3_2.7Gy","rad_5_4.1Gy","rad_20_9.5Gy","rad_30_9.59Gy"))
samples_tab
##### FILTER DATA

## Check the rowMeans
array6.medians <- matrixStats::rowMedians(as.matrix(array6.processed))

par(mfrow=c(1,1))
hist(array6.medians)

## Filter those with median less than 4
array6.filtered <- filter(array6.processed, matrixStats::rowMedians(as.matrix(array6.processed)) > 3)

## Removing multiple mappings

anno_grouped <- group_by(array6.filtered, "PROBEID"=rownames(array6.filtered))
anno_summarized <- 
  dplyr::summarize(anno_grouped, no_of_matches = n_distinct(PROBEID))

sum(anno_summarized$no_of_matches > 1)

dim(array6.filtered)

##### VISUALIZE DATA
table(samples_tab$ID)
## define a variable that will hold colors we want to use for plots. Here we can define colors as numbers, or define them with literal color names
my_cols <- c(1:6)
## alternatively:
#my_cols <- c("green","cyan","blue","red","purple")
names(my_cols) <- c("con_3_0Gy","rad_3_2.7Gy","rad_5_4.1Gy","rad_20_9.5Gy","rad_30_9.59Gy") # assign names to your colors. Here names are sample types
print(my_cols)

RG <- backgroundCorrect(array6.filtered, method="normexp", offset=50)
MA <- normalizeWithinArrays(array6.filtered, weights=NULL)

## plot distribution of counts per sample
par(mfrow=c(1,1))
boxplot(array6.filtered, las=2, cex=.5, cex.axis=.5, main="Raw counts, filtered, log") # we make pseudocounts of the raw data

## Visualize our data as MDS plot - multidimensional scaling plot of distances between gene expression profiles. Remember, MDS is not the same as PCA plot
plotMDS(array6.filtered, las=2, cex=1, col=my_cols[samples_tab$ID])


##### EDGER OBJECT 

## Here we set up an EdgeR object that contains expression data and sample description. We group samples per sample type (CTRL, PT, etc.)
dgeFull <- DGEList(array6.filtered, group=samples_tab$ID)
dgeFull$samples$Individual <- samples_tab$sample
dgeFull$samples

## Here we want to define the order of our sample groups - otherwise they are ordered alphabetically
dgeFull$samples$group <- factor(dgeFull$samples$group, levels=c("con_3_0Gy","rad_3_2.7Gy","rad_5_4.1Gy","rad_20_9.5Gy","rad_30_9.59Gy"))

## PCA

## Here we run PCA analysis to yet again visualize our data, but also to check how it groups based on the biological type. We can also identify here any potential problems, outliers etc.
## Remember - for PCA we need a transposed table of normalized counts - genes as columns and samples as rows
norm_t <- t(array6.filtered)
#norm_t[1:5,1:5]
my.pca <- PCA(norm_t, graph=FALSE)

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
colnames(designMat) <- c("con_3_0Gy","rad_3_2.7Gy","rad_5_4.1Gy","rad_20_9.5Gy","rad_30_9.59Gy")
rownames(designMat) <- rownames(dgeFull$samples)
print(designMat)


contMatrix <- makeContrasts(con_3_0Gy-rad_3_2.7Gy, con_3_0Gy-rad_5_4.1Gy, con_3_0Gy-rad_20_9.5Gy, con_3_0Gy-rad_30_9.59Gy,
                            rad_3_2.7Gy-rad_5_4.1Gy, rad_3_2.7Gy-rad_20_9.5Gy, rad_3_2.7Gy-rad_30_9.59Gy,
                            levels=designMat)
print(contMatrix)
print(designMat)


fit <- lmFit(array6.filtered, designMat)
fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)

topTable(fit2, coef=1, adjust="BH")
results <- decideTests(fit2)

topTable(fit2, number=10, adjust="BH", coef=c(1:5))


array.results <- list("con_3_0Gy-rad_3_2.7Gy"=topTable(fit2, coef=1, adjust="BH", number=8267),
                      "con_3_0Gy-rad_5_4.1Gy"=topTable(fit2, coef=2, adjust="BH", number=8267),
                      "con_3_0Gy-rad_20_9.5Gy"=topTable(fit2, coef=3, adjust="BH", number=8267),
                      "con_3_0Gy-rad_30_9.59Gy"=topTable(fit2, coef=4, adjust="BH", number=8267),
                      "rad_3_2.7Gy-rad_5_4.1Gy"=topTable(fit2, coef=5, adjust="BH", number=8267),
                      "rad_3_2.7Gy-rad_20_9.5Gy"=topTable(fit2, coef=6, adjust="BH", number=8267),
                      "rad_3_2.7Gy-rad_30_9.59Gy"=topTable(fit2, coef=7, adjust="BH", number=8267))

hist(array.results$`con_3_0Gy-rad_3_2.7Gy`$adj.P.Val)
sum(array.results$`con_3_0Gy-rad_3_2.7Gy`$adj.P.Val <= 0.05)

##### CONVERT 
array.results <- lapply(lapply(lapply(lapply(lapply(array.results, function(x) {x$ENTREZID <- rownames(x);x}), function(x) left_join(x, annotTable[,c(3,2)], by="ENTREZID")), function(x) distinct(x)),
                               function(x) filter(x, !SYMBOL == "NA")), function(x) {x[,9] <- unite(x[,c(7:8)], ID, sep="_", remove=T);x})

array.sel.genes <- filter(filter(array.results$`con_3_0Gy-rad_3_2.7Gy`, adj.P.Val <=0.05), !SYMBOL == "NA")

# Remove duplicates by PROBEID and SYMBOL leaving those with the lowest P Val
array.sel.genes.ordered <- array.sel.genes[order(array.sel.genes$adj.P.Val, decreasing = FALSE),]
array.sel.genes.filtered <- array.sel.genes.ordered[!duplicated(array.sel.genes.ordered$ENTREZID),]
array.sel.genes.filtered <- array.sel.genes.filtered[!duplicated(array.sel.genes.filtered$SYMBOL),]


array.heart <- bind_cols(lapply(lapply(array.results, function(x) left_join(array.sel.genes.filtered[,c(1,9)], x, by="ID")), function(x) x[,c(3,7,10)]))[,c(1,2,4,7,10,21)]
head(array.heart)
array.heart.ordered <- array.heart[order(array.heart$adj.P.Val...2, decreasing = FALSE),]
array.heart.filtered <- array.heart.ordered[!duplicated(array.heart.ordered$SYMBOL...21),]
array.heart.filtered <- array.heart.filtered[complete.cases(array.heart.filtered),][,-2]

colnames(array.heart.filtered) <- c("logFC_3","logFC_5","logFC_20","logFC_30", "Gene")
rownames(array.heart.filtered) <- array.heart.filtered$Gene
array.heart.filtered <- array.heart.filtered[,-5]

colnames(array.heart.filtered) <- c("3", "5", "20", "30")

#Figure_S5G
ggplot(melt(filter(array.heart.filtered, array.heart.filtered[,1] > 0)), aes(factor(variable), value)) +
  geom_boxplot(fill = "lightblue", alpha = 0.75, width=0.5) +
  geom_hline(yintercept=0, linetype="dashed", color = "blue", size=1, alpha = 0.5) +
  theme_classic() +
  ylim(-2,2) +
  labs(title="Up-regulated genes\nin blood (Cs 2.7 Gy)", x="Time after exposure (days)", y="logFC", color="Sample") +
  theme(
    axis.title.x = element_text(size = 8, face="bold"),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    axis.title.y = element_text(size = 8, face="bold"),
    title = element_text(size = 8, face = "bold"))

#Figure_S5H
ggplot(melt(filter(array.heart.filtered, array.heart.filtered[,1] < 0)), aes(factor(variable), value)) +
  geom_boxplot(fill = "lightblue", alpha = 0.75, width=0.5) +
  geom_hline(yintercept=0, linetype="dashed", color = "blue", size=1, alpha = 0.5) +
  theme_classic() +
  ylim(-2,2) +
  labs(title="Down-regulated genes\nin blood (Cs 2.7 Gy)", x="Time after exposure (days)", y="logFC", color="Sample") +
  theme(
    axis.title.x = element_text(size = 8, face="bold"),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    axis.title.y = element_text(size = 8, face="bold"),
    title = element_text(size = 8, face = "bold"))

#Trend
anova(lm(colMedians(as.matrix(filter(array.heart.filtered, array.heart.filtered[,1] < 0)))~as.vector(1:4)))[[5]][1]
anova(lm(colMedians(as.matrix(filter(array.heart.filtered, array.heart.filtered[,1] > 0)))~as.vector(1:4)))[[5]][1]

