# 1. Install and load packages ####
library(minfi)
library(limma)
library(data.table)
library(lumi)
library(reshape2)
library(data.table)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(ggplot2)
library(dplyr)
library(ggvenn)
library(karyoploteR)
library(BSgenome.Hsapiens.UCSC.hg19)
library(tibble)
library(readr)
library(kableExtra)
library(rpart)
library(readxl)
library(UpSetR)
library(DMRcate)
library(ExperimentHub)
library(rtracklayer)
library(RColorBrewer)
library(missMethyl)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(bumphunter)
library(rlist)
library(plyr)
library(dplyr)
library(GOSemSim)
library(clusterProfiler)
library(pathview)
library(org.Hs.eg.db)
library(gridExtra)
library(datasets)
library(enrichplot)
library(ggupset) 
library(ggnewscale)
library(annotatr)
library(XML)
library(tidyr)
library(VennDiagram)
library(viridis)
library(GMCM)
library(Kendall)
library(plyranges)
library(mcreplicate)
library(gtools)
library(factoextra)
library(plyranges)

# 2. Import and prepare primary files ####
wpath<- "/mnt/work1/adrian/Nasa/"
setwd(wpath)

#Primary files
data.raw <- data.table::fread("GSE108187_processed_matrix.csv", sep=",", header=TRUE)
data.raw.colnames <- read_excel("GSE108187_colnames.xlsx", sheet="Wszystkie")

cc <- grep("Pval",colnames(data.raw), value = T) # get names of unwanted columns
data.raw <- data.raw[, !cc, with=FALSE]
data.raw <- as.data.frame(data.raw)
rownames(data.raw) <- data.raw$CpG.id # assign row names
data.raw[,c(1:2)] <- NULL # get rid of first two columns
colnames(data.raw) <- paste(data.raw.colnames$Colnames, data.raw.colnames$Title, sep="_")

#Annotations
ann450k<-getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

ucscgenes <- annotateTranscripts(TxDb.Hsapiens.UCSC.hg19.knownGene)

#DMPs.annotations.all.probes.2000 <- matchGenes(ann450k, ucscgenes, promoterDist = 2000) # It takes several hours.

DMPs.annotations.all.probes.2000 <- readRDS("DMPs.annotations.all.probes.2000.rdata")

# 3. Differential methylation analysis ####

# PCA Analysis
m.pca <- lumi::beta2m(data.raw[,c(8:9,48:49,73:75)])
colnames(m.pca) <- c("Fe_1Gy_R1", "Fe_1Gy_R2", "X_1Gy_R1", "X_1Gy_R2", "Si_1Gy_R1", "Si_1Gy_R2", "Si_1Gy_R3")
rownames(m.pca) <- NULL
res.pca <- prcomp(t(na.omit(m.pca)), scale = TRUE)

fviz_eig(res.pca)

#Figure_1B
fviz_pca_ind(res.pca, repel = TRUE)

# Control samples comparison
#Figure_S1A
boxplot(list(Fe=rowMeans(lumi::beta2m(data.raw[,1:3])), X=rowMeans(lumi::beta2m(data.raw[,45:47])), Si=rowMeans(lumi::beta2m(data.raw[,67:69]))), xlab="Radiation Type", ylab="M Value",
        main="Mean Values for 0 Gy samples",pty="s") 

#Pick Samples
f0 <- lumi::beta2m(data.raw[,c(1:3)])
f0.1 <- lumi::beta2m(data.raw[,c(4:5)])
f0.3 <- lumi::beta2m(data.raw[,c(6:7)])
f1 <- lumi::beta2m(data.raw[,c(8:9)])

s0 <- lumi::beta2m(data.raw[,c(67:69)])
s0.3 <- lumi::beta2m(data.raw[,c(70:72)])
s1 <- lumi::beta2m(data.raw[,c(73:75)])

x0 <- lumi::beta2m(data.raw[,c(45:47)])
x1 <- lumi::beta2m(data.raw[,c(48:50)])

#Perform differential methylation analysis for comparisons of interest
get.dmps <- function(a, b){
  x <- na.omit(merge(a, b, by=0))
  rownames(x) <- x$Row.names
  x <- as.matrix(x[,-1])
  
  x <- x[which(apply(as.matrix(x),1,sd)>0),]
  
  groups <- c(rep(as.character(substitute(a)), ncol(a)), rep(as.character(substitute(b)), ncol(b)))
  cbind("sample"=colnames(x), "cat"=as.character(groups))
  
  y <- dmpFinder(as.matrix(x)[,c(which(groups==as.character(substitute(a))), which(groups==as.character(substitute(b))))],
                 as.character(groups[which(groups %in% c(as.character(substitute(a)),as.character(substitute(b))))]),
                 type="categorical")
  
  avDiff <- rowMeans(x[,groups==as.character(substitute(a))], na.rm=TRUE) - rowMeans(x[,groups==as.character(substitute(b))], na.rm=TRUE)
  
  rM1 <- rowMeans(x[,groups==as.character(substitute(a))], na.rm=TRUE)
  
  rM2 <- rowMeans(x[,groups==as.character(substitute(b))], na.rm=TRUE)
  
  Probes <- as.character(names(rM1))
  x.all <- cbind(avDiff, rM1, rM2)
  
  x.all <- data.frame("Probes"=Probes, "avDiff"=avDiff, "rM1"=rM1, "rM2"=rM2)
  
  y$Probes <- rownames(y)
  
  x.all <- left_join(y, x.all, by="Probes")
  
  rownames(x.all) <- c(x.all$Probes)
  
  x.all <- x.all[,-5]
  
  colnames(x.all)[ncol(x.all)-1:0] <- c(as.character(substitute(a)), as.character(substitute(b)))
  
  x.all
}

f0.3s0.3 <- get.dmps(f0.3, s0.3)
f1s1 <- get.dmps(f1, s1)
f1x1 <- get.dmps(f1, x1)

dmpfinder.results <- list("f0.3s0.3"=f0.3s0.3,
                          "f1s1"=f1s1,
                          "f1x1"=f1x1)

# 4. Pick significant DMPs ####

# I condition
dmpfinder.results.filtered <- lapply(dmpfinder.results, function(x) filter(x, qval < 0.05 & abs(avDiff) > 0.58))

# II condition
getprobes <- function(x) {
  x <- left_join(data.frame(x, "cc"=rownames(x)), data.frame("cc"=rownames(data.raw), "control"=lumi::beta2m(rowMeans(data.raw[,c(1:3,45:47,67:69)]))), by="cc")
  x[,10:11] <- x[,6:7] - x[,9]
  y <- filter(x, abs(x[,10]) >= 0.58)
  z <- filter(x, abs(x[,11]) >= 0.58)
  w <- data.frame("cc"=union(y$cc, z$cc))
  x <- left_join(w, x[,c(8,1:7)], by="cc")
  rownames(x) <- x[,1]
  x[,2:8]
}

dmpfinder.results.filtered <- lapply(dmpfinder.results.filtered, function(x) getprobes(x))

# Add methylation levels
probes.for.dmrs <- lapply(lapply(lapply(lapply(dmpfinder.results.filtered, function(x) {x[,8] <- rownames(x);x}), function(x) { z <- left_join(x, cbind(data.frame("V8"=rownames(data.raw), data.raw)), by="V8");
z[,8:ncol(z)]}), function(x) {rownames(x) <- c(x[,1]);x[,-1]}), function(x) lumi::beta2m(x))

probes.for.dmrs.means <- lapply(lapply(probes.for.dmrs, function(x) data.frame("f0"=rowMeans(x[,1:3]), "f0.3"=rowMeans(x[6:7]), "f1"=rowMeans(x[8:9]), "f0.14"=rowMeans(x[,10:12]), "f0.3.14"=rowMeans(x[,16:18]), "f1.14"=rowMeans(x[,19:21]), "f0.22"=rowMeans(x[,22:24]),
                                                                               "f0.3.22"=rowMeans(x[,28:30]), "f1.22"=rowMeans(x[,31:33]), "f0.53"=rowMeans(x[,34:36]), "f0.3.53"=rowMeans(x[,39:41]), "f1.53"=rowMeans(x[,42:44]), "s0"=rowMeans(x[,67:69]), "s0.3"=rowMeans(x[,70:72]), "s1"=rowMeans(x[,73:75]), "s0.13"=rowMeans(x[,76:78]),
                                                                               "s0.3.13"=rowMeans(x[,79:81]), "s1.13"=rowMeans(x[,82:84]), "s0.21"=rowMeans(x[,85:87]), "s0.3.21"=rowMeans(x[,88:90]), "s1.21"=rowMeans(x[,91:93]), "s0.62"=rowMeans(x[,94:96]), "s0.3.62"=rowMeans(x[,97:99]), "s1.62"=rowMeans(x[,100:102]),
                                                                               "x0"=rowMeans(x[,45:47]), "x1"=rowMeans(x[,48:50]), "x0.13"=rowMeans(x[,51:52]), "x1.13"=rowMeans(x[,53:55]), "x0.22"=rowMeans(x[,56:57]), "x1.22"=rowMeans(x[,58:60]), "x0.62"=rowMeans(x[,61:63]), "x1.62"=rowMeans(x[,64:66]))),
                                function(x) {colnames(x) <- gsub("$", ".dmp", colnames(x));x})


# Add Annotations
match.probes <- function(x){
  ann450k.x<-DMPs.annotations.all.probes.2000[match(rownames(x),DMPs.annotations.all.probes.2000$Name),c(1:4,12:ncol(DMPs.annotations.all.probes.2000))]
  ann450k.x <- subset(ann450k.x, !is.na(rownames(ann450k.x)))
  rownames(ann450k.x) <- make.names(rownames(ann450k.x), unique=TRUE)
  ann450k.x <- ann450k.x[,c(1:54,4)]
  y <- rownames(x)
  y <- data.frame("Name.1"=y,x)
  ann450k.x <- left_join(data.frame(ann450k.x), y, by="Name.1")
  ann450k.x <- ann450k.x[complete.cases(ann450k.x[,1:3]),]
}

matched.dmps <- lapply(lapply(lapply(lapply(dmpfinder.results.filtered, function(x) match.probes(x)), function(x) { x <- x[,c(1,2,2,3:ncol(x))]; x}), function(x) {colnames(x)[2:3] <- c("start", "end"); x}),
                       function(x) makeGRangesFromDataFrame(x, keep.extra.columns = T))

# Calculate means for data.raw 
data.raw.M <- apply(data.raw, 2, lumi::beta2m)

data.raw.M.means <- apply(data.frame("f0"=rowMeans(data.raw[,1:3]), "f0.3"=rowMeans(data.raw[6:7]), "f1"=rowMeans(data.raw[8:9]), "f0.14"=rowMeans(data.raw[,10:12]), "f0.3.14"=rowMeans(data.raw[,16:18]), "f1.14"=rowMeans(data.raw[,19:21]), "f0.22"=rowMeans(data.raw[,22:24]),
                                     "f0.3.22"=rowMeans(data.raw[,28:30]), "f1.22"=rowMeans(data.raw[,31:33]), "f0.53"=rowMeans(data.raw[,34:36]), "f0.3.53"=rowMeans(data.raw[,39:41]), "f1.53"=rowMeans(data.raw[,42:44]), "s0"=rowMeans(data.raw[,67:69]), "s0.3"=rowMeans(data.raw[,70:72]), "s1"=rowMeans(data.raw[,73:75]), "s0.13"=rowMeans(data.raw[,76:78]),
                                     "s0.3.13"=rowMeans(data.raw[,79:81]), "s1.13"=rowMeans(data.raw[,82:84]), "s0.21"=rowMeans(data.raw[,85:87]), "s0.3.21"=rowMeans(data.raw[,88:90]), "s1.21"=rowMeans(data.raw[,91:93]), "s0.62"=rowMeans(data.raw[,94:96]), "s0.3.62"=rowMeans(data.raw[,97:99]), "s1.62"=rowMeans(data.raw[,100:102]),
                                     "x0"=rowMeans(data.raw[,45:47]), "x1"=rowMeans(data.raw[,48:50]), "x0.13"=rowMeans(data.raw[,51:52]), "x1.13"=rowMeans(data.raw[,53:55]), "x0.22"=rowMeans(data.raw[,56:57]), "x1.22"=rowMeans(data.raw[,58:60]), "x0.62"=rowMeans(data.raw[,61:63]), "x1.62"=rowMeans(data.raw[,64:66])),
                          2, lumi::beta2m)

# 5. DMRs Building ####

#Building showed for FeSi 1Gy Comparison
f1s1 <- na.omit(probes.for.dmrs$f1s1[,c(8:9, 73:75)])

samps.f1s1 <- list("f1"= colnames(f1s1[,1:2]),
                   "s1"= colnames(f1s1[,3:5]))

design.f1s1 <- matrix(0L, nrow=ncol(f1s1), ncol=length(names(samps.f1s1)))
colnames(design.f1s1) <- names(samps.f1s1)
rownames(design.f1s1) <- colnames(f1s1)
for (e in names(samps.f1s1)){ design.f1s1[samps.f1s1[[e]],e] <- 1 }

contMatrix.f1s1 <- makeContrasts(f1 - s1,
                                 levels=design.f1s1)

myAnnotation.f1s1 <- cpg.annotate(object=as.matrix(f1s1), datatype="array", what="M", analysis.type="differential", design=design.f1s1, contrasts=TRUE, cont.matrix=contMatrix.f1s1,
                                  coef="f1 - s1", annotation=c(array = "IlluminaHumanMethylationEPIC", annotation = "ilm10b4.hg19"))
DMRs.f1s1 <- dmrcate(myAnnotation.f1s1, lambda=1000, C=2)
DMRs.f1s1 <- extractRanges(DMRs.f1s1, genome="hg19")

overlapDMPs <- function(x,y) { #x=GRanges cords; y=matched.dmps$x 
  z <- findOverlaps(x, y)
  
  qHits <- x[queryHits(z)]
  sHits <- y[subjectHits(z)]
  
  qHits$ID <- paste0('E_', 1:length(qHits))
  sHits$ID <- paste0('E_', 1:length(sHits))
  
  qHits <- as.data.frame(qHits)
  sHits <- as.data.frame(sHits)
  sHits <- sHits[,c(6,1:5,7:ncol(sHits))]
  
  a <- left_join(qHits, sHits, by="ID")
  a
}

#Figure_S1I
lapply(Map(overlapDMPs, list("DMRs.f0.3s0.3"=DMRs.f0.3s0.3,
                                  "DMRs.f1s1"=DMRs.f1s1,
                                  "DMRs.f1x1"=DMRs.f1x1), matched.dmps[c(2,3,14)]), function(x) nrow(x))



#Run above code with different lambda value
pct.dmpsINdmrs <- data.frame("Lambda"=c(250,500,1000,2000,3000,4000,5000,6000,7000),
                             "FeSi_0.3Gy"=c(6998,9019,10852,12462,13572,14423,15012,15481,15880),
                             "FeSi_1Gy"=c(10005,12771,14935,16766,17862,18771,19392,19879,20317),
                             "FeX_1Gy"=c(7077,9015,10690,12058,12892,13574,14050,14435,14793))
pct.dmpsINdmrs[,2:4] <- pct.dmpsINdmrs[,2:4]/data.frame(rep(26751, 9), rep(31818,9), rep(24807,9))

ggplot(melt(pct.dmpsINdmrs, id.vars=c(1)), aes(x=Lambda, y=value, colour=variable)) +
  geom_line() +
  geom_vline(xintercept = 1000, linetype="dotted", 
             color = "blue", size=1.5) +
  labs(x="Base Pairs", y="%DMPs in DMRs", colour="") +
  ggtitle("") +
  theme_classic() + 
  theme(
    axis.title.x = element_text(size = 15, face="bold"),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    axis.title.y = element_text(size = 15, face="bold"),
    title = element_text(size = 15, face = "bold"),
    legend.text = element_text(size=15))


# 6. Simulate DMRs building with rows permutation ####

overlapDMPs.test <- function(x,y) { #x=GRanges cords; y=matched.dmps$x 
  z <- findOverlaps(x, y)
  
  qHits <- x[queryHits(z)]
  sHits <- y[subjectHits(z)]
  
  qHits$ID <- paste0('E_', 1:length(qHits))
  sHits$ID <- paste0('E_', 1:length(sHits))
  
  qHits <- as.data.frame(qHits)
  sHits <- as.data.frame(sHits)
  sHits <- sHits[,c(6,1:5,7:ncol(sHits))]
  
  a <- left_join(qHits, sHits, by="ID")
  a
}

#

getprobes.permuted <- function(x) {
  x <- left_join(data.frame(x, "cc"=rownames(x)), data.frame("cc"=rownames(data.raw.permuted), "control"=rowMeans(data.raw.permuted[,c(1:3,45:47,67:69)])), by="cc")
  x[,10:11] <- x[,6:7] - x[,9]
  y <- filter(x, abs(x[,10]) >= 0.58)
  z <- filter(x, abs(x[,11]) >= 0.58)
  w <- data.frame("cc"=union(y$cc, z$cc))
  x <- left_join(w, x[,c(8,1:7)], by="cc")
  rownames(x) <- x[,1]
  x[,2:8]
}

#

random.DM.analysis <- function(x) {
  data.raw.permuted <- x
  rownames(data.raw.permuted) <- permute(rownames(data.raw.permuted))
  data.raw.permuted <- apply(data.raw.permuted, 2, lumi::beta2m)
  
  #Pick Samples
  f0 <- data.raw.permuted[,c(1:3)]
  f0.1 <- data.raw.permuted[,c(4:5)]
  f0.3 <- data.raw.permuted[,c(6:7)]
  f1 <- data.raw.permuted[,c(8:9)]
  
  s0 <- data.raw.permuted[,c(67:69)]
  s0.3 <- data.raw.permuted[,c(70:72)]
  s1 <- data.raw.permuted[,c(73:75)]
  
  x0 <- data.raw.permuted[,c(45:47)]
  x1 <- data.raw.permuted[,c(48:50)]
  
  #DM analysis
  f0.3s0.3 <- get.dmps(f0.3, s0.3)
  f1s1 <- get.dmps(f1, s1)
  f1x1 <- get.dmps(f1, x1)
  
  #Pick DMPs for DMRs
  permuted <- list("f0.3s0.3"=f0.3s0.3,
                   "f1s1"=f1s1,
                   "f1x1"=f1x1)
  # I condition 
  permuted.filtered <- lapply(permuted, function(x) filter(x, qval < 0.05 & abs(avDiff) > 0.58))
  
  # II condition
  permuted.filtered <- lapply(permuted.filtered, function(x) getprobes.permuted(x))
  
  # Add methylation levels
  probes.for.dmrs.permuted <- lapply(lapply(lapply(permuted.filtered, function(x) {x[,8] <- rownames(x);x}), function(x) { z <- left_join(x, cbind(data.frame("V8"=rownames(data.raw.permuted), data.raw.permuted)), by="V8");
  z[,8:ncol(z)]}), function(x) {rownames(x) <- c(x[,1]);x[,-1]})
  
  matched.dmps.permuted <- lapply(lapply(lapply(lapply(permuted.filtered, function(x) match.probes(x)), function(x) { x <- x[,c(1,2,2,3:ncol(x))]; x}), function(x) {colnames(x)[2:3] <- c("start", "end"); x}),
                                  function(x) makeGRangesFromDataFrame(x, keep.extra.columns = T))
  
  #Build DMRs
  f0.3s0.3 <- na.omit(probes.for.dmrs.permuted$f0.3s0.3[,c(6:7,70:72)])
  f1s1 <- na.omit(probes.for.dmrs.permuted$f1s1[,c(8:9, 73:75)])
  
  samps.f0.3s0.3 <- list("f0.3"= colnames(f0.3s0.3[,1:2]),
                         "s0.3"= colnames(f0.3s0.3[,3:5]))
  
  samps.f1s1 <- list("f1"= colnames(f1s1[,1:2]),
                     "s1"= colnames(f1s1[,3:5]))
  
  
  design.f0.3s0.3 <- matrix(0L, nrow=ncol(f0.3s0.3), ncol=length(names(samps.f0.3s0.3)))
  colnames(design.f0.3s0.3) <- names(samps.f0.3s0.3)
  rownames(design.f0.3s0.3) <- colnames(f0.3s0.3)
  for (e in names(samps.f0.3s0.3)){ design.f0.3s0.3[samps.f0.3s0.3[[e]],e] <- 1 }
  
  contMatrix.f0.3s0.3 <- makeContrasts(f0.3 - s0.3,
                                       levels=design.f0.3s0.3)
  
  design.f1s1 <- matrix(0L, nrow=ncol(f1s1), ncol=length(names(samps.f1s1)))
  colnames(design.f1s1) <- names(samps.f1s1)
  rownames(design.f1s1) <- colnames(f1s1)
  for (e in names(samps.f1s1)){ design.f1s1[samps.f1s1[[e]],e] <- 1 }
  
  contMatrix.f1s1 <- makeContrasts(f1 - s1,
                                   levels=design.f1s1)
  
  
  f1x1 <- na.omit(probes.for.dmrs.permuted$f1x1[,c(8:9,48:50)])
  
  samps.f1x1 <- list("f1"= colnames(f1x1[,1:2]),
                     "x1"= colnames(f1x1[,3:5]))
  
  
  design.f1x1 <- matrix(0L, nrow=ncol(f1x1), ncol=length(names(samps.f1x1)))
  colnames(design.f1x1) <- names(samps.f1x1)
  rownames(design.f1x1) <- colnames(f1x1)
  for (e in names(samps.f1x1)){ design.f1x1[samps.f1x1[[e]],e] <- 1 }
  
  contMatrix.f1x1 <- makeContrasts(f1 - x1,
                                   levels=design.f1x1)
  
  
  myAnnotation.f0.3s0.3 <- cpg.annotate(object=as.matrix(f0.3s0.3), datatype="array", what="M", analysis.type="differential", design=design.f0.3s0.3, contrasts=TRUE, cont.matrix=contMatrix.f0.3s0.3,
                                        coef="f0.3 - s0.3", annotation=c(array = "IlluminaHumanMethylationEPIC", annotation = "ilm10b4.hg19"))
  DMRs.f0.3s0.3 <- dmrcate(myAnnotation.f0.3s0.3, lambda=1000, C=2)
  DMRs.f0.3s0.3 <- extractRanges(DMRs.f0.3s0.3, genome="hg19")
  
  myAnnotation.f1s1 <- cpg.annotate(object=as.matrix(f1s1), datatype="array", what="M", analysis.type="differential", design=design.f1s1, contrasts=TRUE, cont.matrix=contMatrix.f1s1,
                                    coef="f1 - s1", annotation=c(array = "IlluminaHumanMethylationEPIC", annotation = "ilm10b4.hg19"))
  DMRs.f1s1 <- dmrcate(myAnnotation.f1s1, lambda=1000, C=2)
  DMRs.f1s1 <- extractRanges(DMRs.f1s1, genome="hg19")
  
  myAnnotation.f1x1 <- cpg.annotate(object=as.matrix(f1x1), datatype="array", what="M", analysis.type="differential", design=design.f1x1, contrasts=TRUE, cont.matrix=contMatrix.f1x1,
                                    coef="f1 - x1", annotation=c(array = "IlluminaHumanMethylationEPIC", annotation = "ilm10b4.hg19"))
  DMRs.f1x1 <- dmrcate(myAnnotation.f1x1, lambda=1000, C=2)
  DMRs.f1x1 <- extractRanges(DMRs.f1x1, genome="hg19")
  
  DMRs.permuted <- list("DMRs.f0.3s0.3"=DMRs.f0.3s0.3,
                        "DMRs.f1s1"=DMRs.f1s1,
                        "DMRs.f1x1"=DMRs.f1x1)
  
  t(data.frame("DMPs"=lapply(matched.dmps.permuted, function(x) length(x)),
               "DMRs"=lapply(DMRs.permuted, function(x) length(x)),
               "DMPsIN"=lapply(Map(overlapDMPs.test, DMRs.permuted, matched.dmps.permuted), function(x) nrow(x))))
}

#DM.permuted.1000 <- mc_replicate(1000, list(random.DM.analysis(data.raw)), refresh = 0.1, mc.cores=30) #Adjust to your dataset.

DM.permuted <- data.frame(bind_cols(DM.permuted))
rownames(DM.permuted) <- c("DMPs_FeSi_0.3Gy", "DMPs_FeSi_1Gy", "DMPs_FeX_1Gy", "DMRs_FeSi_0.3Gy", "DMRs_FeSi_1Gy", "DMRs_FeX_1Gy",
                           "DMPsInDMRs_FeSi_0.3Gy", "DMPsInDMRs_FeSi_1Gy", "DMPsInDMRs_FeX_1Gy")
colnames(DM.permuted) <- NULL
DM.permuted <- t(DM.permuted)

table((DM.permuted[,7]/DM.permuted[,1])<10852/26751) #1
table((DM.permuted[,8]/DM.permuted[,2])<14935/31818) #2
table((DM.permuted[,9]/DM.permuted[,3])<10690/24807) #3

# 7. Load and prepare lists of DMRs ####
DMRs.newest.GRanges <- list.load("DMRs.newest.GRanges.rdata")
DMRs.newest.annotations <- list.load("DMRs.newest.annotations.rdata")
DMRs.newest.annotations <- lapply(DMRs.newest.annotations, function(x) {colnames(x) <- gsub("$", ".dmrs", colnames(x));x})

# Bedtools analysis - below you will find command line code
# singularity shell --bind /mnt/work1/:/mnt/work1 --bind /mnt/data1:/mnt/data1 /opt/container/mab-bio-1.simg
# bedtools nuc -fi GRCh37.p13.genome.fa -bed DMRs21.ranges.gff > filename.gff
# https://www.gencodegenes.org/human/release_19.html
Bed.newest.GRanges <- list.load("Bed.newest.GRanges.rdata")
DMRs.newest.GRanges <- lapply(Map(cbind, Map(cbind, lapply(DMRs.newest.GRanges, function(x) as.data.frame(x)), Bed.newest.GRanges), DMRs.newest.annotations), function(x) makeGRangesFromDataFrame(x, keep.extra.columns = T))

# Overlap picked probes with identified DMRs
overlapDMPs <- function(x,y) { #x=GRanges cords; y=matched.dmps$x 
  z <- findOverlaps(x, y)
  
  qHits <- x[queryHits(z)]
  sHits <- y[subjectHits(z)]
  
  qHits$ID <- paste0('E_', 1:length(qHits))
  sHits$ID <- paste0('E_', 1:length(sHits))
  
  qHits <- as.data.frame(qHits)
  sHits <- as.data.frame(sHits)
  sHits <- sHits[,c(6,1:5,7:ncol(sHits))]
  
  a <- left_join(qHits, sHits, by="ID")
  a$cord <- unite(a[,1:3], cord, remove=T, sep="_")
  
  a <- a[,c(1:106, 104, 105)]
  
  colnames(a)[107:108] <- c("mFe", "mSiX")
  
  a
}

DMRs.newest.GRanges.over <- Map(overlapDMPs, DMRs.newest.GRanges, matched.dmps) #Here are all probes which overlap dmrs 

DMRs.newest.GRanges.over <- lapply(mapply(function(x,y) list(left_join(x, y, by="V109")), x=lapply(DMRs.newest.GRanges.over, function(x) {x[,109] <- x$Name;x}),
                                          y=lapply(probes.for.dmrs.means, function(x) cbind(data.frame("V109"=rownames(x), x)))), function(x) {x$nProbes <- c(1);x})

# Group DMRs with assigned dmps methylation values 
DMRs.newest.GRanges.over.grouped <- lapply(DMRs.newest.GRanges.over, function(x)  x <- x %>%
                                             group_by(cord) %>%
                                             mutate(intercept=mean(intercept)) %>%
                                             mutate(f=mean(f)) %>%
                                             mutate(pval=mean(pval)) %>%
                                             mutate(qval=mean(qval)) %>%
                                             mutate(avDiff=mean(avDiff)) %>%
                                             mutate(mFe=mean(mFe)) %>%
                                             mutate(mSiX=mean(mSiX)) %>%
                                             mutate(f0.dmp=mean(f0.dmp)) %>%
                                             mutate(f0.3.dmp=mean(f0.3.dmp)) %>%
                                             mutate(f1.dmp=mean(f1.dmp)) %>%
                                             mutate(f0.14.dmp=mean(f0.14.dmp)) %>%
                                             mutate(f0.3.14.dmp=mean(f0.3.14.dmp)) %>%
                                             mutate(f1.14.dmp=mean(f1.14.dmp)) %>%
                                             mutate(f0.22.dmp=mean(f0.22.dmp)) %>%
                                             mutate(f0.3.22.dmp=mean(f0.3.22.dmp)) %>%    
                                             mutate(f1.22.dmp=mean(f1.22.dmp)) %>%    
                                             mutate(f0.53.dmp=mean(f0.53.dmp)) %>%    
                                             mutate(f0.3.53.dmp=mean(f0.3.53.dmp)) %>%    
                                             mutate(f1.53.dmp=mean(f1.53.dmp)) %>%    
                                             mutate(s0.dmp=mean(s0.dmp)) %>%   
                                             mutate(s0.3.dmp=mean(s0.3.dmp)) %>%    
                                             mutate(s1.dmp=mean(s1.dmp)) %>%    
                                             mutate(s0.13.dmp=mean(s0.13.dmp)) %>% 
                                             mutate(s0.3.13.dmp=mean(s0.3.13.dmp)) %>%
                                             mutate(s1.13.dmp=mean(s1.13.dmp)) %>%
                                             mutate(s0.21.dmp=mean(s0.21.dmp)) %>%
                                             mutate(s0.3.21.dmp=mean(s0.3.21.dmp)) %>%
                                             mutate(s1.21.dmp=mean(s1.21.dmp)) %>%
                                             mutate(s0.62.dmp=mean(s0.62.dmp)) %>%
                                             mutate(s0.3.62.dmp=mean(s0.3.62.dmp)) %>%
                                             mutate(s1.62.dmp=mean(s1.62.dmp)) %>%    
                                             mutate(x0.dmp=mean(x0.dmp)) %>%    
                                             mutate(x1.dmp=mean(x1.dmp)) %>%    
                                             mutate(x0.13.dmp=mean(x0.13.dmp)) %>%    
                                             mutate(x1.13.dmp=mean(x1.13.dmp)) %>%    
                                             mutate(x0.22.dmp=mean(x0.22.dmp)) %>%    
                                             mutate(x1.22.dmp=mean(x1.22.dmp)) %>%    
                                             mutate(x0.62.dmp=mean(x0.62.dmp)) %>%
                                             mutate(x1.62.dmp=mean(x1.62.dmp)) %>%
                                             mutate(nProbes=sum(nProbes)) %>%
                                             distinct(cord, .keep_all=T) %>%
                                             data.frame())




# 8. Acute methylation change versus control between DMPs ####

dmps.methylation.change <- lapply(matched.dmps, function(x) left_join(as.data.frame(x), data.frame("Name"=rownames(data.raw.M.means), data.raw.M.means), by="Name"))

# FeSi 0.3
control.fesi.03 <- (dmps.methylation.change$f0.3s0.3$f0 + dmps.methylation.change$f0.3s0.3$s0 + dmps.methylation.change$f0.3s0.3$x0)/3
fe03 <- dmps.methylation.change$f0.3s0.3$f0.3.x - control.fesi.03
si03 <- dmps.methylation.change$f0.3s0.3$s0.3.x - control.fesi.03

#Figure_S1B
plot(fe03,si03, col=densCols(fe03, si03), pch = 20, xlab="Fe 0.3 Gy - Control 0 Gy", ylab="Si 0.3 Gy - Control 0 Gy", main="Day2 M Value Change\nDMPs versus Control", pty="s", xlim=c(-2,3.2), ylim=c(-2,1.5),
     abline(v = 0, h = 0, col = c("blue", "blue"), lty = c(2, 2), lwd = c(3, 3)))

# FeSi 1 Gy
control.fesi.1 <- (dmps.methylation.change$f1s1$f0 + dmps.methylation.change$f1s1$s0 + dmps.methylation.change$f1s1$x0)/3
fe1 <- dmps.methylation.change$f1s1$f1.x - control.fesi.1
si1 <- dmps.methylation.change$f1s1$s1.x - control.fesi.1

#Figure_S1C
plot(fe1,si1, col=densCols(fe1, si1), pch = 20, xlab="Fe 1 Gy - Control 0 Gy", ylab="Si 1 Gy - Control 0 Gy", main="Day2 M Value Change\nDMPs versus Control", pty="s", xlim=c(-2,3.2), ylim=c(-2,1.5),
     abline(v = 0, h = 0, col = c("blue", "blue"), lty = c(2, 2), lwd = c(3, 3)))

# FeX 1 Gy
control.fex.1 <- (dmps.methylation.change$f1x1$f0 + dmps.methylation.change$f1x1$x0 + dmps.methylation.change$f1x1$s0)/3
fex1 <- dmps.methylation.change$f1x1$f1.x - control.fex.1
x1 <- dmps.methylation.change$f1x1$x1.x - control.fex.1

#Figure_1C
plot(fex1,x1, col=densCols(fex1, x1), pch = 20, xlab="Fe 1 Gy - Control 0 Gy", ylab="X 1 Gy - Control 0 Gy", main="Day2 M Value Change\nDMPs versus Control", pty="s", xlim=c(-2,3.2), ylim=c(-2,1.5),
     abline(v = 0, h = 0, col = c("blue", "blue"), lty = c(2, 2), lwd = c(3, 3)))

# SiX 1 Gy
six1 <- na.omit(merge(data.frame(dmps.methylation.change$f1s1[,c(65,79,77,89,57)]), data.frame(dmps.methylation.change$f1x1[,c(57,65,77,89,90)]), by="Name.1"))
control.six1 <- rowMeans(six1[,c(2,4,5,6,7,8)])

si1 <- six1[,3] - control.six1
x1 <- six1[,9] - control.six1

#Figure_1D
plot(si1,x1, col=densCols(si1, x1), pch = 20, xlab="Si 1 Gy - Control 0 Gy", ylab="X 1 Gy - Control 0 Gy", main="Day2 M Value Change\nDMPs versus Control", pty="s", xlim=c(-2,3.2), ylim=c(-2,1.5),
     abline(v = 0, h = 0, col = c("blue", "blue"), lty = c(2, 2), lwd = c(3, 3)))

# 9. % of DMPs in DMRs ####
pctprobes <- data.frame("nDMPs"=c(26751,31818,24807),
                        "inDMRs"=c(10852,14935,10690))
rownames(pctprobes) <- c("FeSi_0.3Gy","FeSi_1Gy", "FeX_1Gy")

pctprobes$pct <- round(pctprobes[,2]/pctprobes[,1],2)
pctprobes$pct1 <- round((pctprobes[,1]-pctprobes[,2])/pctprobes[,1],2)

#Figure_S1H
ggplot(pctprobes, aes(x=rownames(pctprobes), y=pct, fill=rownames(pctprobes))) +
  geom_bar(stat="identity", width=0.7) + 
  scale_fill_viridis(discrete=T, name="") +
  labs(x="Sample", y="nDMP (%)") +
  ggtitle("Amount of DMPs building DMRs") +
  theme_classic() + 
  theme(
    axis.title.x = element_text(size = 15, face="bold"),
    axis.text.x = element_text(size = 15, angle = -340, vjust = 1, hjust=1),
    axis.text.y = element_text(size = 15),
    axis.title.y = element_text(size = 15, face="bold"),
    legend.position="none",
    title = element_text(size=15, face="bold")) +
  ylim(0,1)

# 10. DMPs which changed at least of 0.58 or -0.58 M Value versus control ####

DMPs.increased.Fe <- lapply(dmps.methylation.change, function(x) filter(x, (x[,63] - rowMeans(x[,c(65,77,89)])) >= 0.58))
DMPs.decreased.Fe <- lapply(dmps.methylation.change, function(x) filter(x, (x[,63] - rowMeans(x[,c(65,77,89)])) <= -0.58))
DMPs.increased.SiX <- lapply(dmps.methylation.change, function(x) filter(x, (x[,64] - rowMeans(x[,c(65,77,89)])) >= 0.58))
DMPs.decreased.SiX <- lapply(dmps.methylation.change, function(x) filter(x, (x[,64] - rowMeans(x[,c(65,77,89)])) <= -0.58))

Common.notdetailed <- append(append(DMPs.increased.Fe,DMPs.decreased.Fe), append(DMPs.increased.SiX, DMPs.decreased.SiX))

names(Common.notdetailed) <- c("↑Fe_0.3Gy vs. control", "↑Fe_1Gy vs. control", "↑Fe_1Gy vs. control", "↓Fe_0.3Gy vs. control", "↓Fe_1Gy vs. control", "↓Fe_1Gy vs. control",
                               "↑Si_0.3Gy vs. control", "↑Si_1Gy vs. control", "↑X_1Gy vs. control", "↓Si_0.3Gy vs. control", "↓Si_1Gy vs. control", "↓X_1Gy vs. control")

#Figure_S1D
upset.all <- lapply(Common.notdetailed[c(1,4,7,10)], function(x) x$Name.1)

upset(fromList(upset.all), nsets=12, order.by = c("freq"))

#Figure_1E
upset.all <- lapply(Common.notdetailed[c(2,5,8,11)], function(x) x$Name.1)

upset(fromList(upset.all), nsets=12, order.by = c("freq"))

#Figure_1D
upset.all <- lapply(Common.notdetailed[c(3,6,9,12)], function(x) x$Name.1)

upset(fromList(upset.all), nsets=12, order.by = c("freq"))

# 11. Detailed methylation patterns ####

# AA) Fe increased versus control and is negative and the control was negative
# AB) Fe increased versus control and is negative and the control was positive
# M) Fe increased versus control and is positive and the control was negative
# N) Fe increased versus control and is positive and the control was positive
# R) Fe decreased versus control and is positive and the control was negative
# S) Fe decreased versus control and is positive and the control was positive
# W) Fe decreased versus control and is negative and the control was negative
# X) Fe decreased versus control and is negative and the control was positive

# AA) Si or X decreased versus control and is positive and the control was negative
# AB) Si or X decreased versus control and is positive and the control was positive
# O) Si or X decreased versus control and is negative and the control was negative
# P) Si or X decreased versus control and is negative and the control was positive
# T) Si or X increased versus control and is negative and the control was negative
# U) Si or X increased versus control and is negative and the control was positive
# Y) Si or X increased versus control and is positive and the control was negative
# Z) Si or X increased versus control and is positive and the control was positive

DMPS.increased.Fe.AA <- lapply(dmps.methylation.change, function(x) filter(x, ((x[,63] - rowMeans(x[,c(65,77,89)])) >= 0.58 & x[,63] < 0) & rowMeans(x[,c(65,77,89)]) < 0))
DMPS.increased.Fe.AB <- lapply(dmps.methylation.change, function(x) filter(x, ((x[,63] - rowMeans(x[,c(65,77,89)])) >= 0.58 & x[,63] < 0) & rowMeans(x[,c(65,77,89)]) > 0))

DMPS.increased.Fe.M <- lapply(dmps.methylation.change, function(x) filter(x, ((x[,63] - rowMeans(x[,c(65,77,89)])) >= 0.58 & x[,63] > 0) & rowMeans(x[,c(65,77,89)]) < 0))
DMPS.increased.Fe.N <- lapply(dmps.methylation.change, function(x) filter(x, ((x[,63] - rowMeans(x[,c(65,77,89)])) >= 0.58 & x[,63] > 0) & rowMeans(x[,c(65,77,89)]) > 0))
DMPS.decreased.Fe.R <- lapply(dmps.methylation.change, function(x) filter(x, ((x[,63] - rowMeans(x[,c(65,77,89)])) <= -0.58 & x[,63] > 0) & rowMeans(x[,c(65,77,89)]) < 0))
DMPS.decreased.Fe.S <- lapply(dmps.methylation.change, function(x) filter(x, ((x[,63] - rowMeans(x[,c(65,77,89)])) <= -0.58 & x[,63] > 0) & rowMeans(x[,c(65,77,89)]) > 0))
DMPS.decreased.Fe.W <- lapply(dmps.methylation.change, function(x) filter(x, ((x[,63] - rowMeans(x[,c(65,77,89)])) <= -0.58 & x[,63] < 0) & rowMeans(x[,c(65,77,89)]) < 0))
DMPS.decreased.Fe.X <- lapply(dmps.methylation.change, function(x) filter(x, ((x[,63] - rowMeans(x[,c(65,77,89)])) <= -0.58 & x[,63] < 0) & rowMeans(x[,c(65,77,89)]) > 0))

DMPS.decreased.SiX.AA <- lapply(dmps.methylation.change, function(x) filter(x, ((x[,64] - rowMeans(x[,c(65,77,89)])) <= -0.58 & x[,64] > 0) & rowMeans(x[,c(65,77,89)]) < 0))
DMPS.decreased.SiX.AB <- lapply(dmps.methylation.change, function(x) filter(x, ((x[,64] - rowMeans(x[,c(65,77,89)])) <= -0.58 & x[,64] > 0) & rowMeans(x[,c(65,77,89)]) > 0))

DMPS.decreased.SiX.O <- lapply(dmps.methylation.change, function(x) filter(x, ((x[,64] - rowMeans(x[,c(65,77,89)])) <= -0.58 & x[,64] < 0) & rowMeans(x[,c(65,77,89)]) < 0))
DMPS.decreased.SiX.P <- lapply(dmps.methylation.change, function(x) filter(x, ((x[,64] - rowMeans(x[,c(65,77,89)])) <= -0.58 & x[,64] < 0) & rowMeans(x[,c(65,77,89)]) > 0))
DMPS.increased.SiX.T <- lapply(dmps.methylation.change, function(x) filter(x, ((x[,64] - rowMeans(x[,c(65,77,89)])) >= 0.58 & x[,64] < 0) & rowMeans(x[,c(65,77,89)]) < 0))
DMPS.increased.SiX.U <- lapply(dmps.methylation.change, function(x) filter(x, ((x[,64] - rowMeans(x[,c(65,77,89)])) >= 0.58 & x[,64] < 0) & rowMeans(x[,c(65,77,89)]) > 0))
DMPS.increased.SiX.Y <- lapply(dmps.methylation.change, function(x) filter(x, ((x[,64] - rowMeans(x[,c(65,77,89)])) >= 0.58 & x[,64] > 0) & rowMeans(x[,c(65,77,89)]) < 0))
DMPS.increased.SiX.Z <- lapply(dmps.methylation.change, function(x) filter(x, ((x[,64] - rowMeans(x[,c(65,77,89)])) >= 0.58 & x[,64] > 0) & rowMeans(x[,c(65,77,89)]) > 0))

nDMPs.change.increase <- data.frame("nDMPS"=c(as.vector(unlist(lapply(DMPS.increased.Fe.M, function(x) nrow(x)))),
                                              as.vector(unlist(lapply(DMPS.increased.Fe.N, function(x) nrow(x)))),
                                              as.vector(unlist(lapply(DMPS.increased.SiX.T, function(x) nrow(x)))),
                                              as.vector(unlist(lapply(DMPS.increased.SiX.U, function(x) nrow(x)))),
                                              as.vector(unlist(lapply(DMPS.increased.Fe.AA, function(x) nrow(x)))),
                                              as.vector(unlist(lapply(DMPS.increased.Fe.AB, function(x) nrow(x)))),
                                              as.vector(unlist(lapply(DMPS.increased.SiX.Y, function(x) nrow(x)))),
                                              as.vector(unlist(lapply(DMPS.increased.SiX.Z, function(x) nrow(x))))),
                                    "radiation"=rep(c(rep("Fe",6),rep(c("Si","Si","X"),2)),2),
                                    "change"=rep("↑",24),
                                    "dose"=rep(c("0.3Gy","1Gy","1Gy"),8),
                                    "control"=rep(c(rep("Hypomethylated\nat Baseline",3), rep("Hypermethylated\nat Baseline",3)),4),
                                    "methylation"=c(rep("Hypermethylated\nPost-Radiation",6), rep("Hypomethylated\nPost-Radiation",6),
                                                    rep("Hypomethylated\nPost-Radiation",6), rep("Hypermethylated\nPost-Radiation",6)))
                                              
nDMPS.change.decrease <- data.frame("nDMPS"=c(as.vector(unlist(lapply(DMPS.decreased.Fe.W, function(x) nrow(x)))),
                                          as.vector(unlist(lapply(DMPS.decreased.Fe.X, function(x) nrow(x)))),
                                          as.vector(unlist(lapply(DMPS.decreased.SiX.AA, function(x) nrow(x)))),
                                          as.vector(unlist(lapply(DMPS.decreased.SiX.AB, function(x) nrow(x)))),
                                          as.vector(unlist(lapply(DMPS.decreased.Fe.R, function(x) nrow(x)))),
                                          as.vector(unlist(lapply(DMPS.decreased.Fe.S, function(x) nrow(x)))),
                                          as.vector(unlist(lapply(DMPS.decreased.SiX.O, function(x) nrow(x)))),
                                          as.vector(unlist(lapply(DMPS.decreased.SiX.P, function(x) nrow(x))))),
                                "radiation"=rep(c(rep("Fe",6),rep(c("Si","Si","X"),2)),2),
                                "change"=rep("↓",24),
                                "dose"=rep(c("0.3Gy","1Gy","1Gy"),8),
                                "control"=rep(c(rep("Hypomethylated\nat Baseline",3), rep("Hypermethylated\nat Baseline",3)),4),
                                "methylation"=c(rep("Hypomethylated\nPost-Radiation",6), rep("Hypermethylated\nPost-Radiation",6),
                                                rep("Hypermethylated\nPost-Radiation",6), rep("Hypomethylated\nPost-Radiation",6)))

nDMPs.change.increase[,7] <- unite(nDMPs.change.increase[,3:2], fac, sep="", remove = F)
nDMPS.change.decrease[,7] <- unite(nDMPS.change.decrease[,3:2], fac, sep="", remove = F)

nDMPs.change.increase <- rbind(nDMPs.change.increase, nDMPs.change.increase[c(9,12,21,24),])
nDMPs.change.increase[25:28,1] <- 0
nDMPs.change.increase[25:28,4] <- "0.3Gy"

nDMPS.change.decrease <- rbind(nDMPS.change.decrease, nDMPS.change.decrease[c(9,12,21,24),])
nDMPS.change.decrease[25:28,1] <- 0
nDMPS.change.decrease[25:28,4] <- "0.3Gy"

#Figure_S1E
ggplot(nDMPs.change.increase, aes(fill=dose, y=nDMPS, x=fac)) + 
  geom_bar(position="dodge", stat="identity", width=0.7) +
  scale_fill_viridis(discrete=T, name="Dose") +
  labs(x="", y="nDMP", fill="") +
  ggtitle("") +
  theme_classic() + 
  theme(
    axis.title.x = element_text(size = 15, face="bold"),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    axis.title.y = element_text(size = 15, face="bold"),
    title = element_text(size = 15, face = "bold"),
    legend.text = element_text(size=15)) +
  facet_grid(cols=vars(methylation), rows=vars(control))

#Figure_S1F
ggplot(nDMPS.change.decrease, aes(fill=dose, y=nDMPS, x=fac)) + 
  geom_bar(position="dodge", stat="identity", width=0.7) +
  scale_fill_viridis(discrete=T, name="Dose") +
  labs(x="", y="nDMP", fill="") +
  ggtitle("") +
  theme_classic() + 
  theme(
    axis.title.x = element_text(size = 15, face="bold"),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    axis.title.y = element_text(size = 15, face="bold"),
    title = element_text(size = 15, face = "bold"),
    legend.text = element_text(size=15)) +
  facet_grid(cols=vars(methylation), rows=vars(control))


# 12. Select DMRs, DMPs in or out of DMRs and all DMPs for further analysis ####

# DMPs outside of DMRs
cc <- lapply(DMRs.newest.GRanges.over, function(x) data.frame("Name"=x$Name.1))

probes.notINdmrs <- mapply(function(x,y) list(anti_join(x, y, by="Name")), x=lapply(probes.for.dmrs.means, function(x) {x$Name <- rownames(x);x}), y=cc)
probes.notINdmrs <- mapply(function(x,y) list(left_join(x, y, by="Name")), x=probes.notINdmrs, y=lapply(dmpfinder.results.filtered, function(x) {x$Name <- rownames(x);x[,c(8,1:7)]}))

annote.probes.notINdmrs <- lapply(probes.notINdmrs, function(x) {y <- DMPs.annotations.all.probes.2000[match(x$Name, 
                                                                                                             DMPs.annotations.all.probes.2000$Name), c(1:ncol(DMPs.annotations.all.probes.2000))]; cbind(y, x)})
annote.probes.notINdmrs <- lapply(annote.probes.notINdmrs, function(x) {colnames(x)[94] <- c("Name.1");x})

# DMRs
DMRs.increased.Fe.C <- lapply(DMRs.newest.GRanges.over.grouped, function(x) filter(x, (mFe - rowMeans(x[,c(110,122,134)])) >= 0.58))
DMRs.decreased.Fe.C <- lapply(DMRs.newest.GRanges.over.grouped, function(x) filter(x, (mFe - rowMeans(x[,c(110,122,134)])) <= -0.58))
DMRs.increased.SiX.D <- lapply(DMRs.newest.GRanges.over.grouped, function(x) filter(x, (mSiX - rowMeans(x[,c(110,122,134)])) >= 0.58))
DMRs.decreased.SiX.D <- lapply(DMRs.newest.GRanges.over.grouped, function(x) filter(x, (mSiX - rowMeans(x[,c(110,122,134)])) <= -0.58))

# DMPs in DMRs
DMPs.increased.Fe.C <- lapply(DMRs.newest.GRanges.over, function(x) filter(x, (mFe - rowMeans(x[,c(110,122,134)])) >= 0.58))
DMPs.decreased.Fe.C <- lapply(DMRs.newest.GRanges.over, function(x) filter(x, (mFe - rowMeans(x[,c(110,122,134)])) <= -0.58))
DMPs.increased.SiX.D <- lapply(DMRs.newest.GRanges.over, function(x) filter(x, (mSiX - rowMeans(x[,c(110,122,134)])) >= 0.58))
DMPs.decreased.SiX.D <- lapply(DMRs.newest.GRanges.over, function(x) filter(x, (mSiX - rowMeans(x[,c(110,122,134)])) <= -0.58))

# DMPs outside DMRs
notINdmrs.increased.Fe <- lapply(annote.probes.notINdmrs, function(x) filter(x, ((x[,100] - rowMeans(x[,c(62,74,86)])) >= 0.58)))
notINdmrs.decreased.Fe <- lapply(annote.probes.notINdmrs, function(x) filter(x, ((x[,100] - rowMeans(x[,c(62,74,86)])) <= -0.58)))
notINdmrs.increased.SiX <- lapply(annote.probes.notINdmrs, function(x) filter(x, ((x[,101] - rowMeans(x[,c(62,74,86)])) >= 0.58)))
notINdmrs.decreased.SiX <- lapply(annote.probes.notINdmrs, function(x) filter(x, ((x[,101] - rowMeans(x[,c(62,74,86)])) <= -0.58)))

All.samples.not.detailed.1 <- lapply(append(append(Map(rbind, DMRs.increased.Fe.C, DMRs.decreased.Fe.C), Map(rbind, DMRs.increased.SiX.D, DMRs.decreased.SiX.D)),
                                            append(Map(rbind, DMPs.increased.Fe.C, DMPs.decreased.Fe.C), Map(rbind, DMPs.increased.SiX.D, DMPs.decreased.SiX.D))), function(x)
                                            {colnames(x)[1:5] <- gsub("\\.x", "", colnames(x)[1:5]);x})
names(All.samples.not.detailed.1) <- c("DMRs_Fe_0.3Gy", "DMRs_Fe_1Gy", "DMRs_FeX_1Gy", "DMRs_Si_0.3Gy", "DMRs_Si_1Gy", "DMRs_X_1Gy",
                                       "DMPs_Fe_0.3Gy", "DMPs_Fe_1Gy", "DMPs_FeX_1Gy", "DMPs_Si_0.3Gy", "DMPs_Si_1Gy", "DMPs_X_1Gy")

All.samples.not.detailed.1 <- lapply(All.samples.not.detailed.1, function(x) makeGRangesFromDataFrame(x, keep.extra.columns = T))


All.samples.not.detailed.2 <- append(Map(rbind, notINdmrs.increased.Fe, notINdmrs.decreased.Fe), Map(rbind, notINdmrs.increased.SiX, notINdmrs.decreased.SiX))

names(All.samples.not.detailed.2) <- c("notINdmrs_Fe_0.3Gy", "notINdmrs_Fe_1Gy", "notINdmrs_FeX_1Gy", "notINdmrs_Si_0.3Gy", "notINdmrs_Si_1Gy", "notINdmrs_X_1Gy")

All.samples.not.detailed.2 <- lapply(lapply(lapply(lapply(All.samples.not.detailed.2, function(x) x[,c(1,2,2,3:ncol(x))]), function(x) {colnames(x)[2:3] <- c("start", "end");x}), function(x) {x[!is.na(x$start),]}), 
                                     function(x) makeGRangesFromDataFrame(x, keep.extra.columns = T))


# 13. Annotate data.raw and prepare for further normalization ####

# DMRs annotations to Regulatory Feature Groups, Relation to Islands, Genic Regions 
annote.dmps <- DMPs.annotations.all.probes.2000[match(rownames(data.raw), DMPs.annotations.all.probes.2000$Name), c(1:ncol(DMPs.annotations.all.probes.2000))]

# Control - annote.dmps THESE ARE NOT DMPS, BUT DATA.RAW PROBES - REMEMBER NOT TO MISINTERPRETE 
sum.reg <- data.frame("None"=nrow(filter(annote.dmps, Regulatory_Feature_Group == "")),
                      "Gene"=nrow(filter(annote.dmps, Regulatory_Feature_Group == "Gene_Associated" | Regulatory_Feature_Group == "Gene_Associated_Cell_type_specific")),
                      "NonGene"=nrow(filter(annote.dmps, Regulatory_Feature_Group == "NonGene_Associated" | Regulatory_Feature_Group == "NonGene_Associated_Cell_type_specific")),
                      "Promoter"=nrow(filter(annote.dmps, Regulatory_Feature_Group == "Promoter_Associated" | Regulatory_Feature_Group == "Promoter_Associated_Cell_type_specific")),
                      "Unclass"=nrow(filter(annote.dmps, Regulatory_Feature_Group == "Unclassified" | Regulatory_Feature_Group == "Unclassified_Cell_type_specific")))
sum.reg <- rbind(rbind(sum.reg / sum(sum.reg), sum.reg / sum(sum.reg)), sum.reg / sum(sum.reg))

sum.islands <- data.frame("Island"=nrow(filter(annote.dmps, Relation_to_Island == "Island")),
                          "N_Shelf"=nrow(filter(annote.dmps, Relation_to_Island == "N_Shelf")),
                          "N_Shore"=nrow(filter(annote.dmps, Relation_to_Island == "N_Shore")),
                          "OpenSea"=nrow(filter(annote.dmps, Relation_to_Island == "OpenSea")),
                          "S_Shelf"=nrow(filter(annote.dmps, Relation_to_Island == "S_Shelf")),
                          "S_Shore"=nrow(filter(annote.dmps, Relation_to_Island == "S_Shore")))
sum.islands <- rbind(rbind(sum.islands / sum(sum.islands), sum.islands / sum(sum.islands)), sum.islands / sum(sum.islands))


sum.regions <- data.frame(t(data.frame(table(annote.dmps$region))[,2]/sum(table(annote.dmps$region))))[,c(1,2,4,6,7)]
sum.regions <- rbind(rbind(sum.regions / sum(sum.regions), sum.regions / sum(sum.regions)), sum.regions / sum(sum.regions))
colnames(sum.regions) <- c("upstream", "promoter","inside", "close.to.3.", "downstream")

# 14. Annotate Genic Locations ####
All.samples.not.detailed.1 <- lapply(All.samples.not.detailed.1, function(x) as.data.frame(x))
All.samples.not.detailed.2 <- lapply(All.samples.not.detailed.2, function(x) as.data.frame(x))

# Genic Region
INdmrs.regions <- data.frame(rbindlist(lapply(lapply(All.samples.not.detailed.1[7:12], function(x) data.frame("upstream"=nrow(filter(x, region == "upstream")),
                                                                                                              "promoter"=nrow(filter(x, region == "promoter")),
                                                                                                              "inside"=nrow(filter(x, region == "inside")),
                                                                                                              "close to 3'"=nrow(filter(x, region == "close to 3'")),
                                                                                                              "downstream"=nrow(filter(x, region == "downstream")))), function(x) x/sum(x))))

INdmrs.regions <- INdmrs.regions / rbind(sum.regions, sum.regions[rep(rep(1),3),])
INdmrs.regions <- rbind(INdmrs.regions, sum.regions[1,])[c(7,1:6),]
rownames(INdmrs.regions) <- c("allProbes", names(All.samples.not.detailed.1[7:12]))

OUTdmrs.regions <- data.frame(rbindlist(lapply(lapply(All.samples.not.detailed.2, function(x) data.frame("upstream"=nrow(filter(x, region == "upstream")),
                                                                                                         "promoter"=nrow(filter(x, region == "promoter")),
                                                                                                         "inside"=nrow(filter(x, region == "inside")),
                                                                                                         "close to 3'"=nrow(filter(x, region == "close to 3'")),
                                                                                                         "downstream"=nrow(filter(x, region == "downstream")))), function(x) x/sum(x))))

OUTdmrs.regions <- OUTdmrs.regions / rbind(sum.regions, sum.regions[rep(rep(1),3),])
OUTdmrs.regions <- rbind(OUTdmrs.regions, sum.regions[1,])[c(7,1:6),]
rownames(OUTdmrs.regions) <- c("allProbes", names(All.samples.not.detailed.1[7:12]))

regions <- rbind(INdmrs.regions, OUTdmrs.regions[-1,])
rownames(regions) <- c("baseline", "In DMRs Fe 0.3Gy", "In DMRs Fe 1Gy", "In DMRs Fe(X) 1Gy",
                       "In DMRs Si 0.3Gy", "In DMRs Si 1Gy", "In DMRs X 1Gy",
                       
                       "Out DMRs Fe 0.3Gy", "Out DMRs Fe 1Gy", "Out DMRs Fe(X) 1Gy",
                       "Out DMRs Si 0.3Gy", "Out DMRs Si 1Gy", "Out DMRs X 1Gy")

region2 <- regions[c(3,9),]
region2$mod <- rownames(region2)

#Figure_1G
ggplot(data = melt(region2), aes(x = variable, y = value, fill = mod, label=mod)) +
  geom_bar(position="dodge", stat="identity", width=0.7) +
  scale_fill_viridis(discrete=T, name="") +
  labs(x="", y="Enrichment", fill="Genic\nLocation") +
  ggtitle("") +
  theme_classic() + 
  theme(
    axis.title.x = element_text(size = 12, face="bold"),
    axis.text.x = element_text(size = 12, angle = -345, vjust = 1, hjust=1),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 12, face="bold"),
    title = element_text(size = 12, face = "bold")) +
  geom_hline(yintercept=1, linetype="dashed", color = "red")

#Statistics for absolute values
IN.Chisq.regions <- data.frame(rbindlist(lapply(All.samples.not.detailed.1[7:12], function(x) data.frame("upstream"=nrow(filter(x, region == "upstream")),
                                                                                                              "promoter"=nrow(filter(x, region == "promoter")),
                                                                                                              "inside"=nrow(filter(x, region == "inside")),
                                                                                                              "close to 3'"=nrow(filter(x, region == "close to 3'")),
                                                                                                              "downstream"=nrow(filter(x, region == "downstream"))))))

OUT.Chisq.regions <- data.frame(rbindlist(lapply(All.samples.not.detailed.2, function(x) data.frame("upstream"=nrow(filter(x, region == "upstream")),
                                                                                                         "promoter"=nrow(filter(x, region == "promoter")),
                                                                                                         "inside"=nrow(filter(x, region == "inside")),
                                                                                                         "close to 3'"=nrow(filter(x, region == "close to 3'")),
                                                                                                         "downstream"=nrow(filter(x, region == "downstream"))))))

rownames(IN.Chisq.regions) <- names(All.samples.not.detailed.1[7:12])
rownames(OUT.Chisq.regions) <- names(All.samples.not.detailed.2)

chisq.test(cbind(t(IN.Chisq.regions[1,]), t(OUT.Chisq.regions[1,])))
chisq.test(cbind(t(IN.Chisq.regions[2,]), t(OUT.Chisq.regions[2,])))
chisq.test(cbind(t(IN.Chisq.regions[3,]), t(OUT.Chisq.regions[3,])))
chisq.test(cbind(t(IN.Chisq.regions[4,]), t(OUT.Chisq.regions[4,])))
chisq.test(cbind(t(IN.Chisq.regions[5,]), t(OUT.Chisq.regions[5,])))
chisq.test(cbind(t(IN.Chisq.regions[6,]), t(OUT.Chisq.regions[6,])))

# 15. Annotate to CpG Islands ####
All.samples.not.detailed.1.increase.Fe <- lapply(All.samples.not.detailed.1, function(x) filter(x, (mFe - rowMeans(x[,c(110,122,134)])) >= 0.58))
All.samples.not.detailed.1.decrease.Fe <- lapply(All.samples.not.detailed.1, function(x) filter(x, (mFe - rowMeans(x[,c(110,122,134)])) <= -0.58))

All.samples.not.detailed.1.increase.SiX <- lapply(All.samples.not.detailed.1, function(x) filter(x, (mSiX - rowMeans(x[,c(110,122,134)])) >= 0.58))
All.samples.not.detailed.1.decrease.SiX <- lapply(All.samples.not.detailed.1, function(x) filter(x, (mSiX - rowMeans(x[,c(110,122,134)])) <= -0.58))

All.samples.notINdmrs.detailed.1.increase.Fe <- lapply(All.samples.not.detailed.2, function(x) filter(x, (x[,102] - rowMeans(x[,c(64,76,88)])) >= 0.58))
All.samples.notINdmrs.detailed.1.decrease.Fe <- lapply(All.samples.not.detailed.2, function(x) filter(x, (x[,102] - rowMeans(x[,c(64,76,88)])) <= -0.58))

All.samples.notINdmrs.detailed.1.increase.SiX <- lapply(All.samples.not.detailed.2, function(x) filter(x, (x[,103] - rowMeans(x[,c(64,76,88)])) >= 0.58))
All.samples.notINdmrs.detailed.1.decrease.SiX <- lapply(All.samples.not.detailed.2, function(x) filter(x, (x[,103] - rowMeans(x[,c(64,76,88)])) <= -0.58))


INdmrs.increase.Fe.islands <- data.frame(rbindlist(lapply(lapply(All.samples.not.detailed.1.increase.Fe[7:12], function(x) data.frame("Island"=nrow(filter(x, Relation_to_Island == "Island")),
                                                                                                                  "N_Shelf"=nrow(filter(x, Relation_to_Island == "N_Shelf")),
                                                                                                                  "N_Shore"=nrow(filter(x, Relation_to_Island == "N_Shore")),
                                                                                                                  "OpenSea"=nrow(filter(x, Relation_to_Island == "OpenSea")),
                                                                                                                  "S_Shelf"=nrow(filter(x, Relation_to_Island == "S_Shelf")),
                                                                                                                  "S_Shore"=nrow(filter(x, Relation_to_Island == "S_Shore")))), function(x) x/sum(x))))
INdmrs.increase.Fe.islands <- INdmrs.increase.Fe.islands / rbind(sum.islands, sum.islands[rep(rep(1),3), ])

INdmrs.increase.Fe.islands <- rbind(INdmrs.increase.Fe.islands, sum.islands[1,])[c(7,1:6),]
rownames(INdmrs.increase.Fe.islands) <- c("rawProbes", names(All.samples.not.detailed.1[7:12]))


INdmrs.decrease.Fe.islands <- data.frame(rbindlist(lapply(lapply(All.samples.not.detailed.1.decrease.Fe[7:12], function(x) data.frame("Island"=nrow(filter(x, Relation_to_Island == "Island")),
                                                                                                                                      "N_Shelf"=nrow(filter(x, Relation_to_Island == "N_Shelf")),
                                                                                                                                      "N_Shore"=nrow(filter(x, Relation_to_Island == "N_Shore")),
                                                                                                                                      "OpenSea"=nrow(filter(x, Relation_to_Island == "OpenSea")),
                                                                                                                                      "S_Shelf"=nrow(filter(x, Relation_to_Island == "S_Shelf")),
                                                                                                                                      "S_Shore"=nrow(filter(x, Relation_to_Island == "S_Shore")))), function(x) x/sum(x))))
INdmrs.decrease.Fe.islands <- INdmrs.decrease.Fe.islands / rbind(sum.islands, sum.islands[rep(rep(1),3), ])

INdmrs.decrease.Fe.islands <- rbind(INdmrs.decrease.Fe.islands, sum.islands[1,])[c(7,1:6),]
rownames(INdmrs.decrease.Fe.islands) <- c("rawProbes", names(All.samples.not.detailed.1[7:12]))

INdmrs.increase.SiX.islands <- data.frame(rbindlist(lapply(lapply(All.samples.not.detailed.1.increase.SiX[7:12], function(x) data.frame("Island"=nrow(filter(x, Relation_to_Island == "Island")),
                                                                                                                                      "N_Shelf"=nrow(filter(x, Relation_to_Island == "N_Shelf")),
                                                                                                                                      "N_Shore"=nrow(filter(x, Relation_to_Island == "N_Shore")),
                                                                                                                                      "OpenSea"=nrow(filter(x, Relation_to_Island == "OpenSea")),
                                                                                                                                      "S_Shelf"=nrow(filter(x, Relation_to_Island == "S_Shelf")),
                                                                                                                                      "S_Shore"=nrow(filter(x, Relation_to_Island == "S_Shore")))), function(x) x/sum(x))))
INdmrs.increase.SiX.islands <- INdmrs.increase.SiX.islands / rbind(sum.islands, sum.islands[rep(rep(1),3), ])

INdmrs.increase.SiX.islands <- rbind(INdmrs.increase.SiX.islands, sum.islands[1,])[c(7,1:6),]
rownames(INdmrs.increase.SiX.islands) <- c("rawProbes", names(All.samples.not.detailed.1[7:12]))

INdmrs.decrease.SiX.islands <- data.frame(rbindlist(lapply(lapply(All.samples.not.detailed.1.decrease.SiX[7:12], function(x) data.frame("Island"=nrow(filter(x, Relation_to_Island == "Island")),
                                                                                                                                      "N_Shelf"=nrow(filter(x, Relation_to_Island == "N_Shelf")),
                                                                                                                                      "N_Shore"=nrow(filter(x, Relation_to_Island == "N_Shore")),
                                                                                                                                      "OpenSea"=nrow(filter(x, Relation_to_Island == "OpenSea")),
                                                                                                                                      "S_Shelf"=nrow(filter(x, Relation_to_Island == "S_Shelf")),
                                                                                                                                      "S_Shore"=nrow(filter(x, Relation_to_Island == "S_Shore")))), function(x) x/sum(x))))
INdmrs.decrease.SiX.islands <- INdmrs.decrease.SiX.islands / rbind(sum.islands, sum.islands[rep(rep(1),3), ])

INdmrs.decrease.SiX.islands <- rbind(INdmrs.decrease.SiX.islands, sum.islands[1,])[c(7,1:6),]
rownames(INdmrs.decrease.SiX.islands) <- c("rawProbes", names(All.samples.not.detailed.1[7:12]))

OUTdmrs.increase.Fe.islands <- data.frame(rbindlist(lapply(lapply(All.samples.notINdmrs.detailed.1.increase.Fe, function(x) data.frame("Island"=nrow(filter(x, Relation_to_Island == "Island")),
                                                                                                                                      "N_Shelf"=nrow(filter(x, Relation_to_Island == "N_Shelf")),
                                                                                                                                      "N_Shore"=nrow(filter(x, Relation_to_Island == "N_Shore")),
                                                                                                                                      "OpenSea"=nrow(filter(x, Relation_to_Island == "OpenSea")),
                                                                                                                                      "S_Shelf"=nrow(filter(x, Relation_to_Island == "S_Shelf")),
                                                                                                                                      "S_Shore"=nrow(filter(x, Relation_to_Island == "S_Shore")))), function(x) x/sum(x))))
OUTdmrs.increase.Fe.islands <- OUTdmrs.increase.Fe.islands / rbind(sum.islands, sum.islands[rep(rep(1),3), ])

OUTdmrs.increase.Fe.islands <- rbind(OUTdmrs.increase.Fe.islands, sum.islands[1,])[c(7,1:6),]
rownames(OUTdmrs.increase.Fe.islands) <- c("rawProbes", names(All.samples.not.detailed.1[7:12]))


OUTdmrs.decrease.Fe.islands <- data.frame(rbindlist(lapply(lapply(All.samples.notINdmrs.detailed.1.decrease.Fe, function(x) data.frame("Island"=nrow(filter(x, Relation_to_Island == "Island")),
                                                                                                                                      "N_Shelf"=nrow(filter(x, Relation_to_Island == "N_Shelf")),
                                                                                                                                      "N_Shore"=nrow(filter(x, Relation_to_Island == "N_Shore")),
                                                                                                                                      "OpenSea"=nrow(filter(x, Relation_to_Island == "OpenSea")),
                                                                                                                                      "S_Shelf"=nrow(filter(x, Relation_to_Island == "S_Shelf")),
                                                                                                                                      "S_Shore"=nrow(filter(x, Relation_to_Island == "S_Shore")))), function(x) x/sum(x))))
OUTdmrs.decrease.Fe.islands <- OUTdmrs.decrease.Fe.islands / rbind(sum.islands, sum.islands[rep(rep(1),3), ])

OUTdmrs.decrease.Fe.islands <- rbind(OUTdmrs.decrease.Fe.islands, sum.islands[1,])[c(7,1:6),]
rownames(OUTdmrs.decrease.Fe.islands) <- c("rawProbes", names(All.samples.not.detailed.1[7:12]))


OUTdmrs.increase.SiX.islands <- data.frame(rbindlist(lapply(lapply(All.samples.notINdmrs.detailed.1.increase.SiX, function(x) data.frame("Island"=nrow(filter(x, Relation_to_Island == "Island")),
                                                                                                                                        "N_Shelf"=nrow(filter(x, Relation_to_Island == "N_Shelf")),
                                                                                                                                        "N_Shore"=nrow(filter(x, Relation_to_Island == "N_Shore")),
                                                                                                                                        "OpenSea"=nrow(filter(x, Relation_to_Island == "OpenSea")),
                                                                                                                                        "S_Shelf"=nrow(filter(x, Relation_to_Island == "S_Shelf")),
                                                                                                                                        "S_Shore"=nrow(filter(x, Relation_to_Island == "S_Shore")))), function(x) x/sum(x))))
OUTdmrs.increase.SiX.islands <- OUTdmrs.increase.SiX.islands / rbind(sum.islands, sum.islands[rep(rep(1),3), ])

OUTdmrs.increase.SiX.islands <- rbind(OUTdmrs.increase.SiX.islands, sum.islands[1,])[c(7,1:6),]
rownames(OUTdmrs.increase.SiX.islands) <- c("rawProbes", names(All.samples.not.detailed.1[7:12]))


OUTdmrs.decrease.SiX.islands <- data.frame(rbindlist(lapply(lapply(All.samples.notINdmrs.detailed.1.decrease.SiX, function(x) data.frame("Island"=nrow(filter(x, Relation_to_Island == "Island")),
                                                                                                                                        "N_Shelf"=nrow(filter(x, Relation_to_Island == "N_Shelf")),
                                                                                                                                        "N_Shore"=nrow(filter(x, Relation_to_Island == "N_Shore")),
                                                                                                                                        "OpenSea"=nrow(filter(x, Relation_to_Island == "OpenSea")),
                                                                                                                                        "S_Shelf"=nrow(filter(x, Relation_to_Island == "S_Shelf")),
                                                                                                                                        "S_Shore"=nrow(filter(x, Relation_to_Island == "S_Shore")))), function(x) x/sum(x))))
OUTdmrs.decrease.SiX.islands <- OUTdmrs.decrease.SiX.islands / rbind(sum.islands, sum.islands[rep(rep(1),3), ])

OUTdmrs.decrease.SiX.islands <- rbind(OUTdmrs.decrease.SiX.islands, sum.islands[1,])[c(7,1:6),]
rownames(OUTdmrs.decrease.SiX.islands) <- c("rawProbes", names(All.samples.not.detailed.1[7:12]))

islands <- rbind(rbind(rbind(rbind(INdmrs.increase.Fe.islands[1:4,], INdmrs.decrease.Fe.islands[2:4,]), rbind(INdmrs.increase.SiX.islands[5:6,], INdmrs.decrease.SiX.islands[5:6,]),rbind(INdmrs.increase.SiX.islands[7,], INdmrs.decrease.SiX.islands[7,]))),
                 rbind(rbind(rbind(OUTdmrs.increase.Fe.islands[2:4,], OUTdmrs.decrease.Fe.islands[2:4,]), rbind(OUTdmrs.increase.SiX.islands[5:6,], OUTdmrs.decrease.SiX.islands[5:6,]),rbind(OUTdmrs.increase.SiX.islands[7,], OUTdmrs.decrease.SiX.islands[7,]))))

rownames(islands) <- c("baseline", "In DMRs ↑Fe 0.3Gy", "In DMRs ↑Fe 1Gy", "In DMRs ↑Fe(X) 1Gy", "In DMRs ↓Fe 0.3Gy", "In DMRs ↓Fe 1Gy", "In DMRs ↓Fe(X) 1Gy",
                       "In DMRs ↑Si 0.3Gy", "In DMRs ↑Si 1Gy", "In DMRs ↓Si 0.3Gy", "In DMRs ↓Si 1Gy", "In DMRs ↑X 1Gy", "In DMRs ↓X 1Gy",
                       
                       "Out DMRs ↑Fe 0.3Gy", "Out DMRs ↑Fe 1Gy", "Out DMRs ↑Fe(X) 1Gy", "Out DMRs ↓Fe 0.3Gy", "Out DMRs ↓Fe 1Gy", "Out DMRs ↓Fe(X) 1Gy",
                       "Out DMRs ↑Si 0.3Gy", "Out DMRs ↑Si 1Gy", "Out DMRs ↓Si 0.3Gy", "Out DMRs ↓Si 1Gy", "Out DMRs ↑X 1Gy", "Out DMRs ↓X 1Gy")

islands.joined <- t(islands[c(1,3,15,6,18),c(4,2,3,1,5,6)])
islands.joined <- cbind(cbind(islands.joined[,1], rowMeans(islands.joined[,c(2,4)])), rowMeans(islands.joined[,c(3,5)]))
colnames(islands.joined) <- c("baseline", "In DMRs Fe 1Gy", "Out DMRs Fe 1Gy")

#Figure_S1J
ggplot(data = melt(t(islands[c(1,3,6,15,18),c(4,2,3,1,5,6)])), aes(x = Var2, y = value, fill = Var1, label=value)) +
  geom_bar(position="fill", stat="identity", width=0.7) +
  scale_fill_viridis(discrete=T, name="Relation to\nCpG Island") +
  labs(x="", y="Enrichment (%)", fill="") +
  ggtitle("Fe Associated DMPs (In or outside DMRs)") +
  theme_classic() + 
  theme(
    axis.title.x = element_text(size = 12, face="bold"),
    axis.text.x = element_text(size = 12, angle = -330, vjust = 1, hjust=1),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 12, face="bold"),
    title = element_text(size = 12, face = "bold"))

#Figure_S1K
ggplot(data = melt(t(islands[c(1,9,11,21,23),c(4,2,3,1,5,6)])), aes(x = Var2, y = value, fill = Var1, label=value)) +
  geom_bar(position="fill", stat="identity", width=0.7) +
  scale_fill_viridis(discrete=T, name="Relation to\nCpG Island") +
  labs(x="", y="Enrichment (%)", fill="") +
  ggtitle("Si Associated DMPs (In or outside DMRs)") +
  theme_classic() + 
  theme(
    axis.title.x = element_text(size = 12, face="bold"),
    axis.text.x = element_text(size = 12, angle = -330, vjust = 1, hjust=1),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 12, face="bold"),
    title = element_text(size = 12, face = "bold"))

#Figure_S1L
ggplot(data = melt(t(islands[c(1,12:13,24:25),c(4,2,3,1,5,6)])), aes(x = Var2, y = value, fill = Var1, label=value)) +
  geom_bar(position="fill", stat="identity", width=0.7) +
  scale_fill_viridis(discrete=T, name="Relation to\nCpG Island") +
  labs(x="", y="Enrichment (%)", fill="") +
  ggtitle("X ray Associated DMPs (In or outside DMRs)") +
  theme_classic() + 
  theme(
    axis.title.x = element_text(size = 12, face="bold"),
    axis.text.x = element_text(size = 12, angle = -330, vjust = 1, hjust=1),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 12, face="bold"),
    title = element_text(size = 12, face = "bold"))

#Fractions
dmps.joined <- Map(lapply(All.samples.not.detailed.1[7:12], function(x) x[,c(1:5,55)]), 
                   lapply(All.samples.not.detailed.2, function(x) x[,c(1:5,21)]), f=rbind)

test <- data.frame(rbindlist(lapply(lapply(dmps.joined, function(x) data.frame("Island"=nrow(filter(x, Relation_to_Island == "Island")),
                                                                                                                                         "N_Shelf"=nrow(filter(x, Relation_to_Island == "N_Shelf")),
                                                                                                                                         "N_Shore"=nrow(filter(x, Relation_to_Island == "N_Shore")),
                                                                                                                                         "OpenSea"=nrow(filter(x, Relation_to_Island == "OpenSea")),
                                                                                                                                         "S_Shelf"=nrow(filter(x, Relation_to_Island == "S_Shelf")),
                                                                                                                                         "S_Shore"=nrow(filter(x, Relation_to_Island == "S_Shore")))), function(x) x/sum(x))))

test1 <- t(rbind(table(annote.dmps$Relation_to_Island)/sum(table(annote.dmps$Relation_to_Island)), test))
colnames(test1) <- c(names(dmps.joined), "baseline")
colSums(test1)

test2 <- rbind(rbind(test1[c(1,4),], colMeans(test1[2:3,])), colMeans(test1[5:6,]))
rownames(test2) <- c("Island", "OpenSea", "N_Shore_Shelf", "S_Shore_Shelf")
colnames(test2) <- gsub("ALL_DMPs_", "", colnames(test2))
test2 <- test2[c(2,3,1,4),]

#Figure_1H
ggplot(data = melt(test2[,c(2,5,7)]), aes(x = Var1, y = value, fill = Var2, label=Var2)) +
  geom_bar(position="dodge", stat="identity", width=0.7) +
  scale_fill_viridis(discrete=T, name="") +
  labs(x="", y="Fraction (%)", fill="Genic\nLocation") +
  ggtitle("") +
  theme_classic() + 
  theme(
    axis.title.x = element_text(size = 12, face="bold"),
    axis.text.x = element_text(size = 12, angle = -345, vjust = 1, hjust=1),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 12, face="bold"),
    title = element_text(size = 12, face = "bold"))

#Figure_1I
ggplot(data = melt(test2[,c(3,6,7)]), aes(x = Var1, y = value, fill = Var2, label=Var2)) +
  geom_bar(position="dodge", stat="identity", width=0.7) +
  scale_fill_viridis(discrete=T, name="") +
  labs(x="", y="Fraction (%)", fill="Genic\nLocation") +
  ggtitle("") +
  theme_classic() + 
  theme(
    axis.title.x = element_text(size = 12, face="bold"),
    axis.text.x = element_text(size = 12, angle = -345, vjust = 1, hjust=1),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 12, face="bold"),
    title = element_text(size = 12, face = "bold"))

# 16. Chronic methylation change ####
dmps.joined.chronic <- Map(lapply(All.samples.not.detailed.1[7:12], function(x) x[,c(1:5,110:141)]), 
                           lapply(All.samples.not.detailed.2, function(x) x[,c(1:5,64:95)]), f=rbind)

#Fe
chronic.Fe0.3.increase <- filter(dmps.joined.chronic$DMPs_Fe_0.3Gy, (f0.3.dmp - rowMeans(dmps.joined.chronic$DMPs_Fe_0.3Gy[,c(6,18,30)])) >= 0.58)
chronic.Fe0.3.increase <- data.frame(melt(colMeans(na.omit(data.frame("con"=rowMeans(chronic.Fe0.3.increase[,c(6,18,30)]), chronic.Fe0.3.increase[,c(7,10,13,16)])))),
                                     "time"=c("2_control","2","13_14","21_22","53_62"), "sample"=rep("Fe 0.3Gy",5))
chronic.Fe0.3.increase$time <- factor(chronic.Fe0.3.increase$time, levels=c("2_control","2","13_14","21_22","53_62"), labels=c("2_control","2","13_14","21_22","53_62"))

chronic.Fe0.3.decrease <- filter(dmps.joined.chronic$DMPs_Fe_0.3Gy, (f0.3.dmp - rowMeans(dmps.joined.chronic$DMPs_Fe_0.3Gy[,c(6,18,30)])) <= 0.58)
chronic.Fe0.3.decrease <- data.frame(melt(colMeans(na.omit(data.frame("con"=rowMeans(chronic.Fe0.3.decrease[,c(6,18,30)]), chronic.Fe0.3.decrease[,c(7,10,13,16)])))),
                                     "time"=c("2_control","2","13_14","21_22","53_62"), "sample"=rep("Fe 0.3Gy",5))
chronic.Fe0.3.decrease$time <- factor(chronic.Fe0.3.decrease$time, levels=c("2_control","2","13_14","21_22","53_62"), labels=c("2_control","2","13_14","21_22","53_62"))

#
chronic.Fe1.increase <- filter(dmps.joined.chronic$DMPs_Fe_1Gy, (f1.dmp - rowMeans(dmps.joined.chronic$DMPs_Fe_1Gy[,c(6,18,30)])) >= 0.58)
chronic.Fe1.increase <- data.frame(melt(colMeans(na.omit(data.frame("con"=rowMeans(chronic.Fe1.increase[,c(6,18,30)]), chronic.Fe1.increase[,c(8,11,14,17)])))),
                                   "time"=c("2_control","2","13_14","21_22","53_62"), "sample"=rep("Fe 1Gy",5))
chronic.Fe1.increase$time <- factor(chronic.Fe1.increase$time, levels=c("2_control","2","13_14","21_22","53_62"), labels=c("2_control","2","13_14","21_22","53_62"))

chronic.Fe1.decrease <- filter(dmps.joined.chronic$DMPs_Fe_1Gy, (f1.dmp - rowMeans(dmps.joined.chronic$DMPs_Fe_1Gy[,c(6,18,30)])) <= 0.58)
chronic.Fe1.decrease <- data.frame(melt(colMeans(na.omit(data.frame("con"=rowMeans(chronic.Fe1.decrease[,c(6,18,30)]), chronic.Fe1.decrease[,c(8,11,14,17)])))),
                                   "time"=c("2_control","2","13_14","21_22","53_62"), "sample"=rep("Fe 1Gy",5))
chronic.Fe1.decrease$time <- factor(chronic.Fe1.decrease$time, levels=c("2_control","2","13_14","21_22","53_62"), labels=c("2_control","2","13_14","21_22","53_62"))

#
chronic.FeX1.increase <- filter(dmps.joined.chronic$DMPs_FeX_1Gy, (f1.dmp - rowMeans(dmps.joined.chronic$DMPs_FeX_1Gy[,c(6,18,30)])) >= 0.58)
chronic.FeX1.increase <- data.frame(melt(colMeans(na.omit(data.frame("con"=rowMeans(chronic.FeX1.increase[,c(6,18,30)]), chronic.FeX1.increase[,c(8,11,14,17)])))),
                                    "time"=c("2_control","2","13_14","21_22","53_62"), "sample"=rep("FeX 1Gy",5))
chronic.FeX1.increase$time <- factor(chronic.FeX1.increase$time, levels=c("2_control","2","13_14","21_22","53_62"), labels=c("2_control","2","13_14","21_22","53_62"))

chronic.FeX1.decrease <- filter(dmps.joined.chronic$DMPs_FeX_1Gy, (f1.dmp - rowMeans(dmps.joined.chronic$DMPs_FeX_1Gy[,c(6,18,30)])) <= 0.58)
chronic.FeX1.decrease <- data.frame(melt(colMeans(na.omit(data.frame("con"=rowMeans(chronic.FeX1.decrease[,c(6,18,30)]), chronic.FeX1.decrease[,c(8,11,14,17)])))),
                                    "time"=c("2_control","2","13_14","21_22","53_62"), "sample"=rep("FeX 1Gy",5))
chronic.FeX1.decrease$time <- factor(chronic.FeX1.decrease$time, levels=c("2_control","2","13_14","21_22","53_62"), labels=c("2_control","2","13_14","21_22","53_62"))

# Si
chronic.Si0.3.increase <- filter(dmps.joined.chronic$DMPs_Si_0.3Gy, (s0.3.dmp - rowMeans(dmps.joined.chronic$DMPs_Si_0.3Gy[,c(6,18,30)])) >= 0.58)
chronic.Si0.3.increase <- data.frame(melt(colMeans(na.omit(data.frame("con"=rowMeans(chronic.Si0.3.increase[,c(6,18,30)]), chronic.Si0.3.increase[,c(19,22,25,28)])))),
                                     "time"=c("2_control","2","13_14","21_22","53_62"), "sample"=rep("Si 0.3Gy",5))
chronic.Si0.3.increase$time <- factor(chronic.Si0.3.increase$time, levels=c("2_control","2","13_14","21_22","53_62"), labels=c("2_control","2","13_14","21_22","53_62"))

chronic.Si0.3.decrease <- filter(dmps.joined.chronic$DMPs_Si_0.3Gy, (s0.3.dmp - rowMeans(dmps.joined.chronic$DMPs_Si_0.3Gy[,c(6,18,30)])) <= 0.58)
chronic.Si0.3.decrease <- data.frame(melt(colMeans(na.omit(data.frame("con"=rowMeans(chronic.Si0.3.decrease[,c(6,18,30)]), chronic.Si0.3.decrease[,c(19,22,25,28)])))),
                                     "time"=c("2_control","2","13_14","21_22","53_62"), "sample"=rep("Si 0.3Gy",5))
chronic.Si0.3.decrease$time <- factor(chronic.Si0.3.decrease$time, levels=c("2_control","2","13_14","21_22","53_62"), labels=c("2_control","2","13_14","21_22","53_62"))

#
chronic.Si1.increase <- filter(dmps.joined.chronic$DMPs_Si_1Gy, (s1.dmp - rowMeans(dmps.joined.chronic$DMPs_Si_1Gy[,c(6,18,30)])) >= 0.58)
chronic.Si1.increase <- data.frame(melt(colMeans(na.omit(data.frame("con"=rowMeans(chronic.Si1.increase[,c(6,18,30)]), chronic.Si1.increase[,c(20,23,26,29)])))),
                                   "time"=c("2_control","2","13_14","21_22","53_62"), "sample"=rep("Si 1Gy",5))
chronic.Si1.increase$time <- factor(chronic.Si1.increase$time, levels=c("2_control","2","13_14","21_22","53_62"), labels=c("2_control","2","13_14","21_22","53_62"))

chronic.Si1.decrease <- filter(dmps.joined.chronic$DMPs_Si_1Gy, (s1.dmp - rowMeans(dmps.joined.chronic$DMPs_Si_1Gy[,c(6,18,30)])) <= 0.58)
chronic.Si1.decrease <- data.frame(melt(colMeans(na.omit(data.frame("con"=rowMeans(chronic.Si1.decrease[,c(6,18,30)]), chronic.Si1.decrease[,c(20,23,26,29)])))),
                                   "time"=c("2_control","2","13_14","21_22","53_62"), "sample"=rep("Si 1Gy",5))
chronic.Si1.decrease$time <- factor(chronic.Si1.decrease$time, levels=c("2_control","2","13_14","21_22","53_62"), labels=c("2_control","2","13_14","21_22","53_62"))

# X
chronic.X1.increase <- filter(dmps.joined.chronic$DMPs_X_1Gy, (x1.dmp - rowMeans(dmps.joined.chronic$DMPs_X_1Gy[,c(6,18,30)])) >= 0.58)
chronic.X1.increase <- data.frame(melt(colMeans(na.omit(data.frame("con"=rowMeans(chronic.X1.increase[,c(6,18,30)]), chronic.X1.increase[,c(31,33,35,37)])))),
                                  "time"=c("2_control","2","13_14","21_22","53_62"), "sample"=rep("X 1Gy",5))
chronic.X1.increase$time <- factor(chronic.X1.increase$time, levels=c("2_control","2","13_14","21_22","53_62"), labels=c("2_control","2","13_14","21_22","53_62"))

chronic.X1.decrease <- filter(dmps.joined.chronic$DMPs_X_1Gy, (x1.dmp - rowMeans(dmps.joined.chronic$DMPs_X_1Gy[,c(6,18,30)])) <= 0.58)
chronic.X1.decrease <- data.frame(melt(colMeans(na.omit(data.frame("con"=rowMeans(chronic.X1.decrease[,c(6,18,30)]), chronic.X1.decrease[,c(31,33,35,37)])))),
                                  "time"=c("2_control","2","13_14","21_22","53_62"), "sample"=rep("X 1Gy",5))
chronic.X1.decrease$time <- factor(chronic.X1.decrease$time, levels=c("2_control","2","13_14","21_22","53_62"), labels=c("2_control","2","13_14","21_22","53_62"))

fig6 <- do.call("rbind", list(chronic.Fe0.3.increase, chronic.Fe0.3.decrease, chronic.Fe1.increase, chronic.Fe1.decrease, chronic.FeX1.increase, chronic.FeX1.decrease,
                              chronic.Si0.3.increase, chronic.Si0.3.decrease, chronic.Si1.increase, chronic.Si1.decrease,
                              chronic.X1.increase, chronic.X1.decrease))

fig6$direction <- c(rep("↑Fe 0.3Gy", 5), rep("↓Fe 0.3Gy", 5), rep("↑Fe 1Gy", 5), rep("↓Fe 1Gy", 5), rep("↑Fe(X) 1Gy", 5), rep("↓Fe(X) 1Gy", 5), rep("↑Si 0.3Gy", 5), rep("↓Si 0.3Gy", 5), rep("↑Si 1Gy", 5), rep("↓Si 1Gy", 5), rep("↑X 1Gy", 5), rep("↓X 1Gy", 5))
fig6$sample <- c(rep("Fe 0.3Gy", 10), rep("Fe 1Gy", 10), rep("Fe(X) 1Gy", 10), rep("Si 0.3Gy", 10), rep("Si 1Gy", 10), rep("X 1Gy", 10))

#Figure_4A
fig4b<-filter(fig6, sample == "Fe 1Gy")
fig4b$time <- rep(c("2_control", "2", "14", "22", "53"),2)
fig4b$time <- factor(fig4b$time, levels=c("2_control", "2", "14", "22", "53"), labels=c("2_control", "2", "14", "22", "53"))

ggplot(fig4b, aes(x=time, y=value, col=as.factor(direction))) +
  geom_point(size=4) +
  scale_color_viridis(discrete=TRUE) +
  geom_line(aes(group=interaction(direction)), size=2, alpha=0.5) +
  theme_classic() +
  labs(title="", x="Time (Days)", y="Mean Methylation Change versus control (M Value)", color="Radiation") +
  theme(
    axis.title.x = element_text(size = 10, face="bold"),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 10, face="bold"),
    title = element_text(size = 10, face = "bold"))

#Figure_4B
fig4c<-filter(fig6, sample == "Si 1Gy")
fig4c$time <- rep(c("2_control", "2", "13", "21", "62"),2)
fig4c$time <- factor(fig4c$time, levels=c("2_control", "2", "13", "21", "62"), labels=c("2_control", "2", "13", "21", "62"))

ggplot(fig4c, aes(x=time, y=value, col=as.factor(direction))) +
  geom_point(size=4) +
  scale_color_viridis(discrete=TRUE) +
  geom_line(aes(group=interaction(direction)), size=2, alpha=0.5) +
  theme_classic() +
  labs(title="", x="Time (Days)", y="Mean Methylation Change versus control (M Value)", color="Radiation") +
  theme(
    axis.title.x = element_text(size = 10, face="bold"),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 10, face="bold"),
    title = element_text(size = 10, face = "bold"))

#Figure_S4A
ggplot(filter(fig6, !sample == "Fe 1Gy" & !sample == "Si 1Gy"), aes(x=time, y=value, col=as.factor(direction))) +
  geom_point(size=4) +
  scale_color_viridis(discrete=TRUE) +
  geom_line(aes(group=interaction(direction)), size=2, alpha=0.5) +
  theme_classic() +
  labs(title="", x="Time (Days)", y="Mean Methylation Level (M Value)", color="Radiation") +
  theme(
    axis.title.x = element_text(size = 12, face="bold"),
    axis.text.x = element_text(size = 12, angle = -315, vjust = 1, hjust=1),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 12, face="bold"),
    title = element_text(size = 12, face = "bold")) +
  facet_wrap(vars(sample), ncol=6)


# 17. Chromosomes with the highest DMPs frequency ####

frequency <- append(lapply(append(append(Map(rbind, DMRs.increased.Fe.C, DMRs.decreased.Fe.C), Map(rbind, DMRs.increased.SiX.D, DMRs.decreased.SiX.D)),
                                  append(Map(rbind, DMPs.increased.Fe.C, DMPs.decreased.Fe.C), Map(rbind, DMPs.increased.SiX.D, DMPs.decreased.SiX.D))), function(x) as.character(x$seqnames.x)),
                    
                    lapply(append(Map(rbind, notINdmrs.increased.Fe, notINdmrs.decreased.Fe), Map(rbind, notINdmrs.increased.SiX, notINdmrs.decreased.SiX)), function(x) x$chr))
frequency <- append(frequency, mapply(function(x,y) c(x,y), x=frequency[7:12], y=frequency[13:18]))

names(frequency) <- c("DMRs_Fe_0.3Gy", "DMRs_Fe_1Gy", "DMRs_FeX_1Gy", "DMRs_Si_0.3Gy", "DMRs_Si_1Gy", "DMRs_X_1Gy",
                      "DMPs_Fe_0.3Gy", "DMPs_Fe_1Gy", "DMPs_FeX_1Gy", "DMPs_Si_0.3Gy", "DMPs_Si_1Gy", "DMPs_X_1Gy",
                      "notINdmrs_Fe_0.3Gy", "notINdmrs_Fe_1Gy", "notINdmrs_FeX_1Gy", "notINdmrs_Si_0.3Gy", "notINdmrs_Si_1Gy", "notINdmrs_X_1Gy",
                      "ALL_DMPs_Fe_0.3Gy", "ALL_DMPs_Fe_1Gy", "ALL_DMPs_FeX_1Gy", "ALL_DMPs_Si_0.3Gy", "ALL_DMPs_Si_1Gy", "ALL_DMPs_X_1Gy")

frequency <- lapply(lapply(frequency, function(x) data.frame(x)), function(x) {colnames(x) <- c("seqnames.x");x})

chromosome.probes.data.raw <- table(annote.dmps$chr)

freqB <- t(bind_rows(lapply(lapply(frequency, function(x) filter(x, !seqnames.x == "chrY" & !seqnames.x == "chrX")), function(x) table(x$seqnames.x))))
colnames(freqB) <- names(frequency)
freqB <- cbind(freqB, "raw"=chromosome.probes.data.raw[1:22])
freqB[is.na(freqB)] <- 0
chromosome.names <- rownames(freqB)
suma <- t(data.frame(colSums(freqB[,1:24])))
suma2 <- data.frame(freqB[,25]/sum(freqB[,25]))

freqB <- (freqB[,1:24]/ data.frame(suma[rep(1,22),])) / data.frame(suma2[,rep(1,24),])
rownames(freqB) <- chromosome.names
freqB <- cbind(freqB, "chr"=rownames(freqB))
colnames(freqB)
freqB.top23 <- lapply(lapply(lapply(list("DMRs_Fe_0.3Gy"=freqB[,c(1,25)],
                                         "DMRs_Fe_1Gy"=freqB[,c(2,25)],
                                         "DMRs_Si_0.3Gy"=freqB[,c(4,25)],
                                         "DMRs_Si_1Gy"=freqB[,c(5,25)],
                                         "DMRs_X_1Gy"=freqB[,c(6,25)],
                                         "DMPs_Fe_0.3Gy"=freqB[,c(7,25)],
                                         "DMPs_Fe_1Gy"=freqB[,c(8,25)],
                                         "DMPs_Si_0.3Gy"=freqB[,c(10,25)],
                                         "DMPs_Si_1Gy"=freqB[,c(11,25)],
                                         "DMPs_X_1Gy"=freqB[,c(12,25)],
                                         "notINdmrs_Fe_0.3Gy"=freqB[,c(13,25)],
                                         "notINdmrs_Fe_1Gy"=freqB[,c(14,25)],
                                         "notINdmrs_Si_0.3Gy"=freqB[,c(16,25)],
                                         "notINdmrs_Si_1Gy"=freqB[,c(17,25)],
                                         "notINdmrs_X_1Gy"=freqB[,c(18,25)],
                                         "ALL_DMPs_Fe_0.3Gy"=freqB[,c(19,25)],
                                         "ALL_DMPs_Fe_1Gy"=freqB[,c(20,25)],
                                         "ALL_DMPs_Fe(X)_1Gy"=freqB[,c(21,25)],
                                         "ALL_DMPs_Si_0.3Gy"=freqB[,c(22,25)],
                                         "ALL_DMPs_Si_1Gy"=freqB[,c(23,25)],
                                         "ALL_DMPs_X_1Gy"=freqB[,c(24,25)]), function(x) top_n(x, 23, x[,1])), function(x) arrange(x, desc(x[,1]))), function(x) {x$chr <- factor(x$chr, levels=x$chr, labels=x$chr);x})

#Figure_S2B
#Figure_S2C
#Figure_S2D
#Figure_S2E
#Figure_S2F
#Figure_S2G
f1 <- ggplot(melt(freqB.top23$ALL_DMPs_Fe_0.3Gy), aes(x=chr, y=value)) +
  geom_col() +
  scale_color_viridis(discrete=TRUE) +
  theme_classic() +
  labs(title="Normalized frequency of identified DMPs\non chromosomes for Fe 0.3Gy", x="Chromosomes", y="Normalized ratio of DMPs\nto control probes") +
  theme(
    axis.title.x = element_text(size = 12, face="bold"),
    axis.text.x = element_text(size = 12, angle = -315, vjust = 1, hjust=1),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 12, face="bold"),
    title = element_text(size = 12, face = "bold"))

f2 <- ggplot(melt(freqB.top23$ALL_DMPs_Fe_1Gy), aes(x=chr, y=value)) +
  geom_col() +
  scale_color_viridis(discrete=TRUE) +
  theme_classic() +
  labs(title="Normalized frequency of identified DMPs\non chromosomes for Fe 1Gy", x="Chromosomes", y="Normalized ratio of DMPs\nto control probes") +
  theme(
    axis.title.x = element_text(size = 12, face="bold"),
    axis.text.x = element_text(size = 12, angle = -315, vjust = 1, hjust=1),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 12, face="bold"),
    title = element_text(size = 12, face = "bold"))

f3 <- ggplot(melt(freqB.top23$`ALL_DMPs_Fe(X)_1Gy`), aes(x=chr, y=value)) +
  geom_col() +
  scale_color_viridis(discrete=TRUE) +
  theme_classic() +
  labs(title="Normalized frequency of identified DMPs\non chromosomes for Fe(X) 1Gy", x="Chromosomes", y="Normalized ratio of DMPs\nto control probes") +
  theme(
    axis.title.x = element_text(size = 12, face="bold"),
    axis.text.x = element_text(size = 12, angle = -315, vjust = 1, hjust=1),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 12, face="bold"),
    title = element_text(size = 12, face = "bold"))

f4 <- ggplot(melt(freqB.top23$ALL_DMPs_Si_0.3Gy), aes(x=chr, y=value)) +
  geom_col() +
  scale_color_viridis(discrete=TRUE) +
  theme_classic() +
  labs(title="Normalized frequency of identified DMPs\non chromosomes for Si 0.3Gy", x="Chromosomes", y="Normalized ratio of DMPs\nto control probes") +
  theme(
    axis.title.x = element_text(size = 12, face="bold"),
    axis.text.x = element_text(size = 12, angle = -315, vjust = 1, hjust=1),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 12, face="bold"),
    title = element_text(size = 12, face = "bold"))

f5 <- ggplot(melt(freqB.top23$ALL_DMPs_Si_1Gy), aes(x=chr, y=value)) +
  geom_col() +
  scale_color_viridis(discrete=TRUE) +
  theme_classic() +
  labs(title="Normalized frequency of identified DMPs\non chromosomes for Si 1Gy", x="Chromosomes", y="Normalized ratio of DMPs\nto control probes") +
  theme(
    axis.title.x = element_text(size = 12, face="bold"),
    axis.text.x = element_text(size = 12, angle = -315, vjust = 1, hjust=1),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 12, face="bold"),
    title = element_text(size = 12, face = "bold"))

f6 <- ggplot(melt(freqB.top23$ALL_DMPs_X_1Gy), aes(x=chr, y=value)) +
  geom_col() +
  scale_color_viridis(discrete=TRUE) +
  theme_classic() +
  labs(title="Normalized frequency of identified DMPs\non chromosomes for X 1Gy", x="Chromosomes", y="Normalized ratio of DMPs\nto control probes") +
  theme(
    axis.title.x = element_text(size = 12, face="bold"),
    axis.text.x = element_text(size = 12, angle = -315, vjust = 1, hjust=1),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 12, face="bold"),
    title = element_text(size = 12, face = "bold"))


# 18. Assign Hi-C compartments to the data ####
compartments <- read_bigwig("4DNFIHM89EGL.bw")

All.1.change.ov <- lapply(lapply(All.samples.not.detailed.1, function(x) makeGRangesFromDataFrame(x, keep.extra.columns = T)), function(x) findOverlaps(compartments, x))

All.1.change.query <- lapply(All.1.change.ov, function(x) compartments[queryHits(x)])
All.1.change.subject <- mapply(function(x,y) x[subjectHits(y)], x=lapply(All.samples.not.detailed.1, function(x) makeGRangesFromDataFrame(x, keep.extra.columns = T)), y=All.1.change.ov)

lengths(All.1.change.query) == lengths(All.1.change.subject)

# Give IDs and filter
All.1.change.query <- lapply(All.1.change.query, function(x) {x$ID <- paste0('E_', 1:length(x)); x})
All.1.change.subject <- lapply(All.1.change.subject, function(x) {x$ID <- paste0('E_', 1:length(x)); x})

# Annotate Probes
All.1.change.query <- lapply(All.1.change.query, function(x) as.data.frame(x))

All.1.change.subject <- lapply(All.1.change.subject, function(x) as.data.frame(x))

All.1.change <- purrr::map2(All.1.change.query, All.1.change.subject, left_join, by="ID")
All.1.change <- lapply(All.1.change, function(x) { x[,149] <- unite(x[1:3], cord, remove=F, sep="_");x})

#
All.2.change.ov <- lapply(lapply(All.samples.not.detailed.2, function(x) makeGRangesFromDataFrame(x, keep.extra.columns = T)), function(x) findOverlaps(compartments, x))

All.2.change.query <- lapply(All.2.change.ov, function(x) compartments[queryHits(x)])
All.2.change.subject <- mapply(function(x,y) x[subjectHits(y)], x=lapply(All.samples.not.detailed.2, function(x) makeGRangesFromDataFrame(x, keep.extra.columns = T)), y=All.2.change.ov)

lengths(All.2.change.query) == lengths(All.2.change.subject)

# Give IDs and filter
All.2.change.query <- lapply(All.2.change.query, function(x) {x$ID <- paste0('E_', 1:length(x)); x})
All.2.change.subject <- lapply(All.2.change.subject, function(x) {x$ID <- paste0('E_', 1:length(x)); x})

# Annotate Probes
All.2.change.query <- lapply(All.2.change.query, function(x) as.data.frame(x))

All.2.change.subject <- lapply(All.2.change.subject, function(x) as.data.frame(x))

All.2.change <- purrr::map2(All.2.change.query, All.2.change.subject, left_join, by="ID")
All.2.change <- lapply(All.2.change, function(x) { x[,111] <- unite(x[1:3], cord, remove=F, sep="_");x})

# 19. Data preparation for Hi-C analysis ####
annote.dmps.GRanges <- annote.dmps[,c(1,2,2,3,4)]
colnames(annote.dmps.GRanges) <- c("seqnames", "start", "end", "strand", "Name")
annote.dmps.GRanges <- left_join(annote.dmps.GRanges, data.frame("Name"=rownames(data.raw.M.means), data.raw.M.means), by="Name")
annote.dmps.GRanges <- makeGRangesFromDataFrame(annote.dmps.GRanges[complete.cases(annote.dmps.GRanges),], keep.extra.columns = T)
annote.dmps.GRanges.ov <- findOverlaps(compartments, annote.dmps.GRanges)

annote.dmps.GRanges.query <- compartments[queryHits(annote.dmps.GRanges.ov)]
annote.dmps.GRanges.subject <- annote.dmps.GRanges[subjectHits(annote.dmps.GRanges.ov)]

annote.dmps.GRanges.query$ID <- paste0('E_', 1:length(annote.dmps.GRanges.query))
annote.dmps.GRanges.subject$ID <- paste0('E_', 1:length(annote.dmps.GRanges.subject))

annote.dmps.GRanges.query <- as.data.frame(annote.dmps.GRanges.query)
annote.dmps.GRanges.subject <- as.data.frame(annote.dmps.GRanges.subject)

#data.raw comp for normalization
annote.dmps.joined <- left_join(annote.dmps.GRanges.query, annote.dmps.GRanges.subject, by="ID")
annote.dmps.joined <- filter(annote.dmps.joined, !seqnames.x == "chrY")
annote.dmps.joined$control <- rowMeans(annote.dmps.joined[,c(14,26,38)])

#Joined frames
e <- lapply(All.1.change, function(x) filter(x, !seqnames.x == "chrY")) #DMRs, DMPs inside
f <- lapply(lapply(All.2.change, function(x) filter(x, !seqnames.x == "chrY")), function(x) {x$nProbes <- c(1);x}) #DMPs outside

pen.allDMPs.notdetailed <- append(append(append(lapply(e[1:6], function(x) x[,c(1:6,8:12,61,82:83,89:98,105:111,116:147,104,148)]),
                                    lapply(e[7:12], function(x) x[,c(1:6,49:53,61,82:83,89:98,105:111,116:147,104,148)])), 
                             Map(rbind, lapply(e[7:12], function(x) {colnames(x)[49:53] <- gsub("\\.y.y",".y", colnames(x)[49:53]); x[,c(1:6,49:53,61,82:83,89:98,105:111,116:147,104,148)]}), 
                                 lapply(f, function(x) x[,c(1:6,8:12,28,49:50,56:65,104:110,71:102,103,112)]))), lapply(f, function(x) x[,c(1:6,8:12,28,49:50,56:65,104:110,71:102,103,112)]))
names(pen.allDMPs.notdetailed) <- c("DMRs_Fe_0.3Gy", "DMRs_Fe_1Gy", "DMRs_FeX_1Gy", "DMRs_Si_0.3Gy", "DMRs_Si_1Gy", "DMRs_X_1Gy",
                                    "DMPs_Fe_0.3Gy", "DMPs_Fe_1Gy", "DMPs_FeX_1Gy", "DMPs_Si_0.3Gy", "DMPs_Si_1Gy", "DMPs_X_1Gy",
                                    "ALL_DMPs_Fe_0.3Gy", "ALL_DMPs_Fe_1Gy", "ALL_DMPs_FeX_1Gy", "ALL_DMPs_Si_0.3Gy", "ALL_DMPs_Si_1Gy", "ALL_DMPs_X_1Gy",
                                    "notINdmrs_Fe_0.3Gy", "notINdmrs_Fe_1Gy", "notINdmrs_FeX_1Gy", "notINdmrs_Si_0.3Gy", "notINdmrs_Si_1Gy", "notINdmrs_X_1Gy")

# data.raw probes regardless of methylation status
quantile(annote.dmps.joined$score, probs = c(0, 0.2, 0.4, 0.6, 0.8, 1), na.rm=T)

pen.5.annote.dmps.joined <- list("L1"=filter(annote.dmps.joined, score>=-2.33444190 & score<=-0.45549509),
                                  "L2"=filter(annote.dmps.joined, score>-0.45549509 & score<=0.03823841),
                                  "L3"=filter(annote.dmps.joined, score>0.03823841 & score<=0.45039466),
                                  "L4"=filter(annote.dmps.joined, score>0.45039466 & score<=0.77691591),
                                  "L5"=filter(annote.dmps.joined, score>0.77691591 & score<=1.99093556))

pen.5.allDMPs.notdetailed <- lapply(pen.allDMPs.notdetailed, function(x) list("L1"=filter(x, score>=-2.33444190 & score<=-0.45549509),
                                 "L2"=filter(x, score>-0.45549509 & score<=0.03823841),
                                 "L3"=filter(x, score>0.03823841 & score<=0.45039466),
                                 "L4"=filter(x, score>0.45039466 & score<=0.77691591),
                                 "L5"=filter(x, score>0.77691591 & score<=1.99093556)))

# 20. Order and plot data.raw chromosome position ####
annote.dmps.ranges <- annote.dmps[,c(1,2,2,3:61)]
annote.dmps.ranges <- annote.dmps.ranges[!is.na(annote.dmps.ranges$pos),]
colnames(annote.dmps.ranges)[2:3] <- c("start", "end")
annote.dmps.ranges <- makeGRangesFromDataFrame(annote.dmps.ranges, keep.extra.columns = T)
annote.dmps.ranges[1:10,1:6]

# Overlap data.raw probes
data.raw.ov <- findOverlaps(compartments, annote.dmps.ranges)

data.raw.query <- compartments[queryHits(data.raw.ov)]
data.raw.subject <- annote.dmps.ranges[subjectHits(data.raw.ov)]

# Give IDs and filter
data.raw.query$ID <- paste0('E_', 1:length(data.raw.query))
data.raw.subject$ID <- paste0('E_', 1:length(data.raw.subject))

# Annotate DMRs
data.raw.joined <- left_join(as.data.frame(data.raw.query), as.data.frame(data.raw.subject), by="ID")
data.raw.joined$cord <- unite(data.raw.joined[,1:3], cord, remove=T, sep="_")

data.raw.joined <- mutate(data.raw.joined, Compartment =  case_when(
  (score>=-2.33444190 & score<=-0.45549509) ~ "L1",
  (score>-0.45549509 & score<=0.03823841) ~ "L2",
  (score>0.03823841 & score<=0.45039466) ~ "L3",
  (score>0.45039466 & score<=0.77691591) ~ "L4",
  (score>0.77691591 & score<=1.99093556) ~ "L5",
  (is.nan(score) == TRUE) ~ "None"))

test <- filter(data.raw.joined, score <=100)
test <- filter(test, !seqnames.x == "chrY" & !seqnames.x == "chrX")
vec <- names(sort(tapply(test$score, test$seqnames.x, FUN=min), decreasing = F))
test$seqnames.x <- factor(test$seqnames.x, levels=vec, labels=vec)
colnames(test)
test.all <- data.frame("seqnames.x"=c("chr.all"),"value"=test[,6])
test.1001 <- rbind(melt(filter(test[,c(1,6)], seqnames.x == "chr14" | seqnames.x == "chr13" | seqnames.x == "chr10" | seqnames.x == "chr4" | seqnames.x == "chr1"))[,-2], test.all)

#Figure_2B
ggplot(test.1001, aes(x=as.factor(seqnames.x), y=value, fill=seqnames.x)) +
  geom_violin(width=1) +
  scale_fill_viridis(discrete=T, name="") +
  ggtitle("") + 
  theme_classic() +
  theme(legend.position="none",
        axis.title.x = element_text(size = 12, face="bold"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 12, face="bold"),
        title = element_text(size = 12, face = "bold")) +
  xlab("") +
  ylab("Hi-C score (Spatial Position)") #+ #It takes 5-10 minutes to generate this with all layers lines
geom_hline(yintercept=-2.33444190, linetype="dashed", color = "red") +#
  geom_text(aes(0,-2.33444190, label="L1", hjust = -1, vjust= -1.5, size=10)) +
  geom_hline(yintercept=-0.45549509, linetype="dashed", color = "red") +
  geom_text(aes(0,-0.45549509, label="L2", hjust = -1, vjust= -1.5, size=10)) +
  geom_hline(yintercept=0.03823841, linetype="dashed", color = "red") +
  geom_text(aes(0,0.03823841, label="L3", hjust = -1, vjust= -1.5, size=10)) +
  geom_hline(yintercept=0.45039466, linetype="dashed", color = "red") +
  geom_text(aes(0,0.45039466, label="L4", hjust = -1, vjust= -1.5, size=10)) +
  geom_hline(yintercept=0.77691591, linetype="dashed", color = "red") +
  geom_text(aes(0,0.77691591, label="L5", hjust = -1, vjust= -1.5, size=10)) +
  geom_hline(yintercept=1.99093556, linetype="dashed", color = "red")

#Figure_S2A
test.1002 <- rbind(melt(test[,c(1,6)])[,-2], test.all)

ggplot(test.1002, aes(x=as.factor(seqnames.x), y=value, fill=seqnames.x)) +
  geom_violin(width=1) +
  scale_fill_viridis(discrete=T, name="") +
  ggtitle("") + 
  theme_classic() +
  theme(legend.position="none",
        axis.title.x = element_text(size = 12, face="bold"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 12, face="bold"),
        title = element_text(size = 12, face = "bold")) +
  xlab("") +
  ylab("Hi-C score (Spatial Position)")


# 21. Chromosomal frequencies in each nucleus layer ####

#data.raw
chr <- data.raw.joined[,c(1,6,72)] %>%
  filter(!Compartment == "None")

# All not detailed
allDMPs.notdetailed <- lapply(pen.allDMPs.notdetailed, function(x) x[,c(1,6)])

allDMPs.notdetailed <- lapply(lapply(lapply(allDMPs.notdetailed, function(x) mutate(x, Compartment =  case_when(
  (score>=-2.33444190 & score<=-0.45549509) ~ "L1",
  (score>-0.45549509 & score<=0.03823841) ~ "L2",
  (score>0.03823841 & score<=0.45039466) ~ "L3",
  (score>0.45039466 & score<=0.77691591) ~ "L4",
  (score>0.77691591 & score<=1.99093556) ~ "L5",
  (is.nan(score) == TRUE) ~ "None"))), function(x) x[,c(1,2,3)]), function(x) filter(x, !Compartment == "None"))


integrals.all.raw <- t(table(chr[,-2]))[,1:22]
suma.raw <- colSums(integrals.all.raw)

df110 <- integrals.all.raw/data.frame(t(suma.raw)[rep(1,5),])
rownames(df110) <- c("A", "B", "C", "D", "E")
df110$Comp <- rownames(df110)

#Figure_2D
#Figure_S2I
#Figure_S2J

#Fe
ratio.3 <- t(table(allDMPs.notdetailed$ALL_DMPs_Fe_1Gy[,-2])) / t(table(chr[,-2])) #To plot different exposures change a list object
ratio.3 <- filter(filter(data.frame(ratio.3), Freq > 0), !seqnames.x == "chrX")

ratio.4 <- allDMPs.notdetailed$ALL_DMPs_Fe_1Gy %>%  #To plot different exposures change a list object
  filter(!seqnames.x == "chrX") %>%
  group_by(seqnames.x, Compartment) %>%
  summarise(mean=median(score))

test444 <- data.frame(ratio.4, melt(ratio.3), filter(melt(df110), value > 0))[,c(1:3,7,10)]
colnames(test444) <- c("seqnames", "Compartment", "mean", "ratio", "pctChr")
test444$seqnames <- gsub("^chr", "", test444$seqnames)

ggplot(test444, aes(x=mean, y=ratio, colour=as.factor(Compartment), label=as.factor(seqnames))) +
  geom_text(alpha=1, size=5) +
  theme_classic() +
  labs(title="nDMPs / Ratio of nRaw probes in a layer per chr\nto Median Hi-C score per nucleus layers and per chr\nFe 1Gy", x="Median Hi-C score per\nnucleus layers and per chr", 
       y="nDMPs / nRaw probes in a layer per chr", colour="Layer") +
  theme(
    axis.title.x = element_text(size = 12, face="bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 12, face="bold"),
    title = element_text(size = 12, face = "bold"))

#Figure_S2K
#Figure_S2L
#Figure_S2M

#SiX
ratio.3 <- t(table(allDMPs.notdetailed$ALL_DMPs_X_1Gy[,-2])) / t(table(chr[,-2]))  #To plot different exposures change a list object
ratio.3 <- filter(filter(data.frame(ratio.3), Freq >= 0), !seqnames.x == "chrX" & !seqnames.x == "chrY")

ratio.4 <- allDMPs.notdetailed$ALL_DMPs_X_1Gy %>%  #To plot different exposures change a list object
  filter(!seqnames.x == "chrX" & !seqnames.x == "chrY") %>%
  group_by(seqnames.x, Compartment) %>%
  summarise(mean=median(score))

test444 <- data.frame(melt(ratio.3), filter(melt(df110), value > 0))[,c(1:2,4,7)]
colnames(test444) <- c("Compartment", "seqnames", "ratio", "pctChr")
test444$seqnames <- gsub("^chr", "", test444$seqnames)

test444 <- filter(test444, ratio > 0)
test444 <- cbind(test444, ratio.4)[,c(1:4,7)]

ggplot(test444, aes(x=mean, y=ratio, colour=as.factor(Compartment), label=as.factor(seqnames))) +
  geom_text(alpha=1, size=5) +
  theme_classic() +
  labs(title="nDMRs / Ratio of nRaw probes in a layer per chr\nto Median Hi-C score per nucleus layers and per chr", x="Median Hi-C score per\nnucleus layers and per chr", 
       y="nDMRs / nRaw probes in a layer per chr", colour="Layer") +
  theme(
    axis.title.x = element_text(size = 12, face="bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 12, face="bold"),
    title = element_text(size = 12, face = "bold"))

# 22. DMPs frequency and DNA methylation change within nucleus layers ####

#normalization
normal.pen <- data.frame("baseline"=t(data.frame(bind_rows(lapply(pen.5.annote.dmps.joined, function(x) nrow(x))))),
                         t(data.frame(bind_rows(lapply(pen.5.allDMPs.notdetailed, function(x) lapply(x, function(x) nrow(x)))))))

colnames(normal.pen) <- c("baseline", names(pen.5.allDMPs.notdetailed))
normal.pen.pct <- apply(normal.pen,2,function(x){x/sum(x)})

colnames(normal.pen.pct)
df102 <- normal.pen.pct[,c(1,14:19)]
colnames(df102) <- c("control", "Fe 0.3Gy", "Fe 1Gy", "Fe(X) 1Gy", "Si 0.3Gy", "Si 1Gy", "X 1Gy")

#Figure_2C
ggplot(melt(df102[,c(1,3,6)]), aes(x=Var1, y=value, col=as.factor(Var2))) +
  geom_point(size=4) +
  scale_color_viridis(discrete=TRUE) +
  geom_line(aes(group=interaction(Var2)), size=2, alpha=0.5) +
  theme_classic() +
  labs(title="Frequency of DMPs across nucleus layers", x="Nucleus layers", y="nDMPs (%)", color="Sample") +
  theme(
    axis.title.x = element_text(size = 10, face="bold"),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 10, face="bold"),
    title = element_text(size = 10, face = "bold"))

#Figure_S2H
ggplot(melt(df102), aes(x=Var1, y=value, col=as.factor(Var2))) +
  geom_point(size=4) +
  scale_color_viridis(discrete=TRUE) +
  geom_line(aes(group=interaction(Var2)), size=2, alpha=0.5) +
  theme_classic() +
  labs(title="Frequency of DMPs across nucleus layers", x="Nucleus layers", y="nDMPs (%)", color="Sample") +
  theme(
    axis.title.x = element_text(size = 10, face="bold"),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 10, face="bold"),
    title = element_text(size = 10, face = "bold"))

#Abs change
mean.meth.change.pen.5.Fe <- data.frame(t(data.frame(bind_rows(lapply(pen.5.allDMPs.notdetailed[c(1:3,7:9,13:15,19:21)],
                                                                               function(x) lapply(x, function(x) {abs(mean(x[,30] - rowMeans(x[,c(32,44,56)])))}))))), "comp"=rownames(normal.pen.pct))

colnames(mean.meth.change.pen.5.Fe) <- c(names(pen.5.allDMPs.notdetailed[c(1:3,7:9,13:15,19:21)]), "comp")


mean.meth.change.pen.5.SiX <- data.frame(t(data.frame(bind_rows(lapply(pen.5.allDMPs.notdetailed[c(4:6,10:12,16:18,22:24)], 
                                                                                function(x) lapply(x, function(x) {abs(mean(x[,31] - rowMeans(x[,c(32,44,56)])))}))))), "comp"=rownames(normal.pen.pct))

colnames(mean.meth.change.pen.5.SiX) <- c(names(pen.5.allDMPs.notdetailed[c(4:6,10:12,16:18,22:24)]), "comp")

fig1 <- cbind(mean.meth.change.pen.5.Fe[c(7:9)], mean.meth.change.pen.5.SiX[c(7:9,13)])
colnames(fig1) <- c("Fe 0.3Gy", "Fe 1Gy", "FeX 1Gy", "Si 0.3Gy", "Si 1Gy", "X 1Gy", "comp")

#Figure_2E
ggplot(melt(fig1[,c(2,5,7)]), aes(x=comp, y=value, col=as.factor(variable))) +
  geom_point(size=4) +
  scale_color_viridis(discrete=TRUE) +
  geom_line(aes(group=interaction(variable)), size=2, alpha=0.5) +
  theme_classic() +
  labs(title="", x="Nucleus layer", y="Absolute mean methylation change versus control", color="Sample") +
  theme(
    axis.title.x = element_text(size = 12, face="bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 12, face="bold"),
    title = element_text(size = 12, face = "bold"))

#Figure_S2N
ggplot(melt(fig1[,c(1,3,4,6,7)]), aes(x=comp, y=value, col=as.factor(variable))) +
  geom_point(size=4) +
  scale_color_viridis(discrete=TRUE) +
  geom_line(aes(group=interaction(variable)), size=2, alpha=0.5) +
  theme_classic() +
  labs(title="", x="Nucleus layer", y="Absolute mean methylation change versus control", color="Sample") +
  theme(
    axis.title.x = element_text(size = 12, face="bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 12, face="bold"),
    title = element_text(size = 12, face = "bold"))

# 23. Load and prepare histone modifications datasets for methylation analysis ####
# Data.raw probes
dr <- left_join(data.raw.joined, data.frame("Name"=rownames(data.raw.M.means), data.raw.M.means), by="Name")
data.raw.GRanges <- dr[,c(8:12,1:7,13:ncol(dr))]
colnames(data.raw.GRanges)[1:5] <- gsub("\\.y","", colnames(data.raw.GRanges)[1:5])

# Joined probes from DMRs, outside DMRs, both inside and outside and data.raw for increased/decrease 
probes.for.histones.notdetailed <- lapply(lapply(append(lapply(lapply(pen.allDMPs.notdetailed, function(x) x[,c(7:11,1:6,12:24,32:ncol(x))]),
                                                            function(x) {colnames(x)[c(1:5,57)] <- c("seqnames", "start", "end","width", "strand", "Name"); colnames(x)[c(24:56)] <- gsub("\\.dmp", "", colnames(x)[c(24:56)]);x}),
                                                     list(data.raw=data.raw.GRanges[,c(1:5,11,6:10,28,49:50,56:65,73:104,13)])), function(x) makeGRangesFromDataFrame(x, keep.extra.columns = T)),
                                       function(x) {x %>% mutate(Compartment = case_when(
                                         (score>=-2.334442 & score<=-0.45549509) ~ "A",
                                         (score>-0.45549509 & score<=0.03823841) ~ "B",
                                         (score>0.03823841 & score<=0.45039466) ~ "C",
                                         (score>0.45039466 & score<=0.77691591) ~ "D",
                                         (score>0.77691591 & score<=1.99093556) ~ "E",
                                         (is.nan(score) == TRUE) ~ "None"))})

# Load different histone modifications 
H3K9me3 <- as.data.frame(import("/mnt/work1/adrian/Nasa/ChIP/GSM2481428_HBE_H3K9me3_peaks.bed"))
colnames(H3K9me3) <- c("chr", "start", "end")
H3K9me3 <- makeGRangesFromDataFrame(H3K9me3)
H3K9me3 <- filter(as.data.frame(H3K9me3), width <= 2000 & width >= 300)
H3K9me3 <- makeGRangesFromDataFrame(H3K9me3)

H3K27me3 <- as.data.frame(import("/mnt/work1/adrian/Nasa/ChIP/GSM2481424_HBE_H3K27me3_peaks.bed"))[,1:3]
colnames(H3K27me3) <- c("chr", "start", "end")
H3K27me3 <- makeGRangesFromDataFrame(H3K27me3)
H3K27me3 <- filter(as.data.frame(H3K27me3), width <= 2000 & width >= 300)
H3K27me3 <- makeGRangesFromDataFrame(H3K27me3)

H3K4me3 <- as.data.frame(import("/mnt/work1/adrian/Nasa/ChIP/GSM2481422_HBE_H3K4me3_peaks.bed"))[,1:3]
colnames(H3K4me3) <- c("chr", "start", "end")
H3K4me3 <- makeGRangesFromDataFrame(H3K4me3)
H3K4me3 <- filter(as.data.frame(H3K4me3), width <= 2000 & width >= 300)
H3K4me3 <- makeGRangesFromDataFrame(H3K4me3)

#Joined replicates
#H3K27ac1 <- read_bigwig("/mnt/work1/adrian/Nasa/ChIP/GSM4271077_16HBE14o-_H3K27ac_rep2_AHCWR40.bigwig")
#H3K27ac2 <- read_bigwig("/mnt/work1/adrian/Nasa/ChIP/GSM3892733_16HBE14o-_H3K27Ac.bigwig")

#H3K27ac3 <- rbind(as.data.frame(H3K27ac1), as.data.frame(H3K27ac2))

#H3K27ac4 <- filter(H3K27ac3, score >= 10)

#H3K27ac5 <- reduce(makeGRangesFromDataFrame(H3K27ac4))

#H3K27ac6 <- filter(H3K27ac5, width <= 2000 & width >= 300)

#export(data.frame(H3K27ac6),"/mnt/work1/adrian/Nasa/ChIP/H3K27ac_ChIP.bed", "bed")

H3K27ac <- import("/mnt/work1/adrian/Nasa/ChIP/H3K27ac_ChIP.bed")
H3K27ac <- makeGRangesFromDataFrame(H3K27ac)

# Find overlaps between histone modifications
overlap.with.probes <- function(x,y) { #x=hist picks, y=probes to overlap
  z <- findOverlaps(x,y)
  
  x.query <- x[queryHits(z)]
  y.subjectHits <- y[subjectHits(z)]
  
  x.query$ID <- paste0('E_', 1:length(x.query))
  y.subjectHits$ID <- paste0('E_', 1:length(y.subjectHits))
  
  x.query <- as.data.frame(x.query)
  y.subjectHits <- as.data.frame(y.subjectHits)
  
  x.final <- left_join(x.query, y.subjectHits, by="ID")
  
  x.final$cord <- unite(x.final[1:3], cord, sep="_", remove=T)
    
  x.final.grouped <- x.final %>%
                                group_by(cord) %>%
                                mutate(f0=mean(f0)) %>%
                                mutate(f0.3=mean(f0.3)) %>%
                                mutate(f1=mean(f1)) %>%
                                mutate(f0.14=mean(f0.14)) %>%
                                mutate(f0.3.14=mean(f0.3.14)) %>%
                                mutate(f1.14=mean(f1.14)) %>%
                                mutate(f0.22=mean(f0.22)) %>%
                                mutate(f0.3.22=mean(f0.3.22)) %>%    
                                mutate(f1.22=mean(f1.22)) %>%    
                                mutate(f0.53=mean(f0.53)) %>%    
                                mutate(f0.3.53=mean(f0.3.53)) %>%    
                                mutate(f1.53=mean(f1.53)) %>%    
                                mutate(s0=mean(s0)) %>%   
                                mutate(s0.3=mean(s0.3)) %>%    
                                mutate(s1=mean(s1)) %>%    
                                mutate(s0.13=mean(s0.13)) %>% 
                                mutate(s0.3.13=mean(s0.3.13)) %>%
                                mutate(s1.13=mean(s1.13)) %>%
                                mutate(s0.21=mean(s0.21)) %>%
                                mutate(s0.3.21=mean(s0.3.21)) %>%
                                mutate(s1.21=mean(s1.21)) %>%
                                mutate(s0.62=mean(s0.62)) %>%
                                mutate(s0.3.62=mean(s0.3.62)) %>%
                                mutate(s1.62=mean(s1.62)) %>%    
                                mutate(x0=mean(x0)) %>%    
                                mutate(x1=mean(x1)) %>%    
                                mutate(x0.13=mean(x0.13)) %>%    
                                mutate(x1.13=mean(x1.13)) %>%    
                                mutate(x0.22=mean(x0.22)) %>%    
                                mutate(x1.22=mean(x1.22)) %>%    
                                mutate(x0.62=mean(x0.62)) %>%
                                mutate(x1.62=mean(x1.62)) %>%
                                distinct(cord, .keep_all=T) %>%
                                data.frame()
    x.final.grouped
  }

probes.for.histones.notdetailed.H3K9me3 <- lapply(probes.for.histones.notdetailed, function(x) overlap.with.probes(H3K9me3, x))
probes.for.histones.notdetailed.H3K27me3 <- lapply(probes.for.histones.notdetailed, function(x) overlap.with.probes(H3K27me3, x))
probes.for.histones.notdetailed.H3K4me3 <- lapply(probes.for.histones.notdetailed, function(x) overlap.with.probes(H3K4me3, x))
probes.for.histones.notdetailed.H3K27ac <- lapply(probes.for.histones.notdetailed, function(x) overlap.with.probes(H3K27ac, x))


# 24. Assign histone modifications to Hi-C layers ####
#Fe
comp.hist.Fe <- list(H3K9me3=probes.for.histones.notdetailed.H3K9me3[c(1:3,7:9,13:15,19:21,25)],
                     H3K27me3=probes.for.histones.notdetailed.H3K27me3[c(1:3,7:9,13:15,19:21,25)],
                     H3K4me3=probes.for.histones.notdetailed.H3K4me3[c(1:3,7:9,13:15,19:21,25)],
                     H3K27ac=probes.for.histones.notdetailed.H3K27ac[c(1:3,7:9,13:15,19:21,25)])

freq.hist.Fe <- data.frame(bind_rows(lapply(comp.hist.Fe, function(x) lapply(x, function(x) nrow(x)))))
rownames(freq.hist.Fe) <- names(comp.hist.Fe)
freq.hist.Fe$mod <- c("H3K9me3", "H3K27me3", "H3K4me3", "H3K27ac")

#Si
comp.hist.Si <- list(H3K9me3=probes.for.histones.notdetailed.H3K9me3[c(4:6,10:12,16:18,22:24,25)],
                     H3K27me3=probes.for.histones.notdetailed.H3K27me3[c(4:6,10:12,16:18,22:24,25)],
                     H3K4me3=probes.for.histones.notdetailed.H3K4me3[c(4:6,10:12,16:18,22:24,25)],
                     H3K27ac=probes.for.histones.notdetailed.H3K27ac[c(4:6,10:12,16:18,22:24,25)])

freq.hist.Si <- data.frame(bind_rows(lapply(comp.hist.Si, function(x) lapply(x, function(x) nrow(x)))))
rownames(freq.hist.Si) <- names(comp.hist.Si)
freq.hist.Si$mod <- c("H3K9me3", "H3K27me3", "H3K4me3", "H3K27ac")

df117 <- cbind(freq.hist.Fe, freq.hist.Si[,1:12])

#Fe
pen.5.comp.hist.Fe <- lapply(comp.hist.Fe, function(x) lapply(x, function(x) list("L1"=filter(x, score>=-2.33444190 & score<=-0.45549509),
                                                                                   "L2"=filter(x, score>-0.45549509 & score<=0.03823841),
                                                                                   "L3"=filter(x, score>0.03823841 & score<=0.45039466),
                                                                                   "L4"=filter(x, score>0.45039466 & score<=0.77691591),
                                                                                   "L5"=filter(x, score>0.77691591 & score<=1.99093556))))

pen.5.comp.hist.Fe <- lapply(pen.5.comp.hist.Fe, function(x) lapply(x, function(x) lapply(x, function(x) {x[is.na(x)] <- 0;x})))

normal.pen.Fe <- data.frame(data.frame(bind_rows(lapply(pen.5.comp.hist.Fe, function(x) lapply(x, function(x) lapply(x, function(x) nrow(x)))))))
normal.pen.Fe$mod <- c(rep("H3K9me3",5),rep("H3K27me3",5),rep("H3K4me3",5),rep("H3K27ac", 5))
normal.pen.Fe$comp <- c(rep(c("L1","L2","L3","L4","L5"),4))
normal.pen.Fe <- cbind(apply(normal.pen.Fe[,1:13], 2, as.numeric), normal.pen.Fe[,14:15])
normal.pen.Fe[1:5,1:13] <- normal.pen.Fe[1:5,1:13]/ t(data.frame(colSums(normal.pen.Fe[1:5,1:13]))[,rep(1,5)])
normal.pen.Fe[6:10,1:13] <- normal.pen.Fe[6:10,1:13]/ t(data.frame(colSums(normal.pen.Fe[6:10,1:13]))[,rep(1,5)])
normal.pen.Fe[11:15,1:13] <- normal.pen.Fe[11:15,1:13]/ t(data.frame(colSums(normal.pen.Fe[11:15,1:13]))[,rep(1,5)])
normal.pen.Fe[16:20,1:13] <- normal.pen.Fe[16:20,1:13]/ t(data.frame(colSums(normal.pen.Fe[16:20,1:13]))[,rep(1,5)])

#Si
pen.5.comp.hist.Si <- lapply(comp.hist.Si, function(x) lapply(x, function(x) list("L1"=filter(x, score>=-2.33444190 & score<=-0.45549509),
                                                                                  "L2"=filter(x, score>-0.45549509 & score<=0.03823841),
                                                                                  "L3"=filter(x, score>0.03823841 & score<=0.45039466),
                                                                                  "L4"=filter(x, score>0.45039466 & score<=0.77691591),
                                                                                  "L5"=filter(x, score>0.77691591 & score<=1.99093556))))

pen.5.comp.hist.Si <- lapply(pen.5.comp.hist.Si, function(x) lapply(x, function(x) lapply(x, function(x) {x[is.na(x)] <- 0;x})))

normal.pen.Si <- data.frame(data.frame(bind_rows(lapply(pen.5.comp.hist.Si, function(x) lapply(x, function(x) lapply(x, function(x) nrow(x)))))))
normal.pen.Si$mod <- c(rep("H3K9me3",5),rep("H3K27me3",5),rep("H3K4me3",5), rep("H3K27ac",5))
normal.pen.Si$comp <- c(rep(c("L1","L2","L3","L4","L5"),4))
normal.pen.Si <- cbind(apply(normal.pen.Si[,1:13], 2, as.numeric), normal.pen.Si[,14:15])

normal.pen.Si[1:5,1:10] <- normal.pen.Si[1:5,1:10]/ t(data.frame(colSums(normal.pen.Si[1:5,1:10]))[,rep(1,5)])
normal.pen.Si[6:10,1:10] <- normal.pen.Si[6:10,1:10]/ t(data.frame(colSums(normal.pen.Si[6:10,1:10]))[,rep(1,5)])
normal.pen.Si[11:15,1:10] <- normal.pen.Si[11:15,1:10]/ t(data.frame(colSums(normal.pen.Si[11:15,1:10]))[,rep(1,5)])
normal.pen.Si[16:20,1:10] <- normal.pen.Si[16:20,1:10]/ t(data.frame(colSums(normal.pen.Si[16:20,1:10]))[,rep(1,5)])

df120 <- cbind(normal.pen.Fe[,1:13], normal.pen.Si)

pick <- df120[,c(4:6,17:19,27,28)]
colnames(pick) <- c("Fe 0.3Gy", "Fe 1Gy", "Fe(X) 1Gy", "Si 0.3Gy", "Si 1Gy", "X 1Gy", "mod", "comp")

#reg in l/reg
all.reg <- list("H3K27ac"=H3K27ac,
                "H3K27me3"=H3K27me3,
                "H3K4me3"=H3K4me3,
                "H3K9me3"=H3K9me3)
all.reg <- lapply(all.reg, function(x) overlap.reg(x, compartments))

overlap.reg <- function(x,y){
z <- findOverlaps(x, y)

qHits <- x[queryHits(z)]
sHits <- y[subjectHits(z)]

qHits$ID <- paste0('E_', 1:length(qHits))
sHits$ID <- paste0('E_', 1:length(sHits))

qHits <- as.data.frame(qHits)
sHits <- as.data.frame(sHits)

a <- left_join(qHits, sHits, by="ID")
a$cord <- unite(a[,1:3], cord, remove=T, sep="_")
a
}

all.reg.comp <- lapply(all.reg, function(x) list("L1"=filter(x, score>=-2.33444190 & score<=-0.45549509),
                                                             "L2"=filter(x, score>-0.45549509 & score<=0.03823841),
                                                             "L3"=filter(x, score>0.03823841 & score<=0.45039466),
                                                             "L4"=filter(x, score>0.45039466 & score<=0.77691591),
                                                             "L5"=filter(x, score>0.77691591 & score<=1.99093556)))


all.reg.comp.normal <- data.frame(bind_rows(lapply(all.reg.comp, function(x) lapply(x, function(x) nrow(x)))))
rownames(all.reg.comp.normal) <- names(all.reg.comp)
test1003 <- all.reg.comp.normal/rowSums(all.reg.comp.normal)
test1003$mod <- c(rownames(test1003))

#Figure_3B
ggplot(melt(test1003), aes(x=variable, y=value, col=as.factor(mod))) +
  geom_point(size=4) +
  scale_color_viridis(discrete=TRUE) +
  geom_line(aes(group=interaction(mod)), size=2, alpha=0.5) +
  theme_classic() +
  labs(title="Frequency of Histone Modification Peaks", x="Nucleus layers", y="nHistone Modification Picks (%)", color="Histone\nModification") +
  theme(
    axis.title.x = element_text(size = 12, face="bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 12, face="bold"),
    title = element_text(size = 12, face = "bold"))

matrix(as.vector(pick[,c(1,4,7,8)])[,c(1,2,3,4)][,1],byrow=T,ncol=5)/test1003[,-6] #Fe0.3Gy
matrix(as.vector(pick[,c(1,4,7,8)])[,c(1,2,3,4)][,2],byrow=T,ncol=5)/test1003[,-6] #Si0.3Gy
matrix(as.vector(pick[,c(2,5,7,8)])[,c(1,2,3,4)][,1],byrow=T,ncol=5)/test1003[,-6] #Fe1Gy
matrix(as.vector(pick[,c(2,5,7,8)])[,c(1,2,3,4)][,2],byrow=T,ncol=5)/test1003[,-6] #Si1Gy
matrix(as.vector(pick[,c(3,6,7,8)])[,c(1,2,3,4)][,1],byrow=T,ncol=5)/test1003[,-6] #FeX1Gy
matrix(as.vector(pick[,c(3,6,7,8)])[,c(1,2,3,4)][,2],byrow=T,ncol=5)/test1003[,-6] #X1Gy

test1005 <- matrix(as.vector(pick[,c(1,4,7,8)])[,c(1,2,3,4)][,2],byrow=T,ncol=5)/test1003[,-6]
test1005$mod <- rownames(test1005)

#Figure_3C (DMP in Reg/DMP)/(Reg in L/all Reg)
ggplot(melt(test1005), aes(x=variable, y=value, fill=as.factor(mod))) +
  geom_bar(stat="identity", position="dodge") +
  scale_color_viridis(discrete=TRUE) +
  theme_classic() +
  labs(title="", x="Nucleus layers", 
       y="Probability of DMP occurrence in histone mark", fill="Histone\nModification") +
  theme(
    axis.title.x = element_text(size = 12, face="bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 12, face="bold"),
    title = element_text(size = 12, face = "bold")) 

apply(test1005[,-6],1,function(x){#as.numeric(as.vector(test1005[1,1:5]))->x
  anova(lm(x~as.vector(1:5)))[[5]][1]
})

#Mean Absolute Methylation change
mean.meth.change.pen.5.Fe1 <- data.frame(data.frame(bind_rows(lapply(pen.5.comp.hist.Fe, 
                                                                      function(x) lapply(x, function(x) lapply(x, function(x) {mean(abs(x[,33]- rowMeans(x[,c(31,43,55)])))}))))))
mean.meth.change.pen.5.Fe1$mod <- c(rep("H3K9me3",5),rep("H3K27me3",5),rep("H3K4me3",5),rep("H3K27ac",5))
mean.meth.change.pen.5.Fe1$comp <- c(rep(c("L1","L2","L3","L4","L5"),4))
mean.meth.change.pen.5.Fe1 <- cbind(apply(mean.meth.change.pen.5.Fe1[,1:13], 2, as.numeric), mean.meth.change.pen.5.Fe1[,14:15])

#Si1Gy
mean.meth.change.pen.5.Si1 <- data.frame(data.frame(bind_rows(lapply(pen.5.comp.hist.Si, 
                                                                      function(x) lapply(x, function(x) lapply(x, function(x) {mean(abs(x[,45]- rowMeans(x[,c(31,43,55)])))}))))))
mean.meth.change.pen.5.Si1$mod <- c(rep("H3K9me3",5),rep("H3K27me3",5),rep("H3K4me3",5),rep("H3K27ac",5))
mean.meth.change.pen.5.Si1$comp <- c(rep(c("L1","L2","L3","L4","L5"),4))
mean.meth.change.pen.5.Si1 <- cbind(apply(mean.meth.change.pen.5.Si1[,1:13], 2, as.numeric), mean.meth.change.pen.5.Si1[,14:15])

fig4 <- cbind(mean.meth.change.pen.5.Fe1[,c(7:9)], mean.meth.change.pen.5.Si1[,c(7:9,14:15)])
colnames(fig4) <- c("Fe 0.3Gy", "Fe 1Gy", "FeX 1Gy", "Si 0.3Gy", "Si 1Gy", "X 1Gy", "mod", "comp")
fig4$mod <- factor(fig4$mod, levels=c("H3K27me3", "H3K9me3", "H3K4me3","H3K27ac"), labels=c("H3K27me3", "H3K9me3", "H3K4me3", "H3K27ac"))

#Methylation change difference between histones and all DMPs
fig11 <- cbind(fig1, "mod"=rep("All_DMPs",5))
fig5 <- rbind(melt(fig11)[,c(2,1,3,4)], melt(fig4))

#Figure_3D
ggplot(filter(fig5, variable == "Fe 1Gy"), aes(x=comp, y=value, col=as.factor(mod))) +
  geom_point(size=4) +
  scale_color_viridis(discrete=TRUE) +
  geom_line(aes(group=interaction(mod)), size=2, alpha=0.5) +
  theme_classic() +
  labs(title="Fe 1Gy", x="Nucleus layers", y="Mean Absolute Methylation change versus baseline (M Value)", color="") +
  theme(
    axis.title.x = element_text(size = 10, face="bold"),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 10, face="bold"),
    title = element_text(size = 10, face = "bold")) +
  ylim(0,1)

#Figure_3E
ggplot(filter(fig5, variable == "Si 1Gy"), aes(x=comp, y=value, col=as.factor(mod))) +
  geom_point(size=4) +
  scale_color_viridis(discrete=TRUE) +
  geom_line(aes(group=interaction(mod)), size=2, alpha=0.5) +
  theme_classic() +
  labs(title="Si 1Gy", x="Nucleus layers", y="Mean Absolute Methylation change versus baseline (M Value)", color="") +
  theme(
    axis.title.x = element_text(size = 10, face="bold"),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 10, face="bold"),
    title = element_text(size = 10, face = "bold")) +
  ylim(0,1)

#Figure_3F
ggplot(filter(fig5, variable == "X 1Gy"), aes(x=comp, y=value, col=as.factor(mod))) +
  geom_point(size=4) +
  scale_color_viridis(discrete=TRUE) +
  geom_line(aes(group=interaction(mod)), size=2, alpha=0.5) +
  theme_classic() +
  labs(title="X 1Gy", x="Nucleus layers", y="Mean Absolute Methylation change versus baseline (M Value)", color="") +
  theme(
    axis.title.x = element_text(size = 10, face="bold"),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 10, face="bold"),
    title = element_text(size = 10, face = "bold")) +
  ylim(0,1)

#Figure_S3B
ggplot(filter(fig5, variable == "Fe 0.3Gy"), aes(x=comp, y=value, col=as.factor(mod))) +
  geom_point(size=4) +
  scale_color_viridis(discrete=TRUE) +
  geom_line(aes(group=interaction(mod)), size=2, alpha=0.5) +
  theme_classic() +
  labs(title="Fe 0.3Gy", x="Nucleus layers", y="Mean Absolute Methylation change versus baseline (M Value)", color="") +
  theme(
    axis.title.x = element_text(size = 10, face="bold"),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 10, face="bold"),
    title = element_text(size = 10, face = "bold")) +
  ylim(0,1)

#Figure_S3C
ggplot(filter(fig5, variable == "Si 0.3Gy"), aes(x=comp, y=value, col=as.factor(mod))) +
  geom_point(size=4) +
  scale_color_viridis(discrete=TRUE) +
  geom_line(aes(group=interaction(mod)), size=2, alpha=0.5) +
  theme_classic() +
  labs(title="Si 0.3Gy", x="Nucleus layers", y="Mean Absolute Methylation change versus baseline (M Value)", color="") +
  theme(
    axis.title.x = element_text(size = 10, face="bold"),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 10, face="bold"),
    title = element_text(size = 10, face = "bold")) +
  ylim(0,1)

#Figure_S3D
ggplot(filter(fig5, variable == "FeX 1Gy"), aes(x=comp, y=value, col=as.factor(mod))) +
  geom_point(size=4) +
  scale_color_viridis(discrete=TRUE) +
  geom_line(aes(group=interaction(mod)), size=2, alpha=0.5) +
  theme_classic() +
  labs(title="Fe(X) 1Gy", x="Nucleus layers", y="Mean Absolute Methylation change versus baseline (M Value)", color="") +
  theme(
    axis.title.x = element_text(size = 10, face="bold"),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 10, face="bold"),
    title = element_text(size = 10, face = "bold")) +
  ylim(0,1)

# 25. Plot data.raw probes into histone modifications picks and Hi-C layers ####
dr <- left_join(data.raw.joined, data.frame("Name"=rownames(data.raw.M.means), data.raw.M.means), by="Name")
data.raw.GRanges <- dr[,c(8:12,1:7,13:ncol(dr))]
colnames(data.raw.GRanges)[1:5] <- gsub("\\.y","", colnames(data.raw.GRanges)[1:5])
data.raw.GRanges <- makeGRangesFromDataFrame(data.raw.GRanges, keep.extra.columns = T)

# Data.raw
dataraw.H3K9me3.ov <- findOverlaps(H3K9me3, data.raw.GRanges)
dataraw.H3K9me3.queryHits <- H3K9me3[queryHits(dataraw.H3K9me3.ov)]

dataraw.H3K27me3.ov <- findOverlaps(H3K27me3, data.raw.GRanges)
dataraw.H3K27me3.queryHits <- H3K27me3[queryHits(dataraw.H3K27me3.ov)]

dataraw.H3K4me3.ov <- findOverlaps(H3K4me3, data.raw.GRanges)
dataraw.H3K4me3.queryHits <- H3K4me3[queryHits(dataraw.H3K4me3.ov)]

dataraw.H3K27ac.ov <- findOverlaps(H3K27ac, data.raw.GRanges)
dataraw.H3K27ac.queryHits <- H3K27ac[queryHits(dataraw.H3K27ac.ov)]

overlap.dataraw <- function(x,y) {
  z <- findOverlaps(x, y)
  x.query <- x[queryHits(z)]
  y.subjectHits <- y[subjectHits(z)]
  
  x.query$ID <- paste0('E_', 1:length(x.query))
  y.subjectHits$ID <- paste0('E_', 1:length(y.subjectHits)) 
  
  x.final <- left_join(as.data.frame(x.query), as.data.frame(y.subjectHits), by="ID")
  
  x.final$cord <- unite(x.final[1:3], cord, sep="_", remove=T)
  
  x.final.grouped <- x.final %>%
    group_by(cord) %>%
    mutate(f0=mean(f0)) %>%
    mutate(f0.3=mean(f0.3)) %>%
    mutate(f1=mean(f1)) %>%
    mutate(f0.14=mean(f0.14)) %>%
    mutate(f0.3.14=mean(f0.3.14)) %>%
    mutate(f1.14=mean(f1.14)) %>%
    mutate(f0.22=mean(f0.22)) %>%
    mutate(f0.3.22=mean(f0.3.22)) %>%    
    mutate(f1.22=mean(f1.22)) %>%    
    mutate(f0.53=mean(f0.53)) %>%    
    mutate(f0.3.53=mean(f0.3.53)) %>%    
    mutate(f1.53=mean(f1.53)) %>%    
    mutate(s0=mean(s0)) %>%   
    mutate(s0.3=mean(s0.3)) %>%    
    mutate(s1=mean(s1)) %>%    
    mutate(s0.13=mean(s0.13)) %>% 
    mutate(s0.3.13=mean(s0.3.13)) %>%
    mutate(s1.13=mean(s1.13)) %>%
    mutate(s0.21=mean(s0.21)) %>%
    mutate(s0.3.21=mean(s0.3.21)) %>%
    mutate(s1.21=mean(s1.21)) %>%
    mutate(s0.62=mean(s0.62)) %>%
    mutate(s0.3.62=mean(s0.3.62)) %>%
    mutate(s1.62=mean(s1.62)) %>%    
    mutate(x0=mean(x0)) %>%    
    mutate(x1=mean(x1)) %>%    
    mutate(x0.13=mean(x0.13)) %>%    
    mutate(x1.13=mean(x1.13)) %>%    
    mutate(x0.22=mean(x0.22)) %>%    
    mutate(x1.22=mean(x1.22)) %>%    
    mutate(x0.62=mean(x0.62)) %>%
    mutate(x1.62=mean(x1.62)) %>%
    distinct(cord, .keep_all=T) %>%
    data.frame()
  x.final.grouped
}

data.raw.Hist <- lapply(list(H3K9me3 = overlap.dataraw(dataraw.H3K9me3.queryHits, data.raw.GRanges),
                             H3K27me3 = overlap.dataraw(dataraw.H3K27me3.queryHits, data.raw.GRanges),
                             H3K4me3 = overlap.dataraw(dataraw.H3K4me3.queryHits, data.raw.GRanges),
                             H3K27ac = overlap.dataraw(dataraw.H3K27ac.queryHits, data.raw.GRanges)), function(x) filter(x, !seqnames.x == "chrY"))

plot1 <- data.frame("mod"=names(data.raw.Hist), "n"=t(data.frame(bind_rows(lapply(data.raw.Hist, function(x) nrow(x))))))
plot1$mod <- factor(plot1$mod, levels=c("H3K27ac", "H3K4me3", "H3K27me3", "H3K9me3"), labels=c("H3K27ac", "H3K4me3", "H3K27me3", "H3K9me3"))

#Figure_3A
ggplot(plot1, aes(x=mod, y=n, fill=as.factor(mod))) +
  geom_bar(position="dodge", stat="identity", width=0.8) +
  scale_fill_viridis(discrete=T, name="") + 
  theme_classic() +
  labs(title="", x="Epigenetic Mark", y="Frequency in primary datset", fill="") +
  theme(
    legend.position="none",
    axis.title.x = element_text(size = 10, face="bold"),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 10, face="bold"),
    title = element_text(size = 10, face = "bold"))

# Layers
data.raw.Hist.comp <- lapply(data.raw.Hist, function(x) list("L1"=filter(x, score>=-2.33444190 & score<=-0.45549509),
                                                             "L2"=filter(x, score>-0.45549509 & score<=0.03823841),
                                                             "L3"=filter(x, score>0.03823841 & score<=0.45039466),
                                                             "L4"=filter(x, score>0.45039466 & score<=0.77691591),
                                                             "L5"=filter(x, score>0.77691591 & score<=1.99093556)))

# Mean methylation LEVEL among compartments
data.raw.meanM <- data.frame(rbindlist(lapply(data.raw.Hist.comp, function(x) lapply(x, function(x) colMeans(na.omit(x[,c(78,90,102)]))))), "mod"=c(rep(c("H3K9me3"),3),rep(c("H3K27me3"),3),rep(c("H3K4me3"),3),rep(c("H3K27ac"),3)),
           "rad"=rep(c("Fe_1Gy","Si_1Gy","X_1Gy"),4))
data.raw.meanM$mod <- factor(data.raw.meanM$mod, levels=c("H3K27me3", "H3K9me3", "H3K4me3", "H3K27ac"), labels=c("H3K27me3", "H3K9me3", "H3K4me3", "H3K27ac"))

data.raw.meanM <- data.frame(rbindlist(lapply(data.raw.Hist.comp, function(x) lapply(x, function(x) mean(rowMeans(na.omit(x[,c(78,90,102)])))))), "mod"=c("H3K9me3", "H3K27me3", "H3K4me3", "H3K27ac"),
                             "rad"=rep("Baseline",4))
data.raw.meanM$mod <- factor(data.raw.meanM$mod, levels=c("H3K27me3", "H3K9me3", "H3K4me3", "H3K27ac"), labels=c("H3K27me3", "H3K9me3", "H3K4me3", "H3K27ac"))

#Figure_S4A
ggplot(melt(data.raw.meanM), aes(x=variable, y=value, col=as.factor(mod))) +
  geom_point(size=4) +
  scale_color_viridis(discrete=TRUE) +
  geom_line(aes(group=interaction(mod, rad)), size=2, alpha=0.5) +
  theme_classic() +
  labs(title="", x="Nucleus layers", y="Mean Baseline Methylation Level (M Value)", color="Histone\nModification") +
  theme(
    axis.title.x = element_text(size = 10, face="bold"),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 10, face="bold"),
    title = element_text(size = 10, face = "bold"))

# 26. Data preparation for regulatory regions (promoter and enhancers) analysis ####
genecode <- rtracklayer::import("/mnt/data1/marcin/references/GENCODE/gencode.v39/gencode.v39.annotation.gtf")

tss <- filter(data.frame(genecode), type == "gene")

tss.minus <- filter(tss, strand == "-")[,c(1,3,3)]
tss.plus <- filter(tss, strand == "+")[,c(1,2,2)]
colnames(tss.minus) <- c("seqnames", "start", "end")
colnames(tss.plus) <- c("seqnames", "start", "end")

tss.joined <- rbind(tss.minus, tss.plus)

tss.joined <- (filter(tss.joined, duplicated(tss.joined) == FALSE))

tss.joined <- extendRegions(makeGRangesFromDataFrame(tss.joined), extend.start = 2000, extend.end = 2000)

overlapDMPs.prom <- function(x,y) { #x=GRanges cords; y=matched.dmps$x 
  z <- findOverlaps(x, y)
  
  qHits <- x[queryHits(z)]
  sHits <- y[subjectHits(z)]
  
  qHits$ID <- paste0('E_', 1:length(qHits))
  sHits$ID <- paste0('E_', 1:length(sHits))
  
  qHits <- as.data.frame(qHits)
  sHits <- as.data.frame(sHits)
  
  a <- left_join(qHits, sHits, by="ID")
  a
}

#DMPs
promoters <- pen.allDMPs.notdetailed[c(13:18)]
promoters <- lapply(promoters, function(x) x[,c(1:24,30:64)])

#Data.raw
promoters.raw <- data.frame(data.raw.GRanges)[,c(6:11,1:5,28,49:50,56:65,73,73,73:104,13)]

promoters <- append(promoters, list("data.raw"=promoters.raw))
colnames(promoters$data.raw) <- colnames(promoters$ALL_DMPs_Fe_0.3Gy)

promoters <- lapply(lapply(promoters, function(x) x[,c(7:11,1:6,12:ncol(x))]), function(x) {colnames(x)[1:5] <- c("seqnames", "start", "end", "width", "strand");
  makeGRangesFromDataFrame(x, keep.extra.columns = T)})

# Overlap TSS with DMPs
promoters <- lapply(lapply(promoters, function(x) overlapDMPs.prom(x, tss.joined)), function(x) {colnames(x)[1:5] <- c("seqnames", "start", "end", "width", "strand");
  makeGRangesFromDataFrame(x, keep.extra.columns = T)})

# Overlap with Histone Modifications

#H3K4me3
promoters.H3K4me3 <- lapply(promoters, function(x) overlapDMPs.prom(x, H3K4me3))
promoters.H3K4me3 <- lapply(promoters.H3K4me3, function(x) {colnames(x)[1:5] <- c("seqnames", "start", "end", "width", "strand"); makeGRangesFromDataFrame(x, keep.extra.columns = T)})

#H3K27ac
promoters.H3K27ac <- lapply(promoters, function(x) overlapDMPs.prom(x, H3K27ac))
promoters.H3K27ac <- lapply(promoters.H3K27ac, function(x) {colnames(x)[1:5] <- c("seqnames", "start", "end", "width", "strand"); makeGRangesFromDataFrame(x, keep.extra.columns = T)})

#H3K4me3, H3K27ac
promoters.H3K4me3.H3K27ac <- lapply(promoters.H3K4me3, function(x) overlapDMPs.prom(x, H3K27ac))
promoters.H3K4me3.H3K27ac <- lapply(promoters.H3K4me3.H3K27ac, function(x) {colnames(x)[1:5] <- c("seqnames", "start", "end", "width", "strand"); makeGRangesFromDataFrame(x, keep.extra.columns = T)})

#H3K4me3-, H3K27ac-
promoters.df <- lapply(promoters, function(x) data.frame(x))
promoters.H3K4me3.df <- lapply(promoters.H3K4me3, function(x) data.frame(x[,1:60]))
promoters.H3K27ac.df <- lapply(promoters.H3K27ac, function(x) data.frame(x[,1:60]))
promoters.H3K4me3.H3K27ac.df <- lapply(promoters.H3K4me3.H3K27ac, function(x) data.frame(x[,1:60]))

colnames(promoters.df$ALL_DMPs_Fe_0.3Gy) == colnames(promoters.H3K4me3.df$ALL_DMPs_Fe_0.3Gy)

promoters.nohist <- mapply(function(x,y) list(anti_join(x,y,by="Name.1")), x=promoters.df, y=promoters.H3K4me3.df)
promoters.nohist <- mapply(function(x,y) list(anti_join(x,y,by="Name.1")), x=promoters.nohist, y=promoters.H3K27ac.df)
promoters.nohist <- mapply(function(x,y) list(anti_join(x,y,by="Name.1")), x=promoters.nohist, y=promoters.H3K4me3.H3K27ac.df)
promoters.nohist <- lapply(promoters.nohist, function(x) {colnames(x)[1:5] <- c("seqnames", "start", "end", "width", "strand"); makeGRangesFromDataFrame(x, keep.extra.columns = T)})

# Take H3K27ac picks and drop from them TSS windows to get enhancers
enhancers <- pen.allDMPs.notdetailed[c(13:18)]
enhancers <- lapply(enhancers, function(x) x[,c(1:24,30:64)])

#Data.raw
enhancers.raw <- data.frame(data.raw.GRanges)[,c(6:11,1:5,28,49:50,56:65,73,73,73:104,13)]

enhancers <- append(enhancers, list("data.raw"=enhancers.raw))
colnames(enhancers$data.raw) <- colnames(enhancers$ALL_DMPs_Fe_0.3Gy)

enhancers <- lapply(lapply(enhancers, function(x) x[,c(7:11,1:6,12:ncol(x))]), function(x) {colnames(x)[1:5] <- c("seqnames", "start", "end", "width", "strand");
makeGRangesFromDataFrame(x, keep.extra.columns = T)})

enhancers <- lapply(enhancers, function(x) overlapDMPs.prom(x, H3K27ac))

enhancers.df <- lapply(enhancers, function(x) data.frame(x))
enhancers.droped <- mapply(function(x,y) list(anti_join(x,y,by="Name.1")), x=enhancers.df, y=promoters.H3K27ac.df)
enhancers.droped <- lapply(enhancers.droped, function(x) {colnames(x)[1:5] <- c("seqnames", "start", "end", "width", "strand"); makeGRangesFromDataFrame(x, keep.extra.columns = T)})

# 27. Add Hi-C to regulatory regions ####
lapply(enhancers.droped, function(x) table(duplicated(x, by="Name.1")))
colnames(as.data.frame(enhancers.droped$ALL_DMPs_Fe_0.3Gy)) == colnames(as.data.frame(promoters.nohist$ALL_DMPs_Fe_0.3Gy))

pen.combined <- append(append(append(append(promoters.H3K4me3,
                   promoters.H3K27ac),
                   promoters.H3K4me3.H3K27ac),
                   promoters.nohist),
                   enhancers.droped)

names(pen.combined) <- c("p_H3K4me3_Fe_0.3Gy", "p_H3K4me3_Fe_1Gy", "p_H3K4me3_FeX_1Gy", "p_H3K4me3_Si_0.3Gy", "p_H3K4me3_Si_1Gy", "p_H3K4me3_X_1Gy", "p_H3K4me3_raw",
                       "p_H3K27ac_Fe_0.3Gy", "p_H3K27ac_Fe_1Gy", "p_H3K27ac_FeX_1Gy", "p_H3K27ac_Si_0.3Gy", "p_H3K27ac_Si_1Gy", "p_H3K27ac_X_1Gy", "p_H3K27ac_raw",
                       "p_H3K4me3.H3K27ac_Fe_0.3Gy", "p_H3K4me3.H3K27ac_Fe_1Gy", "p_H3K4me3.H3K27ac_FeX_1Gy", "p_H3K4me3.H3K27ac_Si_0.3Gy", "p_H3K4me3.H3K27ac_Si_1Gy", "p_H3K4me3.H3K27ac_X_1Gy", "p_H3K4me3.H3K27ac_raw",
                       "p_noHist_Fe_0.3Gy", "p_noHist_Fe_1Gy", "p_noHist_FeX_1Gy", "p_noHist_Si_0.3Gy", "p_noHist_Si_1Gy", "p_noHist_X_1Gy", "p_noHist_raw",
                       "e_H3K27ac_Fe_0.3Gy", "e_H3K27ac_Fe_1Gy", "e_H3K27ac_FeX_1Gy", "e_H3K27ac_Si_0.3Gy", "e_H3K27ac_Si_1Gy", "e_H3K27ac_X_1Gy", "e_H3K27ac_raw")

#Add_Layers
pen.layers <- lapply(pen.combined, function(x) list("L1"=filter(x, score>=-2.33444190 & score<=-0.45549509),
                                                 "L2"=filter(x, score>-0.45549509 & score<=0.03823841),
                                                 "L3"=filter(x, score>0.03823841 & score<=0.45039466),
                                                 "L4"=filter(x, score>0.45039466 & score<=0.77691591),
                                                 "L5"=filter(x, score>0.77691591 & score<=1.99093556)))

normal.pen.layers <- data.frame(t(data.frame(bind_rows(lapply(pen.layers, function(x) lapply(x, function(x) length(x)))))))

colnames(normal.pen.layers) <- c(names(pen.layers))
normal.pen.layers.pct <- apply(normal.pen.layers,2,function(x){x/sum(x)})

#Figure_S3F
ggplot(melt(normal.pen.layers.pct[,c(7,14,21,28,35)]), aes(x=Var1, y=value, col=as.factor(Var2))) +
  geom_point(size=4) +
  scale_color_viridis(discrete=TRUE) +
  geom_line(aes(group=interaction(Var2)), size=2, alpha=0.5) +
  theme_classic() +
  labs(title="Frequency of primary probes in regulatory\nregions across nucleus layers", x="Nucleus layers", y="Frequency of primary probes (%)", color="Sample") +
  theme(
    axis.title.x = element_text(size = 10, face="bold"),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 10, face="bold"),
    title = element_text(size = 10, face = "bold"))

# 28. DNA Methylation change in regulatory regions ####
pen.layers.Fe.meth.change <- data.frame(t(data.frame(bind_rows(lapply(pen.layers[c(1:3,8:10,15:17,22:24,29:31)],
                                                                      function(x) lapply(x, function(x) {abs(mean(as.data.frame(x)[,25] - rowMeans(as.data.frame(x)[,c(27,39,51)])))}))))), "comp"=rownames(normal.pen.pct))

colnames(pen.layers.Fe.meth.change) <- c(names(pen.layers[c(1:3,8:10,15:17,22:24,29:31)]), "comp")

#Si
pen.layers.SiX.meth.change <- data.frame(t(data.frame(bind_rows(lapply(pen.layers[c(4:6,11:13,18:20,25:27,32:34)], 
                                                                       function(x) lapply(x, function(x) {abs(mean(as.data.frame(x)[,26] - rowMeans(as.data.frame(x)[,c(27,39,51)])))}))))), "comp"=rownames(normal.pen.pct))

colnames(pen.layers.SiX.meth.change) <- c(names(pen.layers[c(4:6,11:13,18:20,25:27,32:34)]), "comp")

#Figure_3G
#Fe 1Gy
ggplot(melt(pen.layers.Fe.meth.change[,c(2,5,8,11,14,16)]), aes(x=comp, y=value, col=as.factor(variable))) +
  geom_point(size=4) +
  scale_color_viridis(discrete=TRUE) +
  geom_line(aes(group=interaction(variable)), size=2, alpha=0.5) +
  theme_classic() +
  labs(title="", x="Nucleus layer", y="Absolute mean methylation change versus control (M Value)", color="Regulatory Regions") +
  theme(
    axis.title.x = element_text(size = 10, face="bold"),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 10, face="bold"),
    title = element_text(size = 10, face = "bold"))

#Figure_S3H
#Fe 0.3Gy
ggplot(melt(pen.layers.Fe.meth.change[,c(1,4,7,10,13,16)]), aes(x=comp, y=value, col=as.factor(variable))) +
  geom_point(size=4) +
  scale_color_viridis(discrete=TRUE) +
  geom_line(aes(group=interaction(variable)), size=2, alpha=0.5) +
  theme_classic() +
  labs(title="", x="Nucleus layer", y="Absolute mean methylation change versus control (M Value)", color="Regulatory Regions") +
  theme(
    axis.title.x = element_text(size = 10, face="bold"),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 10, face="bold"),
    title = element_text(size = 10, face = "bold"))

#Figure_S3I
#FeX 1Gy
ggplot(melt(pen.layers.Fe.meth.change[,c(3,6,9,12,15,16)]), aes(x=comp, y=value, col=as.factor(variable))) +
  geom_point(size=4) +
  scale_color_viridis(discrete=TRUE) +
  geom_line(aes(group=interaction(variable)), size=2, alpha=0.5) +
  theme_classic() +
  labs(title="", x="Nucleus layer", y="Absolute mean methylation change versus control (M Value)", color="Regulatory Regions") +
  theme(
    axis.title.x = element_text(size = 10, face="bold"),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 10, face="bold"),
    title = element_text(size = 10, face = "bold"))

#Figure_S3J
#Si 0.3Gy
ggplot(melt(pen.layers.SiX.meth.change[,c(1,4,7,10,13,16)]), aes(x=comp, y=value, col=as.factor(variable))) +
  geom_point(size=4) +
  scale_color_viridis(discrete=TRUE) +
  geom_line(aes(group=interaction(variable)), size=2, alpha=0.5) +
  theme_classic() +
  labs(title="", x="Nucleus layer", y="Absolute mean methylation change versus control (M Value)", color="Regulatory Regions") +
  theme(
    axis.title.x = element_text(size = 10, face="bold"),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 10, face="bold"),
    title = element_text(size = 10, face = "bold"))

#Figure_S3K
#Si 1Gy
ggplot(melt(pen.layers.SiX.meth.change[,c(2,5,8,11,14,16)]), aes(x=comp, y=value, col=as.factor(variable))) +
  geom_point(size=4) +
  scale_color_viridis(discrete=TRUE) +
  geom_line(aes(group=interaction(variable)), size=2, alpha=0.5) +
  theme_classic() +
  labs(title="", x="Nucleus layer", y="Absolute mean methylation change versus control (M Value)", color="Regulatory Regions") +
  theme(
    axis.title.x = element_text(size = 10, face="bold"),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 10, face="bold"),
    title = element_text(size = 10, face = "bold"))

#Figure_S3L
#X 1Gy
ggplot(melt(pen.layers.SiX.meth.change[,c(3,6,9,12,15,16)]), aes(x=comp, y=value, col=as.factor(variable))) +
  geom_point(size=4) +
  scale_color_viridis(discrete=TRUE) +
  geom_line(aes(group=interaction(variable)), size=2, alpha=0.5) +
  theme_classic() +
  labs(title="", x="Nucleus layer", y="Absolute mean methylation change versus control (M Value)", color="Regulatory Regions") +
  theme(
    axis.title.x = element_text(size = 10, face="bold"),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 10, face="bold"),
    title = element_text(size = 10, face = "bold"))

# 29. Baseline methylation level in regulatory regions ####

baseline.level.mean <- data.frame(t(data.frame(bind_rows(lapply(pen.layers[c(7,14,21,28,35)], function(x) lapply(x, function(x) {mean(colMeans(na.omit(as.data.frame(x)[,c(27,39,51)])))}))))), "comp"=rownames(normal.pen.pct))
colnames(baseline.level.mean) <- c(names(pen.layers[c(7,14,21,28,35)]), "comp")

#Figure_S3G
ggplot(melt(baseline.level.mean), aes(x=comp, y=value, col=as.factor(variable))) +
  geom_point(size=4) +
  scale_color_viridis(discrete=TRUE) +
  geom_line(aes(group=interaction(variable)), size=2, alpha=0.5) +
  theme_classic() +
  labs(title="Baseline mean methylation level in regulatory\nregions overlapped with primary probes", x="Nucleus layer", y="Baseline Mean Methylation Level (M Value)", color="Regulatory Region") +
  theme(
    axis.title.x = element_text(size = 12, face="bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 12, face="bold"),
    title = element_text(size = 12, face = "bold"))

# 30. Add Hi-C to regulatory regions (without DMPs overlapping) ####

# Overlap TSS with Histone
tss.H3K4me3 <- overlapDMPs.prom(tss.joined, H3K4me3)
colnames(tss.H3K4me3)[1:5] <- c("seqnames", "start", "end", "width", "strand")
tss.H3K4me3 <- makeGRangesFromDataFrame(tss.H3K4me3, keep.extra.columns = F)

#H3K27ac
tss.H3K27ac <- overlapDMPs.prom(tss.joined, H3K27ac)
colnames(tss.H3K27ac)[1:5] <- c("seqnames", "start", "end", "width", "strand")
tss.H3K27ac <- makeGRangesFromDataFrame(tss.H3K27ac, keep.extra.columns = F)

#H3K4me3, H3K27ac
tss.H3K4me3.H3K27ac <- overlapDMPs.prom(tss.H3K4me3, H3K27ac)
colnames(tss.H3K4me3.H3K27ac)[1:5] <- c("seqnames", "start", "end", "width", "strand")
tss.H3K4me3.H3K27ac <- makeGRangesFromDataFrame(tss.H3K4me3.H3K27ac, keep.extra.columns = F)

#H3K4me3-, H3K27ac-
tss.nohist <- anti_join(as.data.frame(tss.joined),as.data.frame(tss.H3K4me3))
tss.nohist <- anti_join(as.data.frame(tss.nohist),as.data.frame(tss.H3K27ac))
tss.nohist <- anti_join(as.data.frame(tss.nohist),as.data.frame(tss.H3K4me3.H3K27ac))
tss.nohist <- makeGRangesFromDataFrame(tss.nohist, keep.extra.columns = F)  

#enhancers
enhancers <- overlapDMPs.prom(H3K27ac, tss.H3K27ac)
colnames(enhancers)[1:5] <- c("seqnames", "start", "end", "width", "strand")
enhancers <- makeGRangesFromDataFrame(enhancers, keep.extra.columns = F)
enhancers <- anti_join(as.data.frame(H3K27ac), as.data.frame(enhancers))
enhancers <- makeGRangesFromDataFrame(enhancers, keep.extra.columns = F)

#add layers
layers.pro.enh <- list("p.H3K4me3"=tss.H3K4me3,
                       "p.H3K27ac"=tss.H3K27ac,
                       "p.H3K4me3.H3K27ac"=tss.H3K4me3.H3K27ac,
                       "p.nohist"=tss.nohist,
                       "e.H3K27ac"=enhancers)

layers.pro.enh <- lapply(layers.pro.enh, function(x) overlapDMPs.prom(x, compartments))

layers.pro.enh <- lapply(layers.pro.enh, function(x) list("L1"=filter(x, score>=-2.33444190 & score<=-0.45549509),
                                                          "L2"=filter(x, score>-0.45549509 & score<=0.03823841),
                                                          "L3"=filter(x, score>0.03823841 & score<=0.45039466),
                                                          "L4"=filter(x, score>0.45039466 & score<=0.77691591),
                                                          "L5"=filter(x, score>0.77691591 & score<=1.99093556)))

normal.pen.layers <- data.frame(t(data.frame(bind_rows(lapply(layers.pro.enh, function(x) lapply(x, function(x) nrow(x)))))))

colnames(normal.pen.layers) <- c(names(layers.pro.enh))
normal.pen.layers.pct <- apply(normal.pen.layers,2,function(x){x/sum(x)})

#Figure_S3E
ggplot(melt(normal.pen.layers.pct), aes(x=Var1, y=value, col=as.factor(Var2))) +
  geom_point(size=4) +
  scale_color_viridis(discrete=TRUE) +
  geom_line(aes(group=interaction(Var2)), size=2, alpha=0.5) +
  theme_classic() +
  labs(title="Frequency of regulatory regions across nucleus layers", x="Nucleus layers", y="nRegulatory Regions (%)", color="Location") +
  theme(
    axis.title.x = element_text(size = 10, face="bold"),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 10, face="bold"),
    title = element_text(size = 10, face = "bold"))


