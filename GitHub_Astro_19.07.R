# 1. Load packages ####
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(minfi)
library(data.table)
library(reshape2)
library(ggplot2)
library(dplyr)
library(ggvenn)
library(readxl)
library(UpSetR)
library(RColorBrewer)
library(rlist)
library(plyr)
library(datasets)
library(ggupset) 
library(VennDiagram)
library(plyranges)

# 2. Load files ####
wpath<- "/mnt/work1/adrian/Nasa/"
setwd(wpath)

#https://data.4dnucleome.org/files-processed/4DNFIOZT9QGF/
compartments <- read_bigwig("4DNFIOZT9QGF.bw") #stem cells 

compartments.k562 <- read_bigwig("4DNFI5WH9HQX.bw") #leukemia lymphoblasts https://data.4dnucleome.org/experiment-set-replicates/4DNESU95RUNO/

compartments.hap1 <- read_bigwig("4DNFI2DJ2HNF.bw") #leukemia lymphoblasts https://data.4dnucleome.org/files-processed/4DNFI2DJ2HNF/

compartments.imr90 <- read_bigwig("4DNFIHM89EGL.bw") #lung fibroblasts


ann450k<-getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# 3. Prepare files ####
colnames(ann450k)
genes <- ann450k[,c(1,2,2,3,22)]
rownames(genes) <- NULL
colnames(genes) <- c("chr", "start", "end", "strand", "gene")
genes$gene <- gsub("\\;.*","", genes$gene)
genes <- filter(as.data.frame(genes), !gene == "")

genes <- makeGRangesFromDataFrame(genes, keep.extra.columns = T)

# Overlap
genes.ov <- findOverlaps(compartments, genes)

genes.query <- compartments[queryHits(genes.ov)]
genes.subject <- genes[subjectHits(genes.ov)]

length(genes.query) == length(genes.subject)

# Give IDs and filter
genes.query$ID <- paste0('E_', 1:length(genes.query))
genes.subject$ID <- paste0('E_', 1:length(genes.subject))

# Annotate DMRs
genes.query <- as.data.frame(genes.query)
genes.subject <- as.data.frame(genes.subject)

genes.joined <- left_join(genes.query, genes.subject, by="ID")

genes.joined <- genes.joined %>% 
  group_by(gene) %>%
  distinct(gene, .keep_all=T)

genes.joined <- filter(genes.joined, !score == "NaN")
genes.joined <- filter(genes.joined, !seqnames.x == "chrY"  & !seqnames.x == "chrX")


# 4. Build layers ####
quantile(genes.joined$score, c(0,0.2,0.4,0.6,0.8,1))

genes.layers <- mutate(genes.joined, Layer =  case_when(
  (score>=-1.6999605 & score<=-0.2423967) ~ "L1",
  (score>-0.2423967 & score<=0.1784833) ~ "L2",
  (score>0.1784833 & score<=0.4581482) ~ "L3",
  (score>0.4581482 & score<=0.6208739) ~ "L4",
  (score>0.6208739 & score<=1.2580514) ~ "L5"))

genes.final <- genes.layers[,c(8,9,13,14,6)]
colnames(genes.final) <- c("chr", "cord", "gene", "layer", "score")

# 5. Whole genome - stem cells ####
wpath <- "/mnt/data1/marcin/references/GENCODE"
setwd(wpath)

addnotations <- rtracklayer::import("gencode.v39/gencode.v39lift37.annotation.gtf")

genes <- as.data.frame(filter(addnotations, type == "gene" & gene_type == "protein_coding"))[,c(1,2,3,5,10)]
genes$gene_id <- gsub("\\..*","", genes$gene_id)
rownames(genes) <- NULL
colnames(genes) <- c("chr", "start", "end", "strand", "gene")
genes <- makeGRangesFromDataFrame(genes, keep.extra.columns = T)

# Overlap
genes.stem.ov <- findOverlaps(compartments, genes)

genes.stem.query <- compartments[queryHits(genes.stem.ov)]
genes.stem.subject <- genes[subjectHits(genes.stem.ov)]

length(genes.stem.query) == length(genes.stem.subject)

# Give IDs and filter
genes.stem.query$ID <- paste0('E_', 1:length(genes.stem.query))
genes.stem.subject$ID <- paste0('E_', 1:length(genes.stem.subject))

# Annotate DMRs
genes.stem.query <- as.data.frame(genes.stem.query)
genes.stem.subject <- as.data.frame(genes.stem.subject)

genes.stem.joined <- left_join(genes.stem.query, genes.stem.subject, by="ID")

genes.stem.joined <- filter(genes.stem.joined, !gene == "")

genes.stem.joined <- genes.stem.joined %>% 
  group_by(gene) %>%
  distinct(gene, .keep_all=T)

genes.stem.joined <- filter(genes.stem.joined, !score == "NaN")
genes.stem.joined <- filter(genes.stem.joined, !seqnames.x == "chrY"  & !seqnames.x == "chrX")

#Build layers for the whole genome
quantile(genes.stem.joined$score, c(0,0.2,0.4,0.6,0.8,1))

genes.stem.layers <- mutate(genes.stem.joined, Layer =  case_when(
  (score>=-1.6999605 & score<=-0.2074517) ~ "L1",
  (score>-0.2074517 & score<=0.2107766) ~ "L2",
  (score>0.2107766 & score<=0.4763579) ~ "L3",
  (score>0.4763579 & score<=0.6317974) ~ "L4",
  (score>0.6317974 & score<=1.2580514) ~ "L5"))

genes.stem.final <- genes.stem.layers[,c(8,9,10,11,12,13,14,2,3,6)]

colnames(genes.stem.final) <- c("chr", "start", "end", "width", "strand", "gene", "layer", "layer_start", "layer_end", "score")

# 6. Whole genome - k562 cells ####

# Overlap
genes.k562.ov <- findOverlaps(compartments.k562, genes)

genes.k562.query <- compartments.k562[queryHits(genes.k562.ov)]
genes.k562.subject <- genes[subjectHits(genes.k562.ov)]

length(genes.k562.query) == length(genes.k562.subject)

# Give IDs and filter
genes.k562.query$ID <- paste0('E_', 1:length(genes.k562.query))
genes.k562.subject$ID <- paste0('E_', 1:length(genes.k562.subject))

# Annotate DMRs
genes.k562.query <- as.data.frame(genes.k562.query)
genes.k562.subject <- as.data.frame(genes.k562.subject)

genes.k562.joined <- left_join(genes.k562.query, genes.k562.subject, by="ID")

genes.k562.joined <- filter(genes.k562.joined, !gene == "")

genes.k562.joined <- genes.k562.joined %>% 
  group_by(gene) %>%
  distinct(gene, .keep_all=T)

genes.k562.joined <- filter(genes.k562.joined, !score == "NaN")
genes.k562.joined <- filter(genes.k562.joined, !seqnames.x == "chrY"  & !seqnames.x == "chrX")

#Build layers for the whole genome
quantile(genes.k562.joined$score, c(0,0.2,0.4,0.6,0.8,1))

genes.k562.layers <- mutate(genes.k562.joined, Layer =  case_when(
  (score>=-2.8002751 & score<=-0.1988029) ~ "L1",
  (score>-0.1988029 & score<=0.3726020) ~ "L2",
  (score>0.3726020 & score<=0.6409140) ~ "L3",
  (score>0.6409140 & score<=0.8865776) ~ "L4",
  (score>0.8865776 & score<=1.9667625) ~ "L5"))

genes.k562.final <- genes.k562.layers[,c(8,9,10,11,12,13,14,2,3,6)]

colnames(genes.k562.final) <- c("chr", "start", "end", "width", "strand", "gene", "layer", "layer_start", "layer_end", "score")

# 7. Whole genome - hap1 cells ####

# Overlap
genes.hap1.ov <- findOverlaps(compartments.hap1, genes)

genes.hap1.query <- compartments.hap1[queryHits(genes.hap1.ov)]
genes.hap1.subject <- genes[subjectHits(genes.hap1.ov)]

length(genes.hap1.query) == length(genes.hap1.subject)

# Give IDs and filter
genes.hap1.query$ID <- paste0('E_', 1:length(genes.hap1.query))
genes.hap1.subject$ID <- paste0('E_', 1:length(genes.hap1.subject))

# Annotate DMRs
genes.hap1.query <- as.data.frame(genes.hap1.query)
genes.hap1.subject <- as.data.frame(genes.hap1.subject)

genes.hap1.joined <- left_join(genes.hap1.query, genes.hap1.subject, by="ID")

genes.hap1.joined <- filter(genes.hap1.joined, !gene == "")

genes.hap1.joined <- genes.hap1.joined %>% 
  group_by(gene) %>%
  distinct(gene, .keep_all=T)

genes.hap1.joined <- filter(genes.hap1.joined, !score == "NaN")
genes.hap1.joined <- filter(genes.hap1.joined, !seqnames.x == "chrY"  & !seqnames.x == "chrX")

#Build layers for the whole genome
quantile(genes.hap1.joined$score, c(0,0.2,0.4,0.6,0.8,1))

genes.hap1.layers <- mutate(genes.hap1.joined, Layer =  case_when(
  (score>=-1.7140336 & score<=-0.1877080) ~ "L1",
  (score>-0.1877080 & score<=0.2553671) ~ "L2",
  (score>0.2553671 & score<=0.5731668) ~ "L3",
  (score>0.5731668 & score<=0.8459858) ~ "L4",
  (score>0.8459858 & score<=2.2107592) ~ "L5"))

genes.hap1.final <- genes.hap1.layers[,c(8,9,10,11,12,13,14,2,3,6)]

colnames(genes.hap1.final) <- c("chr", "start", "end", "width", "strand", "gene", "layer", "layer_start", "layer_end", "score")

# 8. Whole genome - imr90 cells ####

# Overlap
genes.imr90.ov <- findOverlaps(compartments.imr90, genes)

genes.imr90.query <- compartments.imr90[queryHits(genes.imr90.ov)]
genes.imr90.subject <- genes[subjectHits(genes.imr90.ov)]

length(genes.imr90.query) == length(genes.imr90.subject)

# Give IDs and filter
genes.imr90.query$ID <- paste0('E_', 1:length(genes.imr90.query))
genes.imr90.subject$ID <- paste0('E_', 1:length(genes.imr90.subject))

# Annotate DMRs
genes.imr90.query <- as.data.frame(genes.imr90.query)
genes.imr90.subject <- as.data.frame(genes.imr90.subject)

genes.imr90.joined <- left_join(genes.imr90.query, genes.imr90.subject, by="ID")

genes.imr90.joined <- filter(genes.imr90.joined, !gene == "")

genes.imr90.joined <- genes.imr90.joined %>% 
  group_by(gene) %>%
  distinct(gene, .keep_all=T)

genes.imr90.joined <- filter(genes.imr90.joined, !score == "NaN")
genes.imr90.joined <- filter(genes.imr90.joined, !seqnames.x == "chrY"  & !seqnames.x == "chrX")

#Build layers for the whole genome
quantile(genes.imr90.joined$score, c(0,0.2,0.4,0.6,0.8,1))

genes.imr90.layers <- mutate(genes.imr90.joined, Layer =  case_when(
  (score>=-2.02147079 & score<=-0.47002986) ~ "L1",
  (score>-0.47002986 & score<=0.02891807) ~ "L2",
  (score>0.02891807 & score<=0.45515209) ~ "L3",
  (score>0.45515209 & score<=0.76953101) ~ "L4",
  (score>0.76953101 & score<=1.99093556) ~ "L5"))

genes.imr90.final <- genes.imr90.layers[,c(8,9,10,11,12,13,14,2,3,6)]

colnames(genes.imr90.final) <- c("chr", "start", "end", "width", "strand", "gene", "layer", "layer_start", "layer_end", "score")


# 9. Compare overlapping genes for further analysis ####

space.genes <- list("All"=Reduce(intersect, list(genes.stem.final$gene,
                                                 genes.k562.final$gene,
                                                 genes.hap1.final$gene)),
                    "L1"=Reduce(intersect, list(filter(genes.stem.final, layer == "L1")$gene,
                                                filter(genes.k562.final, layer == "L1")$gene,
                                                filter(genes.hap1.final, layer == "L1")$gene)),
                    "L2"=Reduce(intersect, list(filter(genes.stem.final, layer == "L2")$gene,
                                                filter(genes.k562.final, layer == "L2")$gene,
                                                filter(genes.hap1.final, layer == "L2")$gene)),
                    "L3"=Reduce(intersect, list(filter(genes.stem.final, layer == "L3")$gene,
                                                filter(genes.k562.final, layer == "L3")$gene,
                                                filter(genes.hap1.final, layer == "L3")$gene)),
                    "L4"=Reduce(intersect, list(filter(genes.stem.final, layer == "L4")$gene,
                                                filter(genes.k562.final, layer == "L4")$gene,
                                                filter(genes.hap1.final, layer == "L4")$gene)),
                    "L5"=Reduce(intersect, list(filter(genes.stem.final, layer == "L5")$gene,
                                                filter(genes.k562.final, layer == "L5")$gene,
                                                filter(genes.hap1.final, layer == "L5")$gene)))
#list.save(space.genes, 'space.genes.rdata')


