---
title: "validation"
author: "Yi"
date: "4/20/2021"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

In this file we will validate the optimized HMM model. The Alexandrov data set will be used in this process and the comparison is between optimized HMM, PCF from Coen, and Alexandrov results by Venn diagram.

```{r}
### Library Loading ---------------------------------------------------

library(BSgenome.Hsapiens.UCSC.hg19)
library(VariantAnnotation)
library(GenomicRanges)
library(plyr)
library(dplyr)
library(stepR)
library(pbapply)
library(rlist)
library(tidyverse)

### Downloading and readying Alexandrov Samples -----------------------------
# use function 'readAli()' to import data
readAli <- function(y){
  data <- readr::read_delim(file=as.character(y), col_names = c('sampleId', 'type', 'chr', 'start', 'end', 'REF', 'ALT', 'origin'), trim_ws = T, col_types = 'cccnnccc', delim = '\t')
  data$tissue <- gsub('_.*', '', gsub('.*/', '', y))
  return(data)
}

# Links to mutations.
alexandrov.links <- readr::read_delim('../Alexandrov/Alexandrov_Datalinks.txt', delim = '\t', col_names = F)

# only use the breast cancer samples
samples.Ali <- lapply(alexandrov.links$X1[4], readAli)
samples.Ali <- do.call(rbind, samples.Ali)

# Only select unique SNVs per sample.
samples.Ali <- samples.Ali %>%
  # Select SNVs.
  dplyr::filter(type == 'subs') %>%
  dplyr::distinct() %>%
 
  # Only take samples with >15 SNVs
  dplyr::group_by(sampleId) %>%
  dplyr::mutate(totalMuts = dplyr::n()) %>%
  dplyr::filter(totalMuts >= 15) %>%
 
  # Only take chromosome with >15 SNVs
  dplyr::group_by(sampleId,chr) %>%
  dplyr::mutate(chrMuts = dplyr::n()) %>%
  dplyr::filter(chrMuts >= 15) %>%

  dplyr::ungroup()

# breast cancer samples with KatRanges details stored in samples.A2.
samples.A2 <- base::lapply(base::split(samples.Ali, samples.Ali$sampleId), function(y){
  m <- VariantAnnotation::VRanges(seqnames = y$chr, ranges = IRanges::IRanges(start = y$start, end = y$end), ref = y$REF, alt = y$ALT, sampleNames = unique(y$sampleId))
  mm <- KatDetect::KatRanges(vr = m, genome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19)
  return(mm)
})


### Perform optimized HMM per chromosome ------------------------------
set.seed(2000)
HMMresults_A <- BiocParallel::bplapply(samples.A2, function(x){
    # Perform HMM per chromosome.
    hmm.PerChr <- lapply(GenomeInfoDb::seqlevelsInUse(x), function(y){
        ParallelLogger::logInfo(sprintf('Working on: %s - %s', levels(sampleNames(x)), y))
        return(performHMM_op(x, chr_num = y))
    })
}, BPPARAM = BiocParallel::MulticoreParam(workers = 1, progressbar = T))

### Conclude the optimized HMM results into two tables, one is on sample level and one is on kataegis level.-------------------------------------

# remove empty rows and samples and get table on kataegis level
HMMresults_A <-BiocParallel::bplapply(HMMresults_A, function(x){
  lapply(x,function(y){
    ifelse(nrow(y)!=0, return(y),return(NULL))})
}, BPPARAM = BiocParallel::MulticoreParam(workers = 1, progressbar = T))

HMMresults_A <- lapply(HMMresults_A, function(x){
  return(Filter(length,x))
})

HMMresults_A <- Filter(length,HMMresults_A)%>%
  lapply(function(x){return(do.call(rbind,x))})

HMMresults_A <- do.call(rbind,HMMresults_A)

# get table on sample level
HMMresults_A2 <- HMMresults_A %>% group_by(sampleNames)%>%
  summarise(samples = n())

##################################################################
### perform PCF() ------------------------------------------
# use the breast cancer samples
samples.Ali <- lapply(alexandrov.links$X1[4], readAli)
samples.Ali <- do.call(rbind, samples.Ali)
#samples.Ali <- as_tibble(samples.Ali)

# Only select unique SNVs per sample.
samples.Ali <- samples.Ali %>%
  # Select SNVs.
  dplyr::filter(type == 'subs') %>%
  dplyr::distinct()

# change the order of chromosomes
samples.Ali$chr = factor(samples.Ali$chr,c(1:22,'X','Y'),ordered = T)
samples.Ali = samples.Ali[with(samples.Ali,order(chr)),]
samples.Ali <-
  na.omit(samples.Ali)
# samples.
samples.A3 <- base::lapply(base::split(samples.Ali, samples.Ali$sampleId), function(y){
  m <- VariantAnnotation::VRanges(seqnames = y$chr, ranges = IRanges::IRanges(start = y$start, end = y$end), ref = y$REF, alt = y$ALT, sampleNames = unique(y$sampleId))
  mm <- KatDetect::KatRanges(vr = m, genome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19)
  return(mm)
})

results_PCF <- lapply(samples.A3,performPCF2,minSegmentSize.Kataegis = 5, maxMeanIMD.Kataegis = 2500, minSegmentSize = 5, usePreComputedMonteCarlo = F, nThreads = 10)

### Conclude the optimized HMM results into two tables, one is on sample level and one is on kataegis level.-------------------------------------
results_PCF <- lapply(results_PCF, function(x){
  n <- data.frame(x[x$couldBeKataegis == TRUE])
  return(n)
})
results_PCF <- do.call(rbind,results_PCF)
# on kataegis level
results_PCF <- dplyr::select(results_PCF,sample,seqnames,start,end,totalVariantsInSegment)
results_PCF <- results_PCF %>% mutate(nMuts= results_PCF$totalVariantsInSegment +1)
# on sample level
results_PCF2 <- results_PCF %>% group_by(sample)%>%
  summarise(samples_num = n())
```

```{r}
### import results of Alexandrov
library(readxl)
# on kataegis level
Alexandrov_KatSamplesDetail <- read_excel("../Alexandrov/Alexandrov_KatSamplesDetail.xls", range = "B4:G877")
Alexanderov_Breast <- Alexandrov_KatSamplesDetail %>%
  filter(Tissue == 'Breast Cancer')
Alexanderov_Breast <- Alexanderov_Breast[order(Alexanderov_Breast$`Sample Name`),]

# on sample level
Alexandrov_SamplesWithKat <- read.delim("~/katdetect/katdetect/Misc/Alexandrov/Alexandrov_SamplesWithKat.txt")
Alexandrov_SamplesWithKat <- Alexandrov_SamplesWithKat %>%
  filter(Tissue == 'Breast Cancer')
```

```{r}
### compare results on sample level----------------------------------------------
# compare optimized_HMM, PCF
setdiff(HMMresults_A$sampleNames,results_PCF$sample)
setdiff(results_PCF$sample,HMMresults_A$sampleNames)
# compare optimized_HMM, Alex
setdiff(Alexandrov_SamplesWithKat$Sample,HMMresults_A$sampleNames)
setdiff(HMMresults_A$sampleNames,Alexandrov_SamplesWithKat$Sample)
# compare Alex, PCF
setdiff(results_PCF$sample,Alexandrov_SamplesWithKat$Sample)
setdiff(Alexandrov_SamplesWithKat$Sample,results_PCF$sample,)
```

```{r}
#### 'Samples Level':
#### make a venn plot to see the overlap
library(devtools)
library(ggvenn)

venn_list = list(
  results_PCF$sample,
  HMMresults_A$sampleNames,
  Alexandrov_SamplesWithKat$Sample
)

names(venn_list) <- c('PCF(72)','HMM(81)',"Alx(67)")

ggvenn(
    venn_list,
    fill_color = c("#0073C2FF", "#EFC000FF",  "#CD534CFF"),
    stroke_size = 0.5, set_name_size = 3
)

```

```{r}
#### 'kataegis Level':
#### make a venn plot to see the overlap
library(devtools)
library(ggvenn)

Alexanderov_Breast <- read.csv("~/katdetect/katdetect/results/Alexanderov_Breast_details", sep="")
Alexanderov_Breast <- Alexanderov_Breast[,-1]
colnames(Alexanderov_Breast) <- c('SampleName','Chr','Start','End','nMut')
Alexanderov_Breast$Chr <- paste0('chr',Alexanderov_Breast$Chr)
Alexanderov_Breast <- mutate(Alexanderov_Breast,Alx = TRUE)

HMM_op_details <- read.csv("~/katdetect/katdetect/results/HMM_op_details", sep="")
colnames(HMM_op_details) <- c('SampleName','Chr','Start','End','nMut')
HMM_op_details <- mutate(HMM_op_details,HMM = TRUE)

PCF_details <- read.csv("~/katdetect/katdetect/results/PCF_details", sep="")
PCF_details <- PCF_details[-5]
colnames(PCF_details) <- c('SampleName','Chr','Start','End','nMut')
PCF_details <- mutate(PCF_details,PCF = TRUE)

k <- full_join(Alexanderov_Breast,HMM_op_details)
k <- full_join(k,PCF_details)
k$rowid = 1:nrow(k)

AlexID <- k[which(k$Alx == TRUE),]$rowid
PCFID <- k[which(k$PCF == TRUE),]$rowid
HMMID <- k[which(k$HMM == TRUE),]$rowid

venn_list2 = list(
  PCFID,
   HMMID,
  AlexID
)

names(venn_list2) <- c('PCF(360)',"HMM(664)",'Alex(456)')

ggvenn(
    venn_list2,
    fill_color = c("#0073C2FF", "#EFC000FF",  "#CD534CFF"),
    stroke_size = 0.5, set_name_size = 3
)

```
