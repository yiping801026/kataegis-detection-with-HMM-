---
title: "Statistics in CPCT-02-DRUP metastatic datasets"
output: html_notebook
---
```{r}
### import data
DRUP_dataset_details <- read.csv("~/katdetect/katdetect/results/DRUP_dataset_details", sep="")

DRUP_dataset_overall <- read.csv("~/katdetect/katdetect/results/DRUP_dataset_overall", sep="")

HMM_optimized_more_details <- read.csv("~/katdetect/katdetect/results/HMM_optimized_more_details", sep="")

########### met nMuts
layout(mat = matrix(c(1,2),2,1, byrow=TRUE),  height = c(1,8))
# Draw the boxplot and the histogram
par(mar=c(0, 3.1, 1.1, 2.1))
boxplot(DRUP_dataset_details$nMuts , horizontal=TRUE , ylim=c(0,350), main="Frequencies of the numbers of variants per detected kataegis event", xaxt="n" ,col= 'orange' , frame=F)
par(mar=c(6, 3.1, 1.1, 2.1))
hist(DRUP_dataset_details$nMuts , breaks=100 , col='orange' , border=T  ,ylab = 'Counts', xlim=c(0,350),ylim = c(0,12000),main = '')



######### met APOBEC
layout(mat = matrix(c(1,2),2,1, byrow=TRUE),  height = c(1,8))
# Draw the boxplot and the histogram
par(mar=c(0, 3.1, 1.1, 2.1))
boxplot(DRUP_dataset_details$p_APO , horizontal=TRUE , ylim=c(0,1), main="Frequencies of different percentages of APOBEC within detected kataegis events", xaxt="n" ,col= 'orange' , frame=F)
par(mar=c(6, 3.1, 1.1, 2.1))
hist(DRUP_dataset_details$p_APO , breaks=40 , col='orange' , border=T  , xlab="Percentage of APOBEC",ylab = 'Counts', xlim=c(0,1),ylim = c(0,12000))
```
