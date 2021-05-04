### R code from vignette source 'Rsubread.Rnw'
#BiocManager::install("Rsubread")
###################################################
### code chunk number 1: Rsubread.Rnw:70-73
###################################################
library(Rsubread)
ref <- system.file("extdata","reference.fa",package="Rsubread")
buildindex(basename="reference_index",reference=ref)

#chose you own file:
#ref <- "insert your file directory here"
#newreads <- "insert your file directory here"

###################################################
### code chunk number 2: Rsubread.Rnw:90-92
###################################################
#import fastq directory
reads <- system.file("extdata","reads.txt.gz",package="Rsubread")

#align reads to the ref
align.stat <- align(index="reference_index",
                      readfile1=reads,
                      output_file="alignResults.BAM",
                      phredOffset=64)

###################################################
### code chunk number 3: Rsubread.Rnw:99-103
###################################################
reads1 <- system.file("extdata","reads1.txt.gz",package="Rsubread")
reads2 <- system.file("extdata","reads2.txt.gz",package="Rsubread")

align.stat2 <- align(index="reference_index",readfile1=reads1,readfile2=reads2,
                     output_file="alignResultsPE.BAM",phredOffset=64)


###################################################
### code chunk number 4: Rsubread.Rnw:125-135
###################################################
ann <- data.frame(
  GeneID=c("gene1","gene1","gene2","gene2"),
  Chr="chr_dummy",
  Start=c(100,1000,3000,5000),
  End=c(500,1800,4000,5500),
  Strand=c("+","+","-","-"),
  stringsAsFactors=FALSE)
ann

fc_SE <- featureCounts("alignResults.BAM",annot.ext=ann)
fc_SE


###################################################
### code chunk number 5: Rsubread.Rnw:142-144
###################################################
fc_PE <- featureCounts("alignResultsPE.BAM",annot.ext=ann,isPairedEnd=TRUE)
fc_PE

####################################################################
#https://hbctraining.github.io/DGE_workshop/lessons/04_DGE_DESeq2_analysis.html

#BiocManager::install("parathyroidSE")
library("parathyroidSE")
data( "parathyroidGenesSE" )

library(DESeq2)

#Creating a Dummry Count Matrix to contrast with exsisting CountMatrix
dummy_cm <- cbind(fc_SE$counts,fc_PE$counts)
colnames(dummy_cm) <- c("Sample1","Sample2")
countmatrix <- assay(parathyroidGenesSE)

#Import ColData and rename Countmatrix
coldata <- colData(parathyroidGenesSE)
experimentalData <- as.data.frame(coldata)

#rename the count matrix and col data
rownames( coldata ) <- coldata$run
experimentalData <- as.data.frame(coldata)

colnames( countmatrix ) <- coldata$run

ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = countmatrix,
  colData = coldata,
  design =  ~time+treatment)

dds <- DESeq(ddsFullCountTable)

###########################################################################

plotDispEsts(dds, ylim = c(1e-6, 1e2) )

res <- results( dds )
res

dataframeresults <- data.frame(res@listData)
rownames(dataframeresults) <- rownames(countmatrix)
  
plotMA( res, ylim = c(-3, 3) )

hist( res$pvalue, breaks=20, col="grey" )
hist( res$padj, breaks=20, col="grey" )

summary(res)

##Using contrast,results(dds, contrast=c("condition","C","B"))
#https://rdrr.io/bioc/DESeq2/man/results.html
res2 <- results( dds, contrast = c("treatment", "OHT", "Control"))
res2
plotMA( res2, ylim = c(-3, 3) )

hist( res2$pvalue, breaks=20, col="grey" )
hist( res2$padj, breaks=20, col="grey" )

summary(res2)

#Performing Shrink
#Approximate posterior estimation for GLM coefficients
resultsNames(dds)
shrink_res2 <- lfcShrink(dds,
                        coef=4,
                        res=res2)

plotMA(shrink_res2 , ylim = c(-3, 3) )

summary(shrink_res2)

######################################################################
res3 <- results(dds, contrast = c("treatment", "DPN", "Control"), alpha=0.998)

plotMA(res3 , ylim = c(-3, 3) )

hist( res3$pvalue, breaks=20, col="grey" )
hist( res3$padj, breaks=20, col="grey" )

#########################################################################
#install.packages("devtools")
#devtools::install_github("stephenturner/annotables")

library(annotables)
library(tidyverse)
grch38 <- grch38

res3_export <- subset(res3, pvalue < 0.05) %>%
                              data.frame() %>%
                              rownames_to_column(var = "geneID")

colnames(res3_export)[1] <-  colnames(grch38)[1]

res3_export3 <- left_join(res3_export,grch38,by="ensgene")

write.csv( as.data.frame(res3_export), file="ensgene" )

