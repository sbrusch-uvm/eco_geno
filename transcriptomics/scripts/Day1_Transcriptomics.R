### code for analyzing RNAseq data using DESeq2 as the package
# load libraries 

BiocManager::install("DESeq2", dependencies=TRUE, force = TRUE)

library(DESeq2) #not working


library(ggplot2)

setwd("~/Projects/eco_geno/transcriptomics/")

#Import counts matrix 

countsTable <- read.table("/gpfs1/cl/pbio3990/Transcriptomics/tonsa_counts.txt",
                          header = TRUE, row.names = 1)
# 21 samples, N new jersey, 1 = 18C, 2= 22C, C=control, sbefore = resequenced

countsTableRound <- round(countsTable)
#DESeq2 doesn't like decimals so we are rounding all the values in the matrix
tail(countsTableRound)
# TRINITY_DN219251 shows the number of counts mapped to gene 

conds <- read.delim("/gpfs1/cl/pbio3990/Transcriptomics/experimental_details.txt",
                    header = TRUE, stringsAsFactors = TRUE, row.names = 1)
#the conditions
head(conds) #heads the data
conds #sees all the data, assigned developmental temp and final temp
## 2 levels Dev 18 and 22
## 3 levels Final Base, A28, A33

dds <- DESeqDataSetFromMatrix(countData= countsTableRound, colData = conds,
                              design = ~ DevTemp + FinalTemp)




