BiocManager::install("DESeq2", dependencies=TRUE, force = TRUE)

library(DESeq2) #differential gene expression analysis
library(ggplot2) # lets us make plots
options(bitmapType = "cairo") # helps make plots
setwd("~/Projects/eco_geno/transcriptomics/") # where the data is located

#Import counts matrix
countsTable <- read.table("/gpfs1/cl/pbio3990/Transcriptomics/tonsa_counts.txt",
                          header = TRUE, row.names = 1)
# 21 samples, N new jersey, 1 = 18C, 2= 22C, C=control, s infront = resequenced
dim(countsTable)

countsTableRound <- round(countsTable) # DESeq2 prefers whole numbers

conds <- read.delim("/gpfs1/cl/pbio3990/Transcriptomics/experimental_details.txt", 
                    header = TRUE, stringsAsFactors = TRUE, row.names = 1) #the conditions
conds #sees all the data 21 samples, assigned developmental temp and final temp, 
## 2 levels Dev 18 and 22
## 3 levels Final Base, A28, A33

# Explore Counts Matrix #
colSums(countsTableRound) #that is summing up all the columns
mean(colSums(countsTableRound)) #this gives the mean of those columns 
# after filtering reads, 18mil reads is excellent

# the average number of counts per gene
rowSums(countsTableRound)
# this transcriptome wasn't made for this experiment which explains all the 0s
mean(rowSums(countsTableRound)) # 3244.739 reads per transcript
median(rowSums(countsTableRound)) # 64 , there are a lot of genes with very high expression and very very low expression
# evidence of over dispersion -> median is very different from the mean

apply(countsTableRound, 2, mean) # 2 is for columns, counts/reads at particular gene 
# one of the first steps DESeq does, gets a sense 
# gives a sense of variation in sequencing effort across samples 

# Analysis Using DESeq2

dds <- DESeqDataSetFromMatrix(countData= countsTableRound, colData = conds,
                              design = ~ DevTemp + FinalTemp)
#design for different expression is the DevTemp + FinalTemp

dim(dds) #dimension of dds

# now filtering
#average number of transcripts per number of reads

dds <- dds[rowSums(counts(dds) >= 10) >= 15, ] #filtering the transcripts 
nrow(dds) # went down to 35,527 transcripts from the original 119million in the countsTable
# = number of transcripts with more than 10 reads in more than or equal to 15 samples

# Run the DESeq model to test for global differential gene expression 
dds <- DESeq(dds) #all the differential gene expression data now exists
resultsNames(dds) # Intercept, DevTemp_D22_vs_D18, FinalTemp_A33_vs_A28, FinalTemp_BASE_vs_A28


library(eulerr)
#start by making groups within DESeq object

dds$group <- factor(paste0(dds$DevTemp, dds$FinalTemp))
design(dds) <- ~ group # groups by the above factors, new factor called group w/ all the new possible levels
dds <- DESeq(dds)
dim(dds) # 35527    21
resultsNames(dds) 

# 1. Compare baseline gene expression D18 between treatment groups A28 and A33
res_D18_BASE_D18_A28 <- results(dds, contrast = c("group", "D18BASE", "D18A28"), alpha = 0.05)
res_D18_BASE_D18_A28 <-  res_D18_BASE_D18_A28[!is.na(res_D18_BASE_D18_A28$padj),]
res_D18_BASE_D18_A28 <- res_D18_BASE_D18_A28[order(res_D18_BASE_D18_A28$padj),]
head(res_D18_BASE_D18_A28)
summary(res_D18_BASE_D18_A28) # 11 genes upregulated, 30 downregulated genes

# make a list of which genes in our comparisons of interest are differential expressed (list of DEGs)
degs_D18_BASE_D18_A28 <- row.names(res_D18_BASE_D18_A28[res_D18_BASE_D18_A28$padj < 0.05,])


plotMA(res_D18_BASE_D18_A28, ylim=c(-4,4))

# 2. compare baseline gene expression between developmental treatment groups
res_D18_BASE_D18_A33 <- results(dds, contrast = c("group", "D18BASE", "D18A33"), alpha = 0.05)
res_D18_BASE_D18_A33 <-  res_D18_BASE_D18_A33[!is.na(res_D18_BASE_D18_A33$padj),]
res_D18_BASE_D18_A33 <- res_D18_BASE_D18_A33[order(res_D18_BASE_D18_A33$padj),]
head(res_D18_BASE_D18_A33)
summary(res_D18_BASE_D18_A33) # 92 genes upregulated, 240 downregulated genes

# make a list of which genes in our comparisons of interest are differential expressed (list of DEGs)
degs_D18_BASE_D18_A33 <- row.names(res_D18_BASE_D18_A33[res_D18_BASE_D18_A33$padj < 0.05,])


plotMA(res_D18_BASE_D18_A33, ylim=c(-4,4))

# 3. compare baseline gene expression between developmental treatment groups
res_D22_BASE_D22_A28 <- results(dds, contrast = c("group", "D22BASE", "D22A28"), alpha = 0.05)
res_D22_BASE_D22_A28 <-  res_D22_BASE_D22_A28[!is.na(res_D22_BASE_D22_A28$padj),]
res_D22_BASE_D22_A28 <- res_D22_BASE_D22_A28[order(res_D22_BASE_D22_A28$padj),]
head(res_D22_BASE_D22_A28)
summary(res_D22_BASE_D22_A28) # 274 genes upregulated, 15 downregulated genes

# make a list of which genes in our comparisons of interest are differential expressed (list of DEGs)
degs_D22_BASE_D22_A28 <- row.names(res_D22_BASE_D22_A28[res_D22_BASE_D22_A28$padj < 0.05,])


plotMA(res_D22_BASE_D22_A28, ylim=c(-4,4))

# 4. compare baseline gene expression between developmental treatment groups
res_D22_BASE_D22_A33 <- results(dds, contrast = c("group", "D22BASE", "D22A33"), alpha = 0.05)
res_D22_BASE_D22_A33 <-  res_D22_BASE_D22_A33[!is.na(res_D22_BASE_D22_A33$padj),]
res_D22_BASE_D22_A33 <- res_D22_BASE_D22_A33[order(res_D22_BASE_D22_A33$padj),]
head(res_D22_BASE_D22_A33)
summary(res_D22_BASE_D22_A33) # 1176 genes upregulated, 388 downregulated genes

# make a list of which genes in our comparisons of interest are differential expressed (list of DEGs)
degs_D22_BASE_D22_A33 <- row.names(res_D22_BASE_D22_A33[res_D22_BASE_D22_A33$padj < 0.05,])


plotMA(res_D22_BASE_D22_A33, ylim=c(-4,4))

#############################

length(degs_D18_BASE_D18_A28) # 41 DEGs between blue
length(degs_D18_BASE_D18_A33) # 332 DEGs between green
length(degs_D22_BASE_D22_A28) # 289 DEGs between orange
length(degs_D22_BASE_D22_A33) # 1564 DEGs between red

length(intersect(degs_D18_BASE_D18_A28, degs_D18_BASE_D18_A33)) # 34 cool
length(intersect(degs_D22_BASE_D22_A28, degs_D22_BASE_D22_A33)) # 144 warm
length(intersect(degs_D18_BASE_D18_A28, degs_D22_BASE_D22_A28)) # 4 blue orange
length(intersect(degs_D18_BASE_D18_A28, degs_D22_BASE_D22_A33)) # 29 blue red
length(intersect(degs_D18_BASE_D18_A33, degs_D22_BASE_D22_A28)) # 14 green orange
length(intersect(degs_D18_BASE_D18_A33, degs_D22_BASE_D22_A33)) # 190 green red


nested_cool <- intersect(degs_D18_BASE_D18_A28, degs_D18_BASE_D18_A33)
nested_warm <- intersect(degs_D22_BASE_D22_A28, degs_D22_BASE_D22_A33)

length(intersect(degs_D22_BASE_D22_A28, nested_cool)) # 4 cool+orange
length(intersect(degs_D22_BASE_D22_A33, nested_cool)) # 29 cool+red
length(intersect(degs_D18_BASE_D18_A28, nested_warm)) # 4 warm + blue
length(intersect(degs_D18_BASE_D18_A33, nested_warm)) # 14 warm + green

length(intersect(nested_cool, nested_warm)) # 4 shared for all



# calculate the number of unique genes in each portion of the Euler plot
41-29-4-29-4-4+4 # blue = -25 -> 0
289-4-14-4-4-14+4 # orange = 253
332-14-190-4-14-29+4 # green = 85
1564-29-190-4-29-14+4 # red = 1302

29-4-29+4 # blue/red = 0
4-4-4+4 # blue/orange = 0
14-4-14+4 # green/orange = 0
190-29-14+4 # green/red = 151

29-4 #cool/red = 25
4-4 # cool/orange = 0
4-4 # warm/blue = 0
14-4 # warm/green = 10


EulerPlot <- euler(c("D18A28"=0, "D22A28"=85, "D18A33"=253, "D22A33"=1302, 
                     "D18A28&D22A28"=0, "D18A33&D22A28"=0, "D18A33&D22A33"=151, "D18A28&D22A33"=0,
                     "D18A28&D22A28&D22A33"=0, "D18A28&D18A33&D22A28"=0, "D18A33&D22A28&D22A33"=10, "D18A28&D18A33&D22A33"=25,
                     "D18A28&D18A33&D22A28&D22A33"=4))
plot(EulerPlot, lty=1:4, quantities=T)


myEuler <- euler(c("BASE"=1807, "A28"=183, "A33"=28, 
                   "BASE&A28"=84, "BASE&A33"=21, "A28&A33"=6,
                   "BASE&A28&A33"=23))
plot(myEuler, lty=1:3, quantities=TRUE)
