## Script for analyzing and visualizing gene correlation networks

library(DESeq2) 
library(ggplot2)
library(WGCNA); options(stringsAsFactors=FALSE);
library(tidyverse)
library(CorLevelPlot)
library(gridExtra)
library(Rmisc)


options(bitmapType = "cairo")

setwd("~/Projects/eco_geno/transcriptomics/")

#Step 1: Importing our counts data

countsTable <- read.table("/gpfs1/cl/pbio3990/Transcriptomics/tonsa_counts.txt",
                          header = TRUE, row.names = 1)
tail(countsTable)
dim(countsTable)

countsTableRound <- round(countsTable)
tail(countsTableRound)

conds <- read.delim("/gpfs1/cl/pbio3990/Transcriptomics/experimental_details.txt",
                    header = TRUE, stringsAsFactors = TRUE, row.names = 1)
head(conds)

traitData <- read.table("/gpfs1/cl/pbio3990/Trait_Data.txt", header = T, row.names = 1)

#filter the matrix to just BASE data (because those are the data for which we have traits measured)
filtered_count_matrix_BASEonly <- countsTable[,conds$FinalTemp == "BASE"]
filtered_sample_metadata_BASEonly <- conds[conds$FinalTemp == "BASE", ]
rounded_filtered_count_matrix <- round(filtered_count_matrix_BASEonly)

#Step 2: Detecting Outliers
#detect outlier genes
gsg <- goodSamplesGenes(t(rounded_filtered_count_matrix))
summary(gsg) # everything according to this assessment seems fine

table(gsg$goodGenes) # 82,203 good genes, 37,235 bad genes
table(gsg$goodSamples) # all 7 samples are true 

#filter our bad genes
data_WGCNA <- rounded_filtered_count_matrix[gsg$goodGenes == TRUE,]
dim(data_WGCNA) #true 82,203 and 7 

# use clustering with a tree dendrogram to identify outlier samples
htree <- hclust(dist(t(data_WGCNA)), method = "average")
plot(htree) #a homework option is to remove the outlier to see how it affects the genes

#PCA - outlier detection method
pca <- prcomp(t(data_WGCNA))
pca_data <- pca$x
#make a dataframe
pca_data <- as.data.frame(pca_data)
pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)

ggplot(pca_data, aes(PC1, PC2))+
  geom_point()+
  geom_text(label = rownames(pca_data))+
  labs(x=paste0("PC1: ", pca.var.percent[1], " %"), 
       y=paste0("PC2: ", pca.var.percent[2], " %"))

# Step 3: Normalization 

colData <- row.names(filtered_sample_metadata_BASEonly)

#run DESeq without any model defined, want WGCNA to cluster without bias
dds_WGCNA <- DESeqDataSetFromMatrix(countData = data_WGCNA, 
                                    colData = filtered_sample_metadata_BASEonly,
                                    design = ~1) #there are no specified groups
dds_WGCNA_75 <- dds_WGCNA[rowSums(counts(dds_WGCNA) >= 15) >=6,]
nrow(dds_WGCNA_75) # filtered down to 29,559 transcripts

dds_norm <- vst(dds_WGCNA_75) # perform variance stabilization

#get and save normalized counts to use below
norm.counts <- assay(dds_norm) %>% 
  t()


# Step 4: Network Construction
#choose a set of soft-thresholding powers

power <- c(c(1:10), seq(from = 12, to = 50, by = 2))

#call the network topology analysis function (takes a couple minutes to run)

sft <- pickSoftThreshold(norm.counts,
                         powerVector = power,
                         networkType = "signed",
                         verbose = 5)
sft.data <- sft$fitIndices
#plot to pick power

a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power))+
  geom_point()+
  geom_text(nudge_y = 0.1)+
  geom_hline(yintercept = 0.8, color = "red")+
  labs(x="Power", y="Scale free topology model fit, signed R^2")+
  theme_classic()



a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power))+
  geom_point()+
  geom_text(nudge_y = 0.1)+
  geom_hline(yintercept = 0.8, color = "red")+
  labs(x="Power", y="Mean Connectivity")+
  theme_classic()

grid.arrange(a1, a2, nrow = 2)

#WGNCA wants us to pick a degree of relatedness
#as connectivity decreases the other increases 








