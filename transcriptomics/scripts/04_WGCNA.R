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

#WGNCA wants us to pick a degree of relatedness as connectivity decreases the other increases 

#scale free topology is the biological connections, the plots we make show strength of correlation between the different gene = power, how do those show a scale free topology 


############ 10/29/2024 ############
# x axis is the strength between the genes, maximizing the top y graph
# second one is so strong all genes are connected
# we will choose 26
# hw is that you could choose different number to see the downstream affects

soft_power <- 26
temp_cor <- cor
cor <- WGCNA::cor # this sets the temp_cor function to use WGCNA's correlation function
#convert our data matrix to numeric

norm.counts[] <- sapply(norm.counts, as.numeric)

# the command below creates the network and identifies modules based on the parameters that we chose
bwnet26 <- blockwiseModules(norm.counts, 
                            maxBlockSize = 30000, 
                            TOMType = "signed",
                            power = soft_power,
                            mergeCutHeight = 0.25,
                            numericLabels = FALSE,
                            randomSeed = 1234,
                            verbose = 3)
cor <- temp_cor #resets the cor function to base R's cor fxn instead of using WGCNA's cor fxn

saveRDS(bwnet26, file = "outputs/bwnet26.rds")

# to load the bwnet file in later use:
## bwnet26 <- readRDS("outputs/bwnet26.rds")


# Step 5) Explore Module Eigengenes 

module_eigengenes <- bwnet26$MEs
head(module_eigengenes)
dim(module_eigengenes)

# Get the number of genes for each module using table funciton
table(bwnet$26$colors)

# Plot the dendrogram and the module colors 
# merge and unmerged, merged based on similarity cutoff (0.25 which we set earlier)

plotDendroAndColors(bwnet26$dendrograms[[1]], cbind(bwnet26$unmergedColors, bwnet26$colors), c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang = 0.03,
                    guideHang = 0.05)


# Step 6) Correlation of modules with traits!

## define the numbers of genes and samples
nSamples <- nrow(norm.counts)
nGenes <- ncol(norm.counts)

## test for a correlation between module eigengenes and trait data
module.trait.corr <- cor(module_eigengenes, traitData, use = "p")

## calculate p valyes for each correlation
module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples)

## utalize module-trait association as a heatmap
heatmap.data <- merge(module_eigengenes, traitData, 
                      by = "row.names")
head(heatmap.data)

### address error of row.names not being numeric
heatmap.data <- heatmap.data %>% 
  column_to_rownames(var = "Row.names")

names(heatmap.data)

## make pretty heatmap of correlations
CorLevelPlot(heatmap.data, 
             x = names(heatmap.data)[46:48], # these values may need to change based on number of eigengenes
             y = names(heatmap.data)[1:45],
             col = c("midnightblue", "deepskyblue", "white", "rosybrown1", "firebrick"))




