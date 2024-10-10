### code for analyzing RNAseq data using DESeq2 as the package
# load libraries 

BiocManager::install("DESeq2", dependencies=TRUE, force = TRUE)

library(DESeq2) #not working
library(ggplot2)
options(bitmapType = "cairo")

setwd("~/Projects/eco_geno/transcriptomics/")


#Import counts matrix 

countsTable <- read.table("/gpfs1/cl/pbio3990/Transcriptomics/tonsa_counts.txt",
                          header = TRUE, row.names = 1)
# 21 samples, N new jersey, 1 = 18C, 2= 22C, C=control, sbefore = resequenced
dim(countsTable)

countsTableRound <- round(countsTable)
#DESeq2 doesn't like decimals so we are rounding all the values in the matrix
tail(countsTableRound)
# TRINITY_DN219251 shows the number of counts mapped to gene 


conds <- read.delim("/gpfs1/cl/pbio3990/Transcriptomics/experimental_details.txt",
                    header = TRUE, stringsAsFactors = TRUE, row.names = 1)
#the conditions
head(conds) #heads the data
conds #sees all the data 21 samples, assigned developmental temp and final temp, 
## 2 levels Dev 18 and 22
## 3 levels Final Base, A28, A33
##################################################
# 
# Explore counts matrix 
# 
##################################################
# let's see how many reads we have from each sample do that with colSums (a base R funciton)

colSums(countsTableRound) #that is summing up all the columns, look at notebook 10/10 page
mean(colSums(countsTableRound)) #this gives the mean of those columns 
# what you typically want in an RNA seq study 20million reads, after filtering reads, 18mil reads is excellent


barplot(colSums(countsTableRound), names.arg = colnames(countsTableRound),
        cex.names = 0.5, las = 2, ylim = c(0,30000000)) # cex.names is so that the names fit 
abline(h=mean(colSums(countsTableRound)), col = "blue4", lwd=2) 
# that gives the average as a line horizontally on the plot

# the average number of counts per gene
rowSums(countsTableRound)
# this transcriptome wasn't made for this experiment which explains all the 0s
mean(rowSums(countsTableRound)) # 3244.739 reads per transcript
median(rowSums(countsTableRound)) # 64 , there are a lot of genes with very high expression and very very low expression
# evidence of over dispersion that median is very different from the mean

apply(countsTableRound, 2, mean) #2 is for columns, counts/reads at particular gene 
# one of the first steps DESeq does, gets a sense 
# gives a sense of variation in sequencing effort across samples 

#############################################################
#
# Start analysis in DESeq2
#
#############################################################

dds <- DESeqDataSetFromMatrix(countData= countsTableRound, colData = conds,
                              design = ~ DevTemp + FinalTemp)
#design for different expression is the DevTemp + FinalTemp

dim(dds) #dimension of dds

# now filtering
#average number of transcripts per number of reads

dds <- dds[rowSums(counts(dds) >= 10) >= 15, ] #filtering the transcripts by doing the rowSums and if you sum across
# if there aren't more than 10 reads in that for at least 15 samples, they we don't want to pay attention to that
# thats 15/21 (75%) samples in that scenario, if there was a gene only majorily expressed in 1 treatment, filter out
# the 000 10000 12000 9000 in physical notebook

nrow(dds) # went down to 35,527 transcripts from the original 119million in the countsTable thingy
# = number of transcripts with more than 10 reads in more than or equal to 15 samples

# Run the DESeq model to test for global differential gene expression 
dds <- DESeq(dds) #all the differential gene expression data now exists

#list the results you've generated with the funtion
resultsNames(dds) # Intercept, DevTemp D22 vs D18, FinalTemp A33 vs A28, FinalTemp BASE vs A28
#those are our results
#now we can look into them
# visualize our global gene expression patterns using PCA
# first we need to transform the data for plotting using variance stabilization 

vsd <- vst(dds, blind=FALSE)

#start building plot
pcaData <- plotPCA(vsd, intgroup=c("DevTemp", "FinalTemp"), returnData=TRUE)
percentVar <- round(100*attr(pcaData, "percentVar")) # percent explained PC1 49, PC2 15

final_temp_colors <- c("BASE" = "grey", "A28"="hotpink", "A33" = "red")
shapes_choose <- c("D18" = 16, "D22" = 18)

p <- ggplot(pcaData, aes(PC1, PC2, color = FinalTemp, shape = DevTemp)) +
  geom_point(size=5) +
  scale_shape_manual(values = shapes_choose) +
  scale_color_manual(values = final_temp_colors)+
  labs(x = paste0('PC1: ', percentVar[1], "%"), 
       y = paste0('PC2: ', percentVar[2], "%")) +
  theme_bw(base_size = 16)
p  

# what do we learn about our data? ->   



