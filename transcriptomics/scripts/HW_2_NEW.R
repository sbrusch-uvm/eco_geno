BiocManager::install("DESeq2", dependencies=TRUE, force = TRUE)

library(DESeq2) 
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

# the average number of counts per gene
rowSums(countsTableRound)
# this transcriptome wasn't made for this experiment which explains all the 0s
mean(rowSums(countsTableRound)) # 3244.739 reads per transcript
median(rowSums(countsTableRound)) # 64 , there are a lot of genes with very high expression and very very low expression
# evidence of over dispersion -> median is very different from the mean

apply(countsTableRound, 2, mean) # 2 is for columns, counts/reads at particular gene 
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
# thats 15/21 (75%) samples in that scenario, if there was a gene only majorly expressed in 1 treatment, filter out
# the 0 0 0 10000 12000 9000 in physical notebook (10/10 page)

nrow(dds) # went down to 35,527 transcripts from the original 119million in the countsTable thingy
# = number of transcripts with more than 10 reads in more than or equal to 15 samples

# Run the DESeq model to test for global differential gene expression 
dds <- DESeq(dds) #all the differential gene expression data now exists

#list the results you've generated with the function
resultsNames(dds) # Intercept, DevTemp_D22_vs_D18, FinalTemp_A33_vs_A28, FinalTemp_BASE_vs_A28


#those are our results
#now we can look into them
# visualize our global gene expression patterns using PCA
# first we need to transform the data for plotting using variance stabilization 

library(pheatmap) # will be making heat maps of gene expression data and how it varies across the samples
resultsNames(dds) # Intercept, DevTemp D22 vs D18, FinalTemp A33 vs A28, FinalTemp BASE vs A28

# pull out the results for developmental temp 22 vs 18
res_D22vsD18 <- results(dds, name = "DevTemp_D22_vs_D18", alpha = 0.05)

#find which ones within that grouping are the most significant
# order by significance 

res_D22vsD18 <- res_D22vsD18[order(res_D22vsD18$padj),]
head(res_D22vsD18) #gives the most sig dif in expression at samples at 18 vs 22
# always a comparison between 2 categories in this case = 22 vs 18 
# whats the gene expression in 22 relative to 18 
# the log2Fold Change, the larger the more sig
# really check the direction of comparison





######## MA plot = logfoldchange vs average gene expression
plotMA(res_D22vsD18, ylim=c(-4,4)) 
# x axis is the counts, y axis is log fold change, 0 is same in both samples
# we see a lot of upregulation at 22 in comparison to 18, showing mounted gene expression
# genes that are super highly expressed are far to the right ex. metabolism, replicating, DNA repair
# visual of over dispersion mean =! variance  



# heatmap is another way to look at change across samples and genes
vsd <- vst(dds, blind = FALSE)

topgenes <- head(rownames(res_D22vsD18), 20)
# there is way to much info to look at it all so we're gonna look at the top most sig. genes
mat <- assay(vsd)[topgenes, ]
df <- as.data.frame(colData(dds)[,c("DevTemp", "FinalTemp")])
pheatmap(mat, annotation_col=df, show_rownames=FALSE, cluster_cols=T, cluster_rows=T)
# columns = samples
# rows = genes
# color refers to expression amount, red band in middle shows high expression of genes across the samples
# some genes show contrast across our samples 
# blue block in middle corresponding with samples from D18, D22 samples have higher expression 


# we started by looking at the result names
# we chose to compare dev temps
# could choose to compare other things

library(eulerr)
#start by making groups within DESeq object

dds$group <- factor(paste0(dds$DevTemp, dds$FinalTemp))
design(dds) <- ~ group # groups by the above factors, new factor called group w/ all the new possible levels
dds <- DESeq(dds)
dim(dds) # 35527    21
resultsNames(dds) # "Intercept"  "group_D18A33_vs_D18A28"  "group_D18BASE_vs_D18A28"
# "group_D22A28_vs_D18A28"  "group_D22A33_vs_D18A28"  "group_D22BASE_vs_D18A28"


#now make separate files for each contrast, so that we have 1 object in our environment
# allows us to compare what is being contrasted in one group or the other group


# 1. compare baseline gene expression between developmental treatment groups
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



#############################

length(degs_D18_BASE_D18_A28) # 41 differentially expressed genes between D18BASE and D22BASE
length(degs_D18_BASE_D18_A33) # 332 DEGs between D18A28 and D22A28

# the higher the acute temperature exposure, the lower differentiated genes 

#############################
# we know overall number of genes, now we will figure out the overlap between 
#look at the overlaps in which genes are differentially expressed in multiple contrasts

length(intersect(degs_D18_BASE_D22_BASE, degs_D18_A28_D22_A28)) # 107
length(intersect(degs_D18_BASE_D22_BASE, degs_D18_A33_D22_A33)) # 44
length(intersect(degs_D18_A28_D22_A28, degs_D18_A33_D22_A33)) # 29
length(intersect(degs_D18_BASE_D22_BASE, 
                 intersect(degs_D18_A28_D22_A28, degs_D18_A33_D22_A33))) # 23


##### Pickup 10/22/24
# calculate the number of unique genes in each portion of the Euler plot
1935-107-44+23 # 1807 genes diff expressed uniquely at Baseline between 18vs22
296-107-29+23 # 183 genes diff expressed when exposed to 28
78-44-29+23 # 28 genes diff expressed when exposed to 33

107-23 # 84 genes unique to BASE and A28
44-23 # 21 genes unique to BASE and A33
29-23 # 6 genes unique to A28 and A33


myEuler <- euler(c("BASE"=1807, "A28"=183, "A33"=28, 
                   "BASE&A28"=84, "BASE&A33"=21, "A28&A33"=6,
                   "BASE&A28&A33"=23))
plot(myEuler, lty=1:2, quantities=TRUE, fill=c("red", "blue", "white"))

########################################################
#Make a scatter plot of responses to A28 when copepods develop at 18vs22

#contrast D18_BASEvsA28

res_D18_BASEvsA28 <- as.data.frame(results(dds, contrast=c("group", "D18BASE", "D18A28"), 
                                           alpha = 0.05))
# contrast D22_BASEvsA28
res_D22_BASEvsA28 <- as.data.frame(results(dds, contrast=c("group", "D22BASE", "D22A28"), 
                                           alpha = 0.05))

# merge dataframes 
res_df28 <- merge(res_D18_BASEvsA28, res_D22_BASEvsA28, by="row.names",
                  suffixes=c(".18", ".22"))
rownames(res_df28) <- res_df28$Row.names
res_df28 <- res_df28[,-1]


library(dplyr)
library(tidyr)
# color based on values in our data frame
# define color mapping logic with the mutate function 

res_df28 <- res_df28 %>% 
  mutate(fill=case_when(
    padj.18 < 0.05 & stat.18 < 0 ~ "lightblue",
    padj.18 < 0.05 & stat.18 > 0 ~ "lightpink",
    padj.22 < 0.05 & stat.22 < 0 ~ "darkblue",
    padj.22 < 0.05 & stat.22 > 0 ~ "magenta4"
  ))

#count the number of points per fill color
color_counts <- res_df28 %>% 
  group_by(fill) %>% 
  summarise(count = n())
  
  
  

label_positions <- data.frame(
  fill=c("darkblue", "lightpink", "magenta4", "lightblue"),
  x_pos=c(1,5,0,-7.5),
  y_pos=c(-5,0,9,3)
)

label_data <- merge(color_counts, label_positions, by="fill")

plot28 <- ggplot(res_df28, aes(x=log2FoldChange.18, y=log2FoldChange.22, color=fill)) +
  geom_point(alpha=0.8)+
  scale_color_identity()+
  geom_text(data=label_data, aes(x=x_pos, y=y_pos, label=count, color=fill),
            size= 5)+
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black")+
  geom_abline(intercept = 0, slope = -1, linetype = "dashed", color = "black")+
  xlim(-10,10)+ylim(-10,10)+
  labs(x="Log2FoldChange 28 vs BASE at 18", 
       y="Log2FoldChange 28 vs BASE at 22",
       title = "How does response to 28C vary by DevTemp?")+
  theme_minimal()
plot28
########################################
#Make a scatter plot of responses to A33 when copepods develop at 18vs22
#contrast D18_BASEvsA33

res_D18_BASEvsA33 <- as.data.frame(results(dds, contrast=c("group", "D18BASE", "D18A33"), 
                                           alpha = 0.05))
# contrast D22_BASEvsA28
res_D22_BASEvsA33 <- as.data.frame(results(dds, contrast=c("group", "D22BASE", "D22A33"), 
                                           alpha = 0.05))

# merge dataframes 
res_df33 <- merge(res_D18_BASEvsA33, res_D22_BASEvsA33, by="row.names",
                  suffixes=c(".18", ".22"))
rownames(res_df33) <- res_df33$Row.names
res_df33 <- res_df33[,-1]

# color based on values in our data frame
# define color mapping logic with the mutate function 

res_df33 <- res_df33 %>% 
  mutate(fill=case_when(
    padj.18 < 0.05 & stat.18 < 0 ~ "lightblue",
    padj.18 < 0.05 & stat.18 > 0 ~ "lightpink",
    padj.22 < 0.05 & stat.22 < 0 ~ "darkblue",
    padj.22 < 0.05 & stat.22 > 0 ~ "magenta4"
  ))

#count the number of points per fill color
color_counts33 <- res_df33 %>% 
  group_by(fill) %>% 
  summarise(count = n())

label_positions33 <- data.frame(
  fill=c("darkblue", "lightpink", "magenta4", "lightblue"),
  x_pos=c(1,5,0,-7.5),
  y_pos=c(-5,0,9,3)
)

label_data33 <- merge(color_counts33, label_positions33, by="fill")


plot33 <- ggplot(res_df33, aes(x=log2FoldChange.18, y=log2FoldChange.22, color=fill)) +
  geom_point(alpha=0.8)+
  scale_color_identity()+
  geom_text(data=label_data33, aes(x=x_pos, y=y_pos, label=count, color=fill),
            size= 5)+
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black")+
  geom_abline(intercept = 0, slope = -1, linetype = "dashed", color = "black")+
  xlim(-10,10)+ylim(-10,10)+
  labs(x="Log2FoldChange 33 vs BASE at 18", 
       y="Log2FoldChange 33 vs BASE at 22",
       title = "How does response to 33C vary by DevTemp?")+
  theme_minimal()
plot33

#put the 2 scatter plots together
library(gridExtra)

combined_scatter_plot <- grid.arrange(plot28, plot33, ncol=2)


ggsave("~/Projects/eco_geno/transcriptomics/figures/combined_scatter_plot.png", 
       combined_scatter_plot, width = 12, height = 6)



