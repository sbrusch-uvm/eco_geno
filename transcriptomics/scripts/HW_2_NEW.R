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

length(degs_D18_BASE_D18_A28) # 41 differentially expressed genes between D18BASE and D18A28
length(degs_D18_BASE_D18_A33) # 332 DEGs between D18BASE and D18A33
length(degs_D22_BASE_D22_A28) # 289 DEGs between D22BASE and D22A28
length(degs_D22_BASE_D22_A33) # 1564 DEGs between D22BASE and D22A33


#############################
# we know overall number of genes, now we will figure out the overlap between 
#look at the overlaps in which genes are differentially expressed in multiple contrasts

length(intersect(degs_D18_BASE_D18_A28, degs_D18_BASE_D18_A33)) # 34 D18
length(intersect(degs_D22_BASE_D22_A28, degs_D22_BASE_D22_A33)) # 144 D22


# calculate the number of unique genes in each portion of the Euler plot
41-34 # 7 genes uniq genes at 18BASE18A28
332-34 # 298 uniq genes at 18BASE18A33

myEuler18 <- euler(c("D18A28"=7, "D18A33"=298, "D18A28&D18A33"=34))


plot(myEuler18, lty=1:2, quantities=TRUE, fill=c("lightblue", "darkseagreen", "aquamarine3"))

289-144 # 145 uniq genes at 22BASE22A28
1564-144 # 1420 

myEuler22 <- euler(c("D22A28"=145, "D22A33"=1420, "D22A28&D22A33"=144))

plot(myEuler22, lty=1:2, quantities=TRUE, fill=c("wheat", "indianred3", "lightcoral"))


########################################################
#Make a scatter plot of responses to A28 when copepods develop at 18vs22

#contrast D18_BASEvsA28

res_D18_BASEvsA28 <- as.data.frame(results(dds, contrast=c("group", "D18BASE", "D18A28"), 
                                           alpha = 0.05))
# contrast D22_BASEvsA28
res_D18_BASEvsA33 <- as.data.frame(results(dds, contrast=c("group", "D18BASE", "D18A33"), 
                                           alpha = 0.05))

# merge dataframes 
res_df18 <- merge(res_D18_BASEvsA28, res_D18_BASEvsA33, by="row.names",
                  suffixes=c(".28", ".33"))
rownames(res_df18) <- res_df18$Row.names
res_df18 <- res_df18[,-1]


library(dplyr)
library(tidyr)
# color based on values in our data frame
# define color mapping logic with the mutate function 

res_df18 <- res_df18 %>% 
  mutate(fill=case_when(
    padj.28 < 0.05 & stat.28 < 0 ~ "lightblue",
    padj.28 < 0.05 & stat.28 > 0 ~ "lightpink",
    padj.33 < 0.05 & stat.33 < 0 ~ "dodgerblue4",
    padj.33 < 0.05 & stat.33 > 0 ~ "deeppink4"
  ))

#count the number of points per fill color
color_counts <- res_df18 %>% 
  group_by(fill) %>% 
  summarise(count = n())
  
  
  

label_positions <- data.frame(
  fill=c("dodgerblue4", "lightpink", "deeppink4", "lightblue"),
  x_pos=c(1,5,0,-7.5),
  y_pos=c(-5,0,9,3)
)

label_data <- merge(color_counts, label_positions, by="fill")

plot18 <- ggplot(res_df18, aes(x=log2FoldChange.28, y=log2FoldChange.33, color=fill)) +
  geom_point(alpha=0.8)+
  scale_color_identity()+
  geom_text(data=label_data, aes(x=x_pos, y=y_pos, label=count, color=fill),
            size= 5)+
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black")+
  geom_abline(intercept = 0, slope = -1, linetype = "dashed", color = "black")+
  xlim(-10,10)+ylim(-10,10)+
  labs(x="Log2FoldChange 28 vs BASE at 18", 
       y="Log2FoldChange 33 vs BASE at 18",
       title = "How does response to 18C vary by FinalTemp?")+
  theme_minimal()
plot18
########################################
#Make a scatter plot of responses to A33 when copepods develop at 18vs22
#contrast D18_BASEvsA33

res_D22_BASEvsA28 <- as.data.frame(results(dds, contrast=c("group", "D22BASE", "D22A28"), 
                                           alpha = 0.05))
# contrast D22_BASEvsA28
res_D22_BASEvsA33 <- as.data.frame(results(dds, contrast=c("group", "D22BASE", "D22A33"), 
                                           alpha = 0.05))

# merge dataframes 
res_df22 <- merge(res_D22_BASEvsA28, res_D22_BASEvsA33, by="row.names",
                  suffixes=c(".28", ".33"))
rownames(res_df22) <- res_df22$Row.names
res_df22 <- res_df22[,-1]

# color based on values in our data frame
# define color mapping logic with the mutate function 

res_df22 <- res_df22 %>% 
  mutate(fill=case_when(
    padj.28 < 0.05 & stat.28 < 0 ~ "lightblue",
    padj.28 < 0.05 & stat.28 > 0 ~ "lightpink",
    padj.33 < 0.05 & stat.33 < 0 ~ "dodgerblue4",
    padj.33 < 0.05 & stat.33 > 0 ~ "deeppink4"
  ))

#count the number of points per fill color
color_counts22 <- res_df22 %>% 
  group_by(fill) %>% 
  summarise(count = n())

label_positions22 <- data.frame(
  fill=c("dodgerblue4", "lightpink", "deeppink4", "lightblue"),
  x_pos=c(1,5,0,-7.5),
  y_pos=c(-5,0,9,3)
)

label_data22 <- merge(color_counts22, label_positions22, by="fill")


plot22 <- ggplot(res_df22, aes(x=log2FoldChange.28, y=log2FoldChange.33, color=fill)) +
  geom_point(alpha=0.8)+
  scale_color_identity()+
  geom_text(data=label_data22, aes(x=x_pos, y=y_pos, label=count, color=fill),
            size= 5)+
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black")+
  geom_abline(intercept = 0, slope = -1, linetype = "dashed", color = "black")+
  xlim(-10,10)+ylim(-10,10)+
  labs(x="Log2FoldChange 28 vs BASE at 22", 
       y="Log2FoldChange 33 vs BASE at 22",
       title = "How does response to 22C vary by FinalTemp?")+
  theme_minimal()
plot22

#put the 2 scatter plots together
library(gridExtra)

combined_plot_18_22 <- grid.arrange(plot18, plot22, ncol=2)


ggsave("~/Projects/eco_geno/transcriptomics/figures/combined_plot_18_22.png", 
       combined_plot_18_22, width = 12, height = 6)



