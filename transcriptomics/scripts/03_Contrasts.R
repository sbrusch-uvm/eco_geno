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
res_D18_BASE_D22_BASE <- results(dds, contrast = c("group", "D18BASE", "D22BASE"), alpha = 0.05)
res_D18_BASE_D22_BASE <-  res_D18_BASE_D22_BASE[!is.na(res_D18_BASE_D22_BASE$padj),]
res_D18_BASE_D22_BASE <- res_D18_BASE_D22_BASE[order(res_D18_BASE_D22_BASE$padj),]
head(res_D18_BASE_D22_BASE)
summary(res_D18_BASE_D22_BASE) # 556 genes upregulated, 1379 downregulated genes

# make a list of which genes in our comparisons of interest are differential expressed (list of DEGs)
degs_D18_BASE_D22_BASE <- row.names(res_D18_BASE_D22_BASE[res_D18_BASE_D22_BASE$padj < 0.05,])


plotMA(res_D18_BASE_D22_BASE, ylim=c(-4,4))

# 2. compare baseline gene expression between developmental treatment groups
res_D18_A28_D22_A28 <- results(dds, contrast = c("group", "D18A28", "D22A28"), alpha = 0.05)
res_D18_A28_D22_A28 <-  res_D18_A28_D22_A28[!is.na(res_D18_A28_D22_A28$padj),]
res_D18_A28_D22_A28 <- res_D18_A28_D22_A28[order(res_D18_A28_D22_A28$padj),]
head(res_D18_A28_D22_A28)
summary(res_D18_A28_D22_A28) # 146 genes upregulated, 150 downregulated genes

# make a list of which genes in our comparisons of interest are differential expressed (list of DEGs)
degs_D18_A28_D22_A28 <- row.names(res_D18_A28_D22_A28[res_D18_A28_D22_A28$padj < 0.05,])


plotMA(res_D18_A28_D22_A28, ylim=c(-4,4))

# 3. compare baseline gene expression between developmental treatment groups
res_D18_A33_D22_A33 <- results(dds, contrast = c("group", "D18A33", "D22A33"), alpha = 0.05)
res_D18_A33_D22_A33 <-  res_D18_A33_D22_A33[!is.na(res_D18_A33_D22_A33$padj),]
res_D18_A33_D22_A33 <- res_D18_A33_D22_A33[order(res_D18_A33_D22_A33$padj),]
head(res_D18_A33_D22_A33)
summary(res_D18_A33_D22_A33) # 47 genes upregulated, 31 downregulated genes

# make a list of which genes in our comparisons of interest are differential expressed (list of DEGs)
degs_D18_A33_D22_A33 <- row.names(res_D18_A33_D22_A33[res_D18_A33_D22_A33$padj < 0.05,])


plotMA(res_D18_A33_D22_A33, ylim=c(-4,4))



#############################

length(degs_D18_BASE_D22_BASE) # 1,935 differentially expressed genes between D18BASE and D22BASE
length(degs_D18_A28_D22_A28) # 296 DEGs between D18A28 and D22A28
length(degs_D18_A33_D22_A33) # 78 DEGs between D18A33 and D22A33

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
plot(myEuler, lty=1:3, quantities=TRUE)

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

ggplot(res_df28, aes(x=log2FoldChange.18, y=log2FoldChange.22, color=fill)) +
  geom_point(alpha=0.8)+
  scale_color_identity()+
  geom_text(data=label_data, aes(x=x_pos, y=y_pos, label=count, color=fill),
            size= 5)+
  labs(x="Log2FoldChange 28 vs BASE at 18", 
       y="Log2FoldChange 28 vs BASE at 22",
       title = "How does response to 28C vary by DevTemp?")+
  theme_minimal()

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


library(dplyr)
library(tidyr)
# color based on values in our data frame
# define color mapping logic with the mutate function 

res_df33 <- res_df33 %>% 
  mutate(fill=case_when(
    padj.18 < 0.05 & stat.18 < 0 ~ "lightblue",
    padj.18 < 0.05 & stat.18 > 0 ~ "lightpink",
    padj.22 < 0.05 & stat.22 < 0 ~ "darkblue",
    padj.22 < 0.05 & stat.22 > 0 ~ "magenta4"
  ))

ggplot(res_df33, aes(x=log2FoldChange.18, y=log2FoldChange.22, color=fill)) +
  geom_point(alpha=0.8)+
  scale_color_identity()+
  labs(x="Log2FoldChange 33 vs BASE at 18", 
       y="Log2FoldChange 33 vs BASE at 22",
       title = "How does response to 33C vary by DevTemp?")+
  theme_minimal()



