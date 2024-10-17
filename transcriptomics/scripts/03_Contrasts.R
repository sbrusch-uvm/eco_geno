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

# 2. compare baseline gene expression between developmental treatment groups
res_D18_A28_D22_A28 <- results(dds, contrast = c("group", "D18BASE", "D22BASE"), alpha = 0.05)
res_D18_A28_D22_A28 <-  res_D18_A28_D22_A28[!is.na(res_D18_A28_D22_A28$padj),]
res_D18_A28_D22_A28 <- res_D18_A28_D22_A28[order(res_D18_A28_D22_A28$padj),]
head(res_D18_A28_D22_A28)
summary(res_D18_A28_D22_A28) # 556 genes upregulated, 1379 downregulated genes

# make a list of which genes in our comparisons of interest are differential expressed (list of DEGs)
degs_D18_A28_D22_A28 <- row.names(res_D18_A28_D22_A28[res_D18_A28_D22_A28$padj < 0.05,])


plotMA(res_D18_A28_D22_A28, ylim=c(-4,4))



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


