#### 10/15/24, load in day 1 transcriptomics script and run DESeq object
options(bitmapType = "cairo")

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


# now look at counts of a specific top gene that we-re interested in to validate that the model is working 


d <- plotCounts(dds, gene = "TRINITY_DN140854_c0_g5_i2", int=(c("DevTemp", "FinalTemp")), returnData = TRUE)
d # for this specific gene, the count numbers per sample

p <- ggplot(d, aes(x=DevTemp, y=count, color=DevTemp, shape=FinalTemp)) +
  theme_minimal() +
  theme(text = element_text(size=20), panel.grid.major=element_line(color="grey")) 
p <- p + geom_point(position = position_jitter(w=0.2, h=0), size=3)
p # reality check of direction and magnitude 
# we would expect expression to be almost 2x lower in 22 than 18 
# great way to make sure table results are similar to visualizing the data 
        
        
        
######## MA plot = logfoldchange vs average gene expression
plotMA(res_D22vsD18, ylim=c(-4,4)) 
# x axis is the counts, y axis is log fold change, 0 is same in both samples
# we see a lot of upregulation at 22 in comparison to 18, showing mounted gene expression
# genes that are super highly expressed are far to the right ex. metabolism, replicating, DNA repair
# visual of over dispersion mean =! variance  

########## Volcano Plot
# convert our DESeq results object into a data frame to plot
res_df <- as.data.frame(res_D22vsD18)

# add a column to data frame to label whether a gene is significantly differentially expressed or not

res_df$Significant <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, "Significant", "Not Significant")

ggplot(res_df, aes(x=log2FoldChange, y=-log10(padj), color = Significant)) +
  geom_point(alpha=0.8)+
  scale_color_manual(values = c("lightblue", "tomato")) +
  labs(x="Log 2 Fold", y="-log10 Adjusted P-Value", title = "Volcano Plot")+
  theme_minimal()+
  theme(legend.position = "top")+
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="orange") +
  geom_vline(xintercept = c(-1,1), linetype = "dashed", color = "orange")

# much more upregulation in 22 in comparison to 18 
# tomato =  significant 
# if upregulated in 22 its downregulated in 22 
# applied additional criteria in the ifelse statement, not only significant but a doubling of gene expression
# each point is a gene


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
        
