library(tidyverse)
library(ggplot2)
library(vcfR)
library(qqman)
library(pcadapt)

vcf <- read.pcadapt("/gpfs1/cl/pbio3990/PopulationGenomics/variants/vcf_final.filtered.vcf",
                    type="vcf")
vcfR <- read.vcfR("/gpfs1/cl/pbio3990/PopulationGenomics/variants/vcf_final.filtered.vcf.gz")
#vcfR <- read.vcfR("~/Projects/eco_geno/population_genomics/outputs/vcf_final.filtered.vcf.gz")
meta <- read.csv("/gpfs1/cl/pbio3990/PopulationGenomics/metadata/meta4vcf.csv")
meta2 <- meta[meta$id %in% colnames(vcfR@gt[,-1]),]

pcadapt.pca <- pcadapt(vcf,
                       K=2,
                       method="componentwise",
                       min.maf = 0.01,
                       LD.clumping = list(size=500, thr=0.2))

#componentwise lets you get separate tests for dif axis
#min.maf = minor allele freq is for filtering out rare alleles
#linkage disequilibrium LD clumping allows for clumping, tests for selection of all the loci
## fits initial PC axis, still gets a test for every single SNP in our dataset 
## represent big peaks in manhattan plots

summary(pcadapt.pca)
# pass number is 9987 out of 15592 that passed our filter
plot(pcadapt.pca, options="scores",
     pop=meta2$region, 
     i=1, j=2,
     K=2)
#basically a manhattan plot w/o loci after filtering out low values
# lower numbers less significant 
# outliers in the high areas


summary(pcadapt.pca)
#the loadings is the vector on graph we drew, twice the amount
#zscores and chi2, test statistic to generate p value
#AF, MAF allele freq and minor allel freq 
#pass how many loci made it though the filter
plot(pcadapt.pca, options="stat.distribution",
     pop=meta2$region, 
     i=1, j=2,
     K=2)

View(head(vcfR@fix))
#all our info at fixed position
# use this data to understand the chromosome and position number from plot

vcfR.fix <- as.data.frame(vcfR@fix[,1:2])
head(vcfR.fix)

chr.main <- unique(vcfR.fix$CHROM)[1:8]
chrnum <- as.data.frame(cbind(chr.main, seq(1,8,1)))
head(chrnum)

#now extract p value
Pvalue <- pcadapt.pca$pvalues

pcadapt.MHplot <- cbind(vcfR.fix, Pvalue)

pcadapt.MHplot <- left_join(chrnum, pcadapt.MHplot, join_by(chr.main==CHROM))

pcadapt.MHplot <- pcadapt.MHplot %>%
  mutate(SNP=paste0(chr.main, "_", POS))
str(pcadapt.MHplot)

pcadapt.MHplot$V2 = as.numeric(pcadapt.MHplot$V2)
pcadapt.MHplot$POS = as.numeric(pcadapt.MHplot$POS)
pcadapt.MHplot$pPC1 = as.numeric(pcadapt.MHplot[,4])
pcadapt.MHplot$pPC2 = as.numeric(pcadapt.MHplot[,5])

pcadapt.MHplot <- pcadapt.MHplot %>% drop_na(pPC1)

manhattan(pcadapt.MHplot,
          chr="V2",
          bp="POS",
          p="pPC1",
          col=c("blue4", "orange"),
          logp = T,
          ylab = "-log10 p-value",
          genomewideline = F,
          main="PCAdapt genome scan for selection (PC1)")

manhattan(pcadapt.MHplot,
          chr="V2",
          bp="POS",
          p="pPC2",
          col=c("blue4", "orange"),
          logp = T,
          ylab = "-log10 p-value",
          genomewideline = F,
          main="PCAdapt genome scan for selection (PC2)")


View(pcadapt.MHplot %>% 
       filter(pPC1<quantile(pcadapt.MHplot$pPC1, 0.001)) %>%
     select(chr.main, POS, pPC1))
#show the most extreme pvalues, very strong evidence with low pvalue
# located near 6,641,802bp
#

# hist(pcadapt.pca$maf) gives histogram of the tail 
