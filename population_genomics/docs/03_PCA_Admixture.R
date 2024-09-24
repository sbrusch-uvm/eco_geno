library(tidyverse)
library(vcfR)
library(SNPfiltR)
library(LEA)

options(bitmapType = "cairo")

setwd("~/Projects/eco_geno/population_genomics/")

vcf <- read.vcfR("outputs/vcf_final.filtered.vcf.gz")

#we need to thin in the SNPs for LD (linkage disequilibrium) before 
#we reun PCA and Admixture analyses to satisfy the assumptions of independence among loci

vcf.thin <- distance_thin(vcf, min.distance = 500)
#new vcf file based on giving it the og vcf file and thin it 500 bp
#we learned most of our SNPs are stacked since 3646/15454 weren't located within 500bp of another SNP

meta <- read.csv("/gpfs1/cl/pbio3990/PopulationGenomics/metadata/meta4vcf.csv")
dim(meta)


meta2 <- meta[meta$id %in% colnames(vcf.thin@gt[,-1]), ]

dim(meta2)

write.vcf(vcf.thin, "outputs/vcf_final.filtered.thinned.vcf.gz")

#now we need to uncompress it
#hide the uncompressed vcf file too big for github outside of repo

system("gunzip -c ~/Projects/eco_geno/population_genomics/outputs/vcf_final.filtered.thinned.vcf.gz > ~/vcf_final.filtered.thinned.vcf")

geno <- vcf2geno(input.file = "/gpfs1/home/s/b/sbrusch/vcf_final.filtered.thinned.vcf", 
                 output.file = "outputs/vcf_final.filtered.thinned.geno")

CentPCA <- LEA::pca("outputs/vcf_final.filtered.thinned.geno", scale=TRUE)

plot(CentPCA$projections, 
     col=as.factor(meta2$region))
legend("bottomright", 
       legend=as.factor(unique(meta2$region)), 
                        fill=as.factor(unique(meta2$region)))

#x axis PC1 and y axis PC2


# 

