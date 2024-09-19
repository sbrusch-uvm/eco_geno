#estimating diversity and genetic differentiation in the filtered centaurea data

library (vcfR)
library(tidyverse)
library(qqman)

#helps solve plotting issues
X11.options(type="cairo")


# read in our VCF file from our repo outputs directory

vcf <- read.vcfR("/gpfs1/cl/pbio3990/PopulationGenomics/variants/Centaurea_filtered.vcf.gz")

# read in our metadata -- info on population of origin

meta <- read.csv("/gpfs1/cl/pbio3990/PopulationGenomics/metadata/meta4vcf.csv")

head(meta) # vcf files has 595 samples
dim(meta) #meta has 629 inds

meta2 <- meta[meta$id %in% colnames(vcf@gt[,-1]),]
#colnames(vcf@gt[,-1]),])[1:5] will give you just the first 5 collum names 

#calculate diversity stats using the genetic_diff fxn in vcfR
vcf.div <- genetic_diff(vcf, pops = as.factor(meta2$region), method = "nei")

