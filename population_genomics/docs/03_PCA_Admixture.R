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


##################### broken off section from 9/24 that we are ignorning?
plot(CentPCA$projections, 
     col=as.factor(meta2$region))
legend("bottomright", 
       legend=as.factor(unique(meta2$region)), 
       fill=as.factor(unique(meta2$region)))

#x axis PC1 and y axis PC2



## edits for 9/26/24 
CentPCA <- load.pcaProject("vcf_final.filtered.thinned.pcaProject")

show(CentPCA)

plot(CentPCA)
#screeplot, = magnitude of PC values decending
#showing the eigenvalue = the axis put through the dots of, PC1 is the top dot
# PCA dots get smaller with each consecuative PCA plot 
#most PCs won't be interesting for us, only focus on the early ones

ggplot(as.data.frame(CentPCA$projections),
       aes(x=V1, y=V2, color=meta2$region, shape=meta2$continent))+
       geom_point(alpha=1) +
        labs(title="Centaurea genetic PCA", x="PC1", y="PC2", color="Region", shape="Continent")
# could be that the allele freq evolved very divergent than the other things
#makes very few assumptions
#first line says what data set to use
#aes assigns the PC lines to the axis with color coordination of region and shapes for continent

ggsave("figures/CentPCA_PC1vsPC2.pdf", width = 6, height = 6, units = "in")
#saved the last plot we made to the figures file

ggplot(as.data.frame(CentPCA$projections),
       aes(x=V2, y=V3, color=meta2$region, shape=meta2$continent))+
  geom_point(alpha=1) +
  labs(title="Centaurea genetic PCA", x="PC1", y="PC2", color="Region", shape="Continent")
#dominantly serperated along PC2 


#CentPCA$eigenvectors[1:5, 1:2] in console gives the head of V1 and V2 the eigenvectors
#each SNP locus has a quantitative value
#value -0.08, going from 0-2 gives a vector of -0.08 and 0.046 so its vector towards neg x axis and pos y axis

#CentPCA$eigenvalues[1:5]
#[1] 33998.20 15892.40 12440.10 10818.00  9887.98
#> sum(CentPCA$eigenvalues)
#[1] 1475410
#> CentPCA$eigenvalues[1]/sum(CentPCA$eigenvalues)
#[1] 0.02304322
#> CentPCA$eigenvalues[2]/sum(CentPCA$eigenvalues)
#[1] 0.01077152


#nice way to look at genetic structure w/o priority
#structure is another way to 


#Admixture
