setwd("~/Projects/eco_geno/population_genomics/docs/")
#this allowed me to save this file to my repo
library(vcfR)
setwd("/gpfs1/cl/pbio3990/PopulationGenomics/")
list.files()

vcf <- read.vcfR("variants/Centaurea_filtered.vcf.gz")
head(vcf)

X11.options(type="cairo")

dna <- ape::read.dna("reference/GCA_030169165.1_ASM3016916v1_genomic.fa.gz", format = "fasta")
gff <- read.table("reference/GCA_030169165.1_ASM3016916v1_genomic_genes.gff.gz", sep = "\t", quote = "")


DP <- extract.gt(vcf, element="DP", as.numeric=T)
DP[1:5, 1:10]

quantile(DP)

DP[DP==0] <- NA

quantile(DP, na.rm=T)

heatmap.bp(DP[1:1000,], rlabels = F, clabels = F)

library(SNPfiltR)



vcf.filt <- hard_filter(vcf, depth=3)


max_depth(vcf.filt)
vcf.filt <- max_depth(vcf.filt, maxdepth=60) 

meta <- read.csv("metadata/meta4vcf.csv", header = T)
meta2 <- meta[,c(1,4)]

names(meta2) <- c("id", "pop")

meta2$id = as.factor(meta2$id)
meta2$pop = as.factor(meta2$pop)

vcf.filt.indMiss <- missing_by_sample(vcf.filt, 
                                      popmap = meta2, 
                                      cutoff = 0.8) #for second graph change value to 0.55


vcf.filt.indMiss <- filter_biallelic(vcf.filt.indMiss)
vcf.filt.indMiss <- min_mac(vcf.filt.indMiss, min.mac = 1)


vcf.filt.indSNPMiss <- missing_by_snp(vcf.filt.indMiss, cutoff = 0.5)

DP2 <- extract.gt(vcf.filt.indSNPMiss, element="DP", as.numeric=T)

heatmap.bp(DP2[1:5000,], 
           rlabels=F, clabels=F)

write.vcf(vcf.filt.indSNPMiss,
          "~/Projects/eco_geno/population_genomics/outputs/vcf_final.filtered0.8.vcf.gz")
# for second filtering save, I changed the name to vcf_final_filtered0.55.vcf.gz


##############################################################################

library(tidyverse)
library(vcfR)
library(SNPfiltR)
library(LEA)

options(bitmapType = "cairo")

setwd("~/Projects/eco_geno/population_genomics/")

vcf <- read.vcfR("outputs/vcf_final.filtered0.8.vcf.gz")

#we need to thin in the SNPs for LD (linkage disequilibrium) before 
#we reun PCA and Admixture analyses to satisfy the assumptions of independence among loci

vcf.thin <- distance_thin(vcf, min.distance = 500)
#new vcf file based on giving it the og vcf file and thin it 500 bp
#we learned most of our SNPs are stacked since 3646/15454 weren't located within 500bp of another SNP

meta <- read.csv("/gpfs1/cl/pbio3990/PopulationGenomics/metadata/meta4vcf.csv")
dim(meta)


meta2 <- meta[meta$id %in% colnames(vcf@gt[,-1]), ]
## edited vcf.thin to vcf 10/1/24
dim(meta2)

write.vcf(vcf.thin, "outputs/vcf_final.filtered.thinned.vcf.gz")

#now we need to uncompress it
#hide the uncompressed vcf file too big for github outside of repo

system("gunzip -c ~/Projects/eco_geno/population_genomics/outputs/vcf_final.filtered.thinned.vcf.gz > ~/vcf_final.filtered.thinned.vcf")

CentPCA <- LEA::pca("outputs/vcf_final.filtered.thinned.geno", scale=TRUE)


plot(CentPCA$projections, 
     col=as.factor(meta2$region))
legend("bottomright", 
       legend=as.factor(unique(meta2$region)), 
       fill=as.factor(unique(meta2$region)))


CentPCA <- load.pcaProject("vcf_final.filtered.thinned.pcaProject")

show(CentPCA)

plot(CentPCA)

ggplot(as.data.frame(CentPCA$projections),
       aes(x=V1, y=V2, color=meta2$region, shape=meta2$continent))+
  geom_point(alpha=1) +
  labs(title="Centaurea genetic PCA", x="PC1", y="PC2", color="Region", shape="Continent")
# cou