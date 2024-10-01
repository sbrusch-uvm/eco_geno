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


meta2 <- meta[meta$id %in% colnames(vcf@gt[,-1]), ]
## edited vcf.thin to vcf 10/1/24
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

## 10/1/24 ##
# loaded in all the important stuff
#Run admixture analyses and create plots 
#for admixture use LEA R package
#function instead LEA is snmf
## runs admixture analyses similar to structure, which runs by basian
## structure takes forever to run on large data sets
## snmf is hella fast, provides very comparable estimates to structure

CentAdmix <- snmf("outputs/vcf_final.filtered.thinned.geno",
                  K=1:10,
                  entropy = T,
                  repetitions = 3,
                  project = "new") # if adding to analysis later, you could choose project="continue", you could add replicates at a different time without having to rewrite
  #somewhere to hold data
# takes the genofile that was produced during the PCA run 
# first have to tell it a K value, give it a range of K
## unlikely that 1 value of k will capture
# make it do the cross validation = mask a portion of data, estimate the admixture model, see how well it predicts the held back data
## entropy= true
## the algorith doesn't get the same value every time
### run repetitions of each value of K = repetitions so that we can choose the best one


par(mfrow=c(2,1))
plot(CentAdmix, col="blue4", main="SNMF") 
# if the dots dont go down, increase K
# K is on x-axis
# look for the elbow in the plot, think of a range of values in the crook for us 3,4,5
# look at multiple K (3,4,5) since a range can be beneficial in understanding
# look for least complicated model to fit your data 
# this plots the cross-entropy score we can use for selecting models that fit our data well
# the other thing we can do is plot PCA
plot(CentPCA$eigenvalues[1:10], ylab="Eigenvalues", xlab="Number of PCs", col="blue4", main="PCA")
#by specifying it, get something more similar to the CentAdmix plot


dev.off()
#resets plotting so that it doesn't assume (2,1)


# to avoid having to change code in a lot of spots, create a placeholder myK
myK=5
#figure out best repeition of K=
CE = cross.entropy(CentAdmix, K=myK)
#put CE into the console
best = which.min(CE)
# put best into the console

myKQ = Q(CentAdmix, K=myK, run=best)
#type View(myKQ) into console to see table 

#we have extracted ancestry score
#now stitch together with 

myKQmeta = cbind(myKQ, meta2)
#cbind pastes columns together, rows must be in same order

my.colors = c("blue2", "gold", "tomato", "lightblue", "pink")

myKQmeta = as_tibble(myKQmeta) %>%
  group_by(continent) %>%
  arrange(region, pop, .by_group = TRUE)

# tibble tidyverse table
# push
  #instead of cluster alphabetically cluster by groups
# group samples by seperatly by contininet
# then by region
# in region by pop


pdf("figures/Admixture_K=5.pdf", width = 10, height=5)
barplot(as.matrix(t(myKQmeta[ , 1:myK])),
        border = NA,
        space=0,
        col = my.colors[1:myK],
        xlab="Geographic Regions",
        ylab="Ancestry Proportions",
        main=paste0("Ancestry Matrix K=", myK)) #pastes 
        
axis(1,
     at=1:length(myKQmeta$region),
     labels=myKQmeta$region,
     tick = F,
     cex.axis=0.5,
     las=3)

dev.off()
#interpretations of plot
##pretty noisy in there
## lots of individuals that have portion of their genome from different clusters
## that makes sense for NE sample = admixed
##PNW stands out as their own cluster, Dark blue and Pink
### that is similar to PCA plot since PNW was clustered 
#### split them across PC1
## CEU cluster defined by yellow chunk
## high admixture in middle also similar to PCA plot


# PCA can give more insight to the Admixture data
# probably 2 different introductions of the PNW species
## similar to PCA again since there is a big cluster of PNW but a couple spread across PC1

