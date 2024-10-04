library(vcfR)
#library(ggplot2)

# set wd to PopulationGenomics folder on gpfs1 class drive

setwd("/gpfs1/cl/pbio3990/PopulationGenomics/")

list.files("variants/")
list.files("reference/")

vcf <- read.vcfR("variants/Centaurea_filtered.vcf.gz")

dna <- ape::read.dna("reference/GCA_030169165.1_ASM3016916v1_genomic.fa.gz", format="fasta")
#dna2 = dna$CM058040.1

gff <- read.table("reference/GCA_030169165.1_ASM3016916v1_genomic_genes.gff.gz", sep="\t", quote="")
#gff2 = grep("CM058040.1",gff)

vcf

head(vcf)

# create the chromR object class
chr1 <- create.chromR(name="Centaurea chr1", vcf=vcf, seq=dna, ann=gff, verbose=F)
#doesn't work: chr2 <- create.chromR(name="Centaurea chr2", vcf=vcf, seq=dna2, ann=gff2, verbose=F)

#head(chr1@var.info)
#chr1
#quantile(chr1@var.info$MQ)
#quantile(chr1@var.info$DP)

# example plot
plot(chr1)
chromoqc(chr1, xlim=c(1e1, 1.1e8))

# example with exporting to pdf -- note path requirements
# also note zoom ability via xlim=c() option
pdf(file="~/Projects/eco_geno/population_genomics/figures/chromoplot111111")
chromoqc(chr1, xlim=c(1e1, 1.1e8))
dev.off()


# Mask poor quality variants
# chr1_masked <- masker(chr1,
#                       min_QUAL=50,
#                       min_DP=1000,
#                       max_DP=10000,
#                       min_MQ=30)
# 
# plot(chr1_masked)
# chromoqc(chr1_masked, xlim=c(1e1, 1.1e8))
# 
# # Now process the chromR object with proc.chromR
# # The default window size = 1000bp, can set with win.size
# chr1_proc <- proc.chromR(chr1_masked, win.size = 1e5)
# plot(chr1_proc)
# chromoqc(chr1_proc, xlim=c(1e1, 1.1e8))

#head(chr1_masked) # see an overview of each list element

####### Left off here last time #######

#How to filter a vcf file for minDP
DP <- extract.gt(vcf, 
                 element="DP", 
                 as.numeric=T, 
                 convertNA=T)
quantile(DP)

DP[DP==0] <- NA #set the 0 depth genos to `NA`

quantile(DP, na.rm=T) # much better

dim(DP) # ensure loci are in rows; samples in columns

#Let's look at mean DP per individual...
# avgDP = colMeans(DP, na.rm=T)
# summary(avgDP)
# hist(avgDP, breaks=50) 
# mean is ~24X, range from 2.8-171X.  
# Pretty good! Ideeally want avg of 15-20X/ind

#What about missingness? We can use the heatmap function:

#pdf(file="~/courses/Ecological_Genomics_24/population_genomics/figures/chromoPlot_chr1.pdf")
heatmap.bp(DP[1:5000,], rlabels=F, clabels=F)
#dev.off()

# set individual genotypes with DP<X to `NA`
# not needed with SNPfiltR??
#vcf@gt[,-1][is.na(DP)==TRUE] <- NA 

vcf # check to see % missing data -- 25.3% if DP==0 <- NA

# Now that we see the data attribtues, let's start filtering 
library(SNPfiltR)

meta <- read.csv("metadata/meta4vcf.csv", header=T)
head(meta)
meta2=meta[,c(1,4)]
names(meta2) = c("id","pop")

# Look at mean depth per ind
hard_filter(vcf)
vcf.filt <- hard_filter(vcf, 
                        depth=3) #What's reasonable while keeping variants?

# Look at allele balance (note these are autotetratploid...0.25, 0.5, 0.75)
vcf.filt <- filter_allele_balance(vcf.filt,
                                  min.ratio = 0.15,
                                  max.ratio=0.85)

vcf.filt <- max_depth(vcf.filt, 
                      maxdepth=60) # generally set filter to 2X mean depth 

# start with no cutoff for exploratory, 
# then add cuttoff=0.8 or 0.75 (N=36 sammples)

vcf.filt.indMiss <- missing_by_sample(vcf.filt, 
                                      popmap = meta,
                                      cutoff=0.75) 
#subset popmap to only include retained individuals
meta <- meta[meta$id %in% colnames(vcf.filt.indMiss@gt),]

# gets rid of monomrophic or multi-allelic sites
vcf.filt.indMiss <- filter_biallelic(vcf.filt.indMiss) 


# use PCA to verify that missingngess is not driving clustering
# maybe try just 1 or 2 threholds at first...note that lower is more permissive
# didn't work for me on the VACC-OOD!
#library(adgenet)
#missPCA <- assess_missing_data_pca(vcfR=vcf.filt.indMiss, 
#popmap = meta, 
#thresholds = 0.5, 
#clustering = FALSE)

#Filter out by SNP missingness -- higher cuttoff is more stringent
vcf.filt.indSNPMiss <- missing_by_snp(vcf.filt.indMiss, 
                                      cutoff=0.5)

vcf.filt.indSNPMiss <- min_mac(vcf.filt.indSNPMiss,
                               min.mac=2)

#assess clustering without MAC cutoff
#miss<-assess_missing_data_tsne(vcf.filt.indSNPMiss, 
# popmap=meta, 
#clustering = FALSE)

DP2 <- extract.gt(vcf.filt.indSNPMiss, 
                  element="DP", 
                  as.numeric=T, 
                  convertNA=T)


heatmap.bp(DP2[1:5000,], rlabels=F, clabels=F)

vcfR::write.vcf(vcf.filt.indSNPMiss, 
                "~/myrepo/outputs/Centaurea_finalfiltered.vcf.gz")

# Can also thin for LD:

vcf.filt.indSNPMiss.thin <- distance_thin(vcf.filt.indSNPMiss,
                                          min.distance=500)

vcfR::write.vcf(vcf.filt.indSNPMiss.thin, 
                "~/myrepo/outputs/Centaurea_finalfiltered_thinned.vcf.gz")

mydiff <- vcfR::genetic_diff(vcf.filt.indSNPMiss,
                             pops=as.factor(meta$region),
                             method="nei")


setwd("~/Projects/eco_geno/population_genomics/docs/")