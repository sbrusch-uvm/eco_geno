setwd("~/Projects/eco_geno/population_genomics/docs/")

library(vcfR)
setwd("/gpfs1/cl/pbio3990/PopulationGenomics/")
list.files()

vcf <- read.vcfR("variants/Centaurea_filtered.vcf.gz")
head(vcf)

X11.options(type="cairo")


dna <- ape::read.dna("reference/GCA_030169165.1_ASM3016916v1_genomic.fa.gz", format = "fasta")
gff <- read.table("reference/GCA_030169165.1_ASM3016916v1_genomic_genes.gff.gz", sep = "\t", quote = "")


chr1 <- create.chromR(name="Chromosome 1", vcf=vcf, seq=dna, ann=gff)

plot(chr1)

pdf(file="~/Projects/eco_geno/population_genomics/figures/Chromoplot.pdf")
chromoqc(chr1, xlim=c(1e1, 1.1e8))
dev.off()


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
                                      cutoff = 0.8)


vcf.filt.indMiss <- filter_biallelic(vcf.filt.indMiss)
vcf.filt.indMiss <- min_mac(vcf.filt.indMiss, min.mac = 1)


vcf.filt.indSNPMiss <- missing_by_snp(vcf.filt.indMiss, cutoff = 0.5)

DP2 <- extract.gt(vcf.filt.indSNPMiss, element="DP", as.numeric=T)

heatmap.bp(DP2[1:5000,], 
           rlabels=F, clabels=F)


write.vcf(vcf.filt.indSNPMiss,
          "~/Projects/eco_geno/population_genomics/outputs/vcf_final.filtered.vcf.gz")

