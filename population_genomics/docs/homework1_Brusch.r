setwd("~/Projects/eco_geno/population_genomics/docs/")
#this allowed me to save this to my repo when first creating this file

library(vcfR)

## anywhere an .8 or .55 is, those are the new thresholds of missingness
## that being said, if the program is running to analyze an 80% cuffoff value, all sections must be specific for 0.8 and vice versa for 0.55

X11.options(type="cairo") # easier for the plots/graphing
setwd("/gpfs1/cl/pbio3990/PopulationGenomics/") #sets the working directory so that I can access the centaurea data
list.files()

vcf <- read.vcfR("variants/Centaurea_filtered.vcf.gz") #accesses the centaurea file within the vairants file
head(vcf)


dna <- ape::read.dna("reference/GCA_030169165.1_ASM3016916v1_genomic.fa.gz", format = "fasta")
gff <- read.table("reference/GCA_030169165.1_ASM3016916v1_genomic_genes.gff.gz", sep = "\t", quote = "")


DP <- extract.gt(vcf, element="DP", as.numeric=T)
DP[1:5, 1:10]
quantile(DP)
DP[DP==0] <- NA
quantile(DP, na.rm=T)

heatmap.bp(DP[1:1000,], rlabels = F, clabels = F) #this heatmap is for all the data

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
#this is where the change was made, cutoff will be either 0.8 or 0.55
#whenever running program, change the labels depending on which filter you are using

vcf.filt.indMiss <- filter_biallelic(vcf.filt.indMiss)
vcf.filt.indMiss <- min_mac(vcf.filt.indMiss, min.mac = 1)
vcf.filt.indSNPMiss <- missing_by_snp(vcf.filt.indMiss, cutoff = 0.5)

DP2 <- extract.gt(vcf.filt.indSNPMiss, element="DP", as.numeric=T)

heatmap.bp(DP2[1:5000,], 
           rlabels=F, clabels=F)

write.vcf(vcf.filt.indSNPMiss,
          "~/Projects/eco_geno/population_genomics/outputs/vcf_final.filtered0.8.vcf.gz")
# for second filtering save, I changed the name to vcf_final_filtered0.55.vcf.gz after rerunning the program

#########################################################################

library(tidyverse)
library(qqman)

options(bitmapType = "cairo")

vcf <- read.vcfR("~/Projects/eco_geno/population_genomics/outputs/vcf_final.filtered0.8.vcf.gz")
#first time with file vcf_final.filtered0.8.vcf.gz then with .filtered0.55.vcf.gz


meta <- read.csv("/gpfs1/cl/pbio3990/PopulationGenomics/metadata/meta4vcf.csv")
head(meta)
dim(meta)
meta2 <- meta[meta$id %in% colnames(vcf@gt[,-1]),]

vcf.div <- genetic_diff(vcf, pops = as.factor(meta2$region), method = "nei")


str(vcf.div)
chr.main <- unique(vcf.div$CHROM)[1:8]
chrnum <- as.data.frame(cbind(chr.main, seq(1, 8, 1)))

vcf.div.MHplot <- left_join(chrnum, vcf.div, join_by(chr.main==CHROM))
vcf.div.MHplot <- vcf.div.MHplot %>% filter(Gst>0) %>% mutate(SNP=paste0(chr.main, "_", POS))
vcf.div.MHplot$V2 = as.numeric(vcf.div.MHplot$V2)
vcf.div.MHplot$POS = as.numeric(vcf.div.MHplot$POS)
names(vcf.div.MHplot)




vcf.div.MHplot %>% 
  as_tibble()%>%
  pivot_longer(c(4:9)) %>%
  ggplot(aes(x=value, fill=name)) + 
  geom_histogram(position = "identity", alpha=0.5, bins=50) +
  labs(title = "Genome-wide expected heterozygosity (Hs)", fill="Regions", 
       x="Gene diversity within Regions", y="Counts of SNPs")
ggsave("Histogram_Genome_Diverity_byRegion0.8.pdf", 
       path="~/Projects/eco_geno/population_genomics/figures/")


vcf.div.MHplot %>% 
  as_tibble()%>%
  pivot_longer(c(4:9)) %>%
  group_by(name) %>% 
  filter(value!=0 & value<0.25) %>% #take out all the SNPs with 0 values, != means does not equal
  summarise(avgHs=mean(value), StdDev_Hs=sd(value), N_Hs=n())

setwd("~/Projects/eco_geno/")

Hs_table_NonZeros <- vcf.div.MHplot %>% 
  as_tibble() %>% 
  pivot_longer(c(4:9)) %>%
  group_by(name) %>%
  filter(value!=0) %>% 
  summarise(avg=mean(value), StdDev_Hs=sd(value), N_Hs=n())

write.csv(Hs_table, "population_genomics/outputs/Hs_table_NonZeros.csv",
          quote=F,
          row.names = F)

Hs_table_Zeros <- vcf.div.MHplot %>% 
  as_tibble() %>% 
  pivot_longer(c(4:9)) %>%
  group_by(name) %>%
  filter(value==0) %>% 
  summarise(avg=mean(value), StdDev_Hs=sd(value), N_Hs=n())

write.csv(Hs_table_Zeros, "population_genomics/outputs/Hs_table_Zeros.csv",
          quote=F,
          row.names = F)





#################################################################################################

library(LEA) 
setwd("~/Projects/eco_geno/population_genomics/") 

vcf <- read.vcfR("outputs/vcf_final.filtered0.8.vcf.gz") 
# first run through is the 0.8 cuttoff and the second runthrough is the 0.55 cutoff

vcf.thin <- distance_thin(vcf, min.distance=500) 

meta <- read.csv("/gpfs1/cl/pbio3990/PopulationGenomics/metadata/meta4vcf.csv")
dim(meta)
meta2 <- meta[meta$id %in% colnames(vcf@gt[, -1]) , ]
dim(meta2)

write.vcf(vcf.thin, "outputs/vcf_final.filtered.thinned.8.vcf.gz")

system("gunzip -c ~/Projects/eco_geno/population_genomics/outputs/vcf_final.filtered.thinned.8.vcf.gz > ~/vcf_final.filtered.thinned.8.vcf")

geno <- vcf2geno(input.file = "/gpfs1/home/s/b/sbrusch/vcf_final.filtered.thinned.8.vcf", 
                 output.file = "outputs/vcf_final.filtered.thinned.8.geno")

CentPCA <- LEA::pca("outputs/vcf_final.filtered.thinned.8.geno", scale=TRUE)

# if you have run the PCA before, you can just open it here without having to do the previous steps
CentPCA <- load.pcaProject("vcf_final.filtered.thinned.8.pcaProject")

show(CentPCA)
plot(CentPCA)


ggplot(as.data.frame(CentPCA$projections),
       aes(x=V1, y=V2, color=meta2$region, shape=meta2$continent)) +
  geom_point(alpha=0.5) +
  labs(title="Centaurea genetic PCA",
       x="PC1",
       y="PC2",
       color="Region",
       shape="Continent") 


ggsave("figures/CentPCA.8_PC1vPC2.pdf", width=6, height=6, units ="in")
