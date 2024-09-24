#estimating diversity and genetic differentiation in the filtered centaurea data

library (vcfR)
library(tidyverse)
library(qqman)

#helps solve plotting issues
X11.options(type="cairo")


# read in our VCF file from our repo outputs directory

vcf <- read.vcfR("~/Projects/eco_geno/population_genomics/outputs/vcf_final.filtered.vcf.gz")

# read in our metadata -- info on population of origin

meta <- read.csv("/gpfs1/cl/pbio3990/PopulationGenomics/metadata/meta4vcf.csv")

head(meta) # vcf files has 595 samples
dim(meta) #meta has 629 inds

meta2 <- meta[meta$id %in% colnames(vcf@gt[,-1]),]
#colnames(vcf@gt[,-1]),])[1:5] will give you just the first 5 colum names 

#calculate diversity stats using the genetic_diff fxn in vcfR
vcf.div <- genetic_diff(vcf, pops = as.factor(meta2$region), method = "nei")


str(vcf.div)
chr.main <- unique(vcf.div$CHROM)[1:8]
chrnum <- as.data.frame(cbind(chr.main, seq(1, 8, 1)))

vcf.div.MHplot <- left_join(chrnum, vcf.div, join_by(chr.main==CHROM))

vcf.div.MHplot <- vcf.div.MHplot %>% filter(Gst>0) %>% mutate(SNP=paste0(chr.main, "_", POS))
#filtered out Gst values greater than 0

vcf.div.MHplot$V2 = as.numeric(vcf.div.MHplot$V2)
vcf.div.MHplot$POS = as.numeric(vcf.div.MHplot$POS)


manhattan(vcf.div.MHplot, 
          chr = "V2", 
          bp = "POS", 
          p="Gst", 
          col = c("blue4", "orange3"), 
          logp = F, 
          ylab = "Fst among regions",
          suggestiveline = quantile(vcf.div.MHplot$Gst, 0.999))
#the manhattan plot dots are SNPs, the Fst is nonrandom mating based on population structure
# large amount of the genome has low Fst 
#the line is the 99.9 percentile (0.999)
# diverence is probably less than 1%

# POS that are multiples of 3 apart are synonymous mutations most likely




## 9/24/24 ##

write.csv(vcf.div.MHplot, "~/Projects/eco_geno/population_genomics/outputs/Genetic_Dif_byRegion.csv",
          quote = F,
          row.names = F)
#this saves the plot to outputs so that you don't have to run it every time

names(vcf.div.MHplot)

options(bitmapType = "cairo")
#now we will make one big column out of the vcf.div file 
#tidy operations %>% means pipe results from one step into the next
vcf.div.MHplot %>% 
  as_tibble()%>%
  pivot_longer(c(4:9)) %>%
  ggplot(aes(x=value, fill=name)) + 
  geom_histogram(position = "identity", alpha=0.5, bins=50) +
  labs(title = "Genome-wide expected heterozygosity (Hs)", fill="Regions", 
       x="Gene diversity within Regions", y="Counts of SNPs")
ggsave("Histogram_Genome_Diverity_byRegion.pdf", 
       path="~/Projects/eco_geno/population_genomics/figures/")

 #make a ggplot filter by all the columns and then 
# Hs=2pq, max value of Hs is 0.5
# most of the allele frequencies are at the extremes, 1 allele is fixed while the other is at a very low freq.
# true for almost any species 
#huge spike at 0 since we saw 0s in out table vcf.div, the pop isnt polymorphic
# the NE group has higher than everyone else freq.
#ggsave saves the last plot you made


vcf.div.MHplot %>% 
  as_tibble()%>%
  pivot_longer(c(4:9)) %>%
  group_by(name) %>% 
  filter(value!=0 & value<0.25) %>% #take out all the SNPs with 0 values, != means does not equal
  summarise(avgHs=mean(value), StdDev_Hs=sd(value), N_Hs=n()) 

# the filter(value!=0) allows you to look at values on the fly
# the & value<0.25 looks at 


