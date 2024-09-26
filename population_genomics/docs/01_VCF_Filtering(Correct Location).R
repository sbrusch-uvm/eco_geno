install.packages("ape")
install.packages("vcfR")
library(vcfR)
setwd("/gpfs1/cl/pbio3990/PopulationGenomics/")
list.files()

vcf <- read.vcfR("variants/Centaurea_filtered.vcf.gz")
head(vcf)

X11.options(type="cairo")
#fixed section contains info of each snit
# genotype section has matrix of rows and columns

# associate reference genome with filtered vcf
dna <- ape::read.dna("reference/GCA_030169165.1_ASM3016916v1_genomic.fa.gz", format = "fasta")
gff <- read.table("reference/GCA_030169165.1_ASM3016916v1_genomic_genes.gff.gz", sep = "\t", quote = "")

#chromosome 1 as example

chr1 <- create.chromR(name="Chromosome 1", vcf=vcf, seq=dna, ann=gff)
#telling it that its vcf use the sequenced dna data and annotated

plot(chr1)

pdf(file="~/Projects/eco_geno/population_genomics/figures/Chromoplot.pdf")
chromoqc(chr1, xlim=c(1e1, 1.1e8))
dev.off()

# x limit for base pair to base pair


##### left off here last time#####

DP <- extract.gt(vcf, element="DP", as.numeric=T)
DP[1:5, 1:10]

quantile(DP)

DP[DP==0] <- NA

quantile(DP, na.rm=T)
#visual the matrix of depth and missingness in our VCF file:

heatmap.bp(DP[1:1000,], rlabels = F, clabels = F)
# white is missing data, varying across individuals
#the bars are the mean depth value
#some loci have some missing data 
#the visual helps us understand and explore the data

library(SNPfiltR)
#lets us visualize the filtering steps and filter, lets us save it out

#you can type ?SNPfiltR into the console and it gives you more info or help


vcf.filt <- hard_filter(vcf, depth=3)
# saying if you dont have at least 3 reads at a SNP position set it to N/A
# this is user customized and defined
# it should be at 35% missing data 
# this takes care of low depth problem
# could explore valyes DP 5, 10, ..... based on your own goals

max_depth(vcf.filt)
#what is the average depth of SNP
#mean value is about 30 reads, doesn't have a huge tail
#rule of thumb, twice the average set at N/A, max depth filter 60

vcf.filt <- max_depth(vcf.filt, maxdepth=60) # filters out genotypes with >60 reads
# didn't lose many SNPs


#the next thing to explore is missingness
meta <- read.csv("metadata/meta4vcf.csv", header = T)
#type View into console and you can get a visual on it 

#suggest we group on region since it will be easier to visualize
#save the 2 columns

meta2 <- meta[,c(1,4)]
#gives all the rows and only columns 1 and 4
#df[rows, columns]
#df[ , c(1,4)]
#type View(meta2) in the console to check it worked 
#type head(meta2) in the console to check it worked

names(meta2) <- c("id", "pop")
#this reassigns the column names since we are calling the regions pops

meta2$id = as.factor(meta2$id)
meta2$pop = as.factor(meta2$pop)
# this makes it factors

vcf.filt.indMiss <- missing_by_sample(vcf.filt, popmap = meta2, cutoff = 0.75)
#individual level of missingness
#the different rows is different options
# cutoff 0 means it doesn't filter anyone out
#this creates violin plots to show distribution top missing data, bottom average depth
#individuals are grouped by regions
#looking at the plots filter ~0.75 os missing data 


vcf.filt.indMiss <- filter_biallelic(vcf.filt.indMiss)
vcf.filt.indMiss <- min_mac(vcf.filt.indMiss, min.mac = 1)
#getting rid of the biallelic SNPs
#min.mac how many times do you have to see it to keep it
#type vcf.filt.indMiss into the console to see how much info there now is after filtering


vcf.filt.indSNPMiss <- missing_by_snp(vcf.filt.indMiss, cutoff = 0.5)
#total SNPs retained, at missing proportion 30%

DP2 <- extract.gt(vcf.filt.indSNPMiss, element="DP", as.numeric=T)

heatmap.bp(DP2[1:5000,], 
           rlabels=F, clabels=F)
#if anything we increased out missing data
#took out a lot of SNPs that had very few reads on them
#patchwork of info is normal RADseq, its cheap 

write.vcf(vcf.filt.indSNPMiss,
          "~/Projects/eco_geno/population_genomics/outputs/vcf_final.filtered.vcf.gz")
#getwd in console, tab to fill in after ~/
# the last part is the name you save it as and the file type 
#write out the most final filtered version of data
