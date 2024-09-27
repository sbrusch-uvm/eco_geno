# Coding and Data notes for the population genomics module

## Author: Sophia Brusch

### 09-10-2024 - Intro to Centaurea GBS data and working with VCF files

We'll be analyzing the GBS from data from 3 regions (Eu, NE, PNW) starting today with variant call formate files (VCF's)

Pushed this new notes file

went to the vacc website to look at the cluster shell access

ls -l to look at list, .gz is a Gzip file --\> zcat allows you to look into the file.gz

zcat WV_9.R1.fq.gz \|head gives you the first 10 lines of data from that data

zcat WV_9.R1.fq.gz \|head -n 4 gives you 4 lines of data

first line is the info for the unique individual CCACA is the barcode

second line is the sequence itself \~150 bp

break line

4th line are other characters = the I's are the quality score as a multi value digit as 1 unit

I is Q-score of 40 which is awesome!!!!!!!!!!! There was a table to reference where the different symbols represented different Q-scores

Some commands we used:

cd change directory

ls -l long list (can be abbreviated as "ll")

zcat view a zipped file

\|head to view first 10 lines of a file

\| "pipe" send output of one fxn to another

### 09-17-2024 - VCF File Stuff

So I was absent last class so this was a lot of catching up to do

Today we filtered the VCF file so that there was more concise data that you could look at this is saved under 01_VCF_Filtering(Correct Location) in my docs folder

the type=cairo is helpful for creating plots, something with the software being weird, should run before trying to make any plots

setwd sets the working directory which makes working with files easier

vcfR allows for manipulation of variant call format (VCR) data which is the centaurea data

the ape:: is a package used for phylogenies and reading DNA info

we made plots of data pulled from the vcf file to make chromosome 1 plots from annotated data

We created a numeric value for depth DP

we made a heatmap to visualize the matrix of depth and missingness of our data the VCF file

We used the SNPfltR to visualize and filter, also lets us save out the info

what is the average depth of SNP

mean value is about 30 reads, doesn't have a huge tail

rule of thumb, twice the average set at N/A, max depth filter 60

we also explored missingness

made a meta variable

we reassigned the column names to something that works better for us

looked at missingness on the individual level

the different rows is different options, cutoff 0 means it doesn't filter anyone out, this creates violin plots to show distribution top missing data, bottom average depth, individuals are grouped by regions, looking at the plots filter \~0.75 is missing data

we got rid of biallelic variants

made another heatmap

saved the datafile into our repo under outputs vcf_final_filtered

### 09-19-2024

lots of notes on 02_Diversity_Differentiation script

used tidyverse and qqman

we pulled our filtered data from last class from our outputs folder and assigned it the name vcf so that it is the updated version.

we looked at the meta data of the file using the head function

we calculate diversity stats using the genetic_diff fxn in vcfR

Today we focused on creating Manhattan plots for our filtered data

str function gives info from the vcf file

\$ allows to extract subdata

made the plots than filtered out Gst values greater than 0 since those are mistakes from stuff

### 9-24-24 Hs averages and stuff

check code for notes, did lots of stuff with Count of SNPs vs Gene diveristy within Regions

saved the manhattan plot to outputs so that we don't have to run it everytime to get to that point

we will make one big column out of the vcf.div file, tidy operations %\>% means pipe results from one step into the next

make a ggplot filter by all the columns and then, Hs=2pq, max value of Hs is 0.5, most of the allele frequencies are at the extremes, 1 allele is fixed while the other is at a very low freq. , true for almost any species , huge spike at 0 since we saw 0s in out table vcf.div, the pop isnt polymorphic , the NE group has higher than everyone else freq. , ggsave saves the last plot you made

filter(value!=0 & value\<0.25) %\>% #take out all the SNPs with 0 values, != means does not equal

the filter(value!=0) allows you to look at values on the fly, the & value\<0.25 looks at

we need to thin in the SNPs for LD (linkage disequilibrium) before

we run PCA and Admixture analyses to satisfy the assumptions of independence among loci

new vcf file based on giving it the og vcf file and thin it 500 bp

we learned most of our SNPs are stacked since 3646/15454 weren't located within 500bp of another SNP

now we need to uncompress it #hide the uncompressed vcf file too big for github outside of repo

### 09-26-2024 - Intro to Admixture and Cross Validation

Okay so, we didn't do too much crazy coding today

made a ggplot of the PCA plot of PC1 vs PC2, we saved the plot to the figures file

We learned how admixture works

-   ex. 2 subpopulations in 1 sample the user makes a k \# of groups

-   assign individuals to 1 of the k groups

-   calc the allele freq in each group

-   calc 2pq and compare that to the observed freq of Het

-   if that previous step wasn't very accurate then repeat the process in different groups

Cross validation splits data into a training set and a testing set

this wont work for self-fertilizing plants

screeplot, = magnitude of PC values decending

showing the eigenvalue = the axis put through the dots of, PC1 is the top dot

PCA dots get smaller with each consecuative PCA plot

most PCs won't be interesting for us, only focus on the early ones

could be that the allele freq evolved very divergent than the other things, makes very few assumptions

first line of ggplot says what data set to use

aes assigns the PC lines to the axis with color coordination of region and shapes for continent

saved the last plot we made to the figures file

looking at PC3 vs PC3 we can say that the Centaurea individuals are dominantly serperated along PC2

CentPCA\$eigenvectors[1:5, 1:2] in console gives the head of V1 and V2 the eigenvectors each SNP locus has a quantitative value

value -0.08, going from 0-2 gives a vector of -0.08 and 0.046 so its vector towards neg x axis and pos y axis

nice way to look at genetic structure without priority
