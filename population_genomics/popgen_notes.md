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

Today we filtered the VCF file so that there was more concise data that you could look at this is saved under 01_VCF_Filtering(Correct Location)

We set the depth

We created heatmaps

### 09-19-2024

lots of notes on 02_Diversity_Differentiation

Today we focused on creating Manhattan plots for our filtered data

we pulled our filtered data from last class from our outputs folder and assigned it the name vcf so that it is the updated version.

we looked at the meta data of the file using the head function

### 9-24-24 Hs averages and stuff

check code for notes, did lots of stuff with Count of SNPs vs Gene diveristy within Regions

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
