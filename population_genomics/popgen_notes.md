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
