## Code for quality control and normalization of MethylEPIC data

Analyses discribed here were performed using R 4.0.3

Raw data include a sample sheet and a total of 192 pairs of .idat files from Illumina MethylationEPIC arrays. These pairs include 4 technical replicates and 1 negative control. 

```
# Install and load required packages

library(minfi)
library(wateRmelon)
library(IlluminaHumanMethylation450kmanifest)
data(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
BiocManager::install("FlowSorted.Blood.450k")
library(FlowSorted.Blood.450k)
BiocManager::install("CpGFilter")
library(CpGFilter)
library(sva)
library(limma)
library(missMethyl)
BiocManager::install("bacon")
library(bacon)
BiocManager::install("qqman")
library(qqman)
BiocManager::install("DMRcate")
library(DMRcate)

# Use minfi function to read sample data from sample sheet
targets <- read.metharray.sheet("~/Methylation/SortedLeukocytes", pattern="methylation_analysis_all_celltypes_2021.csv")

# Read methylation data based on directions in sample sheet
rgSet <- read.metharray.exp(targets=targets, extended = TRUE)

# R object of class RGChannelSetExtended containts 1051815 data points on 192 samples

# Run primary wateRmelon QC 
qcReport(rgSet, sampNames=targets$Sample_Name, sampGroups=targets$Celltype, pdf = "~/Methylation/SortedLeukocytes/qcReport_SortedLeukocytes.pdf")

´´´
Different QC and normalization tools use different R objects as input, with implications for the order of operations. It is appropriate to perform sample filtering before normalization, yet sex checks require mapping to genome, which in the main pipeline is downstream of normalization. An initial non-normalized mapping to the genome is therefore performed for the sake of sex checks only. 

´´´
# Correct annotation
mset.raw <- preprocessRaw(rgSet)
annotation(mset.raw)[2] <- "ilm10b4.hg19"

gmset.raw <- mapToGenome(rgSet)

# Estimate sex based on methylation data
methylSex <- getSex(gmset.raw)

# Plot and compare with database sex
databaseSex <- as.factor(targets$Gender)
databaseSex <- substring(databaseSex, 1, 1)
pdf("~/Methylation/SortedLeukocytes/qcPlot_SexCheck_SortedLeukocytes.pdf")
plot(methylSex$xMed, methylSex$yMed, pch = 19, col = factor(databaseSex))
dev.off()
targets$PassSexCheck <- methylSex$predictedSex == databaseSex

# Rename samples to match methylation data and identify samples failing sex check
targets$PlatePos_ID <- paste(targets$Slide,targets$Array, sep="_")
FailSexCheck <- subset(targets, PassSexCheck == FALSE, select = c(PlatePos_ID))

# Filter out low quality probes and samples:
pfltSet <- pfilter(rgSet)

# 1 samples having 1 % of sites with a detection p-value greater than 0.05 was removed. This was the negative control sample.
# 3089 sites were removed as beadcount <3 in 5 % of samples 
# 5404 sites having 1 % of samples with a detection p-value greater than 0.05 # were removed 

# Remove negative control from target data
targets <- targets[complete.cases(targets$Sample_Group),]
rownames(targets) <- targets$PlatePos_ID

# Identify outliers based on the wateRmelon outlyx function:
pdf("~/Methylation/SortedLeukocytes/qcPlot_Outlyx_SortedLeukocytes.pdf")
outliers <- outlyx(pfltSet, plot = TRUE)
dev.off()

# No outliers identified



´´´
