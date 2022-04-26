## Code for quality control and normalization of MethylEPIC data

Analyses discribed here were performed using R 4.0.3 Data processing was initially performed in NBB (dicovery) first, and an identical pipeline repeated later for the BDR (replication) data.

### Processing of NBB data

```
# Install and load required packages

library(minfi)
library(wateRmelon)
library(IlluminaHumanMethylation450kmanifest)
data(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
BiocManager::install("FlowSorted.DLPFC.450k")
library(FlowSorted.DLPFC.450k)
BiocManager::install("CpGFilter")
library(CpGFilter)

# Read sample data from Sample Sheet (minfi)
targets <- read.metharray.sheet("/EPIC_methylation_analysis", pattern = "NBBall_Samplesheet_corrected_merged_Clinicopath_Clean.csv")

# Read methylation data based on directions in sample sheet
rgSet <- read.metharray.exp(targets=targets, extended = TRUE)

# Run primary wateRmelon QC 
qcReport(rgSet, sampNames=targets$Sample_ID, sampGroups=targets$Neuropath_diagnosis, pdf = "qcReport_NBB.pdf")

```

The cell count estimation tool implemented in minfi requires input data in RGChannelSet format, and must be therefore be performed on raw data, upstream of further QC and filtering.

```
# Use minfy function to estimate cell counts based on DLPFC cell sorted reference data:
CellCounts <- estimateCellCounts(rgSet, compositeCellType = "DLPFC", cellTypes = c("NeuN_neg","NeuN_pos"))

targets$NeuN_neg <- CellCounts[,1]
targets$NeuN_pos <- CellCounts[,2]
targets$NeurProp <- with(targets, NeuN_pos/(NeuN_pos + NeuN_neg))
pairwise.t.test(targets$NeurProp, factor(targets$Neuropath_diagnosis))



```
