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

# Read demographic and phenotype data from Sample Sheet (minfi)
targets <- read.metharray.sheet("/EPIC_methylation_analysis", pattern = "NBBall_Samplesheet_corrected_merged_Clinicopath_Clean.csv")

# Read methylation data based on directions in sample sheet
rgSet <- read.metharray.exp(targets=targets, extended = TRUE)

rgSet
# Class: RGChannelSetExtended 
# dim: 1051815 504

# Run primary wateRmelon QC 
qcReport(rgSet, sampNames=targets$Sample_ID, sampGroups=targets$Neuropath_diagnosis, pdf = "qcReport_NBB.pdf")

```

The cell type composition estimation tool implemented in minfi requires input data in RGChannelSet format, and must be therefore be performed on raw data, upstream of further QC and filtering.

```
# Use algorithm implemented in minfi to estimate cell type composition based on DLPFC cell sorted reference data:
CellCounts <- estimateCellCounts(rgSet, compositeCellType = "DLPFC", cellTypes = c("NeuN_neg","NeuN_pos"))

# Add estimates of NeuN positive and NeuN negative cells to phenotype data and calculate NeuN+ proportion:
targets$NeuN_neg <- CellCounts[,1]
targets$NeuN_pos <- CellCounts[,2]
targets$NeurProp <- with(targets, NeuN_pos/(NeuN_pos + NeuN_neg))

# Assess association between estimated cell type composition and neuropathological diagnoses:
pairwise.t.test(targets$NeurProp, factor(targets$Neuropath_diagnosis))

# Visualize cell composition across groups as boxplot:
boxplot(targets$NeurProp ~ targets$Neuropath_diagnosis, col = "lightblue")

# Boxplot reveals two outliers with low neuronal proportion. Get outlier IDs:
Cell_composition_outliers <- subset(targets, NeurProp <= 0.125, select = c(PlatePos_ID))

```

Different QC and normalization tools use different R objects as input, with implications for the order of operations. It is appropriate to perform sample filtering before normalization, yet sex checks require mapping to genome, which in the main pipeline is downstream of normalization. An initial non-normalized mapping to the genome is therefore performed for the sake of sex checks only. 

```
# Create a MethylSet without normalization (minfi):
mset.raw <- preprocessRaw(rgSet)

# Some trial and error at this point reveals that the MethylSet annotation needs to be updated. Annotate the MethylSet:
annotation(mset.raw)[2] <- "ilm10b4.hg19"

# Map to genome (minfi):
gmset.raw <- mapToGenome(rgSet)

# Estimate sex based on methylation data (minfi)
methylSex <- getSex(gmset.raw)

# Plot and compare with database sex:
databaseSex <- as.factor(targets$Sex)
plot(methylSex$xMed, methylSex$yMed, pch = 19, col = factor(databaseSex))

# Add column with sex check outcome to phenotype data:
targets$PassSexCheck <- methylSex$predictedSex == databaseSex

# Get IDs of 3 sex-check failing samples
targets$PassSexCheck <- methylSex$predictedSex == databaseSex
FailSexCheck <- subset(targets, PassSexCheck == FALSE, select = c(PlatePos_ID))

```

A number of QC tools operate on RGChannelSet, and MethylSet objects.

```
# Filter out low quality probes and samples (wateRmelon):
pfltSet <- pfilter(rgSet)

# No samples and 7930 sites removed based on detection p-value. 9568 sites were removed as beadcount <3 in 5 % of samples.

# Create MethylSet:
mset.pflt <- preprocessRaw(pfltSet)

mset.pflt
# class: MethylSet 
# dim: 849167 504 

# Identify outliers based on the wateRmelon outlyx function. No outliers identified:
outliers <- outlyx(mset.pflt, plot = TRUE)

# Identity low quality samples based on minfi QC:
minfiQC <- getQC(mset.pflt)
plotQC(minfiQC)

# Identify samples failing minfi QC. Plot shows 4 samples failing, falling below the following cuoff:
FailMinfyQC <- subset(minfiQC, mMed < 11 & uMed < 10.05)

# Concatenate lists of all samples failing and passing QC:
Samples_failing_QC <- as.vector(c(Cell_composition_outliers$PlatePos_ID, FailSexCheck$PlatePos_ID, rownames(FailMinfyQC)))
Samples_passing_QC <- targets[!targets$PlatePos_ID %in% Samples_failing_QC,]$PlatePos_ID

# Of samples included in Lewy pathology EWAS, 2 failed sex-check, 2 failed minfyQC and 3 were cell composition outliers.
# Filter out these 10 samples:
mset.pflt.sampleflt <- mset.pflt[,Samples_passing_QC]
targets.flt <- targets[!targets$PlatePos_ID %in% Samples_failing_QC,]

mset.pflt.sampleflt
# class: MethylSet 
# dim: 849167 494
```

A number of different normalization algorithms were tested. Performance was evaluated primarily by assessing the pairwise concordance of beta values for technical duplicates and by inspection of density plots. In order to generate these metrics, a beta matrix had to be generated based on each of the normalization methods, as well as unnormalized betas.

```
# Get raw betas
rawBetas <- getBeta(mset.pflt.sampleflt)

# Make data frame for storing values of mean pairwise difference in probe betas within technical duplicate pairs, of which there are 41:
DupMeanDeltaBeta <- data.frame(raw = 1:41)

# Pheno data includes a column "No_in_DupPair", which identifies duplicate pairs as 1a/1b, 2a/2b etc. 
# Create a for loop that will calculate the mean delta beta across all probes for each duplicate pair, and then add this value to the data frame:  
for (i in c(1:13,15:36,38:41)){deltabeta <- print(mean(na.omit(pmax(rawBetas[,subset(targets, No_in_DupPair == paste(i,"a", sep = ""))$PlatePos_ID]-rawBetas[,subset(targets, No_in_DupPair == paste(i,"b", sep = ""))$PlatePos_ID],-(rawBetas[,subset(targets, No_in_DupPair == paste(i,"a", sep = ""))$PlatePos_ID]-rawBetas[,subset(targets, No_in_DupPair == paste(i,"b", sep = ""))$PlatePos_ID]))))); DupMeanDeltaBeta$raw[[i]] <- deltabeta}
# Note that duplicate 14 and 37 each have one sample filtered out in QC. Fill in missing values for sample 14 and 37
DupMeanDeltaBeta$raw[[14]] <- NA
DupMeanDeltaBeta$raw[[37]] <- NA

# Normalize using the wateRmelon dasen method and get betas:
mset.pflt.sampleflt.dasen <- dasen(mset.pflt.sampleflt)
dasenBetas <- getBeta(mset.pflt.sampleflt.dasen)

# Add dasen column and mean delta beta for each duplicate pair to data frame:
DupMeanDeltaBeta$dasen <- NA
for (i in c(1:13,15:36,38:41)){deltabeta <- print(mean(na.omit(pmax(dasenBetas[,subset(targets, No_in_DupPair == paste(i,"a", sep = ""))$PlatePos_ID]-dasenBetas[,subset(targets, No_in_DupPair == paste(i,"b", sep = ""))$PlatePos_ID],-(dasenBetas[,subset(targets, No_in_DupPair == paste(i,"a", sep = ""))$PlatePos_ID]-dasenBetas[,subset(targets, No_in_DupPair == paste(i,"b", sep = ""))$PlatePos_ID]))))); DupMeanDeltaBeta$dasen[[i]] <- deltabeta}

# Normalize using the minfi preprosessQuantile method and get betas. This normalization method generates a GenomicRatioSet as output:
GRset.pflt.sampleflt.quantile <- preprocessQuantile(mset.pflt.sampleflt, sex = targets$Sex)
quantileBetas <- getBeta(GRset.pflt.sampleflt.quantile)

# Add quantile column and mean delta beta for each duplicate pair to data frame:
DupMeanDeltaBeta$quantile <- NA
for (i in c(1:13,15:36,38:41)){deltabeta <- print(mean(na.omit(pmax(quantileBetas[,subset(targets, No_in_DupPair == paste(i,"a", sep = ""))$PlatePos_ID]-quantileBetas[,subset(targets, No_in_DupPair == paste(i,"b", sep = ""))$PlatePos_ID],-(quantileBetas[,subset(targets, No_in_DupPair == paste(i,"a", sep = ""))$PlatePos_ID]-quantileBetas[,subset(targets, No_in_DupPair == paste(i,"b", sep = ""))$PlatePos_ID]))))); DupMeanDeltaBeta$quantile[[i]] <- deltabeta}

# Normalize using the watRmelon Beta-Mixture Quantile Normalisation method. This output comes as a beta matrix:
mset.pflt.sampleflt.bmiq <- BMIQ(mset.pflt.sampleflt, nfit = 5000)

# Add BMIQ column and mean delta beta for each duplicate pair to data frame:
DupMeanDeltaBeta$bmiq <- NA
for (i in c(1:13,15:36,38:41)){deltabeta <- print(mean(na.omit(pmax(mset.pflt.sampleflt.bmiq[,subset(targets, No_in_DupPair == paste(i,"a", sep = ""))$PlatePos_ID]-mset.pflt.sampleflt.bmiq[,subset(targets, No_in_DupPair == paste(i,"b", sep = ""))$PlatePos_ID],-(mset.pflt.sampleflt.bmiq[,subset(targets, No_in_DupPair == paste(i,"a", sep = ""))$PlatePos_ID]-mset.pflt.sampleflt.bmiq[,subset(targets, No_in_DupPair == paste(i,"b", sep = ""))$PlatePos_ID]))))); DupMeanDeltaBeta$bmiq[[i]] <- deltabeta}

# Summarize duplicate delta beta metrics for each of the normalization methods. The dasen method has the lowest mean value at 0.02614:
summary(DupMeanDeltaBeta)

# Visualize effect of normalization on beta density plots for all samples:
densityPlot(rawBetas, main = "Raw")
densityPlot(dasenBetas, main = "dasen")
densityPlot(quantileBetas, main = "Quantile")
densityPlot(mset.pflt.sampleflt.bmiq, main = "dasen")

# Visualize effect of normalization on beta density by probe type in a random single sample. 
# The GRset and beta matrix lack probe type info, which is extracted from the MethylSet stage and provided for the plot function:
plotBetasByType(mset.pflt.sampleflt[,1], main = "Raw")
plotBetasByType(mset.pflt.sampleflt.dasen[,1], main = "Dasen")
probetypes <- data.frame(Name = 1:848963)
probetypes$Name <- rownames(quantileBetas)
probetypes$Type <- getProbeType(mset.pflt.sampleflt)
plotBetasByType(quantileBetas[,1], probeTypes = probetypes, main = "Quantile")
plotBetasByType(mset.pflt.sampleflt.bmiq[,1], probeTypes = probetypes, main = "BMIQ")

```
