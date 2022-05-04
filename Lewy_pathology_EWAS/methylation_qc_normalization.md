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

The dasen normalization method (WateRmelon) is selected for further analyses. Data are mapped to genome before a number of further filtering steps. 

```
# Map to genome
gmset.pflt.sampleflt.dasen <- mapToGenome(mset.pflt.sampleflt.dasen)

gmset.pflt.sampleflt.dasen
# class: GenomicMethylSet 
# dim: 848963 494 

# Filter out probes with SNPs
gmset.pflt.sampleflt.dasen.dropsnps <- dropLociWithSnps(gmset.pflt.sampleflt.dasen)

gmset.pflt.sampleflt.dasen.dropsnps
# class: GenomicMethylSet 
# dim: 820434 494 

# Filter out probes on sex chromosomes
EPICanno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
keep <- !(featureNames(gmset.pflt.sampleflt.dasen.dropsnps) %in% EPICanno$Name[EPICanno$chr %in% c("chrX", "chrY")])table(keep)
gmset.pflt.sampleflt.dasen.dropsnps.autosomes <- gmset.pflt.sampleflt.dasen.dropsnps[keep,]

gmset.pflt.sampleflt.dasen.dropsnps.autosomes
# class: GenomicMethylSet 
# dim: 801978 494

# Remove cross-reactive probes as reported by Chen et al (https://www.tandfonline.com/doi/full/10.4161/epi.23470):
cross_reactive_probes <- scan("cross_reactive_probes_450k.txt", what = "character")
keep <- !(featureNames(gmset.pflt.sampleflt.dasen.dropsnps.autosomes) %in% cross_reactive_probes)
gmset.pflt.sampleflt.dasen.dropsnps.autosomes.xrflt <- gmset.pflt.sampleflt.dasen.dropsnps.autosomes[keep,]

gmset.pflt.sampleflt.dasen.dropsnps.autosomes.xrflt
# class: GenomicMethylSet 
# dim: 777590 494 

# Get betas: 
betas.pflt.sampleflt.dasen.dropsnps.autosomes.xrflt <- getBeta(gmset.pflt.sampleflt.dasen.dropsnps.autosomes.xrflt)
 
# Filter out probe-wise outliers using the wateRmelon pwod function. Outliers are probably low MAF/SNP heterozygotes
betas.pflt.sampleflt.dasen.dropsnps.autosomes.xrflt.pwoflt <- pwod(betas.pflt.sampleflt.dasen.dropsnps.autosomes.xrflt)
# 242599 probes detected. These are coerced to NA. Note that these are individual data points, not probes across the dataset. 

# Filter out probe-wise outliers using the wateRmelon pwod function. Outliers are probably low MAF/SNP heterozygotes
betas.pflt.sampleflt.dasen.dropsnps.autosomes.xrflt.pwoflt <- pwod(betas.pflt.sampleflt.dasen.dropsnps.autosomes.xrflt)
# 242599 probes detected. 
# These are coerced to NA. Note that these are individual data points, not probes across the dataset.

# Assess technical vs. biological probe variability and plot result:
# First, make sample names where duplicates have identical names, by omitting last characters. 
targets.flt$repdesign <- substr(targets.flt$Sample_Name, start = 1, stop = 5)
# Next, run CpGFilterICC (CpGFilter):
ICCfilter <- CpGFilterICC(betas.pflt.sampleflt.dasen.dropsnps.autosomes.xrflt, targets.flt$repdesign)

# Remove technical duplicates and samples failing genotyping QC:
targets.flt.genoOK <- subset(targets.flt, Fail_NeuroChip_QC == "N")
targets.flt.rmdup <- subset(targets.flt.genoOK, TechRep == "N")
finalBetas <- betas.pflt.sampleflt.dasen.dropsnps.autosomes.xrflt.pwoflt[,targets.flt.rmdup$PlatePos_ID]
dim(finalBetas)
# [1] 777590    442

# Filter out the quartile of probes with the lowest ratio of biological to technical variability:
quantile(ICCfilter)
#        0%       25%       50%       75%      100% 
# 0.0000000 0.1351696 0.3400994 0.5769919 0.9978198 
keep <- names(subset(ICCfilter, ICCfilter > 0.1351696))
finalBetas <- finalBetas[keep,]

dim(finalBetas)
# [1] 583192    442

# The initial dataset included samples from donors with phenotypes not meeting the study inclusion criteria. 
# Extract control and Lewy body spectrum subjects:
aSyn_spctr_pheno <- targets.flt.rmdup[targets.flt.rmdup$Neuropath_GroupDx == "PD" | targets.flt.rmdup$Neuropath_GroupDx == "CONTR" | targets.flt.rmdup$Neuropath_GroupDx == "DLB" | targets.flt.rmdup$Neuropath_GroupDx == "iLBD",]

# Included samples need complete data on relevant covariates:
aSyn_spctr_pheno <- na.omit(subset(aSyn_spctr_pheno, select = c("PlatePos_ID","Sample_Name","Braak_aSyn_stage","Age_death","Sex","PMD_min","NeurProp","Plate","Neuropath_GroupDx")))

# A few samples have an "atypical" Lewy body pattern not classifiable by Braak stage. These are filtered out:
aSyn_spctr_pheno <- aSyn_spctr_pheno[aSyn_spctr_pheno$Braak_aSyn_stage != "Atypical",]

# After removal of "atypical" cases, the Braak stage variable needs to be redifined as numeric:
aSyn_spctr_pheno$Braak_aSyn_stage <- as.numeric(aSyn_spctr_pheno$Braak_aSyn_stage)

# Extract final betas
aSyn_spctr_betas <- finalBetas[,aSyn_spctr_pheno$PlatePos_ID]

# Export NBB betas, pheno and probe data for osca analyses:
write.table(aSyn_spctr_betas, "aSyn_spctr_betas.txt", sep = " ", dec = ".", row.names = TRUE, col.names = TRUE) 
write.table(aSyn_spctr_pheno, "aSyn_spctr_pheno.txt", sep = " ", dec = ".", row.names = FALSE, col.names = TRUE)
EPICannoSubset <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)[match(rownames(finalBetas), getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)$Name), c(1,2,3,22:ncol(getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)))]
write.table(EPICannoSubset[,1:4], "~/Methylation/Brain/NBB_aSyn_spctr_FinalProbes.txt", sep = " ", dec = ".", row.names = TRUE, col.names = FALSE, quote = FALSE)
```

The final output of the QC and normalization pipeline is a matrix of filtered and normalized beta values. 


### Processing of BDR data

Following an identical processing pipeline for the independent replication dataset.

```
# Read pheno data
BDRpheno <- read.delim("BDR_LB_pheno_219.txt", sep = " ", header = TRUE)

# Read methylation data using minfi
rgSet <- read.metharray.exp(base="/BDR/LB/idats", extended = TRUE, force = TRUE)

rgSet
# class: RGChannelSetExtended 
# dim: 1051539 219 

# Sort pheno data frame to match rgSet: 
BDRpheno <- BDRpheno[with(BDRpheno, order(BDRpheno$Basename)),]
colnames(rgSet) == BDRpheno$Basename
# All TRUE

# Run primary wateRmelon QC 
qcReport(rgSet, sampNames=BDRpheno$Basename, sampGroups=BDRpheno$LB_Stage, pdf = "qcReport_BDR.pdf")

# Use minfy function to estimate cell counts based on DLPFC cell sorted reference data:
CellCounts <- estimateCellCounts(rgSet, compositeCellType = "DLPFC", cellTypes = c("NeuN_neg","NeuN_pos"))

BDRpheno$NeuN_neg <- CellCounts[,1]
BDRpheno$NeuN_pos <- CellCounts[,2]
BDRpheno$NeurProp <- with(BDRpheno, NeuN_pos/(NeuN_pos + NeuN_neg))

cor.test(BDRpheno$NeurProp, BDRpheno$LB_stage)
# 	Pearson's product-moment correlation not significant

# Filter out low quality probes and samples (wateRmelon):
pfltSet5 <- pfilter(rgSet, perc = 5)

# 0 samples having 5 % of sites with a detection p-value greater than 0.05 were removed  
# 1052 sites were removed as beadcount <3 in 5 % of samples 
# 45685 sites having 1 % of samples with a detection p-value greater than 0.05 were removed 

# Sex-check has been performed previously on the same samples by Exeter collaborators

mset.pflt <- preprocessRaw(pfltSet5)

mset.pflt
# class: MethylSet 
# dim: 819331 219 

# Identify outliers based on the wateRmelon outlyx function:
outliers <- outlyx(mset.pflt, plot = TRUE)
# No outliers

# Identify low quality samples based on minfi QC:
minfiQC <- getQC(mset.pflt)
pdf("~/Methylation/Brain/BDR/qcPlot_minfiQC_BDR.pdf")
plotQC(minfiQC)
dev.off()

# Identify samples failing minfi QC. Plot shows 19 samples failing, clustering below the following cuoff:
FailMinfiQC <- rownames(subset(minfiQC, mMed < 10.5 & uMed < 10.5))
Passing_QC <- rownames(BDRpheno[!rownames(BDRpheno) %in% FailMinfiQC,])

# Filter out failing samples from mset and pheno data frame: 
mset.pflt.sampleflt <- mset.pflt[,Passing_QC]
BDRpheno.flt <- BDRpheno[Passing_QC,]

mset.pflt.sampleflt
# class: MethylSet 
# dim: 819331 200 

# Map to genome (minfi)
gmset.pflt.sampleflt.dasen <- mapToGenome(mset.pflt.sampleflt.dasen)

gmset.pflt.sampleflt.dasen
# class: GenomicMethylSet 
# dim: 819331 200 

# Filter out probes on sex chromosomes
EPICanno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
keep <- !(featureNames(gmset.pflt.sampleflt.dasen.dropsnps) %in% EPICanno$Name[EPICanno$chr %in% c("chrX", "chrY")])
table(keep)
# keep

gmset.pflt.sampleflt.dasen.dropsnps.autosomes <- gmset.pflt.sampleflt.dasen.dropsnps[keep,]

gmset.pflt.sampleflt.dasen.dropsnps.autosomes
# class: GenomicMethylSet 
# dim: 775240 200  

# Remove cross-reactive probes:
cross_reactive_probes <- scan("~/Methylation/cross_reactive_probes_450k.txt", what = "character")
keep <- !(featureNames(gmset.pflt.sampleflt.dasen.dropsnps.autosomes) %in% cross_reactive_probes)
table(keep)

gmset.pflt.sampleflt.dasen.dropsnps.autosomes.xrflt <- gmset.pflt.sampleflt.dasen.dropsnps.autosomes[keep,]

gmset.pflt.sampleflt.dasen.dropsnps.autosomes.xrflt
# class: GenomicMethylSet 
# dim: 750785 200 

# Get betas: 
BDRbetas <- getBeta(gmset.pflt.sampleflt.dasen.dropsnps.autosomes.xrflt)

# Filter out probe-wise outliers using the wateRmelon pwod function. Outliers are probably low MAF/SNP heterozygotes
BDRbetas.pwoflt <- pwod(BDRbetas)
# 141826 probes detected.

# Export BDR betas, pheno and probe data for osca analyses
write.table(BDRbetas, "BDR_betas.txt", sep = " ", dec = ".", row.names = TRUE, col.names = TRUE) 
write.table(BDRpheno.flt, "BDR_pheno.txt", sep = " ", dec = ".", row.names = FALSE, col.names = TRUE)
EPICannoSubset <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)[match(rownames(BDRbetas), getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)$Name), c(1,2,3,22:ncol(getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)))]
write.table(EPICannoSubset[,1:4], "BDR_FinalProbes.txt", sep = " ", dec = ".", row.names = TRUE, col.names = FALSE, quote = FALSE)

```

Final filtered beta matrices were used for downstream analyses including linear regression in R.


