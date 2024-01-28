## Code for quality control and normalization of whole blood Illumina MethylEPIC data

Analyses discribed here were performed using R 4.0.3

### Oslo dataset

The EPIC array experiment also included a subsert of DLB samples not included in the primary PD EWAS analysis. Raw data include a sample sheet and a total of 192 pairs of .idat files from Illumina MethylationEPIC arrays. These pairs include 4 technical replicates and 1 negative control. 

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
targets <- read.metharray.sheet("/work/users/lassep/EPIC_methylation_WholeBlood", pattern="WholeBlood_Samplesheet_Clean.csv")

# Read methylation data based on directions in sample sheet
rgSet <- read.metharray.exp(targets=targets, extended = TRUE)

# R object of class RGChannelSetExtended containts 1051815 data points on 671 samples

# Use minfy function to estimate cell counts based on cell sorted reference data:
CellCounts <- estimateCellCounts(rgSet, compositeCellType = "Blood", cellTypes = c("CD8T","CD4T","NK","Bcell","Mono","Gran"))
targets <- cbind(targets,CellCounts)

```
Different QC and normalization tools use different R objects as input, with implications for the order of operations. An initial non-normalized mapping to the genome is  performed for the sake of sex checks only. 

```
# Correct annotation
mset.raw <- preprocessRaw(rgSet)
annotation(mset.raw)[2] <- "ilm10b4.hg19"

gmset.raw <- mapToGenome(rgSet)

# Estimate sex based on methylation data
methylSex <- getSex(gmset.raw)

# Plot and compare with database sex
databaseSex <- as.factor(targets$Gender)
databaseSex <- substring(databaseSex, 1, 1)

# Get IDs of samples failing sex-check
targets$PassSexCheck <- methylSex$predictedSex == databaseSex
targets$PlatePos_ID <- paste(targets$Slide,targets$Array, sep="_")
FailSexCheck <- subset(targets, PassSexCheck == FALSE, select = c(PlatePos_ID))

# Samples failing sex-check included 3 PD and 1 control sample

```

Initial QC tools take the raw RGChannelSet as input

```

# Run primary wateRmelon QC 
qcReport(rgSet, sampNames=targets$Sample_ID, sampGroups=targets$Status, pdf = "/work/users/lassep/EPIC_methylation_WholeBlood/qcReport_WholeBlood.pdf")

# Filter out low quality probes and samples:
pfltSet <- pfilter(rgSet)

# 0 samples having 1 % of sites with a detection p-value greater than 0.05 were removed 
# Samples removed:  
# 14621 sites were removed as beadcount <3 in 5 % of samples 
# 4752 sites having 1 % of samples with a detection p-value greater than 0.05 were removed 

# Identify outliers based on the wateRmelon outlyx function:
pdf("~/Methylation/WholeBlood/qcPlot_Outlyx_WholeBlood.pdf")
outliers <- outlyx(pfltSet, plot = TRUE)
dev.off()
FailOutlyx <- subset(outliers, outliers == TRUE)
# 5 outliers identified

# Convert RGChannelSetExtended to MethylSet:
mset.pflt <- preprocessRaw(pfltSet)
mset.pflt

# The unfiltered MethylSet object includes data on 847262 probes and 671 samples (285 controls, 288 PD, 94 DLB and 4 technical replicates)

# Identity low quality samples based on minfi QC:
minfiQC <- getQC(mset.pflt)
pdf("~/Methylation/WholeBlood/qcPlot_minfiQC_WholeBlood.pdf")
plotQC(minfiQC)
dev.off()

# One failing sample (PD):
FailMinfyQC <- subset(minfyQC, mMed <= 10)

# Concatenate list of all samples failing (10 samples in total) and passing QC:

rownames(targets) <- targets$PlatePos_ID
Samples_failing_QC <- as.vector(c(rownames(FailMinfyQC), rownames(FailOutlyx), FailSexCheck$PlatePos_ID))
Samples_passing_QC <- rownames(targets[!rownames(targets) %in% Samples_failing_QC,])

# Filter out samples failing QC for various reasons from both MethylSet and phenotype data object
mset.pflt.sampleflt <- mset.pflt[,Samples_passing_QC]
targets.flt <- targets[Samples_passing_QC,]

# Remaining MethylSet includes data on 847262 probes and 661 individuals (278 Controls, 93 DLB and 281 PD)


```

A wide range of normalization methods are available. We use the dasen method recommended by the wateRmelon developers: 
```

# Analyse imprinting DMRs metric in prenormalized data to evaluate normalization:
raw_dmrse <- dmrse_row(getBeta(mset.pflt.sampleflt))
raw_dmrse
# [1] 0.00118403

# Normalize using the wateRmelon dasen method
mset.pflt.sampleflt.dasen <- dasen(mset.pflt.sampleflt)

# Analyze post-normalization imprinting DMRs (expected improvement shown as lower value):
dasen_dmrse <- dmrse_row(getBeta(mset.pflt.sampleflt.dasen))
dasen_dmrse
# [1] 0.0005411845

# Visualize effect of dasen normalization on Type I vs Type II probe beta density
pdf("/work/users/lassep/EPIC_methylation_WholeBlood/qcPlot_Norm_DensityPlots_WholeBlood.pdf")
par(mfrow=c(1,2))
plotBetasByType(mset.pflt.sampleflt[,1], main = "Raw")
plotBetasByType(mset.pflt.sampleflt.dasen[,1], main = "Dasen")

# Visualize effect of dasen normalization on cross-sample beta density
par(mfrow=c(1,2))
densityPlot(getBeta(mset.pflt.sampleflt), main = "Raw")
densityPlot(getBeta(mset.pflt.sampleflt.dasen), main = "dasen")
dev.off()


```

The normalized data is mapped to the genome, generating a GenomicMethylSet. A number of further filtering and QC steps are performed on this R object

```
## Map to genome
gmset.pflt.sampleflt.dasen <- mapToGenome(mset.pflt.sampleflt.dasen)

# Filter out probes with SNPs
gmset.pflt.sampleflt.dasen.dropsnps <- dropLociWithSnps(gmset.pflt.sampleflt.dasen)

# Filter out probes on sex chromosomes
EPICanno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
keep <- !(featureNames(gmset.pflt.sampleflt.dasen.dropsnps) %in% EPICanno$Name[EPICanno$chr %in% c("chrX", "chrY")])
table(keep)
keep
# FALSE   TRUE 
# 18495 800072
gmset.pflt.sampleflt.dasen.dropsnps.autosomes <- gmset.pflt.sampleflt.dasen.dropsnps[keep,]

# Remove cross-reactive probes:
cross_reactive_probes <- scan("~/Methylation/cross_reactive_probes_450k.txt", what = "character")
keep <- !(featureNames(gmset.pflt.sampleflt.dasen.dropsnps.autosomes) %in% cross_reactive_probes)
table(keep)
keep
# FALSE   TRUE 
# 24277 775795
gmset.pflt.sampleflt.dasen.dropsnps.autosomes.xrflt <- gmset.pflt.sampleflt.dasen.dropsnps.autosomes[keep,]
gmset.pflt.sampleflt.dasen.dropsnps.autosomes.xrflt

# The dataset is now an R object of class GenomicMethylSet with 775795 rows and 661 columns

```
The dataset is converted to methylation beta values. A large number of probes show high measurement errors when replicate samples are used to assess technical vs biological variation. These are filtered out using the CpGFilter package.
Finally, samples lacking age data are removed before beta matrix and phenotype data are exported for association analysis in OSCA.
```

# Get betas: 
betas.pflt.sampleflt.dasen.dropsnps.autosomes.xrflt <- getBeta(gmset.pflt.sampleflt.dasen.dropsnps.autosomes.xrflt)

# Assess technical vs. biological probe variability and plot result:
# First, make sample names where duplicates have identical names, by omitting last characters. 
targets.flt$repdesign <- substr(targets.flt$Sample_Name, start = 1, stop = 7)
# Next, run CpGFilterICC (CpGFilter):
ICCfilter <- CpGFilterICC(betas.pflt.sampleflt.dasen.dropsnps.autosomes.xrflt, targets.flt$repdesign)
# Note that method recommends minimum 8 replicates where we have 4. 

# Filter out low variability probes:
keep <- names(subset(ICCfilter, ICCfilter > 0.1))
all_wb_betas <- betas.pflt.sampleflt.dasen.dropsnps.autosomes.xrflt[keep,]
dim(all_wb_betas)
# [1] 502261    661

# Remove technical duplicates:
targets.flt.rmdup <- subset(targets.flt, TechRep == "N")
all_wb_betas <- betas.pflt.sampleflt.dasen.dropsnps.autosomes.xrflt[,rownames(targets.flt.rmdup)]
dim(all_wb_betas)
# [1] 502261    657

# Remove samples with no Age data
agecomplete <- rownames(subset(targets.flt.rmdup, Age != "NA"))
targets.flt.rmdup <- targets.flt.rmdup[agecomplete,]
all_wb_betas <- all_wb_betas[,agecomplete]
dim(all_wb_betas)
# [1] 502261    652

# Align sample names with genetic data for downstream mQTL analysis: 
rownames(targets.flt.rmdup) <- targets.flt.rmdup$Sample_Name
colnames(all_wb_betas) <- rownames(targets.flt.rmdup)

# Write final beta matrix for osca:
write.table(all_wb_betas, "~/Methylation/WholeBlood/WholeBloodfinalBetas.txt", sep = " ", dec = ".", row.names = TRUE, col.names = TRUE, quote = FALSE)

# Write pheno table: 
write.table(targets.flt.rmdup, "~/Methylation/WholeBlood/WholeBlood_Pheno_allGroups.txt", sep = " ", dec = ".", row.names = FALSE, col.names = TRUE, quote = FALSE)


```
In order to guide what covariates to adjust for in OSCA analyses, potential batch effects are evaluated by association with principal components

```
# Calculate principal components
pca <- prcomp(all_wb_betas[complete.cases(all_wb_betas),])

# This command generates a pca object where the eigenvectors themselves are stored in the variable "rotation"
top_pcas <- pca$rotation[,1:10]

# Make sure categorical covariates are interpreted as factors: 
targets.flt.rmdup$Chip <- as.factor(targets.flt.rmdup$Slide)
targets.flt.rmdup$Plate <- as.factor(targets.flt.rmdup$AMP_Plate)      
targets.flt.rmdup$Position <- as.factor(targets.flt.rmdup$Array)
targets.flt.rmdup$Sex <- as.factor(targets.flt.rmdup$Gender)
targets.flt.rmdup$Status <- as.factor(targets.flt.rmdup$Status)

# Evaluate effect of covariates on PCs

summary(lm(top_pcas[, "PC1"] ~ Plate + Position + Age + Sex + CD8T + CD4T + NK + Bcell + Mono + Status, data = targets.flt.rmdup))

# linear models for PCs 1-5 reveal associations with age, sex and plate, but not for plate position. These covariates are adjusted for in OSCA analysis.
```

### PPMI dataset



```
