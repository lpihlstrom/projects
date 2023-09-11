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
pdf("~/Methylation/SortedLeukocytes/qcPlot_SexCheck_SortedLeukocytes.pdf")
plot(methylSex$xMed, methylSex$yMed, pch = 19, col = factor(databaseSex))
dev.off()
targets$PassSexCheck <- methylSex$predictedSex == databaseSex

# Rename samples to match methylation data and identify 1 sample failing sex check
targets$PlatePos_ID <- paste(targets$Slide,targets$Array, sep="_")
FailSexCheck <- subset(targets, PassSexCheck == FALSE, select = c(PlatePos_ID))

```

Initial QC tools take the raw RGChannelSet as input

```

# Run primary wateRmelon QC 
qcReport(rgSet, sampNames=targets$Sample_Name, sampGroups=targets$Celltype, pdf = "~/Methylation/SortedLeukocytes/qcReport_SortedLeukocytes.pdf")

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

```

A wide range of normalization methods are available. Given our dataset with samples from purified cell populations, we chose functional normalization implemented in the minfi package, as this is recommended for samples with large differences, such as different tissues or cell types: 

```

gmset.funnorm <- preprocessFunnorm(pfltSet, nPCs = 2, ratioConvert = FALSE)

# The output is an R object of class GenomicMethylSet with 857621 data points on 191 samples 

# Visualize effect of funnorm normalization on Type I vs Type II probe beta density in a random sample
mset.raw <- preprocessRaw(pfltSet)

pdf("~/Methylation/SortedLeukocytes/qcPlot_Norm_DensityPlots_SortedLeukocytes.pdf")
par(mfrow=c(1,2))
plotBetasByType(mset.raw[,1], main = "Raw")
funnormBetas <- getBeta(gmset.funnorm)
dim(funnormBetas)
probetypes <- data.frame(Name = 1:857621)
probetypes$Name <- rownames(funnormBetas)
probetypes$Type <- getProbeType(gmset.funnorm)
plotBetasByType(funnormBetas[,1], probeTypes = probetypes, main = "Funnorm")

# Visualize effect of funnorm normalization on cross-sample beta density
par(mfrow=c(1,2))
densityPlot(getBeta(mset.raw), main = "Raw")
densityPlot(funnormBetas, main = "Funnorm")
dev.off()

```

The preprocessFunnorm operation also converted the dataset into a GenomicMethylSet. A number of further filtering and QC steps are performed on this R object

```
# Identify low quality samples based on minfy QC:
minfyQC <- getQC(gmset.funnorm)
pdf("~/Methylation/SortedLeukocytes/qcPlot_minfiQC_SortedLeukocytes.pdf")
plotQC(minfyQC)
dev.off()

# Inspect plot for failing samples: No failing samples

# Filter out samples failing QC. In this project only one sample failing sex check
Samples_passing_QC <- rownames(targets[!rownames(targets) %in% FailSexCheck,])
gmset.funnorm.sampleflt <- gmset.funnorm[,Samples_passing_QC]
targets.flt <- targets[Samples_passing_QC,]

# Filter out probes with SNPs
gmset.funnorm.sampleflt.dropsnps <- dropLociWithSnps(gmset.funnorm.sampleflt)

# Filter out probes on sex chromosomes
EPICanno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
keep <- !(featureNames(gmset.funnorm.sampleflt.dropsnps) %in% EPICanno$Name[EPICanno$chr %in% c("chrX", "chrY")])
table(keep)
keep
#  FALSE   TRUE 
#  18674 810291
gmset.funnorm.sampleflt.dropsnps.autosomes <- gmset.funnorm.sampleflt.dropsnps[keep,]


# Remove cross-reactive probes:
cross_reactive_probes <- scan("~/Methylation/cross_reactive_probes_450k.txt", what = "character")
keep <- !(featureNames(gmset.funnorm.sampleflt.dropsnps.autosomes) %in% cross_reactive_probes)
table(keep)
# keep
#  FALSE   TRUE 
#  24611 785680

gmset.funnorm.sampleflt.dropsnps.autosomes.xrflt <- gmset.funnorm.sampleflt.dropsnps.autosomes[keep,]
gmset.funnorm.sampleflt.dropsnps.autosomes.xrflt

# The dataset is now an R object of class GenomicMethylSet with 785680 rows and 190 columns

```
Evaluation of potential batch effects and surrogate variable analysis requires normalized QC'd beta or M values. Technical replicates are kept for the CpGFilter step, then removed. 

```
# Get betas
betas.funnorm.sampleflt.dropsnps.autosomes.xrflt <- getBeta(gmset.funnorm.sampleflt.dropsnps.autosomes.xrflt)

# Assess technical vs. biological probe variability and plot result:
ICCfilter <- CpGFilterICC(betas.funnorm.sampleflt.dropsnps.autosomes.xrflt, targets.flt$X...Sample_Name)

pdf("~/Methylation/SortedLeukocytes/ICCplot.pdf")
qplot(ICCfilter, geom = "histogram")
dev.off

# Plot reveals ~250000 probes with minimal biological variability - could be filtered out. 

# Filter out low variability probes:
keep <- names(subset(ICCfilter, ICCfilter > 0.2)
finalBetas <- betas.funnorm.sampleflt.dropsnps.autosomes.xrflt[keep,]
dim(finalBetas)
# [1] 495220    190

# Filter out technical replicates:
replicates <- c("203265000010_R03C01","203265000060_R02C01","203293640091_R01C01","203259760075_R04C01")
keep <- rownames(targets.flt[!rownames(targets.flt) %in% replicates,])

finalBetas <- finalBetas[,keep]
targets.flt <- targets.flt[keep,]
dim(finalBetas)
# [1] 495220    186

# Filter out constituitively methylated or unmethylated probes
keep <- rowMeans(finalBetas, na.rm = TRUE) < 0.975
table(keep)
# keep
#  FALSE   TRUE 
#  1806 493414 
 
finalBetas <- finalBetas[keep,]
keep <- rowMeans(finalBetas, na.rm = TRUE) > 0.025
table(keep)
# keep
#  FALSE   TRUE 
#  10944 482470 
finalBetas <- finalBetas[keep,]
dim(finalBetas)
# [1] 482470    186

# Assess association between batch variables and principal components
# Calculate principal components
pca <- prcomp(finalBetas[complete.cases(finalBetas),])

# This command generates a pca object where the eigenvectors themselves are stored in the variable "rotation"
top_pcas <- pca$rotation[,1:10]

# Make sure categorical covariates are interpreted as factors: 
targets.flt$Chip <- as.factor(targets.flt$Slide)
targets.flt$Plate <- as.factor(targets.flt$Sample_Plate)      
targets.flt$Position <- as.factor(targets.flt$Array)
targets.flt$Sex <- as.factor(targets.flt$Gender)

# Evaluate effect of covariates on PCs
summary(lm(top_pcas[, "PC1"] ~ Plate + Position + Age + Sex + Celltype, targets.flt))
summary(lm(top_pcas[, "PC2"] ~ Plate + Position + Age + Sex + Celltype, targets.flt))
summary(lm(top_pcas[, "PC3"] ~ Plate + Position + Age + Sex + Celltype, targets.flt))
summary(lm(top_pcas[, "PC4"] ~ Plate + Position + Age + Sex + Celltype, targets.flt))
summary(lm(top_pcas[, "PC5"] ~ Plate + Position + Age + Sex + Celltype, targets.flt))
# No association between any PC and Plate or Position

# Evaluate potential hidden batch effects using surrogate variable analysis
# Requires comparing the main linear model with a null model
# Main aim of the study is to compare PDvsCO for each celltype. Create a new column in target (samplesheet): 
# Merge Sample_Group and Celltype to a new factor to the dataset that now should be the one to compare (f)
targets.flt$SampleGroup_Celltype <- paste(targets.flt$Sample_Group, targets.flt$Celltype, sep="_")

# Define group of interest:
f <- factor(targets.flt$SampleGroup_Celltype)
# Covariates to adjust for
Age <- as.numeric(targets.flt$Age)
Sex <- factor(targets.flt$Sex)

# Basic linear model without any batch adjustment
design <- model.matrix(~0+f+Sex+Age)
# Null model with just covariates
designNULL <- model.matrix(~0+Sex+Age)

# Assess number of latent factors using sva. We choose M values as the primary methylation measure in the linear regression.
mvalues <- beta2m(finalBetas)
n.sv.m = num.sv(mvalues[is.finite(rowSums(mvalues)),], design, method="leek")
n.sv.m
# [1] 1

# Estimate surrogate variable for M-values using (sva): 
svobj = sva(mvalues[is.finite(rowSums(mvalues)),], design, designNULL, n.sv = 1)

# Add SV1 to pheno object
targets.flt$SV1 <- as.numeric(svobj[[1]])

# add sv to design
SV1 <- as.numeric(targets.flt$SV1)
design <- model.matrix(~0+f+Sex+Age+SV1)
colnames(design) <- make.names(colnames(design))
```
Cell types were mixed on plates and arrays and processed together. Consequently, QC, normalization and evaluation of batch effects and surrogate variables was performed jointly on the combined dataset. We found however, that the wateRmelon pwod function for removal of individual outliers was unsuitable for combined data from multiple cell types. This step was performed separately on each cell type. 

```
IDs_CD14 <- rownames(subset(targets.flt, Celltype == "CD14"))
IDs_noCD14 <- rownames(subset(targets.flt, Celltype != "CD14"))
betas_CD14 <- finalBetas[,IDs_CD14]
betas_noCD14 <- finalBetas[,IDs_noCD14]
betas_CD14_pwoflt <- pwod(betas_CD14)
# 61082 probes detected. These individual outliers are coerced to NA. May reflect rare SNPs. 
betas_CD14_pwoflt <- cbind(betas_CD14_pwoflt, betas_noCD14)
betas_CD14_pwoflt <- betas_CD14_pwoflt[,rownames(targets.flt)]
Mvals_CD14_pwoflt <- beta2m(betas_CD14_pwoflt)
```

