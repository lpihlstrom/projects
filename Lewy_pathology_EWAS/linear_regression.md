## Primary statistical association testing with linear regression

Analyses were conducted in R 4.0.3 directly following quality control and normalization. 

```
# Load additional packages: 
library(sva)
library(limma)
BiocManager::install("bacon")
library(bacon)
install.packages("sailalithabollepalli-EpiSmokEr-54291ea.tar.gz", repos = NULL, type = "source")
library(EpiSmokEr)
```

### Association analysis in NBB data

Surrogate variables were estimated to control for potential unknown confounders

```
# Calculate top 10 surrogate variables for Braak aSyn model:
aSyn_design <- model.matrix(~ 0 + Braak_aSyn_stage + Age_death + Sex + PMD_min + NeurProp + Plate, data = aSyn_spctr_pheno)
aSyn_null_model <- model.matrix(~ 0 + Age_death + Sex + PMD_min + NeurProp + Plate, data = aSyn_spctr_pheno)

svobj_Braak_aSyn <- sva(aSyn_spctr_betas[complete.cases(aSyn_spctr_betas),],aSyn_design,aSyn_null_model, n.sv = 10)
# Number of significant surrogate variables is:  10 
# Iteration (out of 5 ):1  2  3  4  5  > 

# Add surrogate variables to pheno object:
aSyn_spctr_pheno$SV1 <- svobj_Braak_aSyn$sv[,1]
aSyn_spctr_pheno$SV2 <- svobj_Braak_aSyn$sv[,2]
aSyn_spctr_pheno$SV3 <- svobj_Braak_aSyn$sv[,3]
aSyn_spctr_pheno$SV4 <- svobj_Braak_aSyn$sv[,4]
aSyn_spctr_pheno$SV5 <- svobj_Braak_aSyn$sv[,5]
aSyn_spctr_pheno$SV6 <- svobj_Braak_aSyn$sv[,6]
aSyn_spctr_pheno$SV7 <- svobj_Braak_aSyn$sv[,7]
aSyn_spctr_pheno$SV8 <- svobj_Braak_aSyn$sv[,8]
aSyn_spctr_pheno$SV9 <- svobj_Braak_aSyn$sv[,9]
aSyn_spctr_pheno$SV10 <- svobj_Braak_aSyn$sv[,10]

```
Models with different numbers of SVs were tested. We found that including 3 SVs led the inflation measure lambda to fall below 1.2 in the NBB analysis, without further improvement when additional SVs were added.
The same number of SVs were used in alternative models to assess the impact on results using M-values, binarized Braak stage outcome or neuropathological diagnosis as a covariate.

```
# Analyse linear model including 3 SVs:
aSyn_design_3sv <- model.matrix(~ 0 + Braak_aSyn_stage + Age_death + Sex + PMD_min + NeurProp + Plate + SV1 + SV2 + SV3, aSyn_spctr_pheno)
aSyn_fit_3sv <- lmFit(aSyn_spctr_betas, aSyn_design_3sv)
aSyn_fit2_3sv <- eBayes(aSyn_fit_3sv)

# Calculate lambda GC: 
Braak_aSyn_pvals_3sv <- as.vector(aSyn_fit2_3sv$p.value[,"Braak_aSyn_stage"])
qchisq(median(Braak_aSyn_pvals_3sv),1,lower.tail=FALSE)/qchisq(0.5,df=1)
# [1] 1.176968

# Calculate alternative inflation metric using the bacon package: 
bc <- bacon(as.vector(aSyn_fit2_3sv$t[,"Braak_aSyn_stage"]))
inflation(bc)
# sigma.0 
# 1.076503

# Prepare annotation table for results file:
EPICannoSubset <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)[match(rownames(finalBetas), getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)$Name), c(1,2,3,22:ncol(getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)))]

# Write full result table for aSyn with 3SVs
write.table(topTable(aSyn_fit2_3sv, coef = "Braak_aSyn_stage", number = Inf, genelist = EPICannoSubset, confint = TRUE, adjust.method = "BH"), "results_Braak_aSyn_toptable_full.txt", quote = FALSE, sep = " ")


# Test model including 2 SVs:
aSyn_design_2sv <- model.matrix(~ 0 + Braak_aSyn_stage + Age_death + Sex + PMD_min + NeurProp + Plate + SV1 + SV2, aSyn_spctr_pheno)
aSyn_fit_2sv <- lmFit(aSyn_spctr_betas, aSyn_design_2sv)
aSyn_fit2_2sv <- eBayes(aSyn_fit_2sv)

# Calculate lambda GC: 
Braak_aSyn_pvals_2sv <- as.vector(aSyn_fit2_2sv$p.value[,"Braak_aSyn_stage"])
qchisq(median(Braak_aSyn_pvals_2sv),1,lower.tail=FALSE)/qchisq(0.5,df=1)
# [1] 1.812598

# Test model including 4 SVs:
aSyn_design_4sv <- model.matrix(~ 0 + Braak_aSyn_stage + Age_death + Sex + PMD_min + NeurProp + Plate + SV1 + SV2 + SV3 + SV4, aSyn_spctr_pheno)
aSyn_fit_4sv <- lmFit(aSyn_spctr_betas, aSyn_design_4sv)
aSyn_fit2_4sv <- eBayes(aSyn_fit_4sv)

# Calculate lambda GC: 
Braak_aSyn_pvals_4sv <- as.vector(aSyn_fit2_4sv$p.value[,"Braak_aSyn_stage"])
qchisq(median(Braak_aSyn_pvals_4sv),1,lower.tail=FALSE)/qchisq(0.5,df=1)
# [1] 1.201091

# Note minimum lambda observed when including 3SVs - below 1.2 which is in line with AD literature. 

```
Alternative models were explored using M-values, Braak stage as a binary variable and adjusting for neuropathological diagnosis.

```
# Repeat with M values
aSyn_spctr_Mvalues <- beta2m(aSyn_spctr_betas)
aSyn_fit_3sv_Mvals <- lmFit(aSyn_spctr_Mvalues, aSyn_design_3sv)
aSyn_fit2_3sv_Mvals <- eBayes(aSyn_fit_3sv_Mvals)

# Calculate lambda GC: 
Braak_aSyn_pvals_3sv_Mvals <- na.omit(as.vector(aSyn_fit2_3sv_Mvals$p.value[,"Braak_aSyn_stage"]))
qchisq(median(Braak_aSyn_pvals_3sv_Mvals),1,lower.tail=FALSE)/qchisq(0.5,df=1)
# [1] 1.182539

# Similar results with betas and M values
write.table(topTable(aSyn_fit2_3sv_Mvals, coef = "Braak_aSyn_stage", number = Inf, p.value = 0.1, genelist = EPICannoSubset, adjust.method = "BH"), "results_Braak_aSyn_toptable_Mvalues.txt", quote = FALSE, sep = " ")

# Analyze Braak as binary variable

# Generate variable using a cutoff of >= 3
aSyn_spctr_pheno$Braak_binary <- ifelse(aSyn_spctr_pheno$Braak_aSyn_stage>=3, 1, 0)

# Generate model
aSyn_binary_design_3sv <- model.matrix(~ 0 + Braak_binary + Age_death + Sex + PMD_min + NeurProp + Plate + SV1 + SV2 + SV3, aSyn_spctr_pheno)
aSyn_binary_fit_3sv <- lmFit(aSyn_spctr_betas, aSyn_binary_design_3sv)
aSyn_binary_fit2_3sv <- eBayes(aSyn_binary_fit_3sv)

Braak_aSyn_binary_pvals_3sv <- as.vector(aSyn_binary_fit2_3sv$p.value[,"Braak_binary"])

# Calculate lambda GC (updated filtering): 
qchisq(median(Braak_aSyn_binary_pvals_3sv),1,lower.tail=FALSE)/qchisq(0.5,df=1)
# [1] 1.089742

# Calculate alternative inflation metric using the bacon package: 
bc <- bacon(as.vector(aSyn_binary_fit2_3sv$t[,"Braak_binary"]))
inflation(bc)
# sigma.0 
# 1.037972

write.table(topTable(aSyn_binary_fit2_3sv, coef = "Braak_binary", number = Inf, genelist = EPICannoSubset, confint = TRUE, adjust.method = "BH"), "results_Braak_binary_toptable_full.txt", quote = FALSE, sep = " ")

# Analyze using diagnosis as covariate
aSyn_spctr_pheno$Dx <- as.factor(aSyn_spctr_pheno$Neuropath_GroupDx)
aSyn_design_3sv_Dx <- model.matrix(~ 0 + Braak_aSyn_stage + Age_death + Sex + PMD_min + NeurProp + Dx + Plate + SV1 + SV2 + SV3, aSyn_spctr_pheno)
aSyn_fit_3sv_Dx <- lmFit(aSyn_spctr_betas, aSyn_design_3sv_Dx)
aSyn_fit2_3sv_Dx <- eBayes(aSyn_fit_3sv_Dx)

Braak_aSyn_pvals_3sv_Dx <- as.vector(aSyn_fit2_3sv_Dx$p.value[,"Braak_aSyn_stage"])

# Calculate lambda GC (updated filtering): 
qchisq(median(Braak_aSyn_pvals_3sv_Dx),1,lower.tail=FALSE)/qchisq(0.5,df=1)
# [1] 1.084314

# Calculate alternative inflation metric using the bacon package: 
bc <- bacon(as.vector(aSyn_fit2_3sv_Dx$t[,"Braak_aSyn_stage"]))
inflation(bc)
# sigma.0 
# 1.035309 

write.table(limma::topTable(aSyn_fit2_3sv_Dx, coef = "Braak_aSyn_stage", number = Inf, genelist = EPICannoSubset, confint = TRUE, adjust.method = "BH"), "results_Braak_aSyn_adjDx_toptable_full.txt", quote = FALSE, sep = " ")
```
Tools to predict smoking status have been developed based on whole blood data have been shown to be of limited value in brain data. By reviewer's equest we nevertheless ran the smoking prediction algorithm, but could not reproduce an association with Lewy body disorders typically seen in epidemiological studies. 

```
# Estimate predicted smoking and assess association with phenotype

# Predict smoking using EpiSmokEr
result_All <- epismoker(dataset=aSyn_spctr_betas, samplesheet = aSyn_spctr_pheno, method = "all")

# Dataset has 155 of 187 CpGs required for smoking status estimation.
# Dataset has 2 of 4 CpGs required for Zhang et al method
# Dataset has 90 of 121 CpGs required for smoking status estimation.
 
# Inspect results
knitr::kable(result_All)
# Note: All individuals are predicted as former or current smokers. 

# Results come in a text format that needs to be processed into a data frame. Write results to file
write(knitr::kable(result_All), file = "results_EpiSmokEr.txt", sep = "")

# In bash command line:
sed 's/|//g' results_EpiSmokEr.txt | awk '{print $1, $2, $3}' | sed -e '2d' | head

# Read in and sort predicted smoking data:
predSmoke <- read.delim("EpiSmokEr_SSc_MS.txt", header=TRUE, row.names=1, sep = " ")
predSmoke$PlatePos_ID <- rownames(predSmoke)
predSmoke <- predSmoke[match(aSyn_spctr_pheno$PlatePos_ID, predSmoke$PlatePos_ID),]

# Assess correlation with Braak stage
cor.test(aSyn_spctr_pheno$Braak_aSyn_stage, predSmoke$methylationScore)
cor.test(aSyn_spctr_pheno$Braak_aSyn_stage, predSmoke$smokingScore)

# Assess differences across diagnostic groups with pairwise t-test:
pairwise.t.test(predSmoke$methylationScore, factor(aSyn_spctr_pheno$Neuropath_GroupDx))
pairwise.t.test(predSmoke$smokingScore, factor(aSyn_spctr_pheno$Neuropath_GroupDx))

# No significant associations

```


### Association analysis in BDR independent replication dataset

We included the same set of covariates in the replication analysis. 

```
# Ensure correct coding of factors
BDRpheno.flt$Sex <- as.factor(BDRpheno.flt$Gender)
BDRpheno.flt$Plate <- as.factor(BDRpheno.flt$Plate)


# Calculate top 5 surrogate variables for Braak aSyn model:
BDR_design <- model.matrix(~ 0 + LB_stage + Age + Sex + PMI + NeurProp + Plate, BDRpheno.flt)
BDR_null_model <- model.matrix(~ 0 + Age + Sex + PMI + NeurProp + Plate, BDRpheno.flt)

svobj_BDR <- sva(BDRbetas.pwoflt[complete.cases(BDRbetas.pwoflt),],BDR_design,BDR_null_model, n.sv = 5)
# Number of significant surrogate variables is:  5 
# Iteration (out of 5 ):1  2  3  4  5  > 

# Add surrogate variables to pheno object:
BDRpheno.flt$SV1 <- svobj_BDR$sv[,1]
BDRpheno.flt$SV2 <- svobj_BDR$sv[,2]
BDRpheno.flt$SV3 <- svobj_BDR$sv[,3]
BDRpheno.flt$SV4 <- svobj_BDR$sv[,4]
BDRpheno.flt$SV5 <- svobj_BDR$sv[,5]

# Run model including 3 SVs:
BDR_design_3sv <- model.matrix(~ 0 + LB_stage + Age + Sex + PMI + NeurProp + Plate + SV1 + SV2 + SV3, BDRpheno.flt)
BDR_fit_3sv <- lmFit(BDRbetas.pwoflt, BDR_design_3sv)
BDR_fit2_3sv <- eBayes(BDR_fit_3sv)

topTable_BDR_3sv <- topTable(BDR_fit2_3sv, coef = "LB_stage", number = Inf, genelist = EPICannoSubset, confint = TRUE, adjust.method = "BH")

# Calculate lambda GC: 
BDR_pvals_3sv <- as.vector(BDR_fit2_3sv$p.value[,"LB_stage"])

qchisq(median(BDR_pvals_3sv),1,lower.tail=FALSE)/qchisq(0.5,df=1)
# [1] 0.8538896

# Write full result table for BDR 3SVs
write.table(topTable_BDR_3sv, "BDR_N200_3SV_toptable_full.txt", quote = FALSE, sep = " ")

```


