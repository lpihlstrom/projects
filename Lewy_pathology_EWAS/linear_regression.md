## Primary statistical association testing with linear regression

Analyses were conducted in R 4.0.3 directly following quality control and normalization. 

```
# Load additional packages: 
library(sva)
library(limma)
BiocManager::install("bacon")
library(bacon)
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

```


