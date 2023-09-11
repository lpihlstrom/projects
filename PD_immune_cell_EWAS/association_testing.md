### Primary statistical association testing with linear regression

Analyses were conducted in R 4.0.3 following quality control and normalization. 

```
library(limma)

# Model design includinng 1 SV defined after surogate variable analysis: 
design <- model.matrix(~0+f+Sex+Age+SV1)
colnames(design) <- make.names(colnames(design))

# Run linear regression model for CD14 cells using limma
fit <- lmFit(Mvals_CD14_pwoflt, design)
cm_cd14 <- makeContrasts(PDvsCO_cd14 = fPD_CD14 - fcontrol_CD14, levels=design)
fit2_cd14 <- contrasts.fit(fit, cm_cd14)
fit3_cd14 <- eBayes(fit2_cd14)

# Calculate lambda GC: 
CD14_pvals <- as.vector(fit3_cd14$p.value[,"PDvsCO_cd14"])
qchisq(median(na.omit(CD14_pvals)),1,lower.tail=FALSE)/qchisq(0.5,df=1)
# [[1] 1.033282

# Prepare annotation table for results file:
EPICannoSubset <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)[match(rownames(Mvals_CD14_pwoflt), getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)$Name), c(1,2,22:ncol(getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)))]

# Print full DMP results
topTable_PDvsCO_cd14_Mvals_1sv <- limma::topTable(fit3_cd14, coef = "PDvsCO_cd14", number = Inf, genelist = EPICannoSubset, adjust.method = "bonferroni", confint = TRUE, sort.by = "p")
write.table(topTable_PDvsCO_cd14_Mvals_1sv, "~/Methylation/SortedLeukocytes/results_PDvsCO_cd14_Mvals_1sv_toptable_full_Rev.txt", quote = FALSE, sep = " ")

# Testing sex-stratified model
# Create new groups stratified by sex: 
targets.flt$SampleGroup_Celltype_Sex <- paste(targets.flt$SampleGroup_Celltype, targets.flt$Sex, sep="_")

# Define sex stratified design
f_ss <- factor(targets.flt$SampleGroup_Celltype_Sex)
design_ss <- model.matrix(~0+f_ss+Age+SV1)
colnames(design_ss) <- make.names(colnames(design_ss))

# Male specific linear regression CD14:
fit <- lmFit(Mvals_CD14_pwoflt, design_ss)
cm_cd14_Male <- makeContrasts(PDvsCO_cd14_Male = f_ssPD_CD14_Male - f_sscontrol_CD14_Male, levels=design_ss)
fit2_cd14_Male <- contrasts.fit(fit, cm_cd14_Male)
fit3_cd14_Male <- eBayes(fit2_cd14_Male)

# Write result file CD14 male
topTable_PDvsCO_cd14_Mvals_1sv_Male <- limma::topTable(fit3_cd14_Male, coef = "PDvsCO_cd14_Male", number = Inf, genelist = EPICannoSubset, adjust.method = "bonferroni", confint = TRUE, sort.by = "p")
write.table(topTable_PDvsCO_cd14_Mvals_1sv_Male, "~/Methylation/SortedLeukocytes/results_PDvsCO_cd14_Mvals_1sv_Male_toptable_full_Rev.txt", quote = FALSE, sep = " ")

# Female specific linear regression CD14
cm_cd14_Female <- makeContrasts(PDvsCO_cd14_Female = f_ssPD_CD14_Female - f_sscontrol_CD14_Female, levels=design_ss)
fit2_cd14_Female <- contrasts.fit(fit, cm_cd14_Female)
fit3_cd14_Female <- eBayes(fit2_cd14_Female)

# Write result file CD14 female
topTable_PDvsCO_cd14_Mvals_1sv_Female <- limma::topTable(fit3_cd14_Female, coef = "PDvsCO_cd14_Female", number = Inf, genelist = EPICannoSubset, adjust.method = "bonferroni", confint = TRUE, sort.by = "p")
write.table(topTable_PDvsCO_cd14_Mvals_1sv_Female, "~/Methylation/SortedLeukocytes/results_PDvsCO_cd14_Mvals_1sv_Female_toptable_full_Rev.txt", quote = FALSE, sep = " ")

# Repeat for each cell type

```


### Analysis of differentially methylated regions using DMRcate

```
library(DMRcate)

# Analyze differentially methylated regions using DMRcate
myannotation <- cpg.annotate("array", allMvals_CD14_pwoflt, arraytype = "EPIC", analysis.type = "differential", design = design, contrasts = TRUE, cont.matrix = cm_cd14, coef ="PDvsCO_cd14", fdr = 0.05, what = "M")

dmrcateoutput <- dmrcate(myannotation, lambda = 1000, C = 2, min.cpgs = 2)
results.ranges_CD14 <- extractRanges(dmrcateoutput, genome = "hg19")

write.table(results.ranges_CD14, "~/Methylation/SortedLeukocytes/results_DMRcate_CD14_Rev.txt", quote = FALSE, sep = " ")

```

### Pathway analysis

```
library(missMethyl)

# Analyze pathway enrichment using gometh
sign_CpGs_CD14 <- rownames(subset(topTable_PDvsCO_cd14_Mvals_1sv, P.Value < 2.25e-8))
sign_CpGs_CD14_laxthresh <- rownames(subset(topTable_PDvsCO_cd14_Mvals_1sv, P.Value < 0.0001))
all_CpGs <- rownames(topTable_PDvsCO_cd14_Mvals_1sv)

GOanalysis_CD14 <- gometh(sig.cpg = sign_CpGs_CD14, all.cpg = all_CpGs, collection = "GO", array.type = "EPIC", sig.genes = TRUE)
GOanalysis_CD14_KEGG <- gometh(sig.cpg = sign_CpGs_CD14, all.cpg = all_CpGs, collection = "KEGG", array.type = "EPIC", sig.genes = TRUE)

GOanalysis_CD19_laxthresh <- gometh(sig.cpg = sign_CpGs_CD14_laxthresh, all.cpg = all_CpGs, collection = "GO", array.type = "EPIC", sig.genes = TRUE)
GOanalysis_CD19_KEGG_laxthresh <- gometh(sig.cpg = sign_CpGs_CD14_laxthresh, all.cpg = all_CpGs, collection = "KEGG", array.type = "EPIC", sig.genes = TRUE)

# No significant pathways

```

