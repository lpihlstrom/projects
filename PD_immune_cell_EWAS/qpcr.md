## Analyses in R following quantitative PCR of RAB32 and selected reference genes

```
# Read phenotype and top probe data
data <- read.delim("SortedLeukocytesFinalPheno.txt", sep = " ", row.names = 1, header = TRUE)
data$SampleGroup_Celltype <- paste(data$Sample_Group, data$Celltype, sep = "_")
data$SampleGroup_Celltype <- factor(data$SampleGroup_Celltype , levels=c("control_CD4", "PD_CD4", "control_CD8", "PD_CD8", "control_CD19", "PD_CD19", "control_CD14", "PD_CD14"))

# Read data from ViiA7 qpcr experiment
qpcr <- read.delim("Monocytes_qPCR_RAB32_ViiA7export.txt", sep = "\t", row.names = 1, header = TRUE)

# Test correlation between reference genes
cor.test(qpcr$CT_ACTB, qpcr$CT_RPL37A)
cor.test(qpcr$CT_ACTB, qpcr$CT_B2M)
cor.test(qpcr$CT_B2M, qpcr$CT_RPL37A)

# By far the best correlation was seen for B2M and RLP37A (Pearson 0.77, p-value = 1.063e-09 )
# Generate Normalozation factor based on these two: 
qpcr$NF <- sqrt(qpcr$CT_B2M*qpcr$CT_RPL37A)

# Calculate quantity as "delta-CT"
qpcr$quant <- 2^-(qpcr$CT_RAB32-qpcr$NF)

# Calculate quantity relative to the mean across all samples
qpcr$RAB32_rel_quant <- qpcr$quant/mean(na.omit(qpcr$quant))

# Edit rownames to specify monocyte data and prepare merge with data frame from methylation analysis:
rownames(qpcr) <- paste(rownames(qpcr), "CD14", sep = "")

# Merge data
data <- merge(data, qpcr, by = 0, all.x = TRUE, all.y = FALSE)

# Perform linear regression including probes from the RAB32 DMR or disease status as predictors in the model
summary(lm(RAB32_rel_quant ~ cg05420134 + Age + Gender, data = data))
summary(lm(RAB32_rel_quant ~ cg12872357 + Age + Gender, data = data))
summary(lm(RAB32_rel_quant ~ Sample_Group + Age + Gender, data = data))

# Perform logistic regression with disease status as outcome and RAB32 expression as predictor variable
summary(glm(as.factor(Sample_Group) ~ Age + Gender + RAB32_rel_quant, data = data, family = binomial))
# Using Z-score:
data$RAB32_Z <- (data$RAB32_rel_quant-1)/sd(na.omit(data$RAB32_rel_quant))
summary(glm(as.factor(Sample_Group) ~ Age + Gender + data$RAB32_Z, data = data, family = binomial))
summary(lm(RAB32_Z ~ Sample_Group + Age + Gender, data = data))

# Test correlation between methylation and mRNA expression
cor.test(data$RAB32_rel_quant, data$cg05420134)
cor.test(data$RAB32_rel_quant, data$cg12872357)

# Generate methylation-expression scatterplot with regression line
tiff(filename = "~/Work/Methylation_projects/Blood/Sorted_leukocytes/Results/RAB32_scatterplot.tiff", units="in", height = 4, width=6, res=1200)
plot(data$cg05420134, data$RAB32_rel_quant, xlim = c(0.5,0.8), col = ifelse(data$Sample_Group == "PD", "red", "darkgreen"), pch = ifelse(data$Sample_Group == "PD", 17, 19), frame.plot = FALSE, xlab = "Methylation beta", ylab = "RAB32 relative expression")
abline(lm(RAB32_rel_quant.x ~ cg05420134, data = data), col = "blue", lwd = 2)
dev.off()

```
