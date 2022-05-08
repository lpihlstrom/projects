#######################################################################
## Lewy pathology EWAS - R code for meta-analysis and selected plots ##
#######################################################################



## Lasse Pihlstrøm 2021-2022

## Install and load required packages ##

BiocManager::install("meta")
library(meta)
BiocManager::install("qvalue")
library(qvalue)
nBiocManager::install("stats")
library(stats)
library(qqman)



## Load and organize linear regression results from main analysis ##

# Read results from linear regression in NBB data and generate new data frame with key variables:
toptable <- read.delim("~/Work/Methylation_projects/Brain/results_Braak_aSyn_toptable_full.txt", sep = ",", row.names = 1, header = TRUE)
nbb <- as.data.frame(toptable[,"logFC"]*100)
rownames(nbb) <- rownames(toptable)
colnames(nbb)[[1]] <- "Beta_nbb"
nbb$SE_nbb <- (toptable$logFC*100 - toptable$CI.L*100)/1.96
nbb$P_nbb <- toptable$P.Value
nbb$P_nbb_adjust <- toptable$adj.P.Val

# Write full summary stats to file:
write.table(nbb, file = "NBB_linear_regression.txt", quote = FALSE, sep = " ", row.names = TRUE, col.names = TRUE)

# Read results from linear regression in BDR data and generate new data frame with key variables:
toptable_bdr <- read.delim("~/Work/Methylation_projects/Brain/Rep_meta_Exeter/BDR_N200_3SV_toptable_full.txt", sep = " ", row.names = 1, header = TRUE)
bdr <- as.data.frame(toptable_bdr[,"logFC"]*100)
rownames(bdr) <- rownames(toptable_bdr)
colnames(bdr)[[1]] <- "Beta_bdr"
bdr$SE_bdr <- (toptable_bdr$logFC*100 - toptable_bdr$CI.L*100)/1.96
bdr$P_bdr <- toptable_bdr$P.Value
bdr$P_bdr_adjust <- toptable_bdr$adj.P.Val

# Write full smmary stats to file:
write.table(bdr, file = "BDR_linear_regression.txt", quote = FALSE, sep = " ", row.names = TRUE, col.names = TRUE)

# Sort results
nbb <- nbb[sort(rownames(nbb)),]
bdr <- bdr[sort(rownames(bdr)),]

# Merge results
probes <- intersect(rownames(nbb), rownames(bdr))
results_nbb_bdr <- cbind(nbb[probes,],bdr[probes,])



## Perform meta-analysis ##

# For loop to add p-value of meta-analysis to results table
for (n in 1:564387) {results_nbb_bdr$P_meta_fixed[[n]] <- print(metagen(c(results_nbb_bdr$Beta_nbb[[n]], results_nbb_bdr$Beta_bdr[[n]]), c(results_nbb_bdr$SE_nbb[[n]], results_nbb_bdr$SE_bdr[[n]]))$pval.fixed)}

# For loop to add beta values of meta-analysis to results table
for (n in 1:564387) {results_nbb_bdr$Beta_meta_fixed[[n]] <- print(metagen(c(results_nbb_bdr$Beta_nbb[[n]], results_nbb_bdr$Beta_bdr[[n]]), c(results_nbb_bdr$SE_nbb[[n]], results_nbb_bdr$SE_bdr[[n]]))$TE.fixed)}

# For loop to add SE values of meta-analysis to results table
for (n in 1:564387) {results_nbb_bdr$SE_meta_fixed[[n]] <- print(metagen(c(results_nbb_bdr$Beta_nbb[[n]], results_nbb_bdr$Beta_bdr[[n]]), c(results_nbb_bdr$SE_nbb[[n]], results_nbb_bdr$SE_bdr[[n]]))$seTE.fixed)}

# Calculate FDR by the BH method
results_nbb_bdr$P_meta_fixed_BH <- p.adjust(results_nbb_bdr$P_meta_fixed, method = "BH")

# Calculate lambda for meta-analysis
pvals_meta <- unlist(results_nbb_bdr$P_meta_fixed)
qchisq(median(pvals_meta),1,lower.tail = FALSE)/qchisq(0.5,df=1)



## Integrate osca MOA results and prepare tables and figures ##

# Import, organize, and merge osca results
nbb_osca <- read.delim("NBB_aSyn_spctr_adj_moa.moa", sep = "\t", row.names = 2, header = TRUE)
nbb_osca <- nbb_osca[sort(rownames(nbb_osca)),]

colnames(nbb_osca)[[5]] <- "b_nbb_osca"
colnames(nbb_osca)[[6]] <- "se_nbb_osca"
colnames(nbb_osca)[[7]] <- "p_nbb_osca"
nbb_osca$P_nbb_osca_BH_adjust <- p.adjust(nbb_osca$p_nbb_osca, method = "BH")

bdr_osca <- read.delim("BDR_adj_moa.moa", sep = "\t", row.names = 2, header = TRUE)
bdr_osca <- bdr_osca[sort(rownames(bdr_osca)),]

colnames(bdr_osca)[[5]] <- "b_BDR_osca"
colnames(bdr_osca)[[6]] <- "se_BDR_osca"
colnames(bdr_osca)[[7]] <- "p_BDR_osca"

# Calculate lambda for MOA analysis
pvals_meta <- unlist(results_nbb_bdr$P_meta_fixed)
qchisq(median(nbb_osca$p_nbb_osca),1,lower.tail = FALSE)/qchisq(0.5,df=1)

# Merge NBB and BDR osca MOA results
probes_osca <- intersect(rownames(nbb_osca), rownames(bdr_osca))
results_nbb_bdr_osca <- cbind(nbb_osca[probes_osca,5:7],bdr_osca[probes_osca,5:7])

# Merge R linear regression and osca MOA results
probes_R_osca <- intersect(rownames(results_nbb_bdr), rownames(results_nbb_bdr_osca))
results_nbb_bdr_R_osca <- cbind(results_nbb_bdr[probes_R_osca,], results_nbb_bdr_osca[probes_R_osca,])

# Sort by discovery and meta-analysis p-values respectively:
ordered_by_nbb <- results_nbb_bdr_R_osca[order(unlist(results_nbb_bdr_R_osca$P_nbb)),]
ordered_by_meta <- results_nbb_bdr_R_osca[order(unlist(results_nbb_bdr_R_osca$P_meta_fixed)),]

# Extract probes significant in R discovery pipeline
sign_nbb_BH <- subset(ordered_by_nbb, P_nbb_adjust < 0.05)

# Add chromosome and position and update column name: 
sign_nbb_BH <- cbind(sign_nbb_BH, toptable[rownames(sign_nbb_BH),]$chr)
sign_nbb_BH <- cbind(sign_nbb_BH, toptable[rownames(sign_nbb_BH),]$pos)
colnames(sign_nbb_BH)[[19]] <- "chr"
colnames(sign_nbb_BH)[[20]] <- "pos"

# Add gene annotation and update column name: 
sign_nbb_BH <- cbind(sign_nbb_BH, toptable[rownames(sign_nbb_BH),]$UCSC_RefGene_Name)
colnames(sign_nbb_BH)[[21]] <- "UCSC_RefGene_Name"

# Write result to file:
write.table(as.matrix(sign_nbb_BH), file = "NBB_BH_sign.txt", quote = FALSE, sep = " ", row.names = TRUE, col.names = TRUE)

# Generate beta-beta plot and save to file:
tiff(filename = "Beta_beta.tiff", units="in", height = 6, width=6, res=600)
plot(sign_nbb_BH$Beta_nbb,sign_nbb_BH$Beta_bdr, pch = 19, col = "blue", bty = "n", cex = 1.1, cex.axis = 1.2, xlab = "", ylab = "", ylim = c(-0.4,0.6), xlim = c(-0.9,0.9))
segments(sign_nbb_BH$Beta_nbb-sign_nbb_BH$SE_nbb,sign_nbb_BH$Beta_bdr,sign_nbb_BH$Beta_nbb+sign_nbb_BH$SE_nbb,sign_nbb_BH$Beta_bdr, col = "blue")
segments(sign_nbb_BH$Beta_nbb,sign_nbb_BH$Beta_bdr-sign_nbb_BH$SE_bdr,sign_nbb_BH$Beta_nbb,sign_nbb_BH$Beta_bdr+sign_nbb_BH$SE_bdr, col = "blue")
abline(v = 0, col = "red")
abline(h = 0, col = "red")
mtext(side=1, line=3, "Effect size NBB", font=2, cex=1.2)
mtext(side=2, line=3, "Effect size BDR", font=2, cex=1.2)
dev.off()

# Extract CpGs that are both meta-analysis FDR-significant and nominally significant in both NBB and BDR: 
sign_meta_BH_both_nominal <- subset(subset(ordered_by_meta, P_meta_fixed_BH < 0.05), P_bdr < 0.05)

# Add chromosome and position and update column name: 
sign_meta_BH_both_nominal <- cbind(sign_meta_BH_both_nominal, toptable[rownames(sign_meta_BH_both_nominal),]$chr)
sign_meta_BH_both_nominal <- cbind(sign_meta_BH_both_nominal, toptable[rownames(sign_meta_BH_both_nominal),]$pos)
colnames(sign_meta_BH_both_nominal)[[19]] <- "chr"
colnames(sign_meta_BH_both_nominal)[[20]] <- "pos"

# Add gene annotation and update column name: 
sign_meta_BH_both_nominal <- cbind(sign_meta_BH_both_nominal, toptable[rownames(sign_meta_BH_both_nominal),]$UCSC_RefGene_Name)
colnames(sign_meta_BH_both_nominal)[[21]] <- "UCSC_RefGene_Name"

# Write result to file:
write.table(as.matrix(format(sign_meta_BH_both_nominal, digits = 3)), file = "meta_BH_sign_both_nominal_sign.txt", quote = FALSE, sep = " ", row.names = TRUE, col.names = TRUE)


# Generate linear regression vs MOA Z-Z plot and save to file:
tiff(filename = "~/Work/Methylation_projects/Brain/Manuscript/Revision/ZscoreR_ZscoreOSCA.tiff", units="in", height = 6, width=6, res=600)
plot(sign_nbb_BH$Beta_nbb/sign_nbb_BH$SE_nbb,sign_nbb_BH$b_nbb_osca/sign_nbb_BH$se_nbb_osca, pch = 19, col = "blue", bty = "n", cex = 1.1, cex.axis = 1.2, xlab = "", ylab = "")
abline(v = 0, col = "red")
abline(h = 0, col = "red")
mtext(side=1, line=3, "Z score linear regression", font=2, cex=1.2)
mtext(side=2, line=3, "Z score MOA", font=2, cex=1.2)
dev.off()

# Make Manhattan and QQ plot for NBB results

toptable$probe <- rownames(toptable)
toptable$CHR <- as.numeric(substr(toptable$chr, start = 4, stop = 6))
tiff(filename = "~/Work/Methylation_projects/Brain/Manuscript/Revision/NBB_Manhattan.tiff", units="in", height = 3, width=6, res=600)
manhattan(toptable, chr = "CHR", bp = "pos", p = "P.Value", snp = "probe", col = c("lightblue", "darkblue"), chrlabs = c(1:22), suggestiveline = FALSE, genomewideline = FALSE)
dev.off()

# Make Manhattan and QQ plot for BDR results

toptable_bdr$probe <- rownames(toptable_bdr)
toptable_bdr$CHR <- as.numeric(substr(toptable_bdr$chr, start = 4, stop = 6))
tiff(filename = "~/Work/Methylation_projects/Brain/Manuscript/Revision/BDR_Manhattan.tiff", units="in", height = 3, width=6, res=600)
manhattan(toptable_bdr, chr = "CHR", bp = "pos", p = "P.Value", snp = "probe", col = c("lightblue", "darkblue"), chrlabs = c(1:22), suggestiveline = FALSE, genomewideline = FALSE)
dev.off()

tiff(filename = "~/Work/Methylation_projects/Brain/Manuscript/Revision/BDR_QQ.tiff", units="in", height = 6, width=6, res=600)
qq(toptable_bdr$P.Value)
dev.off()

tiff(filename = "~/Work/Methylation_projects/Brain/Manuscript/Revision/NBB_QQ.tiff", units="in", height = 6, width=6, res=600)
qq(toptable$P.Value)
dev.off()

# Make Manhattan and QQ plot for Meta-analysis results

results_nbb_bdr <- cbind(results_nbb_bdr, toptable[rownames(results_nbb_bdr),]$chr)
colnames(results_nbb_bdr)[[13]] <- "CHR"
results_nbb_bdr$CHR <- as.numeric(substr(results_nbb_bdr$CHR, start = 4, stop = 6))

results_nbb_bdr <- cbind(results_nbb_bdr, toptable[rownames(results_nbb_bdr),]$pos)
colnames(results_nbb_bdr)[[14]] <- "pos"

results_nbb_bdr$P_meta_fixed <- as.numeric(results_nbb_bdr$P_meta_fixed)

tiff(filename = "~/Work/Methylation_projects/Brain/Manuscript/Revision/Meta_Manhattan.tiff", units="in", height = 3, width=6, res=600)
manhattan(results_nbb_bdr, chr = "CHR", bp = "pos", p = "P_meta_fixed", col = c("lightblue", "darkblue"), chrlabs = c(1:22), suggestiveline = FALSE, genomewideline = FALSE)
dev.off()

tiff(filename = "~/Work/Methylation_projects/Brain/Manuscript/Revision/Meta_QQ.tiff", units="in", height = 6, width=6, res=600)
qq(results_nbb_bdr$P_meta_fixed)
dev.off()

# Make Manhattan and QQ plot for MOA results

tiff(filename = "~/Work/Methylation_projects/Brain/Manuscript/Revision/MOA_Manhattan.tiff", units="in", height = 3, width=6, res=600)
manhattan(nbb_osca, chr = "Chr", bp = "bp", p = "p_nbb_osca", col = c("lightblue", "darkblue"), chrlabs = c(1:22), suggestiveline = FALSE, genomewideline = FALSE)
dev.off()

tiff(filename = "~/Work/Methylation_projects/Brain/Manuscript/Revision/MOA_QQ.tiff", units="in", height = 6, width=6, res=600)
qq(nbb_osca$p_nbb_osca)
dev.off()



## Load results from alternative linear regression models and generate Z-Z plots ##


# Add M value results
toptable_NBB_Mvals <- read.delim("~/Work/Methylation_projects/Brain/results_Braak_aSyn_toptable_Mvalues.txt", sep = " ", row.names = 1, header = TRUE)
nbb_Mvals <- as.data.frame(toptable_NBB_Mvals[,"logFC"]*100)
rownames(nbb_Mvals) <- rownames(toptable_NBB_Mvals)
colnames(nbb_Mvals)[[1]] <- "b_nbb_Mvals"
nbb_Mvals$SE_nbb_Mvals <- (toptable_NBB_Mvals$logFC*100 - toptable_NBB_Mvals$CI.L*100)/1.96
nbb_Mvals$P_nbb_Mvals <- toptable_NBB_Mvals$P.Value
nbb_Mvals$P_nbb_Mvals_adjust <- toptable_NBB_Mvals$adj.P.Val

# Write full smmary stats to file:
write.table(nbb_Mvals, file = "~/Work/Methylation_projects/Brain/Manuscript/Revision/Data_for_submission/NBB_linear_regression_Mvals.txt", quote = FALSE, sep = " ", row.names = TRUE, col.names = TRUE)

# Add Mvals to results
sign_nbb_BH_all <- cbind(sign_nbb_BH,nbb_Mvals[rownames(sign_nbb_BH),])

#Generate beta vs Mval plot
tiff(filename = "~/Work/Methylation_projects/Brain/Manuscript/Revision/Beta_Mval.tiff", units="in", height = 6, width=6, res=600)
plot(sign_nbb_BH_all$Beta_nbb/sign_nbb_BH_all$SE_nbb,sign_nbb_BH_all$b_nbb_Mvals/sign_nbb_BH_all$SE_nbb_Mvals, pch = 19, col = "blue", bty = "n", cex = 1.1, cex.axis = 1.2, xlab = "", ylab = "")
abline(v = 0, col = "red")
abline(h = 0, col = "red")
mtext(side=1, line=3, "Z score beta values", font=2, cex=1.2)
mtext(side=2, line=3, "Z score M values", font=2, cex=1.2)
dev.off()


# Add results generated using Braak stage as a binary variable
toptable_NBB_Braak_binary <- read.delim("~/Work/Methylation_projects/Brain/results_Braak_binary_toptable_full.txt", sep = " ", row.names = 1, header = TRUE)
nbb_Braak_binary <- as.data.frame(toptable_NBB_Braak_binary[,"logFC"]*100)
rownames(nbb_Braak_binary) <- rownames(toptable_NBB_Braak_binary)
colnames(nbb_Braak_binary)[[1]] <- "b_nbb_Braak_binary"
nbb_Braak_binary$SE_nbb_Braak_binary <- (toptable_NBB_Braak_binary$logFC*100 - toptable_NBB_Braak_binary$CI.L*100)/1.96
nbb_Braak_binary$P_nbb_Braak_binary <- toptable_NBB_Braak_binary$P.Value
nbb_Braak_binary$P_nbb_Braak_binary_adjust <- toptable_NBB_Braak_binary$adj.P.Val

# Write full summary stats to file:
write.table(nbb_Braak_binary, file = "~/Work/Methylation_projects/Brain/Manuscript/Revision/Data_for_submission/NBB_linear_regression_Braak_binary.txt", quote = FALSE, sep = " ", row.names = TRUE, col.names = TRUE)

# Add Braak binary to results
sign_nbb_BH_all <- cbind(sign_nbb_BH_all,nbb_Braak_binary[rownames(sign_nbb_BH),])

#Generate quantitative vs binary plot
tiff(filename = "~/Work/Methylation_projects/Brain/Manuscript/Revision/Braak_quant_binary.tiff", units="in", height = 6, width=6, res=600)
plot(sign_nbb_BH_all$Beta_nbb/sign_nbb_BH_all$SE_nbb,sign_nbb_BH_all$b_nbb_Braak_binary/sign_nbb_BH_all$SE_nbb_Braak_binary, pch = 19, col = "blue", bty = "n", cex = 1.1, cex.axis = 1.2, xlab = "", ylab = "")
abline(v = 0, col = "red")
abline(h = 0, col = "red")
mtext(side=1, line=3, "Z score Braak quantitative", font=2, cex=1.2)
mtext(side=2, line=3, "Z score Braak binary", font=2, cex=1.2)
dev.off()


# Add results adjusting for diagnosis
toptable_NBB_adjDx <- read.delim("~/Work/Methylation_projects/Brain/results_Braak_aSyn_adjDx_toptable_full.txt", sep = " ", row.names = 1, header = TRUE)
nbb_adjDx <- as.data.frame(toptable_NBB_adjDx[,"logFC"]*100)
rownames(nbb_adjDx) <- rownames(toptable_NBB_adjDx)
colnames(nbb_adjDx)[[1]] <- "b_nbb_adjDx"
nbb_adjDx$SE_nbb_adjDx <- (toptable_NBB_adjDx$logFC*100 - toptable_NBB_adjDx$CI.L*100)/1.96
nbb_adjDx$P_nbb_adjDx <- toptable_NBB_adjDx$P.Value
nbb_adjDx$P_nbb_adjDx_adjust <- toptable_NBB_adjDx$adj.P.Val

# Add Dx-adjusted data to results
sign_nbb_BH_all <- cbind(sign_nbb_BH_all,nbb_adjDx[rownames(sign_nbb_BH),])

#Generate Dx-adjusted vs non-adjusted plot
tiff(filename = "~/Work/Methylation_projects/Brain/Manuscript/Revision/Braak_adjDx_ZZplot.tiff", units="in", height = 6, width=6, res=600)
plot(sign_nbb_BH_all$Beta_nbb/sign_nbb_BH_all$SE_nbb,sign_nbb_BH_all$b_nbb_adjDx/sign_nbb_BH_all$SE_nbb_adjDx, pch = 19, col = "blue", bty = "n", cex = 1.1, cex.axis = 1.2, xlab = "", ylab = "")
abline(v = 0, col = "red")
abline(h = 0, col = "red")
mtext(side=1, line=3, "Z score not adjusted for diagnosis", font=2, cex=1.2)
mtext(side=2, line=3, "Z score adjusted for diagnosis", font=2, cex=1.2)
dev.off()



# Write combined results on all NBB FDR-significant probes to file:
write.table(as.matrix(format(sign_nbb_BH_all, digits = 3)), file = "sign_nbb_BH_for_supplementary_info.txt", quote = FALSE, sep = " ", row.names = TRUE, col.names = TRUE)



