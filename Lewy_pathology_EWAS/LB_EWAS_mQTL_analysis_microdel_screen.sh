###########################################################################
## mQTL analysis of tophits from frontal cortex Lewy body pathology EWAS ##
###########################################################################


## Lasse Pihlstrøm 2021 - 2022



# Residuals from CpC methylation linear models with all relevant covariates have been exported from main workflow in R

# Correct difference in ID syntax between genetic and epigenetic data: 
sed 's/YG0/YG/g' residuals_for_eqtl.txt > residuals_for_eqtl_IDsync.txt

# Assessing association between methylation residuals and SNPs with MAF > 0.05 within a 1Mb window surrounding the associated CpG: 


# cg14511218
/Users/lpihlstrom/Work/Tools/plink --bfile /Users/lpihlstrom/Work/Neuropathology_genetics/NBB_genetic_data/NeuroChip_HRCimp_rsq3_filtered_NBBpassQC --keep residuals_for_eqtl_IDsync.txt --pheno residuals_for_eqtl_IDsync.txt --pheno-name cg14511218_residuals --maf 0.05 --linear --ci 0.95 --out cg14511218_eqtl --chr 10 --from-bp 6462843 --to-bp 7462843 --keep-allele-order --allow-no-sex

sort -g -k12 cg14511218_eqtl.assoc.linear

# 2404 SNPs, Bonferroni-corrected threshold 1.9e-5, no significant mQTLs

# cg04011470
/Users/lpihlstrom/Work/Tools/plink --bfile /Users/lpihlstrom/Work/Neuropathology_genetics/NBB_genetic_data/NeuroChip_HRCimp_rsq3_filtered_NBBpassQC --keep residuals_for_eqtl_IDsync.txt --pheno residuals_for_eqtl_IDsync.txt --pheno-name cg04011470_residuals --maf 0.05 --linear --ci 0.95 --out cg04011470_eqtl --chr 8 --from-bp 22579330 --to-bp 21579330 --keep-allele-order --allow-no-sex

sort -g -k12 cg04011470_eqtl.assoc.linear

# 1799 SNPs, threshold 2.7e-5, no significant mQTLs

# cg07107199

/Users/lpihlstrom/Work/Tools/plink --bfile /Users/lpihlstrom/Work/Neuropathology_genetics/NBB_genetic_data/NeuroChip_HRCimp_rsq3_filtered_NBBpassQC --keep residuals_for_eqtl_IDsync.txt --pheno rev_residuals_for_eqtl_IDsync.txt --pheno-name cg07107199_residuals --maf 0.05 --linear --ci 0.95 --out cg07107199_eqtl --chr 1 --from-bp 204715911 --to-bp 205715911 --keep-allele-order --allow-no-sex

sort -g -k12 cg07107199_eqtl.assoc.linear

# 5988 SNPs, threshold 8.4e-6

# cg09985192

/Users/lpihlstrom/Work/Tools/plink --bfile /Users/lpihlstrom/Work/Neuropathology_genetics/NBB_genetic_data/NeuroChip_HRCimp_rsq3_filtered_NBBpassQC --keep residuals_for_eqtl_IDsync.txt --pheno rev_residuals_for_eqtl_IDsync.txt --pheno-name cg09985192_residuals --maf 0.05 --linear --ci 0.95 --out cg09985192_eqtl --chr 14 --from-bp 32297255 --to-bp 33297255 --keep-allele-order --allow-no-sex

sort -g -k12 cg09985192_eqtl.assoc.linear

# 5664 SNPs, threshold 8.8e-6, no significant mQTLs


# Query the mQTL results from Qi et al. Nat Commun 2018. 

smr_Mac --beqtl-summary ~/Downloads/Brain-mMeta/Brain-mMeta --query 1.0e-4 --probe cg14511218 --out cg14511218_mqtl_query 
# Probe not in dataset
smr_Mac --beqtl-summary ~/Downloads/Brain-mMeta/Brain-mMeta --query 1.0e-4 --probe cg04011470 --out cg04011470_mqtl_query
# Probe in dataset, no significant SNPs
smr_Mac --beqtl-summary ~/Downloads/Brain-mMeta/Brain-mMeta --query 1.0e-4 --probe cg07107199 --out cg07107199_mqtl_query
# Probe not in dataset
smr_Mac --beqtl-summary ~/Downloads/Brain-mMeta/Brain-mMeta --query 1.0e-4 --probe cg09985192 --out cg09985192_mqtl_query
# Probe in dataset, no significant SNPs
smr_Mac --beqtl-summary ~/Downloads/Brain-mMeta/Brain-mMeta --query 1.0e-4 --probe cg13986157 --out cg13986157_mqtl_query
# Probe not in dataset



# Analysis of runs of homozygosity on chromosome 22 to screen for 22q11 microdeletion: 
/Users/lpihlstrom/Work/Tools/plink --bfile NeuroChip_HRCimp_rsq3_filtered_NBBpassQC --homozyg --homozyg-kb 2500 --chr 22 --out NeuroChip_NBB_homozyg_chr22
