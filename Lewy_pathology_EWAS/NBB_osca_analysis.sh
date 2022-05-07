#######################################################
## Brain methylation analysis of NBB data using OSCA ##
#######################################################


## Lasse Pihlstrøm - 2021 - 2022 ###



## PREPARE DATA FILES ##


# Osca uses a binary file format, BOD. This is generated from three text files: the beta matrix, a probe info file and a phenotype file. 

# Make oii file from R output: 
tail -n +2 aSyn_spctr_pheno.txt | awk '{print $1,$1,"0","0","NA"}' > NBB_aSyn_spctr.oii

# Make opi file from R output: 
sed -z 's/ \n/ NA\n/g' NBB_aSyn_spctr_FinalProbes.txt | sed 's/chr//g' | awk '{print $2,$1,$3,$5,$4}' > NBB_aSyn_spctr.opi

# Add "IID" to beginning of first line in beta matrix from R: 
sed -i '1s/^/IID /' aSyn_spctr_betas.txt 
 
# Make bod-file:
~/osca_Linux --tefile aSyn_spctr_betas.txt --methylation-beta --make-bod --no-fid --out NBB_aSyn_spctr


# Make discrete covar file for sex and plate: 
tail -n +2 aSyn_spctr_pheno.txt | awk '{print $1,$1,$5,$8}' > NBB_aSyn_spctr.covar

# Make quantitative covar file for age, PMD and estimated cell proportions: 
tail -n +2 aSyn_spctr_pheno.txt | awk '{print $1,$1,$4,$6,$7}' > NBB_aSyn_spctr.qcovar

# Make pheno file:
tail -n +2 aSyn_spctr_pheno.txt | awk '{print $1,$1,$3}' > NBB_aSyn_spctr.phen


# Adjust betas for covariates: 
~/osca_Linux --befile NBB_aSyn_spctr --upper-beta 0.975 --lower-beta 0.025 --covar NBB_aSyn_spctr.covar --qcovar NBB_aSyn_spctr.qcovar --adj-probe --make-bod --out NBB_aSyn_spctr_adj



## Run EWAS analyses ##


# MOA analysis of Braak aSyn stage:
~/osca_Linux --moa --befile NBB_aSyn_spctr_adj --pheno NBB_aSyn_spctr.phen  --out NBB_aSyn_spctr_adj_moa

~/osca_Linux --moment --befile NBB_aSyn_spctr_adj --pheno NBB_aSyn_spctr.phen  --out NBB_aSyn_spctr_adj



## Analysis of variance explained ##

#Generate omics relationship matrix (ORM):
~/osca_Linux --befile NBB_aSyn_spctr --make-orm --out NBB_orm

# Run REML analysis:
~/osca_Linux --reml --orm NBB_orm --pheno NBB_aSyn_spctr.phen --covar NBB_aSyn_spctr.covar --qcovar NBB_aSyn_spctr.qcovar --out NBB_reml

# High variance explained: V(O)/Vp	0.806632	0.182043

# Test if results are sensitive to including surrogate variables as covariates:
# Add SVs to covar file and rerun: 
tail -n +2 aSyn_spctr_pheno.txt | awk '{print $1,$1,$4,$6,$7,$15,$16,$17}' > NBB_aSyn_spctr_3SV.qcovar
~/osca_Linux --reml --orm NBB_orm --pheno NBB_aSyn_spctr.phen --covar NBB_aSyn_spctr.covar --qcovar NBB_aSyn_spctr_3SV.qcovar --out NBB_3SV_reml

# Results essentially equivalent


 