#######################################################
## Brain methylation analysis of BDR data using OSCA ##
#######################################################


## Lasse Pihlstrøm - 2021 - 2022 ###



## PREPARE DATA FILES ##


# Osca uses a binary file format, BOD. This is generated from three text files: the beta matrix, a probe info file and a phenotype file. 

# Make oii file from R output: 
tail -n +2 BDR_pheno.txt | awk '{print $1,$1,"0","0","NA"}' > BDR.oii

# Make opi file from R output: 
sed -z 's/ \n/ NA\n/g' BDR_FinalProbes.txt | sed 's/chr//g' | awk '{print $2,$1,$3,$5,$4}' > BDR.opi

# Add "IID" to beginning of first line in beta matrix from R: 
sed -i '1s/^/IID /' BDR_betas.txt 
 
# Make bod-file:
~/osca_Linux --tefile BDR_betas.txt --methylation-beta --make-bod --no-fid --out BDR


# Make discrete covar file for sex and plate: 
tail -n +2 BDR_pheno.txt | awk '{print $1,$1,$5,$7}' > BDR.covar

# Make quantitative covar file for age, PMD and estimated cell proportions: 
tail -n +2 BDR_pheno.txt | awk '{print $1,$1,$4,$9,$14}' > BDR.qcovar

# Make pheno file:
tail -n +2 BDR_pheno.txt | awk '{print $1,$1,$8}' > BDR.phen


# Adjust betas for covariates: 
~/osca_Linux --befile BDR --upper-beta 0.975 --lower-beta 0.025 --sd-min 0.02 --covar BDR.covar --qcovar BDR.qcovar --adj-probe --make-bod --out BDR_adj



## Run EWAS analyses ##


# MOA analysis of Braak aSyn stage:
~/osca_Linux --moa --befile BDR_adj --pheno BDR.phen  --out BDR_adj_moa

