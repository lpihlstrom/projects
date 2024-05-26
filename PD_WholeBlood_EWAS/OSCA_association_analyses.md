## Code for EWAS analyses using the OSCA software package



### Oslo dataset

Prepare data files. Osca uses a binary file format, BOD. This is generated from three text files: the beta matrix, a probe info file and a phenotype file.
The following files have been prepared from clinical/demographic data: 
WholeBlood_all.covar
WholeBlood_all.qcovar
WholeBlood_all.phen

```
# Make oii file from R output: 
tail -n +2 WholeBlood_Pheno_allGroups.txt | awk '{print $2,$2,"0","0","NA"}' > WholeBlood_all.oii

# Make opi file from R output: 
sed -z 's/ \n/ NA\n/g' WholeBloodFinalProbes.txt | sed 's/chr//g' | awk '{print $2,$1,$3,$5,$4}' > WholeBlood_all.opi

# Add "IID" to beginning of first line in beta matrix from R: 
sed -i '1s/^/IID /' WholeBloodfinalBetas.txt 
 
# Make bod-file:
~/osca_Linux --tefile WholeBloodfinalBetas.txt --methylation-beta --make-bod --no-fid --out WholeBlood_all

```

Filter out constituitively methylated/unmethylated probes and adjust for covariates:

```
~/osca_Linux --befile WholeBlood_all --upper-beta 0.975 --lower-beta 0.025 --covar WholeBlood_all.covar --qcovar WholeBlood_all.qcovar --adj-probe --make-bod --out WholeBlood_all_adj
```

Run MOA EWAS analysis of PD cases vs controls - beta values:

```
~/osca_Linux --moa --befile WholeBlood_all_adj --pheno WholeBlood_all_numPheno.phen --out Oslo_WB_PDvsCO_betas
```

Generate Methylation Profile Score using results from Vallerga et al. 2020 as weights

```
~/osca_Linux --befile WholeBlood_all_adj --score Vallerga2020_NCOMMS_MWAS-meta-analysis-MOA.mlma 2 6 --score-has-header --out OsloByVallerga_methscore
```

Estimate variance explained using OREML:

```
~/osca_Linux --befile WholeBlood_all --make-orm --out WholeBlood_PDvsCO_all_orm
~/osca_Linux --reml --orm WholeBlood_PDvsCO_all_orm --covar WholeBlood_all.covar --qcovar WholeBlood_all.qcovar --pheno WholeBlood_all_numPheno.phen --out WholeBlood_PDvsCO_all_reml
```


### PPMI dataset

Prepare data files.
The following files have been prepared from clinical/demographic data: 
ppmi_sex.covar
ppmi_age_cellprop.qcovar
ppmi.phen

```
# Make oii file from R output: 
cut -d' ' -f6 ppmiSelectedPheno.txt | tail -n +2 | awk '{print $1,$1,"0","0","NA"}' > ppmi.oii

# Add "IID" to beginning of first line in beta matrix from R: 
sed -i '1s/^/IID /' ppmiBetas.txt 

# Make bod-file:
~/osca_Linux --tefile ppmiBetas.txt --methylation-beta --make-bod --no-fid --out ppmi

# Make opi file from R output: 
sed -z 's/ \n/ NA\n/g' ppmiFinalProbes.txt | sed 's/chr//g' | awk '{print $2,$1,$3,$5,$4}' > ppmi.opi

```

Filter out constituitively methylated/unmethylated probes and adjust for covariates:

```
~/osca_Linux --befile ppmi --upper-beta 0.975 --lower-beta 0.025 --covar ppmi_sex.covar --qcovar ppmi_age_cellprop.qcovar --adj-probe --make-bod --out ppmi_adj
```
Run MOA EWAS analysis of PD cases vs controls - beta values:

```
~/osca_Linux --moa --befile ppmi_adj --pheno ppmi.phen  --out PPMI_WB_PDvsCO_betas
```

Generate Methylation Profile Score using results from Vallerga et al. 2020 as weights

```
~/osca_Linux --befile ppmi_adj --score ~/Methylation/WholeBlood/Vallerga2020_NCOMMS_MWAS-meta-analysis-MOA.mlma 2 6 --score-has-header --out ppmiByVallerga_methscore
```




