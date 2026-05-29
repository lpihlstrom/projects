# A genome wide association study of Young Onset Parkinson’s Disease in European ancestry

`GP2 ❤️ Open Science 😍`

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20210180.svg)](https://doi.org/10.5281/zenodo.20210180)

**Last updated:** 07.05.2026  
**Contact info:** ingebhli@uio.no

## Summary

This repository contains key code associated with the manuscript titled **"A genome wide association study of Young Onset Parkinson’s Disease in European ancestry"**.

This is a genome wide association study (GWAS) on a population of young onset Parkinson's Disease cases (YOPD) and controls using both an additive and a recessive GWAS-model, in addition to meta-analysis. 

Please read the methods section of the paper for detailed information regarding upstream quality control and additional details. 

Summary statistics generated from the study will be available for use on the GP2 approved analysis platform (Fill in correct information).


## Data statement
* Data used in the preparation of this article were obtained from the Global Parkinson’s Genetics Program (GP2; https://gp2.org). All GP2 data are hosted in collaboration with the Accelerating Medicines Partnership in Parkinson’s disease, and are available via application on the website (https://amp-pd.org/register-for-amp-pd). For up-to-date information on GP2 data acquisition, access, and policies, visit https://gp2.org/. Tier 1 data can be accessed by completing a form on the Accelerating Medicines Partnership in Parkinson’s Disease (AMP®-PD) website (https://amp-pd.org/register-for-amp-pd). Tier 2 data access requires approval and a Data Use Agreement signed by your institution. In this analysis we used Tier 2 GP2 Release 10 data  ([10.5281/zenodo.15748014](https://doi.org/10.5281/zenodo.15748014)).
* Data from the NeuroGenetics Research Consortium (NGRC) are available from the database of Genotypes and Phenotypes (dbGaP, https://dbgap.ncbi.nlm.nih.gov/),  study accession phs00196.v2.p1. 
* Data from two substudies in the International Parkinson’s Disease Genomics Consortium (IPDGC) are available from dbGaP (study accessions phs000918.v1.p1, phs000089.v3.p2). For other IPDGC subcohorts, raw data have not been made public due to restrictions in privacy regulations, informed consent or ethical approvals. 


### Helpful Links

- [GP2 Website](https://gp2.org/)
  - [GP2 Cohort Dashboard](https://gp2.org/cohort-dashboard-advanced/)
- [Introduction to GP2](https://movementdisorders.onlinelibrary.wiley.com/doi/10.1002/mds.28494)
  - [Other GP2 Manuscripts (PubMed)](https://pubmed.ncbi.nlm.nih.gov/?term=%22global+parkinson%27s+genetics+program%22)

## Repository Orientation
- The `analysis/` directory includes all analyses discussed in the manuscript.

```bash
THIS_REPO/
analysis/
    ├── 00_prepping_data_gp2.html/
    ├── 01_perGroup_GWASes_and_liftover_HG19_to_HG38/
    │   ├── 00_additive_gwas_GP2.html
    │   ├── 01_recessive_gwas_GP2.html
    │   ├── 02_make_gwas_sumstats_additive_model_gp2.html
    │   ├── 03_recessive_gwas_NGRC.Rmd
    │   ├── 04_IPDGC_single_gwas_add_and_rec.txt
    │   └── 05_liftover_single_GWAS_from_hg19_to_hg38.Rmd
    ├── 02_prep_data_METAL/
    │   ├── 00_formatMETAL.Rmd
    │   └── 01_remove_extreme_outliers_recessive_single_gwas.Rmd
    ├── 03_additive_meta_METAL.txt/
    ├── 04_investigate_hits_additive_meta.rmd/
    ├── 05_reformat_meta_for_FUMA/
    │   ├── 00_reformat_meta_for_FUMA.R
    │   ├── 01_liftover_metaGWAS_from_hg38_to_hg19.R
    │   └── 02_liftover_FUMA_risk_loci_from_hg19_to_hg38.R
    ├── 06_comparison_tables.html/
    ├── 07_conditional_analysis_SNCA_locus.html/
    ├── 08_GCTA_GREML.html/
    ├── 09_figures/
    │   ├── 00_meta_figures.Rmd
    │   └── 01_PRS_figures.html
    ├── 10_PRS.html/
    └── 11_ROH_and_WGS_analyses/
        ├── 00_recessive_gwas_wgs_roh_vep.sh
        └── 01_filter_vep.py
```

## Analysis Notebooks
### Languages: Python, bash, and R
| Directory | Notebooks   | Description | 
|-----------|----------------|--------|
|`analyses/`| `00_prepping_data_gp2.html`         | Clean the Gp2 v10 clinical data, add principal components and make the covariate file before running GWAS.|
|`01_perGroup_GWASes_and_liftover_HG19_to_HG38/`| `00_additive_gwas_GP2.html`         |Run additive-model GWAS on GP2 data.|
|`01_perGroup_GWASes_and_liftover_HG19_to_HG38/`| `01_recessive_gwas_GP2.html`         |Run recessive-model GWAS on GP2 data.|
|`01_perGroup_GWASes_and_liftover_HG19_to_HG38/`| `02_make_gwas_sumstats_additive_model_gp2.html`         |Generate merged summary statistics file for GP2 additive GWAS, as the results are originally stored per chromosome.|
|`01_perGroup_GWASes_and_liftover_HG19_to_HG38/`| `03_recessive_gwas_NGRC.Rmd`         |Perform recessive-model GWAS for the NGRC dataset. The additive NGRC GWAS used the same pipeline with the additive model substituted for the recessive model.|
|`01_perGroup_GWASes_and_liftover_HG19_to_HG38/`| `04_IPDGC_single_gwas_add_and_rec.txt`         |Run additive and recessive GWAS for IPDGC-cohorts.|
|`01_perGroup_GWASes_and_liftover_HG19_to_HG38/`| `05_liftover_single_GWAS_from_hg19_to_hg38.Rmd`         |Lift over single-study GWAS summary stats from hg19 to hg38.|
|`02_formatMETAL/`| `00_formatMETAL.Rmd`         |Reformat summary stats into METAL-compatible input|
|`02_formatMETAL/`| `01_remove_extreme_outliers_recessive_single_gwas.Rmd`         |Filter extreme effect sizes / outliers from the recessive GWAS before meta-analysis with METAL.|
|`anayses/`| `03_additive_meta_METAL.txt`         |Run METAL additive-model meta-analysis. The pipeline was similar for recessive-model meta-analysis.|
|`anayses/`| `04_investigate_hits_additive_meta.Rmd`         |Display genome wide significant SNPs in the additive meta-analysis.|
|`05_reformat_meta_for_FUMA/`| `00_reformat_meta_for_FUMA.R`         |Reformat meta-analysis output to FUMA input specifications.|
|`05_reformat_meta_for_FUMA/`| `01_liftover_metaGWAS_from_hg38_to_hg19.R`         |Lift over meta-GWAS summary stats from hg38 to hg19 for FUMA.|
|`05_reformat_meta_for_FUMA/`| `02_liftover_FUMA_risk_loci_from_hg19_to_hg38.R`         |Convert FUMA risk loci coordinates back from hg19 to hg38.|
|`analysis/`| `06_comparison_tables.html`         |Create tables comparing results across meta-analyses (YOPD meta-analysis vs Leonard 2025 and Nalls et al 2019).|
|`analysis/`| `07_conditional_analysis_SNCA_locus.html`         |Run locus-specific conditional analysis in PLINK 2 on GP2 data.|
|`analysis/`| `08_GCTA_GREML.html`         |Estimate h2 using GCTA-GREML on GP2 data. Refdata file is GP2 v10 European imputed and previously QC-ed data, with maf 0.01, no duplicates and chr:bp SNP ids.|
|`09_figures/`| `00_meta_figures.html`         |Script for generation of figures from the meta-analyses (Manhattan plot, QQ-plot, Miami plot, Beta-beta-plots, forest plots.|
|`09_figures/`| `01_PRS_figures.Rmd`         |Script for generation of figures from the PRS data (violin plots and box plots).|
|`analysis/`| `10_PRS.html`         |Calculate polygenic risk scores (PRS) using PRSice-2.|
|`11_ROH_and_WGS_analyses/`| `00_recessive_gwas_wgs_roh_vep.sh`         |Calculates runs of homozygosity (ROHs) around recessive model GWAS hits and outputs potential causative variants from whole genome sequencing (WGS) subset data.|
|`11_ROH_and_WGS_analyses/`| `01_filter_vep.py`         |An auxiliary script to filter a Ensembl Variant Effect Predictor (VEP) annotated VCF for potentially deleterious variants using custom criteria.|

# Software

| Software | Version(s) | Resource URL | RRID | Notes |
|---|---|---|---|---|
| R (RStudio) | 4.4.1 (2024-06-14) and 4.5.2 (2025-10-31)| https://cran.r-project.org/ | RRID:SCR_001905 | Statistical testing and data processing |
| PLINK 2 | plink2_linux_amd_avx2_20250609 and v2.0.0-a.6.9LM 64-bit Intel (29 Jan 2025) | https://www.cog-genomics.org/plink/2.0/ | RRID:SCR_001757 | Genotype QC and association tests |
| GCTA | gcta-1.95.0 | http://cnsgenomics.com/software/gcta/ | N/A | GREML |
| METAL | 2020-05-05 | http://csg.sph.umich.edu/abecasis/Metal/ | RRID:SCR_002013 | Meta-analysis of GWAS summary stats |
| PRSice| PRSice-2 (version:2.3.5) | https://www.prsice.info/ | RRID:SCR_017057| Polygenic risk score calculation|
| Python | 3.10.17 and 3.9.19 | https://www.python.org/ | RRID:SCR_008394 | Scripting and data processing |
| Adobe Photoshop | CS6 Extended 13.0 (2012) | https://helpx.adobe.com/photoshop/ | RRID:SCR_014199 | Image editing for figures |
| liftOver (UCSC) | chain files | https://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/ (hg19→hg38) https://hgdownload.cse.ucsc.edu/goldenPath/hg38/liftOver/ (hg38→hg19) | RRID:SCR_018160 | UCSC chain files for coordinate conversion |
| bedtools | v2.31.1 | https://bedtools.readthedocs.io/ | RRID:SCR_006646 | Genomic interval operations |
| PLINK 1.9 | v1.9.0-b.7.7 64-bit (22 Oct 2024) | https://www.cog-genomics.org/plink/1.9/ | RRID:SCR_001757 | Genotype QC and association tests |
| bcftools | 1.23 | https://samtools.github.io/bcftools/ | RRID:SCR_015687 | VCF/BCF processing |
| VEP (Ensembl Variant Effect Predictor) | 115.2 | https://www.ensembl.org/info/docs/tools/vep/index.html | RRID:SCR_007931 | Variant annotation |
| cyvcf2 (Python package) | 0.31.4 | https://brentp.github.io/cyvcf2/ | RRID:SCR_024000 | Fast VCF parsing and filtering |
