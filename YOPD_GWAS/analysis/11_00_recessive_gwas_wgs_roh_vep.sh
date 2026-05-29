##Project: GP2_EUR_YOPD_metaGWAS##

##Author: MV##

#!/bin/bash

set -xeuo pipefail

# Get genome-wide significant recessive GWAS hits and make BED files
# used for ROH calling around imputed genotypes and for WGS data extraction

awk 'NR == 1 || $10 < 5e-8' \
	METAANALYSIS_YOPD_RECESSIVE.tbl \
	>recessive_gwsig_snps.txt

tail -n+2 recessive_gwsig_snps.txt | cut -f 1 | sed 's/:/ /' |
	awk 'BEGIN{OFS="\t"} {print $1, $2-1, $2}' |
	sort -V >recessive_gwsig_snps.bed

tail -n+2 recessive_gwsig_snps.txt | cut -f 1 | sed 's/:/ /' |
	awk 'BEGIN{OFS="\t"} {print $1, $2-1-5000000, $2+5000000}' |
	sort -V | bedtools merge >recessive_gwsig_snps_5mbflank.bed

tail -n+2 recessive_gwsig_snps.txt | cut -f 1 | sed 's/:/ /' |
	awk 'BEGIN{OFS="\t"} {print $1, $2-1-1000000, $2+1000000}' |
	sort -V | bedtools merge >recessive_gwsig_snps_1mbflank.bed

tail -n+2 recessive_gwsig_snps.txt | cut -f 1 | sed 's/:/ /' |
	awk 'BEGIN{OFS="\t"} {print $1, $2-1-500000, $2+500000}' |
	sort -V | bedtools merge >recessive_gwsig_snps_500kbflank.bed

tail -n+2 recessive_gwsig_snps.txt | cut -f 1 | sed 's/:/ /' |
	awk 'BEGIN{OFS="\t"} {print $1, $2-1-250000, $2+250000}' |
	sort -V | bedtools merge >recessive_gwsig_snps_250kbflank.bed

sed 's/^/chr/' recessive_gwsig_snps_500kbflank.bed \
	>recessive_gwsig_snps_500kbflank_chrprefix.bed

# Get the effect alleles of GWAS hits
awk 'NR > 1 { \
    split($1, pos, ":"); \
    allele_col = ($8 > 0 ? 2 : 3); \
    OFS="\t"; print pos[1], pos[2]-1, pos[2], $1, toupper($allele_col) \
  }' recessive_gwsig_snps.txt |
	sort -V >recessive_gwsig_snps_pd_allele.bed

bedtools intersect -wa \
	-a recessive_gwsig_snps_pd_allele.bed \
	-b YOPD_rec_top_snps.bed \
	>recessive_gwsig_topsnps_pd_allele.bed

# Restrict imputed genotypes to the WGS available subset and GWAS flanking regions.

tail -n+2 v10_individuals_to_extract.txt >yopd_samples.txt
cut -f 1 recessive_gwsig_snps_5mbflank.bed | sort -Vu >flank_chroms.txt

for chrom in $(cat flank_chroms.txt); do
	plink2 \
		--pfile "chr${chrom}_EUR_release10_vwb" \
		--keep yopd_samples.txt \
		--extract bed0 recessive_gwsig_snps_5mbflank.bed \
		--mac 1 \
		--make-pgen \
		--out "recessive_gwsig_flank_subset_${chrom}"
done

sed '1d; s/^/recessive_gwsig_flank_subset_/' flank_chroms.txt >merge_list.txt

plink2 \
	--pfile recessive_gwsig_flank_subset_1 \
	--pmerge-list merge_list.txt \
	--make-pgen \
	--out recessive_gwsig_flank_subset_merged

# Keep common well-spaced SNPs for ROH calling
plink2 \
	--pfile recessive_gwsig_flank_subset_merged \
	--snps-only just-acgt \
	--maf 0.05 \
	--geno 0.01 \
	--bp-space 2000 \
	--make-bed \
	--out recessive_gwsig_flank_subset_merged_easy

plink \
	--bfile recessive_gwsig_flank_subset_merged_easy \
	--homozyg \
	--homozyg-kb 250 \
	--homozyg-het 2 \
	--out roh_recessive_gwsig_flank

gawk 'BEGIN{OFS="\t"} NR > 1 {print "chr"$4, $7-1, $8, $2}' \
	roh_recessive_gwsig_flank.hom | sort -V >rohs.bed

# For each top GWAS hit, find WGS cases/controls homozygous for the effect allele,
# then retain ROHs overlapping that lead variant in those homozygous carriers

while IFS=$'\t' read -r chrom start end name pd_allele; do
	chrom="chr${chrom}"
	echo -e "${chrom}\t${start}\t${end}" >variant_coord.bed

	bedtools intersect -a rohs.bed -b variant_coord.bed -wa >rohs_at_variant.bed

	plink2 \
		--pfile recessive_gwsig_flank_subset_merged \
		--extract bed0 variant_coord.bed \
		--make-bed \
		--out variant_genotypes

	plink2 \
		--bfile variant_genotypes \
		--export A \
		--out variant_genotypes_exp

	allele_1=$(cut -f 7 variant_genotypes_exp.raw | head -n 1 | sed 's/.*://' | cut -d _ -f 1)
	allele_2=$(cut -f 7 variant_genotypes_exp.raw | head -n 1 | sed 's/.*://' | cut -d _ -f 2)

	awk -v "a1=${allele_1}" -v "a2=${allele_2}" 'NR != 1 { \
      OFS="\t"; \
      if ($7 == 0) gt = a1""a1; \
      else if ($7 == 1) gt = a1""a2; \
      else if ($7 == 2) gt = a2""a2; \
      else gt = "X"; \
      print $2, gt \
    }' variant_genotypes_exp.raw >variant_genotypes_exp_formatted.txt

	pd_genotype="${pd_allele}${pd_allele}"

	awk -v "gt=${pd_genotype}" '$2 == gt' variant_genotypes_exp_formatted.txt |
		grep -F -f yopd_wgs_cases.txt | cut -f 1 >"hom_cases_${chrom}_${end}.txt"
	awk -v "gt=${pd_genotype}" '$2 == gt' variant_genotypes_exp_formatted.txt |
		grep -F -f yopd_wgs_controls.txt | cut -f 1 >"hom_controls_${chrom}_${end}.txt"

	grep -F -f "hom_cases_${chrom}_${end}.txt" rohs_at_variant.bed \
		>"hom_cases_rohs_${chrom}_${end}.bed"
	grep -F -f "hom_controls_${chrom}_${end}.txt" rohs_at_variant.bed \
		>"hom_controls_rohs_${chrom}_${end}.bed"
done <recessive_gwsig_topsnps_pd_allele.bed

# Extract WGS variants surrounding the GWAS hits for VEP annotation

for c in {1..22}; do
	bcftools view \
		-S yopd_wgs_cases_controls.txt \
		-R recessive_gwsig_snps_500kbflank_chrprefix.bed \
		-Oz -o "tmp.chr${c}.vcf.gz" \
		"chr${c}.vcf.gz"
done

bcftools concat -Oz -o yopd_gwsig_500kbflank.vcf.gz tmp.chr{1..22}.vcf.gz
bcftools index -t yopd_gwsig_500kbflank.vcf.gz

bcftools annotate \
	--rename-chrs chr_rename.txt \
	-Oz -o yopd_gwsig_500kbflank_chrfix.vcf.gz \
	yopd_gwsig_500kbflank.vcf.gz

# Annotate with VEP
# VEP note: assumes Ensembl VEP GRCh38 cache, VEP_plugins, CADD,
# SpliceAI masked hg38 SNV/indel resources, and dbNSFP installed

vep \
	-i yopd_gwsig_500kbflank_chrfix.vcf.gz \
	-o yopd_gwsig_500kbflank_chrfix_anno.vep.vcf.gz \
	--vcf --compress_output bgzip --force_overwrite \
	--cache --offline --dir_cache /home/rstudio/vep \
	--species homo_sapiens --assembly GRCh38 \
	--dir_plugins /home/rstudio/vep_plugins/VEP_plugins \
	--pick --pick_order mane_select,mane_plus_clinical,canonical,appris,tsl,biotype,ccds,rank,length \
	--mane --symbol --hgvs --variant_class --numbers --canonical --tsl --appris --biotype \
	--regulatory --distance 5000 \
	--af_gnomadg --af_gnomade --max_af \
	--plugin CADD,snv=/home/rstudio/CADD/whole_genome_SNVs.tsv.gz,indels=/home/rstudio/CADD/gnomad.genomes.r4.0.indel.tsv.gz \
	--plugin SpliceAI,snv=/home/rstudio/spliceai/spliceai_scores.masked.snv.hg38.vcf.gz,indel=/home/rstudio/spliceai/spliceai_scores.masked.indel.hg38.vcf.gz \
	--plugin dbNSFP,/home/rstudio/dbnsfp/dbNSFP5.3.1a_grch38.gz,REVEL_score,MetaSVM_pred,MetaLR_pred,ClinPred_score,PrimateAI_score,MutationTaster_pred,PROVEAN_pred \
	--fork 8 \
	--stats_text

bcftools norm \
	-m -any \
	-Oz -o yopd_gwsig_500kbflank_chrfix_anno_split.vep.vcf.gz \
	yopd_gwsig_500kbflank_chrfix_anno.vep.vcf.gz

tabix -p vcf yopd_gwsig_500kbflank_chrfix_anno_split.vep.vcf.gz

# Filter variants of interest using annotations
python3 filter_vep.py \
	--vcf yopd_gwsig_500kbflank_chrfix_anno_split.vep.vcf.gz \
	--cases yopd_wgs_cases.txt \
	--controls yopd_wgs_controls.txt \
	--out-tsv filter_candidates.tsv \
	--min-hom-cases 1 \
	--max-hom-controls 0 \
	--min-cadd 20 \
	--min-spliceai 0.2 \
	--min-revel 0.5 |
	bgzip -c >yopd_gwsig_500kbflank_chrfix_anno_split_filtered.vep.vcf.gz

tabix -p vcf yopd_gwsig_500kbflank_chrfix_anno_split_filtered.vep.vcf.gz
