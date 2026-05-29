# title: "Liftover FUMA results - Riskloci from HG19 to HG38"
#output: html_notebook
  
  
  ## Project: GP2_EUR_YOPD_metaGWAS
  
  ### Author: 
  # IHL
  
  ### Version: 
  # R version 4.5.2 (2025-10-31)
  
  ### Notebook Overview
  # This notebook show how the FUMA meta-analysis results in HG19 were lifted to HG38
  # The same steps were performed on the genomic risk loci files for both the additive and recessive meta analyses 
  
  
  ### Changelog
  # 05.11.25 Notebook made 
  
# Perform liftOver for the coordinates


#This script assumes that the input file is in the default FUMA-format in HG19

# Load libraries
library(rtracklayer)  # Load the rtracklayer package
library(dplyr)        # Load dplyr for data manipulation
library(data.table)

# Read the input file with the genomic risk loci
genomic_riskloci_hg19 <- fread("~/Documents/FUMA/GenomicRiskLoci.txt")

# Recode to dfv for simplicity
df <- genomic_riskloci_hg19

# Prepare the chr column and set up positions
df$chr <- paste("chr", df$chr, sep = "")
df$end <- df$pos
df$start <- df$pos
df$seqnames <- df$chr

# Step 3: Prepare Your SNP List
snp_data <- df %>% select(seqnames, start, end, rsID)  # Include rsID in the selection

# Step 4: Load the LiftOver Chain File
chain <- import.chain("/Users/Documents/USCS_chain/hg19ToHg38.over.chain")

# Step 5: Convert SNP data to GRanges
snp_ranges <- GRanges(
  seqnames = snp_data$seqnames,
  ranges = IRanges(start = snp_data$start, end = snp_data$end)
)

# Step 6: Perform LiftOver
lifted_snps <- liftOver(snp_ranges, chain)

# Step 7: Extract the results into a data frame
lifted_snps_df <- data.frame(
  seqnames_hg38 = seqnames(lifted_snps),  # Extract lifted chromosome names
  start_hg38 = start(lifted_snps),        # Extract lifted start positions
  end_hg38 = end(lifted_snps)             # Extract lifted end positions
)

# Combine old hg19 coordinates with the lifted hg38 coordinates
snp_data_hg19 <- data.frame(
  seqnames_hg19 = snp_data$seqnames,
  start_hg19 = snp_data$start,
  end_hg19 = snp_data$end,
  rsID = snp_data$rsID                         # Keep the rsID column
)

# Combine both data frames (hg19 and hg38)
final_df <- cbind(snp_data_hg19, lifted_snps_df)

# Optionally, view the combined SNPs
print(final_df)

# Step 8: Handle Unlifted SNPs
# Identify unlifted SNPs
# ---------------------------------------------------------------
# Step 6: Identify successfully lifted vs unlifted SNPs
# ---------------------------------------------------------------
is_unlifted <- lengths(lifted_snps) == 0
num_unlifted <- sum(is_unlifted)
num_total <- length(snp_ranges)

cat("LiftOver summary:\n")
cat(num_total - num_unlifted, "SNPs successfully lifted over.\n")
cat(num_unlifted, "SNPs could not be lifted.\n\n")


output <- final_df %>% select("seqnames_hg19", "start_hg19", "rsID", "seqnames_hg38.value",  "start_hg38.value")
head(output)

output_renamed <- output %>% rename(chr_hg19 = seqnames_hg19,
                                    "pos_hg19" = "start_hg19",
                                    "chr_hg38" = "seqnames_hg38.value",
                                    "pos_hg38" = "start_hg38.value")

head(output_renamed)

# Create the MarkerName_hg19 column
output_renamed <- output_renamed %>%
  mutate(
    chr_hg19 = gsub("chr", "", chr_hg19),  # Remove 'chr' from chr_hg19 values
    MarkerName_hg19 = paste0(chr_hg19, ":", pos_hg19)  # Create new column
  )
output_renamed <- output_renamed %>%
  mutate(
    chr_hg38 = gsub("chr", "", chr_hg38),  # Remove 'chr' from chr_hg38 values
    MarkerName_hg38 = paste0(chr_hg38, ":", pos_hg38)  # Create new column
  )
print(output_renamed)

# Step 10: Save the combined data frame 

# Write the data frame to a .txt file
fwrite(output_renamed, "~/Documents/FUMA/GenomicRiskLoci_lifted_from_HG19_to_HG38.txt", sep = "\t")



