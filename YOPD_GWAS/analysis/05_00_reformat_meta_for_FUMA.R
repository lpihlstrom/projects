### title: "Reformat the additive meta-YOPDGWAS sumstats file for FUMA ### 

#The goal of this script is to 
# 1) Format summary statistics for FUMA software 
# - Note that the formatted FUMA-file needs to be converted from HG38 to HG19 to be used in FUMA, see separate script.
# - These commands were also used on the recessive meta-YOPDGWAS sumstatsfile to prepare it for FUMA

#Load libraries
library(tidyverse)
library(data.table)
library(dplyr)
library(tidyr)


#Read in the data
META <- fread('~/Documents/METAANALYSIS_ADDITIVE_YOPDGWAS.tbl')


# Prepare the metaanalysis sumstatfile for FUMA and write to file
# Final format: CHR, BP, A1, A2, P, Beta, SE, N

#change P-value to P for simplicity
META <- META %>%
  rename(P = `P-value`)

# Separate the ID column into CHR and POS
META <- META %>%
  separate(MarkerName, into = c("CHR", "POS"), sep = ":", convert = TRUE, remove = FALSE)

FUMA <- META %>% 
  rename(SE = StdErr,
         BP = POS,
         A1=Allele1,
         A2=Allele2,
         Beta = Effect,
         N = TotalSampleSize)

FUMA <- FUMA %>% select(CHR, BP, A1, A2, P, Beta, SE, N)

# Write the FUMA data table to a specific directory

# Define the directory path
directory_path <- "/Users/Documents"  # Change this to your desired directory

# Define the file name
file_name <- "FUMA_HG38_METAANALYSIS_ADDITIVE_YOPDGWAS.tbl"  

# Combine directory path and file name
full_path <- file.path(directory_path, file_name)

# Write the data.table to the specific directory
write.table(FUMA, file = full_path, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
