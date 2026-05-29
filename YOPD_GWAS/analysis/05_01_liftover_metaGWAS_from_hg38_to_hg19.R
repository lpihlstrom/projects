# Title: "Liftover of meta-analysis results formatted for FUMA from HG38 to HG19"


## Project: GP2_EUR_YOPD_metaGWAS
  
### Author: 
# MV, adapted by IHL

### Version: 
# R version 4.4.1 (2024-06-14)

### Notebook Overview
# This notebook shows how the meta-analysis GWAS in HG38, formatted for FUMA, was lifted from HG38 to HG19 to be able to work with FUMA
# Note that all SNPs are kept in the final file, to give an overview of the SNPs that could not be lifted
# Before uploading the outfile to FUMA, it needs to be cleaned to include only the SNPs lifted HG19 and reformatted to the preferred FUMA-format as displayed in an earlier script.

### Changelog
#06.11.25 Notebook made 

# Perform the liftOver

# Assumptions made of input: no multiallelic coordinates, only autosomal data, no "chr" prefix in chromosome names

#Library these files 
library("BSgenome.Hsapiens.UCSC.hg19")
library("data.table")
library("rtracklayer")
library("GenomicRanges")
library("GenomeInfoDb")
library("Biostrings")

#
infile <- "/Users/Documents/RECESSIVE_METAANALYSIS_FUMA_HG38.tbl"
outfile <- "/Users/Documents/lifted_RECESSIVE_METAANALYSIS_FUMA_from_HG38_to_HG19_ALLROWS_KEPT.tbl"
chain_file <- ("/Users/Documents/USCS_chain/hg38ToHg19.over.chain")

# Read input
dt <- fread(infile, sep="\t", header=TRUE)

# Check column names
stopifnot(all.equal(colnames(dt), c("CHR", "BP", "A1", "A2", "P", "Beta", "SE", "N"))) 

# Create GRanges object from CHR and BP for hg38
gr38 <- GRanges(
  seqnames = paste0("chr", dt$CHR),
  ranges = IRanges(dt$BP, width=1L)
)

# LiftOver (and drop multi-mapping results)
chain <- import.chain(chain_file)
lo <- liftOver(gr38, chain)
# Creating a logical vector to identify which SNPs could be uniquely mapped
uniq <- S4Vectors::elementNROWS(lo) == 1L

# Create a data.frame to hold the output
dt_out <- dt

# Initialize columns for hg19 coordinates
dt_out$Chromosome_19 <- NA_integer_
dt_out$Position_19 <- NA_integer_
dt_out$Allele1_19 <- NA_character_
dt_out$Allele2_19 <- NA_character_

# Fill in the lifted results for those that were uniquely mapped
if (any(uniq)) {
  gr19 <- unlist(lo[uniq])
  dt_out$Chromosome_19[uniq] <- as.integer(sub("^chr", "", as.character(seqnames(gr19))))
  dt_out$Position_19[uniq] <- start(gr19)
  
  a1 <- toupper(dt$A1[uniq])
  a2 <- toupper(dt$A2[uniq])
  
  # Reference matching on hg19 to check alleles
  len1 <- nchar(a1)
  len2 <- nchar(a2)
  
  # Ensure widths of GRanges match alleles
  grA <- gr19
  width(grA) <- len1
  grB <- gr19
  width(grB) <- len2
  
  Hs <- BSgenome.Hsapiens.UCSC.hg19::Hsapiens
  refA <- as.character(getSeq(Hs, grA))
  refB <- as.character(getSeq(Hs, grB))
  
  # Compare alleles with the reference
  ok_direct <- (a1 == refA) | (a2 == refB)
  
  rc1 <- as.character(reverseComplement(DNAStringSet(a1)))
  rc2 <- as.character(reverseComplement(DNAStringSet(a2)))
  ok_rc <- (!ok_direct) & ((rc1 == refA) | (rc2 == refB))
  
  # Use reverse complement if needed
  a1_out <- ifelse(ok_rc & !ok_direct, rc1, a1)
  a2_out <- ifelse(ok_rc & !ok_direct, rc2, a2)
  
  dt_out$Allele1_19[uniq] <- a1_out
  dt_out$Allele2_19[uniq] <- a2_out
}

# Add the hg38 CHR and BP back to the output
dt_out$Chromosome_38 <- dt$CHR
dt_out$Position_38 <- dt$BP

# Rename columns for clarity
setnames(dt_out, old = c("A1", "A2"), new = c("Allele1_38", "Allele2_38"))

# Filter for unique coordinates (if needed)
# dt_out <- dt_out[, if (.N == 1L) .SD, by = .(Chromosome_19, Position_19)]

# Sort output by chr, pos
setorder(dt_out, Chromosome_19, Position_19)

# Write the output
data.table::fwrite(dt_out, outfile, sep="\t", quote=FALSE)
