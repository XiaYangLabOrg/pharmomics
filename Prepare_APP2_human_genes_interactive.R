# =========================================================
# Convert WKDA ENSRNOG IDs to Human Gene Symbols for APP2
# =========================================================
#
# INPUT FILES:
#   wKDA.hubs.txt
#   wKDA.tophits.txt
#
# OUTPUT FILES:
#   APP2_human_genes.txt
#   WKDA_gene_conversion_table.txt
#
# PURPOSE:
#   1. Read WKDA outputs
#   2. Extract ENSRNOG rat Ensembl IDs
#   3. Convert rat IDs -> human ortholog symbols
#   4. Create APP2-ready gene list
#
# =========================================================

# -----------------------------
# Load libraries
# -----------------------------
library(biomaRt)
library(dplyr)
library(readr)

# -----------------------------
# Define input files
# -----------------------------
wkda_hubs_file <- "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/SCING/SCING_output/MSMC_OPCs/kda/wKDA.hubs.txt"
wkda_tophits_file <- "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/SCING/SCING_output/MSMC_OPCs/kda/wKDA.tophits.txt"

# -----------------------------
# Read WKDA files
# -----------------------------
hubs <- read.delim(
  wkda_hubs_file,
  header = TRUE,
  stringsAsFactors = FALSE
)

tophits <- read.delim(
  wkda_tophits_file,
  header = TRUE,
  stringsAsFactors = FALSE
)

# -----------------------------
# Extract ENSRNOG IDs
# -----------------------------
# hubs file
hub_genes <- unique(hubs$HUB)

# tophits file
tophit_genes <- unique(tophits$NODE)

# combine
all_genes <- unique(c(hub_genes, tophit_genes))

# remove missing values
all_genes <- all_genes[!is.na(all_genes)]

# keep only ENSRNOG IDs
all_genes <- all_genes[
  grepl("^ENSRNOG", all_genes)
]

cat("Total ENSRNOG genes found:", length(all_genes), "\n")

# -----------------------------
# Connect to Ensembl
# -----------------------------

rat <- useEnsembl(
  biomart = "genes",
  dataset = "rnorvegicus_gene_ensembl",
  version = 105
)

# -----------------------------
# Convert rat ENSRNOG IDs
# to human gene symbols
# -----------------------------
conversion_table <- getBM(
  attributes = c(
    "ensembl_gene_id",
    "external_gene_name",
    "hsapiens_homolog_associated_gene_name"
  ),
  filters = "ensembl_gene_id",
  values = all_genes,
  mart = rat
)

# -----------------------------
# Clean conversion table
# -----------------------------
colnames(conversion_table) <- c(
  "Rat_Ensembl_ID",
  "Rat_Symbol",
  "Human_Symbol"
)

# remove empty human symbols
conversion_table <- conversion_table[
  conversion_table$Human_Symbol != "",
]

# remove duplicates
conversion_table <- conversion_table %>%
  distinct(Human_Symbol, .keep_all = TRUE)

# -----------------------------
# Create APP2 gene list
# -----------------------------
app2_genes <- unique(
  conversion_table$Human_Symbol
)

# remove blanks
app2_genes <- app2_genes[
  app2_genes != ""
]

# -----------------------------
# Save conversion table
# -----------------------------
write.table(
  conversion_table,
  file = "/u/scratch/v/vturnbil/Pharmomics/WKDA_gene_conversion_table.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# -----------------------------
# Save APP2 gene file
# -----------------------------
app2_output <- data.frame(
  GENE = app2_genes
)

write.table(
  app2_output,
  file = "/u/scratch/v/vturnbil/Pharmomics/APP2_human_genes.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

