library(dplyr)
library(readr)

# =========================================================
# INPUT FILES
# =========================================================

network_file <- "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/SCING/SCING_output/MSMC_OPCs/network_mergeomics.txt"

conversion_file <- "/u/scratch/v/vturnbil/Pharmomics/WKDA_gene_conversion_table.txt"

# =========================================================
# READ NETWORK
# =========================================================

network <- read.delim(
  network_file,
  header = TRUE,
  stringsAsFactors = FALSE
)

cat("Original edges:", nrow(network), "\n")

# =========================================================
# READ CONVERSION TABLE
# =========================================================

conversion <- read.delim(
  conversion_file,
  header = TRUE,
  stringsAsFactors = FALSE
)

# =========================================================
# CREATE LOOKUP TABLE
# =========================================================

lookup <- conversion[, c(
  "Rat_Ensembl_ID",
  "Human_Symbol"
)]

lookup <- unique(lookup)

# =========================================================
# MAP TAIL -> HUMAN
# =========================================================

network2 <- network %>%
  left_join(
    lookup,
    by = c("TAIL" = "Rat_Ensembl_ID")
  )

colnames(network2)[ncol(network2)] <- "TAIL_HUMAN"

# =========================================================
# MAP HEAD -> HUMAN
# =========================================================

network2 <- network2 %>%
  left_join(
    lookup,
    by = c("HEAD" = "Rat_Ensembl_ID")
  )

colnames(network2)[ncol(network2)] <- "HEAD_HUMAN"

# =========================================================
# REMOVE UNMAPPED EDGES
# =========================================================

network2 <- network2 %>%
  filter(
    !is.na(TAIL_HUMAN),
    !is.na(HEAD_HUMAN),
    TAIL_HUMAN != "",
    HEAD_HUMAN != ""
  )

cat("Mapped edges:", nrow(network2), "\n")

# =========================================================
# CREATE FINAL NETWORK
# =========================================================

final_network <- data.frame(
  from = network2$TAIL_HUMAN,
  to = network2$HEAD_HUMAN,
  weight = network2$WEIGHT
)

# remove self loops
final_network <- final_network[
  final_network$from != final_network$to,
]

# remove duplicates
final_network <- distinct(final_network)

# =========================================================
# SAVE APP2 NETWORK
# =========================================================

write.table(
  final_network,
  file = "/u/scratch/v/vturnbil/Pharmomics/SCING_network_human.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

