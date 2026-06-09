args <- commandArgs(TRUE)

if(length(args) == 0){
  
  ##########################################################
  # LOCAL DEBUG MODE
  ##########################################################
  
  sessionID <- "test_session"
  
  cat("Running in local debug mode\n")
  
  ##########################################################
  # RUN ALL 14 CHUNKS LOCALLY
  ##########################################################
  
  chunk_list <- 1:14
  
} else {
  
  ##########################################################
  # COMMAND LINE MODE
  ##########################################################
  
  chunk_list <- as.numeric(args[1])
  sessionID <- args[2]
}

resource_dir <- "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/HP/PharmOmics/app2_resources/"
data_dir <- "/u/scratch/v/vturnbil/Pharmomics/"

library(igraph)
library(readr)
library(ggplot2)
library(matrixStats)

source(
  paste0(
    resource_dir,
    "PharmOmics_app2_seg_utils.R"
  )
)

options(stringsAsFactors = FALSE, warn = -1)

Mouse_symbols2 <- read.delim(
  paste0(
    resource_dir,
    "Human_Mouse_Symbols_Majority_Ensembl_HGNC_supported.txt"
  )
)

load(
  paste0(
    data_dir,
    "WKDA_APP2_Network_server_package.rda"
  )
)

load(
  paste0(
    resource_dir,
    "Final_KDA_frame_Limmav2_12.2020.rda"
  )
)

############################################################
# PARTITION DEFINITIONS
############################################################

splitpoint <- c(
  0,2,5,9,13,18,23,
  32,39,48,57,66,77,88,96
)

############################################################
# RUN THROUGH ALL CHUNKS
############################################################

for(a in chunk_list){
  
  cat("\n=====================================\n")
  cat("RUNNING CHUNK:", a, "\n")
  cat("=====================================\n")
  
  ##########################################################
  # SELECT CURRENT PARTITION
  ##########################################################
  
  currentind <- (
    splitpoint[(a)] * 150 + 1
  ):
    min(
      nrow(finalKDAframe),
      splitpoint[a + 1] * 150
    )
  
  finalKDAframe_subset <- finalKDAframe[currentind,]
  
  ##########################################################
  # MAIN COMPUTATION
  ##########################################################
  
  tableresult <- getResult(
    diseasegenes = diseasegenes,
    full_genes = Genes,
    net = net,
    DOBtable = DOBtable,
    finalKDAframe = finalKDAframe_subset,
    alldistancetable = alldistancetable,
    species = species,
    gene_signatures = gene_signatures
  )
  
  ##########################################################
  # SAVE PARTITION RESULT
  ##########################################################
  
  save(
    tableresult,
    file = paste0(
      resultPartsDir,
      "PART_",
      a,
      "_result.rda"
    )
  )
  
  cat(
    "Finished chunk ",
    a,
    "\n"
  )
}

