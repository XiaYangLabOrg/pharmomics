############################################################
# POSTPROCESSING ONLY
# Assumes PART_1 ... PART_14 already exist
############################################################

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

options(stringsAsFactors = FALSE)

############################################################
# SESSION VARIABLES
############################################################

sessionID <- "WKDA_APP2"

drugNetworksDir <- paste0(
  data_dir,
  sessionID,
  "_drug_networks/"
)

resultPartsDir <- paste0(
  data_dir,
  sessionID,
  "_res_parts/"
)

############################################################
# CREATE OUTPUT DIRECTORY IF NEEDED
############################################################

if(!dir.exists(drugNetworksDir)){
  
  dir.create(
    drugNetworksDir,
    recursive = TRUE
  )
}

############################################################
# LOAD PRECOMPUTED PACKAGE
############################################################

load(
  paste0(
    data_dir,
    sessionID,
    "_Network_server_package.rda"
  )
)

############################################################
# LOAD ADDITIONAL RESOURCES
############################################################

load(
  paste0(
    resource_dir,
    "Final_KDA_frame_Limmav3_w_hepatotoxic_gene_score2022Jan.rda"
  )
)

load(
  paste0(
    resource_dir,
    "Gene_Symbols.rda"
  )
)

ADR_scores <- readRDS(
  paste0(
    resource_dir,
    "ADR_scores.rds"
  )
)

############################################################
# FIND RESULT FILES
############################################################

files <- list.files(
  path = resultPartsDir,
  pattern = "PART_[0-9]+_result.rda"
)

print(files)

cat(
  "\nFound ",
  length(files),
  " result chunks\n"
)

############################################################
# VERIFY CHUNK SIZES
############################################################

chunk_sizes <- c()

for(i in files){
  
  load(file.path(resultPartsDir, i))
  
  cat(
    i,
    ":",
    nrow(tableresult),
    "rows\n"
  )
  
  chunk_sizes <- c(
    chunk_sizes,
    nrow(tableresult)
  )
}

cat(
  "\nTOTAL EXPECTED ROWS:\n"
)

print(sum(chunk_sizes))

############################################################
# CLEAR OLD OBJECTS
############################################################

if(exists("finalresults")){
  
  rm(finalresults)
}

if(exists("tableresult")){
  
  rm(tableresult)
}

gc()

############################################################
# MERGE ALL CHUNKS
############################################################

for(i in files){
  
  load(file.path(resultPartsDir, i))
  
  if(!exists("finalresults")){
    
    finalresults <- tableresult
    
  } else {
    
    finalresults <- rbind.data.frame(
      finalresults,
      tableresult
    )
  }
}

############################################################
# VERIFY MERGED SIZE
############################################################

cat(
  "\nROWS IN FINAL MERGED TABLE:\n"
)

print(nrow(finalresults))

############################################################
# FINAL RANKING
############################################################

tableresult <- finalresults

tableresult <- tableresult[
  order(tableresult$network_result),
]

tableresult$z_scorerank[
  !is.na(tableresult$network_result)
] <-
  rank(
    -tableresult$network_result[
      !is.na(tableresult$network_result)
    ]
  ) /
  length(
    tableresult$network_result[
      !is.na(tableresult$network_result)
    ]
  )

tableresult$z_scorepvalue <- pnorm(
  tableresult$network_result,
  lower.tail = TRUE
)

############################################################
# NETWORK CONSTRUCTION
############################################################

subnet <- make_ego_graph(
  net,
  order = 1,
  nodes = diseasegenes,
  mode = "all",
  mindist = 0
)

finalnet <- subnet[[1]]

for(i in 2:length(subnet)){
  
  finalnet <- igraph::union(
    finalnet,
    subnet[[i]]
  )
}

############################################################
# NETWORK NAMING
############################################################

colnames(tableresult) <- gsub(
  "all",
  "",
  colnames(tableresult)
)

tableresult$Drug_network_name <- paste(
  tableresult$drugs,
  tableresult$species,
  tableresult$tissues,
  tableresult$status,
  tableresult$time,
  tableresult$`dose with unit`,
  tableresult$dataset,
  sep = "_"
)

############################################################
# CLEAN FILENAMES
############################################################

tableresult$Drug_network_name <- gsub(
  "Rattus norvegicus",
  "rat",
  tableresult$Drug_network_name
)

tableresult$Drug_network_name <- gsub(
  "Homo sapiens",
  "human",
  tableresult$Drug_network_name
)

tableresult$Drug_network_name <- gsub(
  "limma combined",
  "",
  tableresult$Drug_network_name
)

tableresult$Drug_network_name <- gsub(
  " mg/kg",
  "mgperkg",
  tableresult$Drug_network_name
)

tableresult$Drug_network_name <- gsub(
  " ug/kg",
  "ugperkg",
  tableresult$Drug_network_name
)

tableresult$Drug_network_name <- gsub(
  " ng/mL",
  "ngpermL",
  tableresult$Drug_network_name
)

tableresult$Drug_network_name <- gsub(
  " ug/mL",
  "ugpermL",
  tableresult$Drug_network_name
)

tableresult$Drug_network_name <- gsub(
  " mg/mL",
  "mgpermL",
  tableresult$Drug_network_name
)

tableresult$Drug_network_name <- gsub(
  " ng/ml",
  "ngpermL",
  tableresult$Drug_network_name
)

tableresult$Drug_network_name <- gsub(
  "/",
  "_",
  tableresult$Drug_network_name
)

############################################################
# GENERATE TOP NETWORK PLOTS
############################################################

for(i in 1:50){
  
  drug <- tableresult$Drug_network_name[i]
  
  cat(
    "\nPlotting:",
    drug,
    "\n"
  )
  
  druggene <- unique(c(
    unlist(
      strsplit(
        tableresult$signatures_up[i],
        ","
      )
    ),
    unlist(
      strsplit(
        tableresult$signatures_down[i],
        ","
      )
    )
  ))
  
  if(length(druggene) > 0){
    
    prepareplot(
      diseasegenes = diseasegenes,
      druggene = druggene,
      net = finalnet,
      filelocation = drugNetworksDir,
      drugname = drug,
      species = species
    )
  }
}

############################################################
# FINAL OUTPUT TABLE
############################################################

write.table(
  tableresult,
  paste0(
    data_dir,
    sessionID,
    "_app2result.txt"
  ),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t"
)

cat("\nFinished postprocessing\n")

