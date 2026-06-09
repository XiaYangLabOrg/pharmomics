# run_app2.php writes species, tissue (if sample), network_file (if custom), and sessionID variables here

sessionID <- "WKDA_APP2"

species <- "Human"

gene_signatures <- "top500"

network_file <- "/u/scratch/v/vturnbil/Pharmomics/SCING_network_human.txt"


resource_dir <- "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/HP/PharmOmics/app2_resources/"
data_dir <- "/u/scratch/v/vturnbil/Pharmomics/"

library(igraph)
library(readr)
library(ggplot2)
library(matrixStats)

source(paste0(resource_dir,"PharmOmics_app2_seg_utils.R"))
options(stringsAsFactors = FALSE)

# this location is where the network files (visualization) is output
drugNetworksDir <- paste0(data_dir,sessionID,"_drug_networks/")
# this location is where the 32 parts will go
resultPartsDir <- paste0(data_dir,sessionID,"_res_parts/")
dir.create(path = drugNetworksDir)
dir.create(path = resultPartsDir)

Mouse_symbols2 <- read.delim(paste0(resource_dir,"Human_Mouse_Symbols_Majority_Ensembl_HGNC_supported.txt"))
#load(paste0(resource_dir,"Final_KDA_frame_Limmav2_1.2021.rda"))
load(paste0(resource_dir,"Final_KDA_frame_Limmav3_w_hepatotoxic_gene_score2022Jan.rda"))
load(paste0(resource_dir, "Gene_Symbols.rda"))
ADR_scores <- readRDS(paste0(resource_dir, "ADR_scores.rds"))

Genes <- read.delim(paste0("/u/scratch/v/vturnbil/Pharmomics/APP2_human_genes.txt"))
Genes <- Genes$GENE

if(length(Genes) > 500){Genes <- Genes[1:500]}
if(species %in% "Human"){
  if(mean(Genes %in% Mouse_symbols2$human_symbol) < 0.05){Genes <- unique(Mouse_symbols2$human_symbol[Mouse_symbols2$mouse_symbol %in% Genes])}
}else if(species %in% "Mouse"){
  if(mean(Genes %in% Mouse_symbols2$mouse_symbol) < 0.05){Genes <- unique(Mouse_symbols2$mouse_symbol[Mouse_symbols2$human_symbol %in% Genes])}
}else{stop("Something unexpected went wrong, please contact developer")}
if(length(Genes) < 20){stop("Error- please check your genes are either human or mouse genes")}

if(species %in% "Human"){
  geneSymbols <- human_genes
}else if(species %in% "Mouse"){
  geneSymbols <- mouse_genes
}else{stop("Something unexpected went wrong, please contact developer")}

if(!exists("network_file")){ # user selected sample
  load(paste0("/u/scratch/v/vturnbil/Pharmomics/SCING_network_human.txt"))
} else{ # if user uploaded own network
  # networks <- read_delim(network_file, 
  #                        "\t", escape_double = FALSE, trim_ws = TRUE)
  network_file <- unlist(strsplit(network_file,"/"))[length(unlist(strsplit(network_file,"/")))]
  networks <- read.delim(paste0(data_dir,network_file))
  
  net <- graph_from_data_frame(d=networks,  directed=T)
  
  if(length(V(net)[V(net)$name %in% geneSymbols])>12500){
    # submit to hoffman2
    save(net, geneSymbols, file = paste0(sessionID,"_App2_Network.rda"))
    if(length(V(net)[V(net)$name %in% geneSymbols])<=18000){
      mem = "8G"
    } else if(length(V(net)[V(net)$name %in% geneSymbols])<=25000){
      mem = "16G"
    } else if(length(V(net)[V(net)$name %in% geneSymbols])<=30000){
      mem = "24G"
    } else{ # user uploaded network should not exceed 35000
      mem = "32G"
    }
    system(paste0("qsub -cwd -V -m bea -l h_data=",mem,",h_rt=1:00:00,highp ", resource_dir,"buildDistMat.sh ", sessionID," seg"))
    # wait for dist mat to build and cp over
    Sys.sleep(60)
    files <- list.files(data_dir)
    clock = 0
    while(!(paste0(sessionID,"_App2_DistMat.rda") %in% files)){
      Sys.sleep(20)
      clock = clock + 20
      files <- list.files(data_dir)
      if(clock==3600){
        cat("Distance matrix building did not finish\n")
        break()
      }
    }
    Sys.sleep(60) # wait for file to be fully written
    notLoaded = TRUE
    while(notLoaded){
      tryCatch({load(paste0(data_dir,sessionID,"_App2_DistMat.rda"))})
      if(exists("alldistancetable")) notLoaded = FALSE
      Sys.sleep(30)
    }
  } else{
    alldistancetable <- distances(net, 
                                  v = V(net)[V(net)$name %in% geneSymbols], 
                                  to = V(net)[V(net)$name %in% geneSymbols], 
                                  mode = "all", weights = NULL,"unweighted")
    alldistancetable[is.infinite(alldistancetable)] <- max(alldistancetable[!is.infinite(alldistancetable)])
  }
  
  test <- degree(net, v = V(net), mode = "all",
                 loops = TRUE, normalized = FALSE)
  DOBtable <- data.frame(gene = names(test),degree = test)
  DOBtable <- DOBtable[DOBtable$gene %in% geneSymbols,]
  degtable <- as.data.frame(table(DOBtable$degree))
  lessthan50genes <- degtable[degtable$Freq < 50,]
  cutoff <- max(as.numeric(degtable$Var1[degtable$Freq > 50]))
  newbinnumber <- round(sum(degtable$Freq[as.numeric(degtable$Var1) > cutoff])/50)
  DOBtable2 <- DOBtable[DOBtable$degree > cutoff,]
  DOBtable <- DOBtable[DOBtable$degree <= cutoff,]
  DOBtable$wt <- DOBtable$degree
  DOBtable2$wt <- tryCatch({DOBtable2$wt <- as.numeric(cut_number(DOBtable2$degree,newbinnumber))+cutoff},
                           error=function(e){
                             print(paste0("Known error: ", e))
                             cat("\nChanging bin number...\n")
                             newbinnumber = newbinnumber - 1
                             weight <- as.numeric(cut_number(DOBtable2$degree,newbinnumber))+cutoff
                             return(weight)
                           })
  DOBtable <- rbind.data.frame(DOBtable,DOBtable2)
  DOBtable$gene <- as.character(DOBtable$gene)
}

diseasegenes <- Genes[Genes %in% V(net)$name]
# disease genes need to be part of DOBtable as well
diseasegenes <- diseasegenes[diseasegenes %in% geneSymbols]

save(diseasegenes, Genes, net, DOBtable, 
     alldistancetable,species, drugNetworksDir, 
     resultPartsDir, gene_signatures,
     file = paste0(data_dir,sessionID, "_Network_server_package.rda"))

