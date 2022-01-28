# run_app2.php writes species, tissue (if sample), network_file (if custom), and sessionID variables here

resource_dir <- "/u/home/m/mergeome/PharmOmics_resource/"
data_dir <- "/u/scratch/m/mergeome/app2seg/"

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

Genes <- read.delim(paste0(data_dir,sessionID,"_genes.txt"))
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
  load(paste0(resource_dir,"Network_repos_pkg_", species,"_",tissue,".rda"))
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

if(gene_signatures=="top500"){
  system(paste0("qsub -cwd -V -N ", "RUN_",sessionID,"_APP2 ", "-l h_data=6G,h_rt=4:00:00,highp -q mergeome_pod.q -t 1-14:1 ",resource_dir,"wrapper_app2.sh ",sessionID))
} else {
  system(paste0("qsub -cwd -V -N ", "RUN_",sessionID,"_APP2 ", "-l h_data=6G,h_rt=6:00:00,highp -q mergeome_pod.q -t 1-14:1 ",resource_dir,"wrapper_app2.sh ",sessionID))
}
current <- system("myjobs",intern = T)
ID <- current[grep(substr(sessionID, 1, 6), current)][1]
ID <- strsplit(ID,"\\s+")[[1]][2]
save(ID, file = paste0(data_dir, sessionID, "_first_round_submission_ID.rda"))

in_myjobs = TRUE
while(in_myjobs){ # break out of loop when job ID no longer in myjobs output meaning all jobs have finished
  Sys.sleep(100)
  current <- system("myjobs",intern = T)
  # could track with job id but better to do sessionID since jobs may be resubmitted
  #if(sum(grepl(ID, current))==0) in_myjobs = FALSE
  if(sum(grepl(substr(sessionID, 1, 6), current))==0) in_myjobs = FALSE
  files <- list.files(path = resultPartsDir,
                      pattern = "PART_[0-9]+_result.rda")
  if(length(files)==14) in_myjobs = FALSE
}
cat("No more jobs running\n")

files <- list.files(path = resultPartsDir,
                    pattern = "PART_[0-9]+_result.rda")

# each job also has a check to rerun failed jobs but if the last job fails or something else goes wrong,
# this script will resubmit the remaining failed job(s)
if(length(files)<14){
  missingpart <- setdiff(1:14, as.numeric(sapply(strsplit(files,"_"), `[[`, 2)))
  # assuming that job failed, not that it ran out of time (otherwise this would not help)
  for(i in missingpart){
    if(gene_signatures=="top500"){
      system(paste0("qsub -cwd -V -N ", "RUN_",sessionID,"_APP2 -M jading@mail ", "-l h_data=6G,h_rt=4:00:00,highp ",resource_dir,"wrapper_app2_redo.sh ",i, " ",sessionID))
    } else {
      system(paste0("qsub -cwd -V -N ", "RUN_",sessionID,"_APP2 -M jading@mail ", "-l h_data=6G,h_rt=6:00:00,highp ",resource_dir,"wrapper_app2_redo.sh ",i, " ",sessionID))
    }
  }
  clock = 0
  while(length(files)<14){
    #check every 100 seconds, capped at 15 hours (job is currently just 12 h though)
    files <- list.files(path = resultPartsDir,
                        pattern = "PART_[0-9]+_result.rda")
    clock = clock + 100
    Sys.sleep(100)
    if(clock > 54000){
      files <- list.files(path = resultPartsDir,
                          pattern = "PART_[0-9]+_result.rda")
      if(length(files) < 14){cat("Something wrong in part of the computation, will go anyway")}
      break()
    }
  }
}
cat("Obtained all results for ", sessionID, "\n")

setwd(resultPartsDir)
for(i in files){
  load(i)
  if(!exists("finalresults")){finalresults <- tableresult
  }else{finalresults <- rbind.data.frame(finalresults,tableresult)}
}

tableresult <- finalresults
tableresult <- tableresult[order(tableresult$network_result),]
tableresult$z_scorerank[!is.na(tableresult$network_result)] <- rank(-tableresult$network_result[!is.na(tableresult$network_result)])/length(tableresult$network_result[!is.na(tableresult$network_result)])
tableresult$z_scorepvalue <- pnorm(tableresult$network_result, lower.tail = T)

subnet <- make_ego_graph(net, order = 1, nodes = diseasegenes, mode = "all",
                         mindist = 0)
finalnet <- subnet[[1]]
for(i in 2:length(subnet)){
  finalnet <- igraph::union(finalnet,subnet[[i]])
}

# take top 50 to create networks
colnames(tableresult) <- gsub("all","",colnames(tableresult))
tableresult$Drug_network_name <- paste(tableresult$drugs, tableresult$species,tableresult$tissues,tableresult$status,
                                       tableresult$time, tableresult$`dose with unit`, tableresult$dataset, sep = "_")
tableresult$Drug_network_name <- gsub("Rattus norvegicus", "rat", tableresult$Drug_network_name)
tableresult$Drug_network_name <- gsub("Homo sapiens", "human", tableresult$Drug_network_name)
tableresult$Drug_network_name <- gsub("limma combined", "", tableresult$Drug_network_name)
tableresult$Drug_network_name <- gsub(" mg/kg", "mgperkg", tableresult$Drug_network_name)
tableresult$Drug_network_name <- gsub(" ug/kg", "ugperkg", tableresult$Drug_network_name)
tableresult$Drug_network_name <- gsub(" ng/mL", "ngpermL", tableresult$Drug_network_name)
tableresult$Drug_network_name <- gsub(" ug/mL", "ugpermL", tableresult$Drug_network_name)
tableresult$Drug_network_name <- gsub(" mg/mL", "mgpermL", tableresult$Drug_network_name)
tableresult$Drug_network_name <- gsub(" ng/ml", "ngpermL", tableresult$Drug_network_name) # should be ng/mL 
tableresult$Drug_network_name <- gsub("  uM", "uM", tableresult$Drug_network_name)

hepatotoxicity_table <- tableresult[grepl("Hepatotoxicity|CTD_Chemical", tableresult$Drug),]
tableresult <- tableresult[!grepl("Hepatotoxicity|CTD_Chemical", tableresult$Drug),]

for(i in 1:50){
  drug <- tableresult$Drug_network_name[i]
  druggene <- unique(c(unlist(strsplit(tableresult$signatures_up[i],",")),unlist(strsplit(tableresult$signatures_down[i],","))))
  if(length(druggene)>0){
    prepareplot(diseasegenes = diseasegenes, druggene = druggene,
                net = finalnet, filelocation = drugNetworksDir, 
                drugname=drug, species=species)
  }
}
tableresult$signatures_top <- NULL
sig_adr_info <- tableresult[,grep("signatures|^ADR|^hepatotox|SIDER|DILI",colnames(tableresult), value = TRUE)]
tableresult <- tableresult[,!grepl("signatures|^ADR|^hepatotox|SIDER|DILI",colnames(tableresult))]
tableresult$KEGGname <- NULL
tableresult$KEGGlisted <- NULL
tableresult$controlsample <- NULL
tableresult$treatmentsample <- NULL

colnames(tableresult) <- c("Drug","Species","Tissue","Study","Time","Dose",
                           "Control sample size","Treatment sample size",
                           "Database","Method",
                           # added on from analysis
                           "Jaccard_score","Z_score",
                           "All_drug_gene_input_gene_overlap",
                           "Drug_gene_input_gene_overlap_in_network",
                           "Drug_genes_directly_connected_to_input_gene",
                           "Input_genes_directly_connected_to_drug_gene",
                           "Z_score_rank","Pvalue", "Drug_network_name")
tableresult <- cbind(tableresult, sig_adr_info)
# SIDER link is 83rd column
# display rank and p-value

hepatotoxicity_table <- hepatotoxicity_table[,c("drugs","Jaccard_result","network_result",
                                                "All_drug_gene_input_gene_overlap",
                                                "Network_drug_gene_input_gene_overlap",
                                                "Network_drug_gene_input_gene_first_neighbor_overlap",
                                                "Input_genes_directly_connected_to_drug_gene")]

stats <- t(apply(hepatotoxicity_table, 1, function(x){
  if(is.na(x['network_result'])){
    return('no score')
  } else {
    mod <- gsub("Hepatotoxicity_","",x['drugs'])
    if(mod=="CTD_Chemicalinducedliverinjury"){
      col <- "hepatotox_complete_signature"
    } else {
      col <- paste0("ADRscore_",mod)
    }
    scores <- ADR_scores[,col]
    scores <- as.numeric(c(x['network_result'], scores))
    ranks <- rank(-scores)/length(scores)
    pvals <- pnorm(scores, lower.tail = T)
    return(c(ranks[1],pvals[1]))
  }
}))

hepatotoxicity_table <- cbind(hepatotoxicity_table, stats)
colnames(hepatotoxicity_table) <- c("Hepatotoxicity Submodule","Jaccard_score","Z_score",
                                    "Hepatotoxicity_gene_input_gene_overlap",
                                    "Hepatotoxicity_gene_input_gene_overlap_in_network",
                                    "Hepatotoxicity_genes_directly_connected_to_input_gene",
                                    "Input_genes_directly_connected_to_Hepatotoxicity_gene", 
                                    "Z_score_rank","P_value")

# hepatotoxicity_table$`Jaccard_score_rank` <- "none"
# scores <- ADR_scores$hepatotox_complete_signature_Jaccard
# scores <- c(hepatotoxicity_table$Jaccard_score[hepatotoxicity_table$`Hepatotoxicity Submodule`=="CTD_Chemicalinducedliverinjury"],
#             scores)
# ranks <- rank(scores)/length(scores)
# pvals <- pnorm(scores, lower.tail = F)
# hepatotoxicity_table$`Jaccard_score_rank`[hepatotoxicity_table$`Hepatotoxicity Submodule`=="CTD_Chemicalinducedliverinjury"] <-
#   ranks[1]
# hepatotoxicity_table$`Jaccard_score_rank` <- substr(hepatotoxicity_table$`Jaccard_score_rank`, start = 1, stop = 5)
# hepatotoxicity_table$`Jaccard_score_p_value` <- "none"
# hepatotoxicity_table$`Jaccard_score_p_value`[hepatotoxicity_table$`Hepatotoxicity Submodule`=="CTD_Chemicalinducedliverinjury"] <-
#   pvals[1]
# hepatotoxicity_table$`Jaccard_score_p_value` <- substr(hepatotoxicity_table$`Jaccard_score_p_value`, start = 1, stop = 5)

hepatotoxicity_table <- hepatotoxicity_table[,c("Hepatotoxicity Submodule",
                                                "Z_score",
                                                "Z_score_rank",
                                                "P_value",
                                                #"Jaccard_score",
                                                #"Jaccard_score_rank",
                                                #"Jaccard_score_p_value",
                                                "Hepatotoxicity_gene_input_gene_overlap",
                                                "Hepatotoxicity_gene_input_gene_overlap_in_network",
                                                "Hepatotoxicity_genes_directly_connected_to_input_gene",
                                                "Input_genes_directly_connected_to_Hepatotoxicity_gene" 
                                                )]

hepatotoxicity_table$`Hepatotoxicity Submodule` <- gsub("Hepatotoxicity_mod","Submodule ",
                                                        hepatotoxicity_table$`Hepatotoxicity Submodule`)
hepatotoxicity_table$`Hepatotoxicity Submodule` <- gsub("Chemicalinducedliverinjury",
                                                        "Chemical induced liver injury full signature",
                                                        hepatotoxicity_table$`Hepatotoxicity Submodule`)

hepatotoxicity_table$Hepatotoxicity_gene_input_gene_overlap <- 
  gsub(",",", ",hepatotoxicity_table$Hepatotoxicity_gene_input_gene_overlap)

tableresult$hepatotox_complete_signature_rank[is.na(tableresult$hepatotox_complete_signature_rank)] <-
  "none"
tableresult$hepatotox_complete_signature_rank <- substr(tableresult$hepatotox_complete_signature_rank, start = 1, stop = 6)

tableresult$Jaccard_score <- NULL

write.table(hepatotoxicity_table, paste0(data_dir, sessionID, "_app2result_hepatotox.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(tableresult, paste0(data_dir, sessionID, "_app2result.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")

folder <- unlist(strsplit(drugNetworksDir,"/"))[length(unlist(strsplit(drugNetworksDir,"/")))]
system(paste0("~/bin/sshpass -p \"mergeomics729@\" scp -r ", drugNetworksDir," smha118@164.67.172.118:/home/www/abhatta3-webserver/Data/Pipeline/Resources/shinyapp2_temp/", folder,"/"))
system(paste0("~/bin/sshpass -p \"mergeomics729@\" scp ", data_dir, sessionID,"_app2result.txt smha118@164.67.172.118:/home/www/abhatta3-webserver/Data/Pipeline/Results/shinyapp2/"))
system(paste0("~/bin/sshpass -p \"mergeomics729@\" scp ", data_dir, sessionID,"_app2result_hepatotox.txt smha118@164.67.172.118:/home/www/abhatta3-webserver/Data/Pipeline/Results/shinyapp2/"))
system(paste0("touch ", sessionID, "_is_done"))
system(paste0("~/bin/sshpass -p \"mergeomics729@\" scp ", sessionID, "_is_done smha118@164.67.172.118:/home/www/abhatta3-webserver/Data/Pipeline/Resources/shinyapp2_temp/"))
