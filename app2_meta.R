# run_app2.php writes species, tissue (if sample), network_file (if custom), and sessionID variables here

.libPaths(c("/home/www/abhatta3-webserver/Rnew2/lib64/R/library/",.libPaths()))
library(igraph)
library(readr)
library(ggplot2)
library(matrixStats)
source("/home/www/abhatta3-webserver/include/pharmomics/PharmOmics_app2_utils.R")
options(stringsAsFactors = FALSE)
load("/home/www/abhatta3-webserver/include/pharmomics/Final_KDA_frame_meta_Jan2022.rda")
load("/home/www/abhatta3-webserver/include/pharmomics/Gene_Symbols.rda")
ADR_scores <- readRDS("/home/www/abhatta3-webserver/include/pharmomics/ADR_scores.rds")
load("/home/www/abhatta3-webserver/include/pharmomics/hepatotox_genes.rda")

app2_resource = "/home/www/abhatta3-webserver/Data/Pipeline/Resources/shinyapp2_temp/"
pharm_resource = "/home/www/abhatta3-webserver/Data/Pipeline/Resources/app2/"

if(!exists("Genes")){
  Genes <- read.delim(paste0(app2_resource, sessionID,"_genes.txt"))
}
Genes <- unique(Genes$GENE)
drugNetworksDir <- paste0(app2_resource,sessionID,"_drug_networks/")
dir.create(path = drugNetworksDir)

Mouse_symbols2 <- read_delim(paste0(pharm_resource,"Human_Mouse_Symbols_Majority_Ensembl_HGNC_supported.txt"), 
                             "\t", escape_double = FALSE, trim_ws = TRUE)
if(length(Genes) > 500){Genes <- Genes[1:500]}
if(species %in% "Human"){
  if(mean(Genes %in% Mouse_symbols2$human_symbol) < 0.05){Genes <- unique(Mouse_symbols2$human_symbol[Mouse_symbols2$mouse_symbol %in% Genes])}
}else if(species %in% "Mouse"){
  if(mean(Genes %in% Mouse_symbols2$mouse_symbol) < 0.05){Genes <- unique(Mouse_symbols2$mouse_symbol[Mouse_symbols2$human_symbol %in% Genes])}
}else{stop("Something unexpected went wrong, please contact developer")}
if(length(Genes) < 1){stop("Error- please check your genes are either human or mouse genes")}

if(species %in% "Human"){
  geneSymbols <- human_genes
}else if(species %in% "Mouse"){
  geneSymbols <- mouse_genes
}else{stop("Something unexpected went wrong, please contact developer")}

if(!exists("network_file")){ # user selected sample
  load(paste0(pharm_resource,"Network_repos_pkg_", species,"_",tissue,".rda"))
} else{ # if user uploaded own network
  
  networks <- read_delim(network_file, "\t", escape_double = FALSE, trim_ws = TRUE)
  
  net <- graph_from_data_frame(d=networks,  directed=T)
  
  if(length(V(net)[V(net)$name %in% geneSymbols])>12500){
    # submit to hoffman2
    resource_dir <- "/u/home/m/mergeome/PharmOmics_resource/"
    data_dir <- "/u/scratch/m/mergeome/app2seg/"
    login = "sshpass -p \"pharmomics129@\" "
    save(net, geneSymbols, file = paste0(sessionID,"_Meta_App2_Network.rda"))
    cpnet = paste0("scp ",paste0(sessionID,"_Meta_App2_Network.rda")," mergeome@192.154.2.204:",data_dir)
    system(paste0(login,cpnet))
    
    mem = "8G"
    
    if(length(V(net)[V(net)$name %in% geneSymbols])<=25000){
      mem = "16G"
    } else if(length(V(net)[V(net)$name %in% geneSymbols])<=30000){
      mem = "24G"
    } else{ # user uploaded network should not exceed 35000
      mem = "32G"
    }
    cmds = paste0("'source /etc/profile;module load R;cd ",data_dir,
                  ";qsub -cwd -V -m bea -l h_data=",mem,",h_rt=1:00:00,highp ", resource_dir,"buildDistMat.sh ", sessionID," meta'")
    system(paste0(login,"ssh mergeome@192.154.2.201 ",cmds))
    # wait for dist mat to build and cp over
    Sys.sleep(120)
    files <- list.files(app2_resource)
    clock = 0
    while(!(paste0(sessionID,"_Meta_App2_DistMat.rda") %in% files)){
      Sys.sleep(20)
      clock = clock + 20
      files <- list.files(app2_resource)
      if(clock==3600){
        cat("Distance matrix building did not finish\n")
        break()
      }
    }
    load(paste0(app2_resource,sessionID,"_Meta_App2_DistMat.rda"))
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

if(species=="Human"){
  finalKDAframe <- finalKDAframe[finalKDAframe$Species=="Homo sapiens",]
} else{
  finalKDAframe <- finalKDAframe[finalKDAframe$Species=="Mus musculus" |
                                   finalKDAframe$Species=="Rattus norvegicus",]
}

tableresult <- getResult(diseasegenes = diseasegenes, full_genes = Genes, net = net, 
                         DOBtable = DOBtable, finalKDAframe = finalKDAframe, 
                         alldistancetable = alldistancetable, species = species, 
                         drugNetworksDir = drugNetworksDir)
if(species=="Human"){
  hepatotox_genes <- unique(Mouse_symbols2$human_symbol[Mouse_symbols2$mouse_symbol %in% hepatotox_genes])
}
hepatotox_jaccard <- Jaccard2(hepatotox_genes,Genes)
scores <- ADR_scores$hepatotox_complete_signature_Jaccard
scores <- c(hepatotox_jaccard,scores)
ranks <- rank(scores)/length(scores)
pvals <- pnorm(scores, lower.tail = F)

hepatotoxicity_table <- data.frame("Adverse drug reaction"="CTD_Chemical induced liver injury full signature",
                                   "Jaccard score"=hepatotox_jaccard,
                                   "Jaccard score rank"=ranks[1],
                                   "Jaccard p value"=pvals[1],
                                   "ADR genes and input gene overlap"=do.call("paste",c(intersect(hepatotox_genes, Genes),
                                                                                        list("sep"=", "))))

if(nrow(tableresult)==0 | !exists("tableresult")){
  tableresult = data.frame("Result"="Analysis gave no results.")
}

write.table(tableresult, paste0("/home/www/abhatta3-webserver/Data/Pipeline/Results/shinyapp2/",sessionID, "_app2result.txt"), 
            row.names=FALSE, quote = FALSE, sep ="\t")
write.table(hepatotoxicity_table, paste0("/home/www/abhatta3-webserver/Data/Pipeline/Results/shinyapp2/",sessionID, "_app2result_hepatotox.txt"), 
            row.names=FALSE, quote = FALSE, sep ="\t")
