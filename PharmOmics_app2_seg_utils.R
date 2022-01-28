Jaccard2 <- function(set1,set2){
  I <- length(intersect(set1,set2))
  return(I/(length(set1)+length(set2)-I))
}

concatenate=function(myvect, mysep="")
{
  if(length(myvect)==0) return(myvect)
  if(length(myvect)==1) return(myvect)
  string = ""
  for(item in myvect){
    string = paste(string, item, sep = mysep)
  }
  string = substring(string, first=(nchar(mysep)+1))
  return(string)
}

preparenetpermutation <- function(genes,net,DOBtable){
  genes <- genes[genes %in% V(net)$name]
  Diseasedeg <- degree(net, v = genes, mode = "all",
                       loops = TRUE, normalized = FALSE)
  DiseaseDEGtable <- as.data.frame(table(Diseasedeg), stringsAsFactors = FALSE)
  DiseaseDEGtable$Diseasedeg <- as.numeric(DiseaseDEGtable$Diseasedeg)
  DiseaseDEGtable$Freq <- as.numeric(DiseaseDEGtable$Freq)
  #rm(finalsample)
  for(i in 1:nrow(DiseaseDEGtable)){
    genenumber <- DiseaseDEGtable$Freq[i]
    DEGnumber <- DiseaseDEGtable$Diseasedeg[i]
    
    geneweight <- unique(DOBtable$wt[DOBtable$degree %in% DEGnumber])
    genepools <- DOBtable$gene[DOBtable$wt %in% geneweight]
    
    sampledgenes <- replicate(1000,sample(genepools,genenumber,replace = TRUE))
    if(!exists("finalsample")){
      finalsample <- sampledgenes
    }else{
      finalsample <- rbind(finalsample,sampledgenes)
    }
    
  }
  
  #finalsample <- split(finalsample, rep(1:ncol(finalsample), each = nrow(finalsample)))
  if(class(finalsample)[1] %in% "character"){
    finalsample <- as.list(finalsample)
  }else{
    finalsample <- split(finalsample, rep(1:ncol(finalsample), each = nrow(finalsample)))
  }
  return(finalsample)
}

calculatenetworkscore <- function(drugDEG, diseasegene, permutedrugDEG,permutediseaseDEG,alldistancetable, net){
  finalresult <- NULL
  allscores <- NULL
  for(k in 1:1000){
    if(sum(drugDEG %in% V(net)$name)==1){ # JD change
      allscores[k] <- mean(alldistancetable[permutedrugDEG[[k]],permutediseaseDEG[[k]]], na.rm = T)
    }else{
      #allscores[k] <- mean(apply( alldistancetable[permutedrugDEG[[k]],permutediseaseDEG[[k]]], 1, min, na.rm = T ), na.rm = T)
      allscores[k] <- mean(rowMins( alldistancetable[permutedrugDEG[[k]],permutediseaseDEG[[k]]], na.rm = T ), na.rm = T)
    }
  }
  
  if(sum(drugDEG %in% V(net)$name)==1 | sum(rownames(alldistancetable) %in% drugDEG)==1){ # JD change
    res <- mean(alldistancetable[rownames(alldistancetable) %in% drugDEG,colnames(alldistancetable) %in% diseasegene], na.rm = T)
  }else{
    res <- mean(apply( alldistancetable[rownames(alldistancetable) %in% drugDEG, colnames(alldistancetable) %in% diseasegene], 1,min, na.rm = T ), na.rm = T)
    #res <- mean(apply( alldistancetable[drugDEG,diseasegene], 1,min, na.rm = T ), na.rm = T)
  }
  #permutemean <- mean(allscores)
  #resd[j] <- res
  finalresult["network"] <- (res-mean(allscores, na.rm = T))/sd(allscores,na.rm = T)
  #finalresult["Jaccard"] <- Jaccard2(drugDEG,diseasegene)
  return(finalresult)
}

prepareplot <- function(diseasegenes,druggene,net,filelocation,drugname,species){
  #Prepare disease subnetwork
  if(species %in% "Human"){  # switch druggene to appropriate for network
    if(mean(druggene %in% Mouse_symbols2$human_symbol) < 0.3){druggene <- unique(Mouse_symbols2$human_symbol[Mouse_symbols2$mouse_symbol %in% druggene])}
  }else if(species %in% "Mouse"){
    if(mean(druggene %in% Mouse_symbols2$mouse_symbol) < 0.3){druggene <- unique(Mouse_symbols2$mouse_symbol[Mouse_symbols2$human_symbol %in% druggene])}
  }else{stop("Something unexpected went wrong, please contact developer")}
  
  neighbor1disease <- names(V(net))
  neighbor1disease <- setdiff(neighbor1disease, diseasegenes)
  
  currentgene <- druggene
  currentgene_net <- intersect(currentgene,diseasegenes)
  currentgene_netneighbor <- intersect(currentgene,neighbor1disease)
  
  V(net)$label.cex = 0.1
  V(net)$color = "#E0DEDE" 
  V(net)$size = 1
  
  pre_url = "http://chart.apis.google.com/chart?cht=p&chs=200x200&chf=bg,s,00000000&chd="
  
  #Non-overlapped disease genes
  V(net)[diseasegenes[diseasegenes %in% V(net)$name]]$color <- "#F56A79"
  V(net)[diseasegenes[diseasegenes %in% V(net)$name]]$url <- paste0(pre_url, "t:1&chco=F56A79")
  V(net)[diseasegenes[diseasegenes %in% V(net)$name]]$label.cex <- 0
  V(net)[diseasegenes[diseasegenes %in% V(net)$name]]$size <- 10
  V(net)[diseasegenes[diseasegenes %in% V(net)$name]]$status <- "Input_gene"
  
  #1-neighr overlapped disease genes
  V(net)[currentgene_netneighbor[currentgene_netneighbor %in% V(net)$name]]$color <- "#51ADCF"
  V(net)[currentgene_netneighbor[currentgene_netneighbor %in% V(net)$name]]$url <- paste0(pre_url, "t:1&chco=51ADCF")
  V(net)[currentgene_netneighbor[currentgene_netneighbor %in% V(net)$name]]$label.cex <- 0
  V(net)[currentgene_netneighbor[currentgene_netneighbor %in% V(net)$name]]$size <- 5
  V(net)[currentgene_netneighbor[currentgene_netneighbor %in% V(net)$name]]$status <- "First_neighbor_drug_gene"
  
  #direct overlapping disease genes
  V(net)[currentgene_net[currentgene_net %in% V(net)$name]]$color <- "#FFFFFF"
  V(net)[currentgene_net[currentgene_net %in% V(net)$name]]$url <- paste0(pre_url, "t:1,1&chco=F56A79|51ADCF")
  V(net)[currentgene_net[currentgene_net %in% V(net)$name]]$label.cex <- 0
  V(net)[currentgene_net[currentgene_net %in% V(net)$name]]$size = 10
  V(net)[currentgene_net[currentgene_net %in% V(net)$name]]$status = "Overlap_drug_input_gene"
  
  othernodes <- setdiff(names(V(net)),union(c(currentgene_net,currentgene_netneighbor),diseasegenes))
  finalnet2 <- delete_vertices(net, othernodes)
  
  edges <- as.data.frame(get.edgelist(finalnet2))
  colnames(edges) <- c("SOURCE","TARGET")
  edges$WEIGHT = 1
  write.table(edges, paste0(filelocation,drugname,"_cytoscape_edges.txt"), row.names = FALSE, sep = "\t", quote = FALSE)
  nodes <- data.frame("DRUG"=drugname,
                      "NODE"=vertex_attr(finalnet2)$name, 
                      "COLOR"=vertex_attr(finalnet2)$color, 
                      "URL"=vertex_attr(finalnet2)$url,"SIZE"="50","SHAPE"="Ellipse",
                      "STATUS"=vertex_attr(finalnet2)$status,
                      stringsAsFactors = FALSE)
  write.table(nodes, paste0(filelocation, drugname,"_cytoscape_nodes.txt"), row.names = FALSE, sep = "\t", quote = FALSE)
  
  #write_graph(finalnet2, filelocation, "graphml")
  
}

getResult <- function(diseasegenes, full_genes, net, DOBtable, finalKDAframe, 
                      alldistancetable,species = "Human", 
                      gene_signatures = ""){
  randomgenes_disease <- preparenetpermutation(genes = diseasegenes,
                                               net = net, DOBtable = DOBtable)
  
  resd <- NULL
  Jaccard_resd <- NULL
  direct_overlap <- NULL
  drug_input_gene_network_overlap <- NULL
  one_neighbor_overlap <- NULL
  input_neighbor_to_drug_gene <- NULL
  
  subnet <- make_ego_graph(net, order = 1, nodes = diseasegenes, mode = "all",
                           mindist = 0)
  finalnet <- subnet[[1]]
  for(i in 2:length(subnet)){
    finalnet <- igraph::union(finalnet,subnet[[i]])
  }
  
  breaks <- levels(cut(1:nrow(finalKDAframe), breaks = 100))
  upperbounds <- sapply(breaks, function(x){return(unlist(strsplit(x, ","))[2])})
  upperbounds <- as.numeric(gsub("]","",upperbounds))
  names(upperbounds) <- 1:100
  changegenes_list <- list()
  
  percent=0
  for(j in 1:nrow(finalKDAframe)){
    if(j>upperbounds[percent+1]){
      percent = percent + 1
      cat(percent,"%","\n")
    }
    
    if(gene_signatures=="top500"){
      changegenes <- unique(unlist(strsplit(finalKDAframe$allsignatures_top[j],",")))
    }
    else{
      changegenes <- unique(c(unlist(strsplit(finalKDAframe$allsignatures_up[j],",")),
                              unlist(strsplit(finalKDAframe$allsignatures_down[j],","))))
    }
    
    if(length(changegenes)==0){
      Jaccard_resd[j] <- NA
      direct_overlap[j] <- "none"
      one_neighbor_overlap[j] <- "none"
      drug_input_gene_network_overlap[j] <- "none"
      resd[j] <- NA
      input_neighbor_to_drug_gene[j] <- "none"
      next 
    }
    
    if(species %in% "Human"){  # JD change 0.05 to 0.3
      if(mean(changegenes %in% Mouse_symbols2$human_symbol) < 0.2){changegenes <- unique(Mouse_symbols2$human_symbol[Mouse_symbols2$mouse_symbol %in% changegenes])}
    }else if(species %in% "Mouse"){
      if(mean(changegenes %in% Mouse_symbols2$mouse_symbol) < 0.2){changegenes <- unique(Mouse_symbols2$mouse_symbol[Mouse_symbols2$human_symbol %in% changegenes])}
    }else{stop("Something unexpected went wrong, please contact developer")}
    
    if(length(changegenes)==0 | sum(changegenes %in% V(net)$name)==0){
      #missing[j] <- "none after species symbols"
      Jaccard_resd[j] <- NA
      direct_overlap[j] <- "none"
      one_neighbor_overlap[j] <- "none"
      drug_input_gene_network_overlap[j] <- "none"
      resd[j] <- NA
      input_neighbor_to_drug_gene[j] <- "none"
      next
    }  
    
    if(sum(changegenes %in% V(net)$name)==0){
      #missing[j] <- "none in net"
      Jaccard_resd[j] <- ifelse(length(intersect(changegenes,full_genes))>0,
                                Jaccard2(changegenes,full_genes),
                                NA)
      resd[j] <- NA
      direct_overlap[j] <- ifelse(length(intersect(changegenes,full_genes))>0,
                                  concatenate(intersect(changegenes,full_genes), ","),
                                  "none")
      one_neighbor_overlap[j] <- "none"
      drug_input_gene_network_overlap[j] <- "none"
      input_neighbor_to_drug_gene[j] <- "none"
      next # add Jan 2021 not sure why wasn't here
    } 
  
    Jaccard_resd[j] <- Jaccard2(changegenes,full_genes)
    direct_overlap[j] <- ifelse(length(intersect(changegenes,full_genes))>0,
                                concatenate(intersect(changegenes,full_genes), ","),
                                "none")
    
    # drug genes whose first neighbors are input genes
    neighbor1disease <- names(V(finalnet))
    neighbor1disease <- setdiff(neighbor1disease, diseasegenes)
    currentgene_netneighbor <- intersect(changegenes,neighbor1disease)
    one_neighbor_overlap[j] <- ifelse(length(currentgene_netneighbor)>0,
                                      concatenate(currentgene_netneighbor,","),
                                      "none")
    
    currentgene_net <- intersect(changegenes,diseasegenes)
    currentgene_net <- currentgene_net[currentgene_net %in% V(net)$name]
    drug_input_gene_network_overlap[j] <- ifelse(length(currentgene_net)>0,
                                                 concatenate(currentgene_net,","),
                                                 "none")
    
    # drug genes whose first neighbors are input genes
    edges <- as.data.frame(get.edgelist(finalnet))
    colnames(edges) <- c("SOURCE","TARGET")
    edges <- edges[(edges$SOURCE %in% changegenes) | (edges$TARGET %in% changegenes),]
    nodes <- unique(c(edges$SOURCE, edges$TARGET))
    input_neighbor_to_drug_gene[j] <- ifelse(length(diseasegenes[diseasegenes %in% nodes])>0,
                                             concatenate(diseasegenes[diseasegenes %in% nodes], ","),
                                             "none")
    

    randomgenes_drug <- preparenetpermutation(genes = changegenes,
                                              net = net,DOBtable = DOBtable)
    resd[j] <- calculatenetworkscore(drugDEG = changegenes, 
                                     diseasegene = diseasegenes, 
                                     permutedrugDEG = randomgenes_drug,
                                     permutediseaseDEG = randomgenes_disease,
                                     alldistancetable = alldistancetable, 
                                     net=net) # JD added net

  }
  
  #z_scorerank <- rank(-resd)/length(resd)
  #z_scorep_pvalue <- pnorm(resd, lower.tail = T)
  
  tableresult <- data.frame(Jaccard_result = Jaccard_resd, network_result = resd, 
                            All_drug_gene_input_gene_overlap=direct_overlap, 
                            Network_drug_gene_input_gene_overlap=drug_input_gene_network_overlap,
                            Network_drug_gene_input_gene_first_neighbor_overlap=one_neighbor_overlap,
                            Input_genes_directly_connected_to_drug_gene=input_neighbor_to_drug_gene)
  #save(tableresult, missing, file = "tableresult_debug.rda")
  tableresult <- cbind.data.frame(finalKDAframe, tableresult)
  #tableresult <- tableresult[order(tableresult$pvalue),]
  
  return(tableresult)
}