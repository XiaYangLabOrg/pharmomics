# commented out portions are added as a function of the user input and random_string

# Genes_up <- read.delim("./Data/Pipeline/Resources/shinyapp3_temp/sessionID_up_genes.txt", stringsAsFactors = FALSE)
# Genes_up <- Genes_up$GENE
# Genes_down <- read.delim("./Data/Pipeline/Resources/shinyapp3_temp/sessionID_down_genes.txt", stringsAsFactors = FALSE)
# Genes_down <- Genes_down$GENE
.libPaths(c("/home/www/abhatta3-webserver/Rnew2/lib64/R/library/",.libPaths()))
library(GeneOverlap)
Jaccard2 <- function(set1,set2){
  I <- length(intersect(set1,set2))
  return(I/(length(set1)+length(set2)-I))
}
Jaccard3 <- function(set1,set2,genomesize){
  go.obj <- newGeneOverlap(set1,
                           set2,
                           genome.size=genomesize)
  go.obj <- testGeneOverlap(go.obj)
  return(list(go.obj@Jaccard, go.obj@odds.ratio, go.obj@pval))
}
load("/home/www/abhatta3-webserver/include/pharmomics/Jaccard_app_databasev3_No_GeoDE.rda")
ADR_scores <- readRDS("/home/www/abhatta3-webserver/include/pharmomics/ADR_scores.rds")
load("/home/www/abhatta3-webserver/include/pharmomics/hepatotox_genes.rda")
cat("5%\n")

count <- 5
total <- nrow(rbind(Humanframe, Ratframe, Mouseframe))
breaks <- levels(cut(1:total, breaks = 95))
upperbounds <- sapply(breaks, function(x){return(unlist(strsplit(x, ","))[2])})
upperbounds <- as.numeric(gsub("]","",upperbounds))
names(upperbounds) <- 6:100
percent <- 5

Genes_up <- unique(Genes_up)
if(length(Genes_down)>0){
  Genes_down <- unique(Genes_down)
}
bothind <- ifelse(length(Genes_down) != 0, T,F)
diseasegenes_up <- Genes_up[Genes_up %in% HUGO_symbols2$`Approved Symbol`]
if(bothind){diseasegenes_down <- Genes_down[Genes_down %in% HUGO_symbols2$`Approved Symbol`]}

species <- "Human"

if(length(diseasegenes_up)/length(Genes_up) < 0.05){
  species <- "non-Human"
  diseasegenes_up_Rat <- Genes_up[Genes_up %in% RAT_symbols2$rat_symbol]
  if(bothind){diseasegenes_down_Rat <- Genes_down[Genes_down %in% RAT_symbols2$rat_symbol]}
  
  diseasegenes_up <-  unique(RAT_symbols2$human_symbol[RAT_symbols2$rat_symbol %in% Genes_up])
  if(bothind){diseasegenes_down <-  unique(RAT_symbols2$human_symbol[RAT_symbols2$rat_symbol %in% Genes_down])}
  
  diseasegenes_up_Mouse <-  Genes_up[Genes_up %in% Mouse_symbols2$mouse_symbol]
  if(bothind){diseasegenes_down_Mouse <-  Genes_up[Genes_up %in% Mouse_symbols2$mouse_symbol]}
}else{
  diseasegenes_up_Rat <- unique(RAT_symbols2$rat_symbol[RAT_symbols2$human_symbol %in% Genes_up])
  if(bothind){diseasegenes_down_Rat <- unique(RAT_symbols2$rat_symbol[RAT_symbols2$human_symbol %in% Genes_down])}
  
  diseasegenes_up_Mouse <-  unique(Mouse_symbols2$mouse_symbol[Mouse_symbols2$human_symbol %in% Genes_up])
  if(bothind){diseasegenes_down_Mouse <-  unique(Mouse_symbols2$mouse_symbol[Mouse_symbols2$human_symbol %in% Genes_down])}
}



Jaccardscore <- NULL
OR <- NULL
pval <- NULL
overlap <- NULL
nOverlap <- NULL
for(r in 1:nrow(Ratframe)){
  count = count + 1
  if(count>upperbounds[percent+1]){
    percent = percent + 1
    cat(percent,"%","\n")
  }
  
  druggenes_up <- Ratgenesup[[r]]
  druggenes_down <- Ratgenesdown[[r]]
  if(bothind){
    Jaccardscore[r] <- Jaccard2(druggenes_up, diseasegenes_up_Rat) + Jaccard2(druggenes_down, diseasegenes_down_Rat)-Jaccard2(druggenes_up,diseasegenes_down_Rat)-Jaccard2(druggenes_down,diseasegenes_up_Rat)
    OR[r] <- Jaccard3(unique(c(druggenes_up,druggenes_down)),unique(c(diseasegenes_up_Rat,diseasegenes_down_Rat)),nrow(RAT_symbols2))[[2]]
    pval[r] <- Jaccard3(unique(c(druggenes_up,druggenes_down)),unique(c(diseasegenes_up_Rat,diseasegenes_down_Rat)),nrow(RAT_symbols2))[[3]]
    intersected <- intersect(unique(c(druggenes_up,druggenes_down)), unique(c(diseasegenes_up_Rat,diseasegenes_down_Rat)))
  }else{
    Jaccardscore[r] <- Jaccard2(unique(c(druggenes_up,druggenes_down)),diseasegenes_up_Rat)
    OR[r] <- Jaccard3(unique(c(druggenes_up,druggenes_down)),diseasegenes_up_Rat,nrow(RAT_symbols2))[[2]]
    pval[r] <- Jaccard3(unique(c(druggenes_up,druggenes_down)),diseasegenes_up_Rat,nrow(RAT_symbols2))[[3]]
    intersected <- intersect(unique(c(druggenes_up,druggenes_down)),diseasegenes_up_Rat)
  }
  if(length(intersected)==0){
    overlap[r] = "None"
    nOverlap[r] <- 0
  } else{
    overlap[r] <- do.call("paste", c(intersected, list(sep=", ")))
    nOverlap[r] <- length(intersected)
  }
}
Ratframe$Jaccardscore <- Jaccardscore
Ratframe$OR <- OR
Ratframe$pval <- pval
if(bothind){
  rank_within_species <- rep(NA,nrow(Ratframe))
  rank_within_species[which(Ratframe$Jaccardscore > 0)] <- rank( Ratframe$Jaccardscore[Ratframe$Jaccardscore > 0])/nrow( Ratframe[Ratframe$Jaccardscore > 0,])
  rank_within_species[which(Ratframe$Jaccardscore < 0)] <- -rank( -Ratframe$Jaccardscore[Ratframe$Jaccardscore < 0])/nrow( Ratframe[Ratframe$Jaccardscore < 0,])
  Ratframe$rank_within_species <- rank_within_species
}else{
  Ratframe$rank_within_species <- rank( Ratframe$Jaccardscore)/nrow( Ratframe)
} 
#cat("30%\n")
Ratframe$Overlap <- overlap
Ratframe$nOverlap <- nOverlap

Jaccardscore <- NULL
OR <- NULL
pval <- NULL
overlap <- NULL
nOverlap <- NULL
for(r in 1:nrow(Mouseframe)){
  count = count + 1
  if(count>upperbounds[percent+1]){
    percent = percent + 1
    cat(percent,"%","\n")
  }
  
  druggenes_up <- Mousegenesup[[r]]
  druggenes_down <- Mousegenesdown[[r]]
  if(bothind){
    Jaccardscore[r] <- Jaccard2(druggenes_up, diseasegenes_up_Mouse) + Jaccard2(druggenes_down, diseasegenes_down_Mouse)-Jaccard2(druggenes_up,diseasegenes_down_Mouse)-Jaccard2(druggenes_down,diseasegenes_up_Mouse)
    OR[r] <- Jaccard3(unique(c(druggenes_up,druggenes_down)),unique(c(diseasegenes_up_Mouse,diseasegenes_down_Mouse)),nrow(Mouse_symbols2))[[2]]
    pval[r] <- Jaccard3(unique(c(druggenes_up,druggenes_down)),unique(c(diseasegenes_up_Mouse,diseasegenes_down_Mouse)),nrow(Mouse_symbols2))[[3]]
    intersected <- intersect(unique(c(druggenes_up,druggenes_down)),unique(c(diseasegenes_up_Mouse,diseasegenes_down_Mouse)))
  }else{
    Jaccardscore[r] <- Jaccard2(unique(c(druggenes_up,druggenes_down)),diseasegenes_up_Mouse)
    OR[r] <- Jaccard3(unique(c(druggenes_up,druggenes_down)),diseasegenes_up_Mouse,nrow(Mouse_symbols2))[[2]]
    pval[r] <- Jaccard3(unique(c(druggenes_up,druggenes_down)),diseasegenes_up_Mouse,nrow(Mouse_symbols2))[[3]]
    intersected <- intersect(unique(c(druggenes_up,druggenes_down)),diseasegenes_up_Mouse)
  }
  if(length(intersected)==0){
    overlap[r] = "None"
    nOverlap[r] <- 0
  } else{
    overlap[r] <- do.call("paste", c(intersected, list(sep=", ")))
    nOverlap[r] <- length(intersected)
  }
}
Mouseframe$Jaccardscore <- Jaccardscore
Mouseframe$OR <- OR
Mouseframe$pval <- pval
if(bothind){
  rank_within_species <- rep(NA,nrow(Mouseframe))
  rank_within_species[which(Mouseframe$Jaccardscore > 0)] <- rank( Mouseframe$Jaccardscore[Mouseframe$Jaccardscore > 0])/nrow( Mouseframe[Mouseframe$Jaccardscore > 0,])
  rank_within_species[which(Mouseframe$Jaccardscore < 0)] <- -rank( -Mouseframe$Jaccardscore[Mouseframe$Jaccardscore < 0])/nrow( Mouseframe[Mouseframe$Jaccardscore < 0,])
  Mouseframe$rank_within_species <- rank_within_species
}else{
  Mouseframe$rank_within_species <- rank( Mouseframe$Jaccardscore)/nrow( Mouseframe)
}
Mouseframe$Overlap <- overlap
Mouseframe$nOverlap <- nOverlap

#cat("50%\n")

Jaccardscore <- NULL
OR <- NULL
pval <- NULL
overlap <- NULL
nOverlap <- NULL
for(r in 1:nrow(Humanframe)){
  count = count + 1
  if(count>upperbounds[percent+1]){
    percent = percent + 1
    cat(percent,"%","\n")
  }
  
  druggenes_up <- Humangenesup[[r]]
  druggenes_down <- Humangenesup[[r]]
  if(bothind){
    Jaccardscore[r] <- Jaccard2(druggenes_up, diseasegenes_up) + Jaccard2(druggenes_down, diseasegenes_down)-Jaccard2(druggenes_up,diseasegenes_down)-Jaccard2(druggenes_down,diseasegenes_up)
    OR[r] <- Jaccard3(unique(c(druggenes_up,druggenes_down)),unique(c(diseasegenes_up,diseasegenes_down)),nrow(HUGO_symbols2))[[2]]
    pval[r] <- Jaccard3(unique(c(druggenes_up,druggenes_down)),unique(c(diseasegenes_up,diseasegenes_down)),nrow(HUGO_symbols2))[[3]]
    intersected <- intersect(unique(c(druggenes_up,druggenes_down)),unique(c(diseasegenes_up,diseasegenes_down)))
  }else{
    Jaccardscore[r] <- Jaccard2(unique(c(druggenes_up,druggenes_down)),diseasegenes_up)
    OR[r] <- Jaccard3(unique(c(druggenes_up,druggenes_down)),diseasegenes_up,nrow(HUGO_symbols2))[[2]]
    pval[r] <- Jaccard3(unique(c(druggenes_up,druggenes_down)),diseasegenes_up,nrow(HUGO_symbols2))[[3]]
    intersected <- intersect(unique(c(druggenes_up,druggenes_down)),diseasegenes_up)
  }
  if(length(intersected)==0){
    overlap[r] = "None"
    nOverlap[r] <- 0
  } else{
    overlap[r] <- do.call("paste", c(intersected, list(sep=", ")))
    nOverlap[r] <- length(intersected)
  }
}

Humanframe$Jaccardscore <- Jaccardscore
Humanframe$OR <- OR
Humanframe$pval <- pval
if(bothind){
  rank_within_species <- rep(NA,nrow(Humanframe))
  rank_within_species[which(Humanframe$Jaccardscore > 0)] <- rank( Humanframe$Jaccardscore[Humanframe$Jaccardscore > 0])/nrow( Humanframe[Humanframe$Jaccardscore > 0,])
  rank_within_species[which(Humanframe$Jaccardscore < 0)] <- -rank( -Humanframe$Jaccardscore[Humanframe$Jaccardscore < 0])/nrow( Humanframe[Humanframe$Jaccardscore < 0,])
  Humanframe$rank_within_species <- rank_within_species
}else{
  Humanframe$rank_within_species <- rank( Humanframe$Jaccardscore)/nrow( Humanframe)
}
Humanframe$Overlap <- overlap
Humanframe$nOverlap <- nOverlap

#cat("70%\n")

result <- rbind.data.frame(Humanframe,Mouseframe,Ratframe)
result = result[order(result$Jaccardscore, decreasing = TRUE),]

result <- result[,c("dataset","method","alldrugs","allspecies","alltissues","allstatus",
                    "dose with unit","alltime","Jaccardscore","OR","pval","rank_within_species","Overlap","nOverlap",
                    "allsignatures_up","allsignatures_down","allsignatures_combined","SIDER_link")]

result$dataset <- gsub(",",", ",result$dataset) 

new_names = c("Database","Method","Drug","Species","Tissue or Cell Line","Study","Dose",
              "Time","Jaccard Score","Odds Ratio","P value","Within Species Rank","Overlap","nOverlap",
              "allsignatures_up","allsignatures_down","allsignatures_combined","SIDER_link")
colnames(result) <- new_names

# ADR rank calculation ----
Genes <- unique(Genes_up)
if(length(Genes_down)>0){
  Genes <- unique(c(Genes, Genes_down))
}

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
if(exists("sessionID")){
  write.table(result, paste0("/home/www/abhatta3-webserver/Data/Pipeline/Results/shinyapp3/",sessionID, "_app3result.txt"), 
              row.names=FALSE, quote = FALSE, sep ="\t")
  write.table(hepatotoxicity_table, paste0("/home/www/abhatta3-webserver/Data/Pipeline/Results/shinyapp3/",sessionID, "_app3result_hepatotox.txt"), 
              row.names=FALSE, quote = FALSE, sep ="\t")
}

#write.table(result, "app3result.txt", row.names=FALSE, quote = FALSE, sep = "\t")



