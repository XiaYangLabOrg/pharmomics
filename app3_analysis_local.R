# commented out portions are added as a function of the user input and random_string
# Parse command line arguments
# args[1] = consider_direction (mandatory: TRUE/FALSE or T/F or 1/0)
# args[2] = GENECOUNT_THRESHOLD (optional, default: 40)
args <- commandArgs(trailingOnly = TRUE)

if(length(args) < 1) {
  stop("Error: consider_direction is required as the first command line argument (TRUE/FALSE or T/F or 1/0)")
}

# Parse consider_direction (mandatory)
consider_direction_arg <- toupper(trimws(args[1]))
if(consider_direction_arg %in% c("TRUE", "T", "1")) {
  consider_direction <- TRUE
} else if(consider_direction_arg %in% c("FALSE", "F", "0")) {
  consider_direction <- FALSE
} else {
  stop("Error: consider_direction must be TRUE/FALSE, T/F, or 1/0, got: ", args[1])
}

# Parse GENECOUNT_THRESHOLD (optional, default 40)
GENECOUNT_THRESHOLD <- if(length(args) > 1 && !is.na(as.numeric(args[2]))) {
  as.numeric(args[2])
} else {
  40
}

cat("consider_direction =", consider_direction, "\n")
cat("GENECOUNT_THRESHOLD =", GENECOUNT_THRESHOLD, "\n")

# Load gene files based on consider_direction setting
if(consider_direction) {
  # When considering direction, load from separate up and down files
  if(!file.exists("up_genes.txt")) {
    stop("Error: 'up_genes.txt' file not found. When consider_direction=TRUE, the script requires 'up_genes.txt' in the working directory.")
  }
  Genes_up <- read.delim("up_genes.txt", stringsAsFactors = FALSE)
  if(!"GENE" %in% colnames(Genes_up)) {
    stop("Error: 'up_genes.txt' must contain a column named 'GENE'")
  }
  Genes_up <- Genes_up$GENE
  
  # When considering direction, down_genes.txt is optional but must have GENE column if it exists
  if(file.exists("down_genes.txt")) {
    Genes_down_df <- read.delim("down_genes.txt", stringsAsFactors = FALSE)
    if(!"GENE" %in% colnames(Genes_down_df)) {
      stop("Error: 'down_genes.txt' must contain a column named 'GENE'")
    }
    Genes_down <- Genes_down_df$GENE
  } else {
    # If down_genes.txt doesn't exist, use empty character vector
    Genes_down <- character(0)
  }
} else {
  # When not considering direction, load from single 'genes.txt' file
  if(!file.exists("genes.txt")) {
    stop("Error: 'genes.txt' file not found. When consider_direction=FALSE, the script requires 'genes.txt' in the working directory.")
  }
  Genes_all <- read.delim("genes.txt", stringsAsFactors = FALSE)
  if(!"GENE" %in% colnames(Genes_all)) {
    stop("Error: 'genes.txt' must contain a column named 'GENE'")
  }
  Genes_up <- Genes_all$GENE
  Genes_down <- Genes_all$GENE  # Use same genes for both when not considering direction
}

library(GeneOverlap)
# Count comma-delimited gene names (handles NA/blank/extra spaces)
count_genes <- function(x) {
  ifelse(
    is.na(x) | trimws(x) == "",
    0L,
    vapply(strsplit(x, ",", fixed = TRUE), function(v) {
      sum(nzchar(trimws(v)))
    }, integer(1))
  )
}

j2_dbg <- function(a, b) {
  a <- unique(a[!is.na(a)])
  b <- unique(b[!is.na(b)])
  u <- union(a, b)
  i <- intersect(a, b)
  ulen <- length(u); ilen <- length(i)
  val <- if (ulen == 0) NaN else ilen / ulen
  list(
    val   = val,
    nA    = length(a),
    nB    = length(b),
    nInt  = ilen,
    nUnion= ulen,
    emptyA= (length(a) == 0),
    emptyB= (length(b) == 0)
  )
}

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

# Helper function to calculate rank within species
calculate_rank_within_species <- function(frame, consider_direction) {
  if(consider_direction){
    rank_within_species <- rep(NA, nrow(frame))
    pos_idx <- which(frame$Jaccardscore > 0)
    neg_idx <- which(frame$Jaccardscore < 0)
    if(length(pos_idx) > 0) {
      rank_within_species[pos_idx] <- rank(frame$Jaccardscore[pos_idx]) / length(pos_idx)
    }
    if(length(neg_idx) > 0) {
      rank_within_species[neg_idx] <- -rank(-frame$Jaccardscore[neg_idx]) / length(neg_idx)
    }
    return(rank_within_species)
  } else {
    return(rank(frame$Jaccardscore) / nrow(frame))
  }
}

# Helper function to print statistics
print_statistics <- function(frame, species_name, include_detailed = FALSE) {
  js <- frame$Jaccardscore
  pos <- which(js > 0)
  neg <- which(js < 0)
  zer <- which(js == 0)
  cat(species_name, ": nrow=", nrow(frame), 
      " len(js)=", length(js),
      " len(pos)=", length(pos),
      " len(neg)=", length(neg),
      " len(zer)=", length(zer),
      " len(rank(js[neg]))=", length(rank(js[neg])),
      " len(rank(js[pos]))=", length(rank(js[pos])), "\n")
  
  if(include_detailed) {
    only_na   <- sum(is.na(js) & !is.nan(js))
    only_nan  <- sum(is.nan(js))
    pos_inf   <- sum(js ==  Inf, na.rm = TRUE)
    neg_inf   <- sum(js == -Inf, na.rm = TRUE)
    finite    <- sum(is.finite(js))
    
    cat("Only NA: ", only_na,  "\n",
        "Only NaN:", only_nan, "\n",
        "+Inf:     ", pos_inf, "\n",
        "-Inf:     ", neg_inf, "\n",
        "Finite:   ", finite,  "\n", sep = "")
  }
}

# Helper function to process a species frame
process_species_frame <- function(frame_thresh, idx_orig, genesup_list, genesdown_list,
                                 diseasegenes_up, diseasegenes_down, consider_direction,
                                 symbol_table, species_name, count, percent, upperbounds,
                                 debug_mode = FALSE) {
  n <- nrow(frame_thresh)
  Jaccardscore <- numeric(n)
  OR <- numeric(n)
  pval <- numeric(n)
  overlap <- character(n)
  nOverlap <- integer(n)
  
  for(r in seq_len(n)) {
    count <- count + 1
    if(count > upperbounds[percent + 1]) {
      percent <- percent + 1
      cat(percent, "%", "\n")
    }
    
    r_orig <- idx_orig[r]
    druggenes_up <- genesup_list[[r_orig]]
    druggenes_down <- genesdown_list[[r_orig]]
    
    if(consider_direction) {
      if(debug_mode && species_name == "Rat") {
        # Use debug version for Rat
        p1 <- j2_dbg(druggenes_up,   diseasegenes_up)    # up vs up
        p2 <- j2_dbg(druggenes_down, diseasegenes_down)  # down vs down
        p3 <- j2_dbg(druggenes_up,   diseasegenes_down)  # up vs down
        p4 <- j2_dbg(druggenes_down, diseasegenes_up)    # down vs up
        
        val <- p1$val + p2$val - p3$val - p4$val
        Jaccardscore[r] <- val
        if (is.nan(Jaccardscore[r])) {
          cat(sprintf(
            "NaN @ r=%d | up-up: val=%s (A=%d B=%d ∩=%d ∪=%d emptyA=%s emptyB=%s) | down-down: val=%s (A=%d B=%d ∩=%d ∪=%d emptyA=%s emptyB=%s) | up-down: val=%s (A=%d B=%d ∩=%d ∪=%d) | down-up: val=%s (A=%d B=%d ∩=%d ∪=%d)\n",
            r,
            ifelse(is.nan(p1$val), "NaN", sprintf("%.4f", p1$val)), p1$nA, p1$nB, p1$nInt, p1$nUnion, p1$emptyA, p1$emptyB,
            ifelse(is.nan(p2$val), "NaN", sprintf("%.4f", p2$val)), p2$nA, p2$nB, p2$nInt, p2$nUnion, p2$emptyA, p2$emptyB,
            ifelse(is.nan(p3$val), "NaN", sprintf("%.4f", p3$val)), p3$nA, p3$nB, p3$nInt, p3$nUnion,
            ifelse(is.nan(p4$val), "NaN", sprintf("%.4f", p4$val)), p4$nA, p4$nB, p4$nInt, p4$nUnion
          ))
        }
      } else {
        # Standard calculation for Mouse and Human
        Jaccardscore[r] <- Jaccard2(druggenes_up, diseasegenes_up) + 
                          Jaccard2(druggenes_down, diseasegenes_down) - 
                          Jaccard2(druggenes_up, diseasegenes_down) - 
                          Jaccard2(druggenes_down, diseasegenes_up)
      }
      
      drug_combined <- unique(c(druggenes_up, druggenes_down))
      disease_combined <- unique(c(diseasegenes_up, diseasegenes_down))
      j3_result <- Jaccard3(drug_combined, disease_combined, nrow(symbol_table))
      OR[r] <- j3_result[[2]]
      pval[r] <- j3_result[[3]]
      intersected <- intersect(drug_combined, disease_combined)
    } else {
      drug_combined <- unique(c(druggenes_up, druggenes_down))
      Jaccardscore[r] <- Jaccard2(drug_combined, diseasegenes_up)
      j3_result <- Jaccard3(drug_combined, diseasegenes_up, nrow(symbol_table))
      OR[r] <- j3_result[[2]]
      pval[r] <- j3_result[[3]]
      intersected <- intersect(drug_combined, diseasegenes_up)
    }
    
    if(length(intersected) == 0) {
      overlap[r] <- "None"
      nOverlap[r] <- 0
    } else {
      overlap[r] <- paste(intersected, collapse = ", ")
      nOverlap[r] <- length(intersected)
    }
  }
  
  frame_thresh$Jaccardscore <- Jaccardscore
  frame_thresh$OR <- OR
  frame_thresh$pval <- pval
  frame_thresh$Overlap <- overlap
  frame_thresh$nOverlap <- nOverlap
  
  return(list(frame = frame_thresh, count = count, percent = percent))
}
# Check and load required database files
if(!file.exists("Jaccard_app_databasev3_No_GeoDE.rda")) {
  stop("Error: 'Jaccard_app_databasev3_No_GeoDE.rda' file not found in the working directory.")
}
load("Jaccard_app_databasev3_No_GeoDE.rda")

if(!file.exists("ADR_scores.rds")) {
  stop("Error: 'ADR_scores.rds' file not found in the working directory.")
}
ADR_scores <- readRDS("ADR_scores.rds")

if(!file.exists("hepatotox_genes.rda")) {
  stop("Error: 'hepatotox_genes.rda' file not found in the working directory.")
}
load("hepatotox_genes.rda")
cat("5%\n")

cat("Before filtering:\n")
cat("  Humanframe:", nrow(Humanframe), "\n")
cat("  Ratframe:  ", nrow(Ratframe), "\n")
cat("  Mouseframe:", nrow(Mouseframe), "\n")

Humanframe_thresh <- Humanframe[count_genes(Humanframe$allsignatures_combined) >= GENECOUNT_THRESHOLD, , drop = FALSE]
idx_human <- which(count_genes(Humanframe$allsignatures_combined) >= GENECOUNT_THRESHOLD)
Ratframe_thresh   <- Ratframe[count_genes(Ratframe$allsignatures_combined) >= GENECOUNT_THRESHOLD,   , drop = FALSE]
idx_rat <- which(count_genes(Ratframe$allsignatures_combined) >= GENECOUNT_THRESHOLD)
Mouseframe_thresh <- Mouseframe[count_genes(Mouseframe$allsignatures_combined) >= GENECOUNT_THRESHOLD, , drop = FALSE]
idx_mouse <- which(count_genes(Mouseframe$allsignatures_combined) >= GENECOUNT_THRESHOLD)

cat("\nAfter filtering:\n")
cat("  Humanframe:", nrow(Humanframe_thresh), "\n")
cat("  Ratframe:  ", nrow(Ratframe_thresh), "\n")
cat("  Mouseframe:", nrow(Mouseframe_thresh), "\n")

count <- 5
total <- nrow(rbind(Humanframe_thresh, Ratframe_thresh, Mouseframe_thresh))
breaks <- levels(cut(1:total, breaks = 95))
upperbounds <- sapply(breaks, function(x){return(unlist(strsplit(x, ","))[2])})
upperbounds <- as.numeric(gsub("]","",upperbounds))
names(upperbounds) <- 6:100
percent <- 5

Genes_up <- unique(Genes_up)
if(length(Genes_down)>0){
  Genes_down <- unique(Genes_down)
} else {
  Genes_down <- character(0)  # Ensure Genes_down is defined even if empty
}
diseasegenes_up <- Genes_up[Genes_up %in% HUGO_symbols2$`Approved Symbol`]
if(consider_direction){
  diseasegenes_down <- Genes_down[Genes_down %in% HUGO_symbols2$`Approved Symbol`]
}

species <- "Human"

if(length(diseasegenes_up)/length(Genes_up) < 0.05){
  species <- "non-Human"
  diseasegenes_up_Rat <- Genes_up[Genes_up %in% RAT_symbols2$rat_symbol]
  if(consider_direction){diseasegenes_down_Rat <- Genes_down[Genes_down %in% RAT_symbols2$rat_symbol]}
  
  diseasegenes_up <-  unique(RAT_symbols2$human_symbol[RAT_symbols2$rat_symbol %in% Genes_up])
  if(consider_direction){diseasegenes_down <-  unique(RAT_symbols2$human_symbol[RAT_symbols2$rat_symbol %in% Genes_down])}
  
  diseasegenes_up_Mouse <-  Genes_up[Genes_up %in% Mouse_symbols2$mouse_symbol]
  if(consider_direction){diseasegenes_down_Mouse <-  Genes_up[Genes_up %in% Mouse_symbols2$mouse_symbol]}
}else{
  diseasegenes_up_Rat <- unique(RAT_symbols2$rat_symbol[RAT_symbols2$human_symbol %in% Genes_up])
  if(consider_direction){diseasegenes_down_Rat <- unique(RAT_symbols2$rat_symbol[RAT_symbols2$human_symbol %in% Genes_down])}
  
  diseasegenes_up_Mouse <-  unique(Mouse_symbols2$mouse_symbol[Mouse_symbols2$human_symbol %in% Genes_up])
  if(consider_direction){diseasegenes_down_Mouse <-  unique(Mouse_symbols2$mouse_symbol[Mouse_symbols2$human_symbol %in% Genes_down])}
}

# Process Rat frame
result_rat <- process_species_frame(
  Ratframe_thresh, idx_rat, Ratgenesup, Ratgenesdown,
  diseasegenes_up_Rat, diseasegenes_down_Rat, consider_direction,
  RAT_symbols2, "Rat", count, percent, upperbounds, debug_mode = TRUE
)
Ratframe_thresh <- result_rat$frame
count <- result_rat$count
percent <- result_rat$percent

print_statistics(Ratframe_thresh, "Rat", include_detailed = TRUE)
Ratframe_thresh$rank_within_species <- calculate_rank_within_species(Ratframe_thresh, consider_direction)

# Process Mouse frame
result_mouse <- process_species_frame(
  Mouseframe_thresh, idx_mouse, Mousegenesup, Mousegenesdown,
  diseasegenes_up_Mouse, diseasegenes_down_Mouse, consider_direction,
  Mouse_symbols2, "Mouse", count, percent, upperbounds, debug_mode = FALSE
)
Mouseframe_thresh <- result_mouse$frame
count <- result_mouse$count
percent <- result_mouse$percent

print_statistics(Mouseframe_thresh, "Mouse", include_detailed = FALSE)
Mouseframe_thresh$rank_within_species <- calculate_rank_within_species(Mouseframe_thresh, consider_direction)

# Process Human frame
result_human <- process_species_frame(
  Humanframe_thresh, idx_human, Humangenesup, Humangenesdown,
  diseasegenes_up, diseasegenes_down, consider_direction,
  HUGO_symbols2, "Human", count, percent, upperbounds, debug_mode = FALSE
)
Humanframe_thresh <- result_human$frame
count <- result_human$count
percent <- result_human$percent

print_statistics(Humanframe_thresh, "Human", include_detailed = FALSE)
Humanframe_thresh$rank_within_species <- calculate_rank_within_species(Humanframe_thresh, consider_direction)

result <- rbind.data.frame(Humanframe_thresh,Mouseframe_thresh,Ratframe_thresh)
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

input_hepatotox <- intersect(hepatotox_genes, Genes)

if(length(input_hepatotox)==0){
  input_hepatotox <- ""
}

hepatotoxicity_table <- data.frame("Adverse drug reaction"="CTD_Chemical induced liver injury full signature",
                                   "Jaccard score"=hepatotox_jaccard,
                                   "Jaccard score rank"=ranks[1],
                                   "Jaccard p value"=pvals[1],
                                   "ADR genes and input gene overlap"=do.call("paste",c(input_hepatotox,
                                                                                        list("sep"=", "))))
#if(exists("sessionID")){
  write.table(result, paste0("_app3result.txt"), 
              row.names=FALSE, quote = FALSE, sep ="\t")
  write.table(hepatotoxicity_table, paste0("_app3result_hepatotox.txt"), 
              row.names=FALSE, quote = FALSE, sep ="\t")
#}

#write.table(result, "app3result.txt", row.names=FALSE, quote = FALSE, sep = "\t")



