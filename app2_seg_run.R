args=(commandArgs(TRUE))
if(length(args)==0){
  stop("No arguments supplied.")
}else{
  a = as.numeric(args[1])
  sessionID = args[2]
}

resource_dir <- "/u/home/m/mergeome/PharmOmics_resource/"
data_dir <- "/u/scratch/m/mergeome/app2seg/"

library(igraph)
library(readr)
library(ggplot2)
library(matrixStats)
source(paste0(resource_dir,"PharmOmics_app2_seg_utils.R"))

options(stringsAsFactors = FALSE,warn = -1) 

Mouse_symbols2 <- read.delim(paste0(resource_dir,"Human_Mouse_Symbols_Majority_Ensembl_HGNC_supported.txt"))
load(paste0(data_dir,sessionID, "_Network_server_package.rda"))

load(paste0(resource_dir,"Final_KDA_frame_Limmav2_12.2020.rda"))
splitpoint <- c(0,2,5,9, 13,18 ,23 ,32 ,39, 48, 57 ,66 ,77 ,88, 96)
currentind <- (splitpoint[(a)]*150+1):min(nrow(finalKDAframe),splitpoint[a+1]*150)
finalKDAframe <- finalKDAframe[currentind,]

tableresult <- getResult(diseasegenes = diseasegenes, net = net, 
                         DOBtable = DOBtable, finalKDAframe = finalKDAframe, 
                         alldistancetable = alldistancetable,
                         species = species, gene_signatures = gene_signatures)

save(tableresult, file = paste0(resultPartsDir, "PART_",a,"_result.rda"))

files <- list.files(path = resultPartsDir,
                    pattern = "PART_[0-9]+_result.rda")
if(length(files)!=14){
  # check if any jobs failed from first round
  # load ID from first round of submission (job array)
  load(paste0(data_dir, sessionID, "_first_round_submission_ID.rda"))
  # app2_seg_prepare.R will check the 'redo' jobs that failed
  
  current <- system("myjobs",intern = T)
  currently_running <- current[grep(ID, current)]
  if(length(currently_running)>1){
    inprogress_jobs <- sapply(currently_running, function(x){
      if(strsplit(x,"\\s+")[[1]][6]=="r"){
        return(strsplit(x,"\\s+")[[1]][11])
      } else {
        arr <- strsplit(x,"\\s+")[[1]][10]
        arr <- gsub(":1","",arr)
        return(strsplit(arr, split = "-")[[1]][1]:strsplit(arr, split = "-")[[1]][2])
      }
    })
  }else{
    if(strsplit(currently_running,"\\s+")[[1]][6]=="r"){
      inprogress_jobs <- strsplit(currently_running,"\\s+")[[1]][11]
    } else{
      arr <- strsplit(currently_running,"\\s+")[[1]][10]
      arr <- gsub(":1","",arr)
      inprogress_jobs <- strsplit(arr, split = "-")[[1]]:strsplit(arr, split = "-")[[2]]
    }
  }
  
  if(class(inprogress_jobs)=="list"){
    inprogress_jobs <- unlist(inprogress_jobs)
  }
  
  done_jobs <- setdiff(1:14, inprogress_jobs)
  # check that files were output
  if(length(done_jobs)>0){
    done_files <- sprintf("PART_%s_result.rda",as.character(done_jobs))
    current_files <- list.files(resultPartsDir)
    failed_jobs <- setdiff(done_files, current_files)
    if(length(failed_jobs)>0){
      failed_jobs <- gsub("PART_|_result.rda","",failed_jobs)
      for(f in failed_jobs){
        if(file.exists(paste0(data_dir, sessionID, "_redone.txt"))){
          redone_jobs <- read.delim(paste0(data_dir, sessionID, "_redone.txt"), header = FALSE)
          if(f %in% redone_jobs$V1){
            next
          }
          else{
            system(paste0("echo ", f," >> ",data_dir, sessionID, "_redone.txt"))
          }
        }
        else{
          system(paste0("echo ", f," > ",data_dir, sessionID, "_redone.txt"))
        }
        if(gene_signatures=="top500"){
          system(paste0("qsub -cwd -V -N ", "RUN_",sessionID,"_APP2 ", "-l h_data=6G,h_rt=4:00:00,highp ",resource_dir,"wrapper_app2_redo.sh ",f, " ",sessionID))
        } else {
          system(paste0("qsub -cwd -V -N ", "RUN_",sessionID,"_APP2 ", "-l h_data=6G,h_rt=10:00:00,highp ",resource_dir,"wrapper_app2_redo.sh ",f, " ",sessionID))
        }
      }
    }
  }
}
