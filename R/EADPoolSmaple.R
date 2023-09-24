#############################################################################################
################# Expression Aware Demultiplexing per sample pool ###########################
#############################################################################################

#'@import dplyr
#'@import data.table
#'@import stringr
#'@import tidyverse
#'@export
#'@name EADDonorPool
#'@title Expression Aware Demultiplexing per sample pool
#'@description .
#'@author {Isar Nassiri}
#'@param InputDir
#'A folder including matrix of cell-gene expression after Gene Expression Recovery in the SAVER/ folder with 'AllCells.txt' extensions, output of genetic demultiplexing of the sample pool using vireo (donor_ids.tsv), and output of IDCAvis function with extension of "merged.txt" in "IDCA_Analysis/IDCA_Plots/" folder.
#'@return You can find the results in the IDCA_Analysis/Expression_Aware_Cell_Assignment/ folder called "Expression_Aware_Cell_Assignment.txt".
#'@examples
#'library("scDIV")
#'InputDir = system.file("extdata", package = "scDIV")
#'EADDonorPool( InputDir )
#'@export

#----------------------------------
# library(data.table)
# library(stringr)
# library(tidyverse)
# library(dplyr)
# InputDir=paste0('/Users/isarnassiri/Documents/OGC/DCOX_scRNAseq_saver/Package_scripts_inputs/Rpackage/Inputs/')
# EADDonorPool(InputDir)
#----------------------------------

EADDonorPool <- NULL
EADDonorPool <- function(InputDir)
{
  path = paste0(InputDir,
                '/IDCA_Analysis/Expression_Aware_Cell_Assignment/')
  
  #---------------------------------- find the input file
  listfiles = list.files(path, pattern = '_Expression_Aware_Cell_Assignment.txt')
  
  if (length(grep(
    '_Results_Expression_Aware_Cell_Assignment.txt',
    listfiles
  )) > 0)
  {
    listfiles = listfiles[-grep('_Results_Expression_Aware_Cell_Assignment.txt',
                                listfiles)]
  }
  #----------------------------------
  
  #---------------------------------- list of donors
  sp = str_split_fixed(listfiles, "_", 4)
  FCs = unique(sp[, 1])
  
  Donors = unique(c(sp[, 2], sp[, 3]))[order(unique(c(sp[, 2], sp[, 3])))]
  
  cat(paste0("\033[0;", 47, "m", "List of Detected Donors:: ", "\033[0m", "\n"))
  print(gsub(" ", "','", capture.output(cat(Donors))))
  
  #----------------------------------
  
  setwd(path)
  
  for (i in 1:length(FCs))
  {
    FC = FCs[i]
    
    #---------------------------------- merge files
    listfiles_subet = listfiles[grep(FC, listfiles)]
    
    for (t in 1:length(listfiles_subet))
    {
      temp = fread(
        listfiles_subet[t],
        select = c(
          'BARCODE',
          'predicted_clusters',
          'donor_id',
          'best_singlet',
          'best_doublet',
          'Donor1',
          'Donor2',
          'FC'
        ),
        stringsAsFactors = F,
        header = T
      )
      if (t == 1) {
        merged_Pair = temp
      } else{
        merged_Pair = rbind(merged_Pair, temp)
      }
    }
    merged_Pair = as.data.frame(merged_Pair)
    #----------------------------------
    
    iter = 1
    
    
    for (n in Donors)
    {
      #---------------------------------- cells assigned to the donor n
      merged_Pair_subset = merged_Pair[which(merged_Pair$Donor1 == n |
                                               merged_Pair$Donor2 == n), ]
      #----------------------------------
      
      #---------------------------------- remove duplicated cells
      merged_Pair2 = merged_Pair_subset %>% group_by(BARCODE, sort = T) %>% distinct() %>% count(BARCODE, sort = T)
      merged_Pair2 = as.data.frame(merged_Pair2)
      #----------------------------------
      
      #---------------------------------- count number of times which a cell has assigned to a donor
      merged_Pair2 = merged_Pair2[which(merged_Pair2$n > 1), ]
      
      merged_Pair_subset = merge(merged_Pair_subset, merged_Pair2[, c('BARCODE', 'n')], by = 'BARCODE')
      colnames(merged_Pair_subset)[which(colnames(merged_Pair_subset) ==
                                           'n')] = 'NumDonorPairs'
      
      merged_Pair2 = merged_Pair_subset %>% group_by(BARCODE, sort = T) %>% distinct() %>% count(predicted_clusters, sort = T)
      merged_Pair2 = as.data.frame(merged_Pair2)
      
      merged_Pair_subset = merge(merged_Pair_subset, merged_Pair2[, c('BARCODE', 'n')], by = 'BARCODE')
      colnames(merged_Pair_subset)[which(colnames(merged_Pair_subset) ==
                                           'n')] = 'NumAssignedtoDonor'
      
      merged_Pair_subset = merged_Pair_subset[order(merged_Pair_subset$NumAssignedtoDonor, decreasing = T), ]
      merged_Pair_subset = merged_Pair_subset[!duplicated(merged_Pair_subset$BARCODE), ]
      #----------------------------------
      
      #---------------------------------- We confirm that the cell belongs to the donor if we successfully assign it to the donor for an equal or greater number of all pairs of donors minus 1.
      merged_Pair_subset$EA_Assignment = rep('NULL', dim(merged_Pair_subset)[1])
      merged_Pair_subset$EA_Assignment[which(merged_Pair_subset$NumAssignedtoDonor >= (unique(merged_Pair_subset$NumDonorPairs) -
                                                                                         1))] = 'Assigned'
      merged_Pair_subset$EA_Assignment[which(merged_Pair_subset$NumAssignedtoDonor < (unique(merged_Pair_subset$NumDonorPairs) -
                                                                                        1))] = 'unAssigned'
      
      #---------------------------------- save results
      if (iter == 1) {
        RESULT = merged_Pair_subset
        iter = iter + 1
      } else{
        RESULT = rbind(merged_Pair_subset, RESULT)
      }
      
      #----------------------------------
    }
    
    RESULT = RESULT[, which(
      colnames(RESULT) %in% c(
        'BARCODE',
        'predicted_clusters',
        'FC',
        'NumDonorPairs',
        'NumAssignedtoDonor',
        'EA_Assignment'
      )
    )]
    
    #---------------------------------- print outputs for all donors
    print('----------')
    print(paste0('All Cells in ', FC, ': ', dim(RESULT)[1]))
    print('----- Demultiplexing Summary -----')
    RESULT_subset = RESULT[which(RESULT$predicted_clusters == RESULT$best_singlet), ]
    print(paste0('% of assigned cell: ', round(
      dim(RESULT_subset)[1] / dim(RESULT)[1], digits = 3
    ) * 100))
    print(paste0('% of unassigned cell: ', round(1 - (
      dim(RESULT_subset)[1] / dim(RESULT)[1]
    ), digits = 3) * 100))
    print(c(dim(RESULT_subset)[1], dim(RESULT)[1]))
    
    cat(
      paste0('All Cells in ', FC, ':', dim(RESULT)[1]),
      file = paste0(path, FC, '_Summary.txt'),
      append = TRUE
    )
    cat("\n",
        file = paste0(path, FC, '_Summary.txt'),
        append = TRUE)
    cat(
      '----- Mixed -----',
      file = paste0(path, FC, '_Summary.txt'),
      append = TRUE
    )
    cat("\n",
        file = paste0(path, FC, '_Summary.txt'),
        append = TRUE)
    cat(
      paste0('% of assigned cell: ', round(
        dim(RESULT_subset)[1] / dim(RESULT)[1], digits = 3
      ) * 100),
      file = paste0(path, FC, '_Summary.txt'),
      append = TRUE
    )
    cat("\n",
        file = paste0(path, FC, '_Summary.txt'),
        append = TRUE)
    cat(
      paste0('% of unassigned cell: ', round(1 - (
        dim(RESULT_subset)[1] / dim(RESULT)[1]
      ), digits = 3) * 100),
      file = paste0(path, FC, '_Summary.txt'),
      append = TRUE
    )
    cat("\n",
        file = paste0(path, FC, '_Summary.txt'),
        append = TRUE)
    cat('----------',
        file = paste0(path, FC, '_Summary.txt'),
        append = TRUE)
    cat("\n",
        file = paste0(path, FC, '_Summary.txt'),
        append = TRUE)
    cat("\n",
        file = paste0(path, FC, '_Summary.txt'),
        append = TRUE)
    #----------------------------------
    
    #---------------------------------- save outputs for all donors
    fwrite(
      RESULT,
      paste0(
        path,
        FC,
        '_Results_Expression_Aware_Cell_Assignment.txt'
      ),
      quote = F,
      row.names = F,
      sep = '\t'
    )
    #----------------------------------
  }
  
  cat(
    paste0(
      "\033[0;",
      47,
      "m",
      "You can find the results in: ",
      "\033[0m",
      "\n",
      InputDir,
      "/IDCA_Analysis/Expression_Aware_Cell_Assignment/"
    ),
    paste0(
      ": ",
      FC,
      '_Summary.txt & ',
      paste0(FC, '_Results_Expression_Aware_Cell_Assignment.txt')
    )
  )
  
}