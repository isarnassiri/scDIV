#############################################################################################
################# Expression Aware Demultiplexing per Donor Pair ############################
#############################################################################################

#'@import dplyr
#'@import data.table
#'@import stringr
#'@import flexmix
#'@export 
#'@name EADDonorPair
#'@title Expression Aware Demultiplexing per Donor Pair
#'@description .
#'@author {Isar Nassiri}
#'@param InputDir
#'#'A folder including matrix of cell-gene expression after Gene Expression Recovery in the SAVER/ folder with 'AllCells.txt' extensions, output of genetic demultiplexing of the sample pool using vireo (donor_ids.tsv), and output of IDCAvis function with extension of "merged.txt" in "IDCA_Analysis/IDCA_Plots/" folder.
#'@return You can find the results in the IDCA_Analysis/Expression_Aware_Cell_Assignment/ folder called "Expression_Aware_Cell_Assignment.txt".
#'@examples
#'library("scDIV")
#'InputDir = system.file("extdata", package = "scDIV")
#'EADDonorPair( InputDir )
#'@export

#----------------------------------
# library(data.table)
# library(dplyr)
# library(stringr)
# library(flexmix)
# InputDir=paste0('/Users/isar/Documents/Expression_aware_demultiplexing/Inputs/')
# EADDonorPair(InputDir)
#----------------------------------

EADDonorPair <- NULL
EADDonorPair <- function( InputDir )
{

#---------------------------------- read in input files
listfiles = list.files(paste0(InputDir, '/SAVER/'), pattern = '_AllCells.txt')
dir.create(paste0('IDCA_Analysis/Expression_Aware_Cell_Assignment'), recursive = T)

for (S in 1:length(listfiles))
{

  #---------------------------------- generate name for objects
  S0 = gsub('_AllCells.txt','',listfiles[S])
  
  sp = str_split_fixed(S0, "_", 3)
 
  S1=sp[1,1]
  S2=sp[1,2]
  FC=sp[1,3]
  #----------------------------------
  
  #---------------------------------- read in the transformed gene expression profile for an indicated pair of donors
  estimate = fread(paste0(paste0(InputDir, '/SAVER/'), listfiles[S]), stringsAsFactors = F, header = T)
  estimate = as.data.frame(estimate)
  
  row.names(estimate) = make.names(estimate$GeneID, unique = T)
  estimate = estimate[,-c(1)]
  
  colnames(estimate) = gsub('.1', '-1', colnames(estimate))
  #----------------------------------
  
  #---------------------------------- if previous steps provided results for S0
  if(file.exists(paste0(InputDir, 'IDCA_Analysis/IDCA_Plots/',S0,'_merged.txt')))
  {
    merged_Pair = fread( paste0(InputDir, '/IDCA_Analysis/IDCA_Plots/',S0,'_merged.txt'), stringsAsFactors = F, header = T)
    merged_Pair = merged_Pair[order(-as.numeric(merged_Pair$TPrate), merged_Pair$pValDiff_adj),]

    #---------------------------------- read in results of genetic demultiplexing
    vireo = fread(paste0(InputDir,'/donor_ids.tsv'), stringsAsFactors = F, header = T)
    vireo = as.data.frame(vireo)
    colnames(vireo)[1] = 'BARCODE'
    
    #---------------------------------- select top first inter-individual differential correlated gene pair
    Gene1=merged_Pair$Gene2[1] # in RESULT_merged, Gene2 is the original Gene1
    Gene2=merged_Pair$Gene1[1]
    #----------------------------------
    
    #---------------------------------- subset gene expression profile
    estimate_Gene1 = estimate[which(row.names(estimate) == Gene1),]
    estimate_Gene2 = estimate[which(row.names(estimate) == Gene2),]
    #----------------------------------
    
    #----------------------------------
    input_mixture_model = t(rbind(as.matrix(estimate_Gene1), as.matrix(estimate_Gene2)))
    input_mixture_model = as.data.frame(input_mixture_model, stringsAsFactors = F)
    
    mixed_model = flexmix(cbind(as.numeric(input_mixture_model[,1]),as.numeric(input_mixture_model[,2]))~1,
    k=2, data=input_mixture_model,
    model = FLXMCmvnorm(diag = T), #diag = T to get fix results FLXMCmvpois(), FLXMCmvnorm
    control = list(tolerance = 1e-15, iter.max = 1000))
    
    #--------------------------- assignments comparison
    input_mixture_model$predicted_clusters = clusters(mixed_model)
    identical(row.names(input_mixture_model), colnames(estimate))
    
    input_mixture_model$BARCODE = row.names(input_mixture_model)
    merged = merge(input_mixture_model, vireo, by = 'BARCODE')
    #----------------------------------
    
    #----------------------------------
    if(length(unique((input_mixture_model$predicted_clusters))) > 1)
    {
      freq_c1 = data.frame(table(merged[which(merged$predicted_clusters == 1), "best_singlet"]), stringsAsFactors = F)
      freq_c1 = freq_c1[order(freq_c1$Freq, decreasing = T),]
 
      merged[which(merged$predicted_clusters == 1),"predicted_clusters"] = as.character(freq_c1$Var1[1])
      
      freq_c2 = data.frame(table(merged[which(merged$predicted_clusters == 2),"best_singlet"]), stringsAsFactors = F)
      freq_c2 = freq_c2[order(freq_c2$Freq, decreasing = T),]
 
      merged[which(merged$predicted_clusters == 2),"predicted_clusters"] = as.character(freq_c2$Var1[1])
      merged$Donor1=rep(S1, dim(merged)[1])
      merged$Donor2=rep(S2, dim(merged)[1])
      merged$FC=rep(FC, dim(merged)[1])
      
      #----------------------------------
      fwrite(merged, paste0('IDCA_Analysis/Expression_Aware_Cell_Assignment/', FC, '_', S1, '_', S2, '_Expression_Aware_Cell_Assignment.txt'), row.names = F, quote = F, sep = '\t')
      #----------------------------------
 
    }else{print(paste0(S," had one predicted cluster."))}
  }
}

cat(paste0("\033[0;", 47, "m", "You can find the results in: ", "\033[0m","\n", InputDir, "/IDCA_Analysis/Expression_Aware_Cell_Assignment/"))

}
