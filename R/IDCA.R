#############################################################################################
############# Inter-individual Differential gene Correlation Analysis (IDCA) ################
#############################################################################################

#'@import glmnet
#'@import DGCA
#'@import grDevices
#'@import grid
#'@import stringr
#'@import data.table
#'@import flexmix
#'@import mvtnorm
#'@export 
#'@name IDCA
#'@title Inter-individual Differential gene Correlation Analysis (IDCA)
#'@description The functions uses correlation coefficients and performed Inter-individual Differential gene Correlation Analysis (IDCA) for two donors (D1 and D2) and genes (G1 and G2).
#'@author {Isar Nassiri}
#'@param InputDir
#' A folder including gene-cell count matrix (read_count.csv), filtered feature-barcode matrices from 10X CellRanger count (barcodes.tsv.gz), and output of genetic demultiplexing of the sample pool using vireo (donor_ids.tsv).
#'@param Donors
#' A pair of selected donors (e.g., donor6_donor2)
#'@param FC
#' A name for flow cell (e.g., FAI5649A17).
#'@param ERP
#'Matrix of cell-gene expression after Gene Expression Recovery in the SAVER/ folder with 'AssignedCells.txt' and 'AllCells.txt' extensions.
#'@param TEST
#'use TEST = TRUE to get sample outputs.
#'@return You can find the results in the IDCA_Analysis/ folder.
#'@examples
#'library("scDIV")
#'ERP = "donor6_donor2_FAI5649A17_AssignedCells.txt"
#'Donors='donor6_donor2'
#'FC='FAI5649A17'
#'InputDir = system.file("extdata", package = "scDIV")
#'IDCA( InputDir, Donors, FC, ERP, TEST = T )
#'@export

#----------------------------------
# library(DGCA)
# library(grDevices)
# library(grid)
# library(stringr)
# library(data.table)
# library(flexmix)
# library(mvtnorm)
# library(glmnet)
# 
# ERP = "donor6_donor2_FAI5649A17_AssignedCells.txt"
# Donors='donor6_donor2'
# FC='FAI5649A17'
# InputDir = '/Users/isar/Documents/Expression_aware_demultiplexing/Inputs/'
# IDCA( InputDir, Donors, FC, ERP, TEST = T )
#----------------------------------

IDCA <- NULL
IDCA <- function( InputDir, Donors, FC, ERP, TEST = FALSE )
{
  
#---------------
setwd(InputDir)
dir.create('IDCA_Analysis/')
OutPutPath <- paste0('IDCA_Analysis/')
#---------------

#---------------------------------- read in list of barcodes
estimate <- fread(paste0('SAVER/', ERP), stringsAsFactors = FALSE, header = TRUE)
estimate <- as.data.frame(estimate)

row.names(estimate) <- make.names(estimate$GeneID, unique = TRUE)
estimate <- estimate[,-c(1)]
#----------------------------------

#---------------------------------- read in results of genetic based demultiplexing (vireo)
vireo <- fread('donor_ids.tsv', stringsAsFactors = FALSE, header = TRUE)
vireo <- as.data.frame(vireo)
colnames(vireo)[1] <- 'BARCODE'
#----------------------------------

#---------------------------------- list of donors
sp <- str_split_fixed(Donors, "_", 2)

D1<-sp[1,1]
D2<-sp[1,2]

#----------------------------------

#---------------------------------- subset gene expression profiles per selected pair of donors
S1_BARCODE <- vireo[which(vireo$donor_id == D1),'BARCODE']
S2_BARCODE <- vireo[which(vireo$donor_id == D2),'BARCODE']

S1_BARCODE <- gsub('-', '.', S1_BARCODE)
S2_BARCODE <- gsub('-', '.', S2_BARCODE)

estimate <- estimate[,which(colnames(estimate) %in% c(S1_BARCODE, S2_BARCODE))]
#----------------------------------

#---------------------------------- inter-individual differential gene correlation analysis (IDCA) for two donors and two genes

considered_GenePairs <- vector()

if(!TEST)
{
for (i in 1:dim(estimate)[1]) {

SelectedGene <- row.names(estimate)[i] # e.g., 'RPL31'

AllGenes<-row.names(estimate)[-which(row.names(estimate) == SelectedGene)]

#---------------------------------- obtain the most representative genes subset using least absolute shrinkage and selection operator (LASSO) 
x <- t(estimate[-which(row.names(estimate) == SelectedGene),])
y <- as.numeric(estimate[which(row.names(estimate) == SelectedGene),])

set.seed(123)
train<-sample(seq(dim(x)[1]),(dim(x)[1]/2), replace=FALSE)

lasso.tr<-glmnet(x[train,],y[train], standardize= FALSE)
pred<-predict(lasso.tr,x[-train,])
rmse<- sqrt(apply((y[-train]-pred)^2,2,mean))
lam.best<-lasso.tr$lambda[order(rmse)[1]]

lasso.tr<-glmnet(x ,y, standardize= FALSE)

Lasso_coefficient <- coef(lasso.tr, s=lam.best)
inds<-which(Lasso_coefficient!=0)
variables<-row.names(Lasso_coefficient)[inds]
variables<-variables[!(variables %in% '(Intercept)')];

length(which(0 != Lasso_coefficient[,1]))
results <- as.data.frame(Lasso_coefficient[which(0 != Lasso_coefficient[,1]),1])
colnames(results) <- "coef"
results <- as.data.frame(results[-1,,FALSE]) 

AllGenes <- row.names(results)
#----------------------------------

RESULT_merged <- data.frame()

if(length(AllGenes) > 0)
{
  for(j in 1:length(AllGenes))
  {
    Gene1<- SelectedGene; 
    Gene2<-AllGenes[j];
    
    #---------------------------------- consider each possible gene pair one time
    if( length(considered_GenePairs[which(considered_GenePairs %in% c(paste(Gene1, Gene2, sep = '_'), paste(Gene2, Gene1, sep = '_')))]) > 0)
    {
      print('already is considered.')
    }else{
    
    considered_GenePairs[length(considered_GenePairs)+1] <- paste(Gene1, Gene2, sep = '_')
    considered_GenePairs[length(considered_GenePairs)+1] <- paste(Gene2, Gene1, sep = '_')
    #----------------------------------
    
    estimate_Gene1 <- estimate[which(row.names(estimate) == Gene1),]
    estimate_Gene2 <- estimate[which(row.names(estimate) == Gene2),]

    #---------------------------------- Mixture Modeling
    input_mixture_model <- t(rbind(as.matrix(estimate_Gene1), as.matrix(estimate_Gene2)))
    input_mixture_model <- as.data.frame(input_mixture_model, stringsAsFactors = FALSE)
    
    mixed_model <- flexmix(cbind(as.numeric(input_mixture_model[,1]),as.numeric(input_mixture_model[,2]))~1, k=2, data=input_mixture_model, model = FLXMCmvnorm(diag = TRUE), #diag = TRUE to get fix results FLXMCmvpois(), FLXMCmvnorm
    control <- list(tolerance = 1e-15, iter.max = 1000))
    
    if(length(unique(clusters(mixed_model)))==2)
    {
    #---------------------------------- IDCA
    design_matrix <- cbind(colnames(estimate), colnames(estimate))
    design_matrix[-which(design_matrix[,1]%in%S1_BARCODE),1]=0
    design_matrix[which(design_matrix[,1]%in%S1_BARCODE),1]=1
    design_matrix[-which(design_matrix[,2]%in%S2_BARCODE),2]=0
    design_matrix[which(design_matrix[,2]%in%S2_BARCODE),2]=1
    colnames(design_matrix) <- c("G0", "G1")
    
    design_matrix <- as.data.frame(design_matrix)
    design_matrix[,1] <- as.numeric(design_matrix[,1])
    design_matrix[,2] <- as.numeric(design_matrix[,2])
    design_matrix <- as.matrix(design_matrix)
    ddcor_G1_G0 <- ddcorAll(inputMat = estimate[which(row.names(estimate) %in% c(Gene1,Gene2)),], design = design_matrix,compare = c("G1","G0"),adjust = "fdr", nPerm = 0, corrType = "spearman", splitSet = Gene1, sigThresh = 0.05, sortBy = "pValDiff_adj", verbose = TRUE)
    
    #---------------------------------- compare the results of mixture model and genetic based demultiplexing
    input_mixture_model$real_clusters <- colnames(estimate)
    input_mixture_model$real_clusters[which(input_mixture_model$real_clusters%in%S1_BARCODE)]<-1
    input_mixture_model$real_clusters[which(input_mixture_model$real_clusters%in%S2_BARCODE)]<-2
    
    input_mixture_model$predicted_clusters <- clusters(mixed_model)
    #----------------------------------
    
    #---------------------------------- revise of cluster names
    if(input_mixture_model$predicted_clusters[1] == 2)
    {
      input_mixture_model$predicted_clusters[which(input_mixture_model$predicted_clusters == 1)]<-0
      input_mixture_model$predicted_clusters[which(input_mixture_model$predicted_clusters == 2)]<-1
      input_mixture_model$predicted_clusters[which(input_mixture_model$predicted_clusters == 0)]<-2
    }
    #----------------------------------
    
    #---------------------------------- save the results
    table <- table(input_mixture_model$real_clusters,input_mixture_model$predicted_clusters)
    TPrate <- (table[1,1] + table[2,2])/length(clusters(mixed_model))
    SM <- summary(mixed_model)
    
    RESULT <- data.frame()
    RESULT <- ddcor_G1_G0
    RESULT$AIC <- SM@AIC
    RESULT$TPrate <- TPrate
    RESULT <- RESULT[which(RESULT$pValDiff_adj < 1e-10 ),]# & RESULT$G0_pVal != 0.0
    
    if(dim(RESULT)[1]>0)
    {
      if(j==1){RESULT_merged = RESULT}
      if(j!=1 & dim(RESULT_merged)[1] != 0){RESULT_merged <- rbind(RESULT_merged, RESULT)}
      if(j!=1 & dim(RESULT_merged)[1] == 0){RESULT_merged <- RESULT}
      rm(list = c('RESULT', 'TPrate', 'ddcor_G1_G0', 'SM'))
    }
    print(RESULT_merged)
    
    }  
  }  
}  

#---------------------------------- save the results
if(dim(RESULT_merged)[1] != 0 )
{
  S0 <- gsub('_AssignedCells.txt','', ERP)
  write.table(RESULT_merged, paste0(OutPutPath, '/', S0, '-', Gene1, '_IDCA.txt'), quote = FALSE, row.names = FALSE, sep = '\t')
} 
#----------------------------------

}
}}else{
i=1; 
SelectedGene <- 'RPL31';
AllGenes<-row.names(estimate)[-which(row.names(estimate) == SelectedGene)]

#---------------------------------- obtain the most representative genes subset using least absolute shrinkage and selection operator (LASSO) 
x <- t(estimate[-which(row.names(estimate) == SelectedGene),])
y <- as.numeric(estimate[which(row.names(estimate) == SelectedGene),])

set.seed(123)
train<-sample(seq(dim(x)[1]),(dim(x)[1]/2), replace=FALSE)

lasso.tr<-glmnet(x[train,],y[train], standardize= FALSE)
pred<-predict(lasso.tr,x[-train,])
rmse<- sqrt(apply((y[-train]-pred)^2,2,mean))
lam.best<-lasso.tr$lambda[order(rmse)[1]]

lasso.tr<-glmnet(x ,y, standardize= FALSE)

Lasso_coefficient <- coef(lasso.tr, s=lam.best)
inds<-which(Lasso_coefficient!=0)
variables<-row.names(Lasso_coefficient)[inds]
variables<-variables[!(variables %in% '(Intercept)')];

length(which(0 != Lasso_coefficient[,1]))
results <- as.data.frame(Lasso_coefficient[which(0 != Lasso_coefficient[,1]),1])
colnames(results) <- "coef"
results <- as.data.frame(results[-1,,FALSE]) 

AllGenes <- row.names(results)
#----------------------------------

RESULT_merged = data.frame()

if(length(AllGenes) > 0)
{
  for(j in 1:length(AllGenes))
  {
    Gene1<- SelectedGene; 
    Gene2<-AllGenes[j];
    
    #---------------------------------- consider each possible gene pair one time
    if( length(considered_GenePairs[which(considered_GenePairs %in% c(paste(Gene1, Gene2, sep = '_'), paste(Gene2, Gene1, sep = '_')))]) > 0)
    {
      print('already is considered.')
    }else{
      
      considered_GenePairs[length(considered_GenePairs)+1] <- paste(Gene1, Gene2, sep = '_')
      considered_GenePairs[length(considered_GenePairs)+1] <- paste(Gene2, Gene1, sep = '_')
      #----------------------------------
      
      estimate_Gene1 <- estimate[which(row.names(estimate) == Gene1),]
      estimate_Gene2 <- estimate[which(row.names(estimate) == Gene2),]
      
      #---------------------------------- Mixture Modeling
      input_mixture_model <- t(rbind(as.matrix(estimate_Gene1), as.matrix(estimate_Gene2)))
      input_mixture_model <- as.data.frame(input_mixture_model, stringsAsFactors = FALSE)
      
      mixed_model <- flexmix(cbind(as.numeric(input_mixture_model[,1]),as.numeric(input_mixture_model[,2]))~1, 
                            k=2, data=input_mixture_model,
                            model = FLXMCmvnorm(diag = TRUE), #diag = TRUE to get fix results FLXMCmvpois(), FLXMCmvnorm
                            control = list(tolerance = 1e-15, iter.max = 1000))
      
      if(length(unique(clusters(mixed_model)))==2)
      {
        #---------------------------------- IDCA
        design_matrix = cbind(colnames(estimate), colnames(estimate))
        design_matrix[-which(design_matrix[,1]%in%S1_BARCODE),1]<-0
        design_matrix[which(design_matrix[,1]%in%S1_BARCODE),1]<-1
        design_matrix[-which(design_matrix[,2]%in%S2_BARCODE),2]<-0
        design_matrix[which(design_matrix[,2]%in%S2_BARCODE),2]<-1
        colnames(design_matrix) <- c("G0", "G1")
        
        design_matrix <- as.data.frame(design_matrix)
        design_matrix[,1] <- as.numeric(design_matrix[,1])
        design_matrix[,2] <- as.numeric(design_matrix[,2])
        design_matrix <- as.matrix(design_matrix)
        ddcor_G1_G0 <- ddcorAll(inputMat = estimate[which(row.names(estimate) %in% c(Gene1,Gene2)),], design = design_matrix,compare = c("G1","G0"),adjust = "fdr", nPerm = 0, corrType = "spearman", splitSet = Gene1, sigThresh = 0.05, sortBy = "pValDiff_adj", verbose = TRUE)
        
        #---------------------------------- compare the results of mixture model and genetic based demultiplexing
        input_mixture_model$real_clusters <- colnames(estimate)
        input_mixture_model$real_clusters[which(input_mixture_model$real_clusters%in%S1_BARCODE)]<-1
        input_mixture_model$real_clusters[which(input_mixture_model$real_clusters%in%S2_BARCODE)]<-2
        
        input_mixture_model$predicted_clusters <- clusters(mixed_model)
        #----------------------------------
        
        #---------------------------------- revise of cluster names
        if(input_mixture_model$predicted_clusters[1] == 2)
        {
          input_mixture_model$predicted_clusters[which(input_mixture_model$predicted_clusters == 1)]<-0
          input_mixture_model$predicted_clusters[which(input_mixture_model$predicted_clusters == 2)]<-1
          input_mixture_model$predicted_clusters[which(input_mixture_model$predicted_clusters == 0)]<-2
        }
        #----------------------------------
        
        #---------------------------------- save the results
        table <- table(input_mixture_model$real_clusters,input_mixture_model$predicted_clusters)
        TPrate <- (table[1,1] + table[2,2])/length(clusters(mixed_model))
        SM <- summary(mixed_model)
        
        RESULT <- data.frame()
        RESULT <- ddcor_G1_G0
        RESULT$AIC <- SM@AIC
        RESULT$TPrate <- TPrate
        RESULT <- RESULT[which(RESULT$pValDiff_adj < 1e-10 ),]# & RESULT$G0_pVal != 0.0
        
        if(dim(RESULT)[1]>0)
        {
          if(j==1){RESULT_merged = RESULT}
          if(j!=1 & dim(RESULT_merged)[1] != 0){RESULT_merged <- rbind(RESULT_merged, RESULT)}
          if(j!=1 & dim(RESULT_merged)[1] == 0){RESULT_merged <- RESULT}
          rm(list = c('RESULT', 'TPrate', 'ddcor_G1_G0', 'SM'))
        }
        print(RESULT_merged)
        
      }  
    }  
  }  
  
  #---------------------------------- save the results
  if(dim(RESULT_merged)[1] != 0 )
  {
    S0 <- gsub('_AssignedCells.txt','', ERP)
    write.table(RESULT_merged, paste0(OutPutPath, '/', S0, '-', Gene1, '_IDCA.txt'), quote = FALSE, row.names = FALSE, sep = '\t')
  } 
  #----------------------------------
  
}

}

cat(paste0("\033[0;", 47, "m", "You can find the results in: ", "\033[0m","\n", InputDir, "/IDCA_Analysis/"))

}
