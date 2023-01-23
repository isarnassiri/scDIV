#############################################################################################
########################## Gene Expression Recovery #################################
#############################################################################################

#'@import SAVER
#'@import stringr
#'@importFrom data.table
#'@export
#'@name GeneExpressionRecovery
#'@title Gene Expression Recovery
#'@description This function uses SAVER (single-cell analysis via expression recovery), an expression recovery method for unique molecule index (UMI)-based scRNA-seq data to provide accurate expression estimates for all genes in a scRNA-seq profile.
#'@author {Isar Nassiri}
#'@param InputDir
#'#'A folder including gene-cell count matrix (read_count.csv), filtered feature-barcode matrices from 10X CellRanger count (barcodes.tsv.gz), list of QC passed barcodes (Cells_passed_QC_noDoublet.txt) in QC/ folder, and output of genetic demultiplexing of the sample pool using vireo (donor_ids.tsv).
#'@param Donors
#' A pair of selected donors (e.g., donor6_donor2)
#'@param FC
#' A name for flow cell (e.g., FAI5649A17).
#'#'@return You can find the results in the SAVER/ folder with 'AssignedCells.txt' and 'AllCells.txt' extensions.
#'@examples
#'library("scDIV")
#'csQCEAdir <- system.file("extdata", package = "scDIV")
#'Donors='donor6_donor2'
#'FC='FAI5649A17'
#'GeneExpressionRecovery( InputDir = csQCEAdir, Donors = Donors, FC = FC )
#'@export

#----------------------------------
# library(data.table)
# library(SAVER)
# library(stringr)
# csQCEAdir = '/Users/isar/Documents/Expression_aware_demultiplexing/Inputs'
# Donors='donor6_donor2'
# FC='FAI5649A17'
# GeneExpressionRecovery( InputDir = csQCEAdir, Donors = Donors, FC = FC )
#----------------------------------

GeneExpressionRecovery <- NULL
GeneExpressionRecovery <- function( InputDir, Donors, FC )
{

#---------------------------------- Gene expression recovery for single-cell RNA sequencing
OutPutPath = 'SAVER'
dir.create(OutPutPath, recursive = T)

#---------------------------------- read in gene expression profile
read_count = fread('read_count.csv', stringsAsFactors = F, header = T)
read_count = as.data.frame(read_count)
colnames(read_count)[1] = 'GeneID'

#---------------------------------- add gene symbols
GTF = fread('features.tsv.gz', stringsAsFactors = F, header = F)
GTF = as.data.frame(GTF)
colnames(GTF) = c('GeneID','gene_name','Annotation')

estimate = merge(GTF,read_count,by='GeneID')
row.names(estimate) = make.names(estimate$gene_name, unique = T)
estimate = estimate[,-c(1:3)]

#---------------------------------- keep expressed genes
estimate = estimate[apply(estimate, 1, function(x) !all(x==0)),]
#----------------------------------

#---------------------------------- subsed gene expression profile using the list of QC passed cells
Cells_passed_QC = fread('QC/Cells_passed_QC_noDoublet.txt', stringsAsFactors = F, header = F)
Cells_passed_QC = as.data.frame(Cells_passed_QC)
colnames(Cells_passed_QC)[1] = 'BARCODE'
estimate = estimate[,which(colnames(estimate) %in% Cells_passed_QC$BARCODE)]
#----------------------------------

#---------------------------------- read in the results of genetic based demultiplexing (vireo)
vireo = fread('donor_ids.tsv', stringsAsFactors = F, header = T)
vireo = as.data.frame(vireo)
colnames(vireo)[1] = 'BARCODE'
#----------------------------------

#---------------------------------- list of donors
sp = str_split_fixed(Donors, "_", 2)

D1=sp[1,1]
D2=sp[1,2]

#----------------------------------

#---------------------------------- gene expression profile of assigned cells by genetic demultiplexing per pair of donors
S1 = vireo[which(vireo$donor_id == D1),]
S2 = vireo[which(vireo$donor_id == D2),]

read_count_S1S2 = estimate[,which(colnames(estimate) %in% c(S1$BARCODE,S2$BARCODE))]
read_count_S1S2 = read_count_S1S2[rowSums(read_count_S1S2[])>as.numeric(summary(rowSums(read_count_S1S2[]))[4]),] # > Median [3], mean [4]
#----------------------------------

#---------------------------------- run SAVER (Single-cell Analysis Via Expression Recovery)
set.seed(123)
estimate_S1S2 <- saver(read_count_S1S2, ncores = 6, estimates.only = TRUE)
estimate_S1S2 = data.frame(GeneID=row.names(estimate_S1S2), estimate_S1S2)
#----------------------------------

#---------------------------------- save the transformed gene expression profile
fwrite(estimate_S1S2, paste0(OutPutPath, '/', D1, '_', D2, '_', FC, '_', 'AssignedCells.txt'), row.names = F, quote = F, sep = '\t')
#----------------------------------

#---------------------------------- gene expression profile of all cells per pair of donors
vireo_subset = vireo[which(vireo$donor_id %in% c('doublet', 'unassigned')),]
S3 = vireo_subset[which(vireo_subset$best_singlet %in% c(D1, D2)),] # sometime best_singlet does not exist in best_doublet
read_count_S3 = estimate[,which(colnames(estimate) %in% S3$BARCODE)]

read_count_S3 = read_count_S3[which(row.names(read_count_S3) %in% row.names(read_count_S1S2)),]
# identical(row.names(read_count_S3), row.names(read_count_S1S2))

read_count_S1S2S3 = cbind(read_count_S1S2, read_count_S3)
#----------------------------------

#---------------------------------- run SAVER (Single-cell Analysis Via Expression Recovery)
estimate_S1S2S3 <- saver(read_count_S1S2S3, ncores = 6, estimates.only = TRUE)
estimate_S1S2S3 = data.frame(GeneID=row.names(estimate_S1S2S3), estimate_S1S2S3)

#---------------------------------- save the transformed gene expression profile
fwrite(estimate_S1S2S3, paste0(OutPutPath, '/', D1, '_', D2, '_', FC, '_', 'AllCells.txt'), row.names = F, quote = F, sep = '\t')
#----------------------------------

cat(paste0("\033[0;", 47, "m", "You can find the results in: ", "\033[0m","\n", csQCEAdir, "/SAVER/"))

}

