
#############################################################################################
########################## Single Cell data Quality Control #################################
#############################################################################################

#'@import dplyr
#'@import SingleCellExperiment
#'@import scater
#'@import scDblFinder
#'@importFrom data.table fread 
#'@export
#'@name scQC
#'@title Single Cell data Quality Control
#'@description Identify potential doublet cells based on the local density of simulated doublet expression profiles using scDblFinder package.
#'@author {Isar Nassiri}
#'@param InputDir
#'#'A folder including gene-cell count matrix (read_count.csv) and filtered feature-barcode matrices from 10X CellRanger count (barcodes.tsv.gz)
#'@return You can find the results in the QC/ folder under the title of 'Cells_passed_QC_noDoublet.txt'.
#'@examples
#'library("scDIV")
#'csQCEAdir <- system.file("extdata", package = "scDIV")
#'scQC( InputDir = csQCEAdir ) 
#'@export

#----------------------------------
# library(data.table)
# library(dplyr)
# library(SingleCellExperiment)
# library(scater)
# library(scDblFinder)
# InputDir = '/Users/isar/Documents/Expression_aware_demultiplexing/Inputs/'
# scQC(InputDir)
#----------------------------------

scQC <- NULL
scQC <- function( InputDir )
{

setwd(InputDir)
dir.create('QC/')

#---------------------------------- read in gene expression profile
eislet = fread('read_count.csv', stringsAsFactors = F,header = T)
eislet = as.data.frame(eislet)
row.names(eislet) = eislet$V1
eislet = eislet[,-1]
#----------------------------------

#---------------------------------- keep expressed genes
eislet = eislet[apply(eislet, 1, function(x) !all(x==0)),]
counts.fname = as.matrix(eislet)
#----------------------------------

#---------------------------------- add gene symbols
genesID_names = fread(paste0('features.tsv.gz'), stringsAsFactors = F, header = F)
genesID_names = as.data.frame(genesID_names)
colnames(genesID_names) = c('Ensembl','Symbol','Annotation')

temp = data.frame(Ensembl = row.names(counts.fname))
genesID_names = left_join(temp, genesID_names, by='Ensembl')
#----------------------------------

#---------------------------------- convert matrix to SingleCellExperiment object
sce <- SingleCellExperiment(
list(counts=as(counts.fname, "dgCMatrix")), 
rowData=rownames(counts.fname), 
colData=DataFrame(Barcode=colnames(counts.fname))
)

rowData(sce)$Ensembl = genesID_names$Ensembl
rowData(sce)$Symbol = genesID_names$Symbol

rownames(sce) <- uniquifyFeatureNames(rowData(sce)$Ensembl, rowData(sce)$Symbol)
colnames(sce) <- sce$Barcode
#----------------------------------

#---------------------------------- Doublet detection by simulation
set.seed(100)

dbl.dens <- computeDoubletDensity(sce) 
sce$DoubletScore <- dbl.dens
sce = sce[,dbl.dens < 2]
#----------------------------------

#---------------------------------- Save results
write.table(colnames((assay(sce))), 'QC/Cells_passed_QC_noDoublet.txt', quote = F, row.names = F, col.names = F)
#----------------------------------

cat(paste0("\033[0;", 47, "m", "You can find the results in: ", "\033[0m","\n", InputDir, "/QC/Cells_passed_QC_noDoublet.txt"))

}

