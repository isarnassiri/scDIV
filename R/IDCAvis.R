#############################################################################################
############################### visualization of IDCA outputs ###############################
#############################################################################################

#'@import stringr
#'@import data.table
#'@import ggplot2
#'@export
#'@name IDCAvis
#'@title Visualization of IDCA outputs of Inter-individual Differential gene Correlation Analysis (IDCA) results
#'@description The function visualizes the outputs of IDC analysis.
#'@author {Isar Nassiri}
#'@param InputDir
#'A folder including matrix of cell-gene expression after Gene Expression Recovery in the SAVER/ folder with 'AssignedCells.txt' extension, output of genetic demultiplexing of the sample pool using vireo (donor_ids.tsv), and output of IDCA function with extension of "_IDCA.txt" in IDCA_Analysis/ folder.
#'@return You can find the results in the IDCA_Analysis/IDCA_Plots/ as pdf file(s).
#'@examples
#'library("scDIV")
#'InputDir = system.file("extdata", package = "scDIV")
#'IDCAvis( InputDir )
#'@export

#----------------------------------
# library(data.table)
# library(stringr)
# library(ggplot2)
# InputDir='/Users/isarnassiri/Documents/OGC/DCOX_scRNAseq_saver/Package_scripts_inputs/Rpackage/Inputs/'
# IDCAvis( InputDir )
#----------------------------------

IDCAvis <- NULL
IDCAvis <- function(InputDir)
{
  #---------------------------------- find and read input file(s)
  setwd(InputDir)
  
  listfiles <-
    list.files(paste0(InputDir, '/IDCA_Analysis/'), pattern = '_IDCA.txt')
  Pairs <- unique(gsub('-.*', '', listfiles))
  
  for (i in 1:length(Pairs))
  {
    listfiles_Pair <-
      list.files(paste0(InputDir, '/IDCA_Analysis/'), pattern = "Pairs[i]|.txt")
    
    setwd(paste0(InputDir, '/IDCA_Analysis/'))
    merged_Pair <- rbindlist(lapply(listfiles_Pair, fread))
    merged_Pair <- as.data.frame(merged_Pair)
    
    merged_Pair <-
      merged_Pair[order(-as.numeric(merged_Pair$TPrate),
                        merged_Pair$pValDiff_adj), ]
    
    dir.create('IDCA_Plots/')
    fwrite(
      merged_Pair,
      paste0('IDCA_Plots/', Pairs[i], '_merged.txt'),
      quote = FALSE,
      row.names = FALSE,
      sep = '\t'
    )
    
    #---------------------------------- read in gene expression profile
    estimate <-
      fread(
        paste0(InputDir, '/SAVER/', Pairs[i], '_AssignedCells.txt'),
        stringsAsFactors = FALSE,
        header = TRUE
      )
    estimate <- as.data.frame(estimate)
    row.names(estimate) <-
      make.names(estimate$GeneID, unique = TRUE)
    estimate <- estimate[, -c(1)]
    #----------------------------------
    
    #---------------------------------- generate names for outputs
    S0 <- gsub('_AssignedCells.txt', '', Pairs[i])
    sp <- str_split_fixed(S0, "_", 3)
    
    S1 <- sp[1, 1]
    S2 <- sp[1, 2]
    FC <- sp[1, 3]
    #----------------------------------
    
    #---------------------------------- read in the results of genetic demultiplexing (vireo)
    vireo <-
      fread(
        paste0(InputDir, '/donor_ids.tsv'),
        stringsAsFactors = FALSE,
        header = TRUE
      )
    vireo <- as.data.frame(vireo)
    colnames(vireo)[1] <- 'BARCODE'
    
    S1_BARCODE <- vireo[which(vireo$donor_id == S1), 'BARCODE']
    S2_BARCODE <- vireo[which(vireo$donor_id == S2), 'BARCODE']
    
    S1_BARCODE <- gsub('-', '.', S1_BARCODE)
    S2_BARCODE <- gsub('-', '.', S2_BARCODE)
    #----------------------------------
    
    #---------------------------------- subset the gene expression profile
    estimate <-
      estimate[, which(colnames(estimate) %in% c(S1_BARCODE, S2_BARCODE))]
    #----------------------------------
    
    #---------------------------------- visualization
    Gene1 <-
      merged_Pair$Gene2[1] # in RESULT_merged, Gene2 is the original Gene1
    Gene2 <- merged_Pair$Gene1[1]
    
    estimate_Gene1 <- estimate[which(row.names(estimate) == Gene1), ]
    estimate_Gene2 <- estimate[which(row.names(estimate) == Gene2), ]
    
    input <-
      t(rbind(as.matrix(estimate_Gene1), as.matrix(estimate_Gene2)))
    input <- as.data.frame(input, stringsAsFactors = FALSE)
    
    input$real_clusters <- colnames(estimate)
    input$real_clusters[which(input$real_clusters %in% S1_BARCODE)] <-
      1
    input$real_clusters[which(input$real_clusters %in% S2_BARCODE)] <-
      2
    
    S1 <- unique(input$real_clusters)[1]
    S2 <- unique(input$real_clusters)[2]
    
    expr0 <-
      data.frame(Gene1 = as.numeric(input[which(input[, 'real_clusters'] == S1), 1]), Gene2 =
                   as.numeric(input[which(input[, 'real_clusters'] == S1), 2]))
    expr1 <-
      data.frame(Gene1 = as.numeric(input[which(input[, 'real_clusters'] == S2), 1]), Gene2 =
                   as.numeric(input[which(input[, 'real_clusters'] == S2), 2]))
    
    expr0$Sample <- rep(paste0('Sample-0'), dim(expr0)[1])
    expr1$Sample <- rep(paste0('Sample-1'), dim(expr1)[1])
    
    input.plot <- do.call(rbind, list(expr0, expr1))
    input.plot$Sample <-
      factor(input.plot$Sample, levels = unique(input.plot$Sample))
    
    theme_set(theme_bw())
    g <-
      ggplot(input.plot, aes(Gene1, Gene2)) + labs(
        subtitle = "Sample-specific correlation for expression values",
        title = "",
        x = Gene1,
        y = Gene2
      )
    p <-
      g + geom_jitter(aes(col = Sample, size = Gene1)) + geom_smooth(aes(col =
                                                                           Sample), method = "lm", se = FALSE) + scale_color_manual(values = c("#999999", "#E69F00", "#56B4E9")) +
      geom_rug(aes(color = Sample)) + theme(
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.subtitle = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18)
      )
    
    #---------------------------------- save the image
    pdf(
      file = paste0(
        paste0(InputDir, '/IDCA_Analysis/IDCA_Plots/'),
        Gene1,
        '_',
        Gene2,
        '-',
        Pairs[i],
        "_RealCluster_correlationplot.pdf"
      ),
      width = 15,
      height = 10,
      useDingbats = FALSE
    )
    print(p)
    dev.off()
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
      "/IDCA_Analysis/IDCA_Plots/"
    )
  )
  
}