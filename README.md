# scDIV
scDIV: Single Cell RNA Sequencing Data Demultiplexing using Interindividual Variations

## **scDIV**

### Introduction 
This documentation gives an introduction and usage manual of scDIV (acronym of the Single Cell RNA sequencing data Demultiplexing using Interindividual Variations) an R package to use inter-individual differential co-expression patterns for demultiplexing the pooled samples without any extra experimental steps. 
<br />

### Infer genetic variants from scRNA-seq data 
cellsnp-lite is used to pileup the expressed alleles in single-cell data, which can be directly used for donor deconvolution in multiplexed single-cell RNA-seq data, which assigns cells to donors without genotyping reference [(LINK)](https://github.com/single-cell-genetics/cellsnp-lite). 

cellsnp-lite gets bam file and list of barcodes as variable inputs, a Variant Call Format (vcf) file listing all candidate SNPs (regionsVCF) as backend input variable:

```{r,eval=FALSE}
cellsnp-lite -s possorted_genome_bam.bam -b barcodes.tsv.gz -O FOLDER-NAME -R regionsVCF -p 22 --minMAF 0.05 --minCOUNT 10 --gzip 
```

cellsnp-lite generates a vcf file including called genetic variants as follows:

| ![Figure 1](/cellsnp-lite.png) | 
|:--:| 
| *Figure 1. Example of vcf file generated by cellsnp-lite.* |

### Demultiplex pooled samples 
We use Vireo (Variational Inference for Reconstructing Ensemble Origin) for donor deconvolution using expressed SNPs in multiplexed scRNA-seq data [(LINK)](https://vireosnp.readthedocs.io/en/latest/).

Vireo gets s variants info file provided by cellsnp-lite as an input: 

```{r,eval=FALSE}
vireo -c input-vcf-file -o output-folder --randSeed 2 -N Number-of-donors -t GP  
```

We use "donor_ids.tsv" file from outputs of Vireo for downstream analysis:

| ![Figure 2](/vireo.png) | 
|:--:| 
| *Figure 2. Example of output from Vireo.* |


### Generate gene-cell count matrices for all possible pairs of individuals
Raw count data from 10X CellRanger (outs/read_count.csv) or other single-cell experiments has the gene as a row (the gene name should be the human or mouse Ensembl gene ID) and the cell as a column. You can convert an HDF5 Feature-Barcode Matrix [(LINK)](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/h5_matrices) to a gene-cell count matrix using the cellranger mat2csv [(LINK)](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/matrices#mat2csv) command provided by 10Xgenomics. The cells in the `read_count.csv` file are from the filtered feature-barcode matrix generated by cell ranger.

A filtered feature-barcode matrix generated by the cell ranger can be converted from HDF5 feature-barcode matrix to a gene-cell count matrix using the cellranger mat2csv (command provided by 10Xgenomics) as follows:

```{r,eval=FALSE}
cellranger mat2csv filtered_feature_bc_matrix.h5 read_count.csv
```





