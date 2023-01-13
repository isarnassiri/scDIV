# scDIV
scDIV: Single Cell RNA Sequencing Data Demultiplexing using Interindividual Variations

## **scDIV**

### Introduction 
This documentation gives an introduction and usage manual of scDIV (acronym of the Single Cell RNA sequencing data Demultiplexing using Interindividual Variations) an R package to use inter-individual differential co-expression patterns for demultiplexing the pooled samples without any extra experimental steps. 
<br />
 
Raw count data from 10X CellRanger (outs/read_count.csv) or other single-cell experiments has the gene as a row (the gene name should be the human or mouse Ensembl gene ID) and the cell as a column. You can convert an HDF5 Feature-Barcode Matrix [(LINK)](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/h5_matrices) to a gene-cell count matrix using the cellranger mat2csv [(LINK)](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/matrices#mat2csv) command provided by 10Xgenomics. The cells in the `read_count.csv` file are from the filtered feature-barcode matrix generated by cell ranger.

