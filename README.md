# lenhancer
Enhancer identification using [LASSO](https://en.wikipedia.org/wiki/Lasso_(statistics)). 

Combining a novel strategy to identify communities of related control elements with a penalized regression approach, we developed individual gene-by-gene models to identify the potential control elements predictive of the expression of the gene. 

The package requires the following data:
* expression: gene expression data in say *n* cell types
* regions: A set of putative enhancer regions with a given distance of the genes
* regulation signals: regulation signals such as active chromatin marks, DNaseI-seq or ATAC-seq in the *n* cell types. You can normalise the way you prefer for library sizes, and region sizes.
* tfbs: Transcription factor binding signals in these regions. You can use as many 

The method consists of the following steps:


