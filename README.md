# lenhancer
Enhancer identification using [LASSO](https://en.wikipedia.org/wiki/Lasso_(statistics)). 

Combining a novel strategy to identify communities of related control elements with a penalized regression approach, we developed individual gene-by-gene models to identify the potential control elements predictive of the expression of the gene. 

The package requires the following data:
* expression (E): gene expression data in say *n* cell types
* regions (R): A set of putative enhancer CREs with a given distance of the genes
* regulation signals (S): regulation signals such as active chromatin marks, DNaseI-seq or ATAC-seq in the *n* cell types. You can normalise the way you prefer for library sizes, and region sizes.
* tfbs (T): Transcription factor binding signals in these regions. You can use as many TF binding data as you want, but TF ChIP-seq in the *n* cell types is highly recommended.

For each gene, the method consists of the following steps:
* *m* = consider peaks that are with 100kB (can be changed)
* Using correlation of *T* and *S* between these *m* regions, a community clustering algorithm is used to group the CREs into clusters called coCREs.
* singleton CREs within 10Kb and coCREs are then considered as predictors,  and the gene's expression across cell types as response variable.
* Penalised regression (L1 regularisation, LASSO) along with cross-validation is used to identify enhancers for the gene.
* [covariance testing](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4285373/) is used to compute the p-value. 

