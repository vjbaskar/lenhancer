# lenhancer: [LASSO](https://en.wikipedia.org/wiki/Lasso_(statistics)) based Enhancer identification 

## Introduction
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

## Installation

```R
install.packages("devtools")
devtools::install_github("hadley/devtools")
devtools::install_github(c("vjbaskar/covTest", "vjbaskar/lenhancer"))
```

## Example run

The below example run is for gene "Runx1". Please not the format of the four input data below and keep them the same for your input as well. 


```R
#### Input variables
gene = "Runx1"
coCRE_cutoff = 100000 # 100Kb upper limit for considering regions as being mapped to the gene
singleton_cutoff = 20000 # 20Kb upper limit for considering regions that are not coCREs in the model
coCRE_corr_cutoff = 0.5 # The higher the value the tighter the community of CREs
min_tfb_events = 2 # lower limit for considering a region for walktrap community clustering
alphaVal=1 # 1 is the lasso penalty, and 0 the ridge penalty.
scale.predictors=TRUE # scale predictors True or False
family = "gaussian" # "gaussian","binomial","poisson","multinomial","cox","mgaussian"
nfold = 10 # nfold cross-validation: >= 3 and can be as large as the sample (# of cell types) size (leave-one-out CV) 


data("expression")
data("regulation_signal")
data("regulation_tfbs")
data("region_gene_mapping")
gene_preds = findEnhancer(gene, expression, regulation_signal, regulation_tfbs, region_gene_mapping, min_tfb_events = 2, coCRE_corr_cutoff = 0.5, singleton_cutoff = 20000, alphaVal = 1, scale.predictors = TRUE, family = "gaussian", nfoldxval = 10)
```
