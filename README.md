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

## Input data types

### expression

*expression* data consists of gene expression values (such as *log<sub>2</sub>(FPKM)*). gene as rows, cell types as columns.
```R
> data(expression)
> head(expression)
             ESC      CESC       MES      CMES        CM        CP        HB        HE        HP
Sergef  1.137991  1.044118  1.895715  2.217683  1.666260  0.646738  2.008090  1.244708  1.802633
Bcl7a   2.206694  2.550776  2.844359  3.496080  3.116237  3.478903  2.129927  1.403524  2.419382
Lnx2    2.053856  2.500270  1.863102  1.802951  1.930653  1.332876  2.289595  2.377188  2.195265
Ppia   10.074714 10.342072  9.826837 10.160430  9.599557  9.509954  1.916361  1.533936  2.086379
Gkn1   -9.965784 -9.965784 -9.965784 -9.965784 -9.965784 -9.965784 -9.965784 -9.965784 -9.965784
Lrat   -6.122783 -6.090128 -6.229214 -7.160255 -6.348240 -7.412639 -8.432611 -5.641031 -6.950677
              MAC
Sergef  2.6563687
Bcl7a   0.2551161
Lnx2    2.8721921
Ppia    0.8198570
Gkn1   -9.9657843
Lrat   -5.6209541
```

### regions
*region_gene_mapping* data frame consists of all the regulatory regions that are mapped to the genes in *expression*.
There are many software such as [ChIPseeker](http://bioconductor.org/packages/release/bioc/html/ChIPseeker.html), [ChIPpeakAnno](http://bioconductor.org/packages/release/bioc/html/ChIPpeakAnno.html) and [HOMER](http://homer.ucsd.edu/homer/index.html) that can map regions to genes. If you map peaks to genes, and the peak lies within the gene you can mark their *distance* as any value less than *singleton_cutoff* (default, 20Kb)

```R
> data("region_gene_mapping")
> head(region_gene_mapping)
       chr    start      end         peak  gene distance
7786 chr16 92514276 92515568 merge.105629 Runx1    85898
7792 chr16 92515907 92516390 merge.105630 Runx1    85076
7797 chr16 92522863 92525592 merge.105632 Runx1    75874
7815 chr16 92572528 92574113 merge.105643 Runx1    27353
7820 chr16 92579265 92579735 merge.105644 Runx1    21731
7824 chr16 92605464 92606424 merge.105650 Runx1     3000
```
*NOTE*: Please keep the colnames the same. *findEnhancer* function looks for `c("peak", "gene", "distance")` columns specifically in its codes and is case-senstive.


### regulation signal

*regulation_signal* data consists of regions as rows and chromatin activity signal (such as *log<sub>2</sub>(FPKM)*) for each cell type as columns. 

*NOTE* that the colnames for *expression* and *regulation_signal* are the same.

```R
> data("regulation_signal")
> head(regulation_signal)
                    ESC       CESC       MES      CMES        CM        CP        HB         HE
merge.105629 2.56820271  2.5806279 0.9168354 1.5256741 0.4129932 0.5330959 1.6784691 -0.1319497
merge.105630 2.44097491  1.6429856 0.5797611 2.1639965 1.5515023 1.7488139 3.1067132  0.7631582
merge.105632 1.30810596  0.7371804 0.9214790 1.5322401 1.1576157 1.1091539 1.2070316  0.7889895
merge.105643 0.02400262 -0.3086307 0.6754163 0.2950902 1.3850548 0.7898872 0.7038305  3.2788514
merge.105644 1.24260338  1.2721741 1.9149712 2.2436335 2.2543188 2.6906041 3.7203829  2.8612361
merge.105650 2.54822747  3.0395300 3.8857761 3.0925726 2.1348828 1.5124602 4.0995385  2.2658194
                    HP        MAC
merge.105629 0.9904566  0.8144609
merge.105630 1.6498018  2.3616716
merge.105632 1.2059928  4.0246107
merge.105643 1.5171622  5.4947875
merge.105644 3.3901037  0.6582151
merge.105650 1.3444846 -3.3652031
```
