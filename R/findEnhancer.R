
#' \code{findEnhancer} is the main function that identifies enhancers for a given gene
#' @param gene gene of interest eg. Runx1 (case sensitive)
#' @param expression data.frame/matrix of expression with cell types as cols and genes as rows
#' @param regulation_signal data.frame/matrix of enhancer signals (active histone marks/DNaseI-seq/ATAC-seq etc) with cell types as cols and putative enhancer regions as rows.
#' @param regulation_tfbs data.frame/matrix with TF ChIP-seq data as cols and putative enhancer regions as rows. This is a binary matrix with 0 = no binding and 1 = binding
#' @param region_gene_mapping data.frame with region to gene mapping. Please refer to data(region_gene_mapping) for the format.
#' @param min_tfb_events the minimum number of TFB events per site to be consider for community CRE calculation. default = 2
#' @param coCRE_corr_cutoff the cutoff above which two putative enhancer regions (CREs) are considered to be correlated. Default = 0.5
#' @param coCRE_cutoff the cutoff beyond which the regions will not be considered for coCRE calculation. Default = 10000 (100kB)
#' @param singleton_cutoff the cutoff for considering singleton CRE. Default = 20000 (20Kb)
#' @param alphaVal alpha value for glmnet. Read glmnet manual
#' @param scale.predictors (TRUE/FALSE). Default = TRUE
#' @param family response family. Default = gaussian
#' @param nfoldxval n-fold cross validation for lmin/l1se. Default = 10 (for leave out cross validation this is equal to total number of cols in expression data)
#' @return list of regions and p-values
#' @export
#' @import glmnet
#' @import igraph
#' @import arules
#' @import covTest
#' @import selectiveInference
#' @import lars

findEnhancer <- function(gene, expression, regulation_signal, regulation_tfbs, region_gene_mapping, min_tfb_events = 2, coCRE_corr_cutoff = 0.5, coCRE_cutoff = 100000, singleton_cutoff = 20000, alphaVal = 1, scale.predictors = TRUE, family = "gaussian", nfoldxval = 10, ...){
    proceed_fwd = 1
    cat("Performing basic checks\n")
    celltypes <- colnames(expression)
    temp <- celltypes [ ! celltypes %in% colnames(regulation_signal) ]
    if (length(temp) > 0){
        message("There are extra cell types in your regulation signal than your expression signal. Please equat them")
        proceed_fwd <- 0
        return(NULL)
    }
    cat("Getting gene expression","\n")
    gene_expression <- t(expression [ "Runx1", ])

    cat("Getting regions within mapping distance of gene","\n")
    gene_meta <- getPeaks(region_gene_mapping, gene, coCRE_cutoff)
    cat("Total regions found = ", nrow(gene_meta),"\n")
    if(nrow(gene_meta) <= 1){
        proceed_fwd = 0
        message("Total regions should be at least 2. Aborting ...")
        return(NULL)
    }
    gene_region_ids <- unique(as.character(gene_meta$peak))
    cat("Getting the regions' regulation signal data","\n")
    gene_reg_signal <- getRegulation(regulation_signal, gene_region_ids)

    if(nrow(gene_reg_signal) != length(gene_region_ids)){
        message("Looks like some of your regions do not have signal data. Double check")
        proceed_fwd = 0
        return(NULL)
    }

    cat("Getting the regions' TFBS data","\n")
    gene_reg_tfbs <- getRegulation(regulation_tfbs, gene_region_ids)
    if(nrow(gene_reg_tfbs) != length(gene_region_ids)){
        message("Looks like some of your regions do not have TF binding data. If none of your tfs bind to some regions, add in them as a all zero row. You can use min_tfb_events to filter off low occupancy sites")
        proceed_fwd = 0
        return(NULL)
    }

    cat("Computing coCREs","\n")
    gene_coCRE_signal <- get_coCRE(gene_reg_tfbs, gene_reg_signal, gene_meta, min_tfb_events, coCRE_corr_cutoff, singleton_cutoff)
    if(nrow(gene_coCRE_signal) <= 1){
        message("After performing coCRE it seems like you are left with only one predictor. Cannot use LASSO here. Aborting")
        proceed_fwd = 0
        return(gene_coCRE_signal)
    }
    cat("Running lasso","\n")
    lasso_outpt <- run_lasso(gene_coCRE_signal, gene_expression, ltype= "l1se", alphaVal = alphaVal, scale.predictors = scale.predictors, family = family, nfold = nfoldxval)
    cat("Collating data\n")
    outpt = list()
    outpt[[gene]][["expression"]] <- gene_expression
    outpt[[gene]][["meta"]] <- gene_meta
    outpt[[gene]][["reg_signal"]] <- gene_reg_signal
    outpt[[gene]][["reg_tfbs"]] <- gene_reg_tfbs
    outpt[[gene]][["lambda_1se"]] <- lasso_outpt[["l1se_pred"]]
    outpt[[gene]][["lambda_min"]] <- lasso_outpt[["lmin_pred"]]
    outpt[[gene]][["p"]] <- lasso_outpt[["p"]]
    return(outpt)
}


#' Expression data:
#'
#' rownames should be genes
#'
#' colnames should be cell types
#'
#' the gene name and cell types are both case-sensitive
#'
#' they should be used consistently across other input data as well.
#' @name expression
#' @docType data
#' @keywords expression
NULL

#' Dataframe that maps regions to genes:
#'
#' All the values are case-sensitive
#'
#' Colnames are case sensitive and should be preserved as shown.
#'
#' chr: chromosome id
#'
#' start: start of region
#'
#' end: end of region
#'
#' peak: name of the region.
#'
#' gene: name of gene.
#'
#' distance: shortest dist between region and peak. If the peak is within a gene, you can use any value that is less than "singleton_cutoff"
#'
#' colnames should be cell types
#'
#' the gene name and cell types are both case-sensitive
#'
#' they should be used consistently across other input data as well.
#' @name region_gene_mapping
#' @docType data
#' @keywords metadata
NULL

#' Regulation Signal
#'
#' DataFrame with row.names as regions and col.names as cell types
#'
#' Region names should be preserved in region_gene_mapping and regulation_tfbs
#'
#' Cell types names should be preserved in regulation_tfbs and expression DFs
#'
#' @name regulation_signal
#' @docType data
NULL

#' Regulation TFBS
#'
#' DataFrame with row.names as regions and col.names as TF ChIP-seq.
#'
#' A binary data, with 0 denoting absence of binding and 1 denoting binding
#'
#' Region names should be preserved in region_gene_mapping and regulation_signal
#'
#' Cell types names should be preserved in regulation_signal and expression DFs
#'
#' @name regulation_tfbs
#' @docType data
NULL




