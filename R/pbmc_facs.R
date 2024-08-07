#' @name pbmc_facs
#'
#' @title Mixture of 10 FACS-purified PBMC Single-Cell RNA-seq data
#'
#' @docType data
#' 
#' @description These data are a selection of the reference
#'   transcriptome profiles generated via single-cell RNA sequencing
#'   (RNA-seq) of 10 bead-enriched subpopulations of PBMCs (Donor A),
#'   described in Zheng \emph{et al} (2017). The data are unique
#'   molecular identifier (UMI) counts for 16,791 genes in 3,774 cells.
#'   (Genes with no expression in any of the cells were removed.) Since
#'   the majority of the UMI counts are zero, they are efficiently
#'   stored as a 3,774 x 16,791 sparse matrix. These data are used in
#'   the vignette illustrating how 'fastTopics' can be used to analyze to
#'   single-cell RNA-seq data. Data for a separate set of 1,000 cells is
#'   provided as a \dQuote{test set} to evaluate out-of-sample predictions.
#'
#' @format \code{pbmc_facs} is a list with the following elements:
#' 
#' \describe{
#'
#'   \item{counts}{3,774 x 16,791 sparse matrix of UMI counts, with
#'      rows corresponding to samples (cells) and columns corresponding to
#'      genes. It is an object of class \code{"dgCMatrix"}).}
#'
#'   \item{counts_test}{UMI counts for an additional test set of 100
#'     cells.}
#'
#'   \item{samples}{Data frame containing information about the
#'     samples, including cell barcode and source FACS population
#'     (\dQuote{celltype} and \dQuote{facs_subpop}).}
#'
#'   \item{samples_test}{Sample information for the additional test
#'      set of 100 cells.}
#' 
#'   \item{genes}{Data frame containing information and the genes,
#'     including gene symbol and Ensembl identifier.}
#'
#'   \item{fit}{Poisson non-negative matrix factorization (NMF) fitted
#'     to the UMI count data \code{counts}, with rank \code{k = 6}. See
#'     the vignette how the Poisson NMF model fitting was performed.}}
#'
#' \url{https://www.10xgenomics.com/resources/datasets}
#' 
#' @references
#' G. X. Y. Zheng \emph{et al} (2017). Massively parallel digital
#' transcriptional profiling of single cells. \emph{Nature Communications}
#' \bold{8}, 14049. \doi{10.1038/ncomms14049}
#' 
#' @keywords data
#'
#' @examples
#' library(Matrix)
#' data(pbmc_facs)
#' cat(sprintf("Number of cells: %d\n",nrow(pbmc_facs$counts)))
#' cat(sprintf("Number of genes: %d\n",ncol(pbmc_facs$counts)))
#' cat(sprintf("Proportion of counts that are non-zero: %0.1f%%.\n",
#'             100*mean(pbmc_facs$counts > 0)))
#' 
NULL
