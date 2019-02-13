#' Calculate percentage of transcripts of gene list.
#'
#' Get percentage of transcripts of gene list compared to all transcripts per cell.
#' @param object Seurat object.
#' @param genes List(s) of genes.
#' @keywords seurat cerebro
#' @export
#' @import Seurat
#' @examples
#' calculatePercentGenes(object = seurat, genes = gene_list)
calculatePercentGenes <- function(
  object,
  genes
) {
  ##--------------------------------------------------------------------------##
  ## get for every supplied gene list, get the genes that are present in the
  ## data set and calculate the percentage of transcripts that they account for
  ##--------------------------------------------------------------------------##
  result <- pbapply::pblapply(
    genes,
    function(x) {
      genes_here <- intersect(x, rownames(object@raw.data))
      Matrix::colSums(object@raw.data[genes_here,]) / Matrix::colSums(object@raw.data)
    }
  )
  ##--------------------------------------------------------------------------##
  ## return list with results
  ##--------------------------------------------------------------------------##
  return(result)
}