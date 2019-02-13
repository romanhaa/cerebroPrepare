#' Calculate percentage of transcripts of gene list.
#'
#' Get percentage of transcripts of gene list compared to all transcripts per cell.
#' @param object Seurat object.
#' @param genes List(s) of genes.
#' @keywords seurat cerebro
#' @export
#' @examples
#' calculatePercentGenes(object = seurat, genes = gene_list)
calculatePercentGenes <- function(
  object,
  genes
) {
  require("Matrix")
  require("Seurat")
  require("pbapply")
  result <- pblapply(
    genes,
    function(x) {
      genes_here <- intersect(x, rownames(object@raw.data))
      Matrix::colSums(object@raw.data[genes_here,]) / Matrix::colSums(object@raw.data)
    }
  )
  return(result)
}
