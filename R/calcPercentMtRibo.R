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

  function(
  # - get genes from list and check which are present in data set
  # - get colSums of present genes
    Matrix::colSums(object@raw.data[genes_mt_here,]) / 
    Matrix::colSums(object@raw.data)
  )

  result <- lapply(
    function()
  )

  percent_mt <- Matrix::colSums(seurat@raw.data[genes_mt_here,]) / 
    Matrix::colSums(seurat@raw.data)

  percent_ribo <- Matrix::colSums(seurat@raw.data[genes_ribo_here,]) /
    Matrix::colSums(seurat@raw.data)

  return(
    list()
  )
}

#' Calculate percentage of transcripts of gene list.
#'
#' Get percentage of transcripts of gene list compared to all transcripts per cell.
#' @param object Seurat object.
#' @param genes List(s) of genes.
#' @keywords seurat cerebro
#' @export
#' @examples
#' calculatePercentGenes(object = seurat, genes = gene_list)
addPercentMtRibo <- function(
  object
) {
  require("dplyr")
  require("Seurat")
  require("readr")

  genes_mt <- c()
  genes_ribo <- c()


  if ( organism %in% c("hg","mm") ) {

    genes_mt <- read_tsv(
        paste0(cerebro_prefix, "/resources/", organism, "/genes_mt_names.txt"),
        col_names = FALSE
      ) %>%
      select(1) %>%
      t() %>%
      as.vector()

    genes_ribo <- read_tsv(
      paste0(cerebro_prefix, "/resources/", organism, "/genes_ribo_names.txt"),
      col_names = FALSE
    ) %>%
    select(1) %>%
    t() %>%
    as.vector()


    genes_mt_here <- intersect(genes_mt, rownames(object@raw.data))
    genes_ribo_here <- intersect(genes_ribo, rownames(object@raw.data))

    object@misc$gene_lists$mitochondrial_genes <- genes_mt_here
    object@misc$gene_lists$ribosomal_genes <- genes_ribo_here

    values <- calculatePercentGenes(
      object,
      list(
        "genes_mt" = genes_mt,
        "genes_ribo" = genes_ribo
      )
    )
    object <- AddMetaData(
      object,
      data.frame(
        row.names(object@cell.names),
        "percent_mt" = values[[genes_mt]],
        "percent_ribo" = values[[genes_ribo]]
      )
    )
    return(object)
  } else {
    stop("Organism-specific list of mitochondrial and ribosomal genes not available.")
  }
}









