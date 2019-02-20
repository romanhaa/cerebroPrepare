#' Add percentage of mitochondrial and ribosomal transcripts.
#' @title Add percentage of mitochondrial and ribosomal transcripts.
#' @description Get percentage of transcripts of gene list compared to all transcripts per cell.
#' @param object Seurat object.
#' @param organism Organism, can be either human ("hg") or mouse ("mm"). Genes need to annotated as gene symbol, e.g. MKI67 (human) / Mki67 (mouse).
#' @keywords seurat cerebro
#' @export
#' @import dplyr
#' @import Seurat
#' @examples
#' calculatePercentGenes(object = seurat, genes = gene_list)
addPercentMtRibo <- function(
  object,
  organism
) {
  ##--------------------------------------------------------------------------##
  ## check if organism is supported
  ##--------------------------------------------------------------------------##
  if ( !(organism %in% c("hg","mm")) ) {
    stop("Organism-specific list of mitochondrial and ribosomal genes not available.")
  }
  ##--------------------------------------------------------------------------##
  ## load mitochondrial and ribosomal gene lists from extdata
  ##--------------------------------------------------------------------------##
  genes_mt <- read_tsv(
      # "/hpcnfs/scratch/PGP/rhillje/cerebroPrepare/inst/extdata/genes_mt_mm.txt",
      system.file(
        "extdata",
        paste0("genes_mt_", organism, "_names.txt"), package = "cerebroPrepare"
      ),
      col_types = cols(),
      col_names = FALSE
    ) %>%
    dplyr::select(1) %>%
    t() %>%
    as.vector()
  genes_ribo <- read_tsv(
      # "/hpcnfs/scratch/PGP/rhillje/cerebroPrepare/inst/extdata/genes_ribo_mm.txt",
      system.file(
        "extdata",
        paste0("genes_ribo_", organism, "_names.txt"), package = "cerebroPrepare"
      ),
      col_types = cols(),
      col_names = FALSE
    ) %>%
    dplyr::select(1) %>%
    t() %>%
    as.vector()
  ##--------------------------------------------------------------------------##
  ## keep only genes that are present in data set
  ##--------------------------------------------------------------------------##
  genes_mt_here <- intersect(genes_mt, rownames(object@raw.data))
  genes_ribo_here <- intersect(genes_ribo, rownames(object@raw.data))
  ##--------------------------------------------------------------------------##
  ## save gene lists in Seurat object and create place if not existing yet
  ##--------------------------------------------------------------------------##
  if ( is.null(object@misc$gene_lists) ) {
    object@misc$gene_lists <- list()
  }
  object@misc$gene_lists$mitochondrial_genes <- genes_mt_here
  object@misc$gene_lists$ribosomal_genes <- genes_ribo_here
  ##--------------------------------------------------------------------------##
  ## calculate percentage of transcripts for mitochondrial and ribosomal genes
  ##--------------------------------------------------------------------------##
  message("Calculate percentage of mitochondrial and ribosomal transcripts...")
  values <- calculatePercentGenes(
    object,
    list(
      "genes_mt" = genes_mt_here,
      "genes_ribo" = genes_ribo_here
    )
  )
  ##--------------------------------------------------------------------------##
  ## add results to Seurat object
  ##--------------------------------------------------------------------------##
  object <- Seurat::AddMetaData(
    object,
    data.frame(
      colnames(object@raw.data),
      "percent_mt" = values[["genes_mt"]],
      "percent_ribo" = values[["genes_ribo"]]
    )
  )
  ##--------------------------------------------------------------------------##
  ##
  ##--------------------------------------------------------------------------##
  return(object)
}









