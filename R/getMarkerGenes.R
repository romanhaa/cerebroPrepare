#' Get marker genes for every sample and cluster in Seurat object.
#' @title Get marker genes for every sample and cluster in Seurat object.
#' @description This function gets marker genes for every sample and cluster of the Seurat object.
#' @param object Seurat object.
#' @param organism Organism information for pulling info about presence of marker genes of cell surface; can be omitted if already saved in Seurat object; defaults to NULL.
#' @param column_sample Column in object@meta.data that contains information about sample; defaults to "sample".
#' @param column_cluster Column in object@meta.data that contains information about cluster; defaults to "cluster".
#' @param only.pos Identify only over-expressed genes; defaults to TRUE.
#' @param min.pct Only keep genes that are expressed in at least n \% of current group of cells, defaults to 0.70 (70\%).
#' @param thresh.use Only keep genes that have an FDR of less than or equal to n, defaults to 0.25.
#' @param test.use Statistical test used, defaults to "t" (t-test).
#' @param print.bar Print progress bar; defaults to TRUE.
#' @param ... Further parameters can be passed to control Seurat::FindAllMakers().
#' @keywords seurat cerebro
#' @export
#' @import dplyr
#' @import Seurat
#' @examples
#' getMarkerGenes(object = seurat)

getMarkerGenes <- function(
  object,
  organism = NULL,
  column_sample = "sample",
  column_cluster = "cluster",
  only.pos = TRUE,
  min.pct = 0.70,
  thresh.use = 0.25,
  test.use = "t",
  print.bar = TRUE,
  ...
) {
  ##--------------------------------------------------------------------------##
  ## Get list of genes in cell surface through gene ontology term GO:0009986.
  ##--------------------------------------------------------------------------##
  if ( organism == "hg" || organism == "human" ) {
    genes_on_cell_surface <- biomaRt::getBM(
      attributes = "hgnc_symbol",
      filters = "go",
      values = "GO:0009986",
      mart = biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    )[,1]
  } else if ( organism == "mm" || organism == "mouse" ) {
    genes_on_cell_surface <- biomaRt::getBM(
      attributes = "external_gene_name",
      filters = "go",
      values = "GO:0009986",
      mart = biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")
    )[,1]
  } else {
    message("No information about genes on cell surface because organism is either not specified or not human/mouse.")
  }
  ##--------------------------------------------------------------------------##
  ## make copy of Seurat object
  ##--------------------------------------------------------------------------##
  temp_seurat <- object
  ##--------------------------------------------------------------------------##
  ## samples
  ## - check if column_sample is provided and exists in meta data
  ## - get sample names
  ## - check if more than 1 sample exists
  ## - run Seurat::FindAllMarkers()
  ## - check if any marker genes were found
  ## - sort and rename columns
  ## - add column with cell surface genes if present
  ##--------------------------------------------------------------------------##
  #
  if ( !is.null(column_sample) && (column_sample %in% names(temp_seurat@meta.data)) ) {
    #
    if ( is.factor(temp_seurat@meta.data[[column_sample]]) ) {
      sample_names <- levels(temp_seurat@meta.data[[column_sample]])
    } else {
      sample_names <- unique(temp_seurat@meta.data[[column_sample]])
    }
    #
    if ( length(sample_names) > 1 ) {
      temp_seurat <- SetAllIdent(temp_seurat, id = column_sample)
      temp_seurat@ident <- factor(temp_seurat@ident, levels = sample_names)
      message("Get marker genes by sample...")
      markers_by_sample <- Seurat::FindAllMarkers(
          temp_seurat,
          only.pos = only.pos,
          min.pct = min.pct,
          thresh.use = thresh.use,
          test.use = test.use,
          print.bar = print.bar,
          ...
        )
      #
      if ( nrow(markers_by_sample) > 0 ) {
        markers_by_sample <- markers_by_sample %>%
          dplyr::select(c("cluster", "gene", "p_val", "avg_logFC", "pct.1", "pct.2", "p_val_adj")
          ) %>%
          dplyr::rename(
            sample = cluster
          )
        #
        if ( exists("genes_surface") ) {
          markers_by_sample <- markers_by_sample %>%
            mutate(on_cell_surface = gene %in% genes_surface)
        }
      } else {
        message("No marker genes found for any of the samples.")
        markers_by_sample <- "no_markers_found"
      }
    } else {
      message("Sample column provided but only 1 sample found.")
    }
  } else {
    warning("Provided column name with sample information cannot be found.")
  }
  ##--------------------------------------------------------------------------##
  ## clusters
  ## - check if column_cluster is provided and exists in meta data
  ## - get cluster names
  ## - check if more than 1 cluster exists
  ## - run Seurat::FindAllMarkers()
  ## - check if any marker genes were found
  ## - sort and rename columns
  ## - add column with cell surface genes if present
  ##--------------------------------------------------------------------------## 
  #
  if ( !is.null(column_cluster) & column_cluster %in% names(temp_seurat@meta.data) ) {
    if ( is.factor(temp_seurat@meta.data[[column_cluster]]) ) {
      cluster_names <- levels(temp_seurat@meta.data[[column_cluster]])
    } else {
      cluster_names <- sort(unique(temp_seurat@meta.data[[column_cluster]]))
    }
    #
    if ( length(cluster_names) > 1 ) {
      temp_seurat <- SetAllIdent(temp_seurat, id = column_cluster)
      temp_seurat@ident <- factor(temp_seurat@ident, levels = cluster_names)
      message("Get marker genes by cluster...")
      markers_by_cluster <- Seurat::FindAllMarkers(
          temp_seurat,
          only.pos = only.pos,
          min.pct = min.pct,
          thresh.use = thresh.use,
          test.use = test.use,
          print.bar = print.bar,
          ...
        )
      #
      if ( nrow(markers_by_cluster) > 0 ) {
        markers_by_cluster <- markers_by_cluster %>%
          dplyr::select(c("cluster", "gene", "p_val", "avg_logFC", "pct.1", "pct.2", "p_val_adj")
          )
        #
        if ( exists("genes_surface") ) {
          markers_by_cluster <- markers_by_cluster %>%
            mutate(on_cell_surface = gene %in% genes_surface)
        }
      } else {
        message("No marker genes found for any of the clusters.")
        markers_by_cluster <- "no_markers_found"
      }
    } else {
      message("Cluster column provided but only 1 cluster found.")
    }
  } else {
    warning("Provided column name with cluster information cannot be found.")
  }
  ##--------------------------------------------------------------------------##
  ## store results in Seurat object
  ##--------------------------------------------------------------------------##
  if ( is.null(object@misc$marker_genes) ) {
    temp_seurat@misc$marker_genes <- list()
  }
  temp_seurat@misc$marker_genes$by_sample <- markers_by_sample
  temp_seurat@misc$marker_genes$by_cluster <- markers_by_cluster
  ##--------------------------------------------------------------------------##
  ## return Seurat object
  ##--------------------------------------------------------------------------##
  return(temp_seurat)
}










