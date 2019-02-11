#' Annotate marker genes with EnrichR.
#'
#' This function uses the enrichR API to look for enriched pathways in marker gene sets of samples and clusters.
#' @param object Seurat object.
#' @param column_sample Column in object@meta.data that contains information about sample; defaults to "sample".
#' @param column_cluster Column in object@meta.data that contains information about cluster; defaults to "cluster".
#' @param adj_p_cutoff Cut-off for adjusted p-value of enriched pathways; defaults to 0.01,
#' @param max_terms Save only first n entries of each database; defaults to 100.
#' @keywords seurat cerebro
#' @export
#' @examples
#' annotateMarkerGenes(object = seurat)

annotateMarkerGenes <- function(
  object,
  column_sample = "sample",
  column_cluster = "cluster",
  adj_p_cutoff = 0.01,
  max_terms = 100
) {
  # try to load enrichR package and complain if it's not available
  if (!requireNamespace("enrichR", quietly = TRUE)) {
    stop(
      "Package 'enrichR' needed for this function to work. Please install it.",
      call. = FALSE
    )
  }
  #
  enrichr_dbs <- c(
      "GO_Biological_Process_2018",
      "GO_Cellular_Component_2018",
      "GO_Molecular_Function_2018",
      "KEGG_2016",
      "WikiPathways_2016",
      "Reactome_2016",
      "Panther_2016",
      "Human_Gene_Atlas",
      "Mouse_Gene_Atlas"
    )
  #
  seurat_object <- object
  # check if marker genes are present
  if ( is.null(seurat_object@misc$marker_genes) ) {
    stop("Please run 'getMarkerGenes()' first.", call. = FALSE)
  }
  # check if marker genes by sample are available
  if ( !is.null(seurat_object@misc$marker_genes$by_sample) ) {
    # if sample column is already a factor, take the levels from there
    if ( is.factor(seurat_object@meta.data[column_sample]) ) {
      sample_names <- levels(seurat_object@meta.data[column_sample])
    } else {
      sample_names <- unique(seurat_object@meta.data[column_sample])
    }
    markers_by_sample <- seurat_object@misc$marker_genes$by_sample
    # annotate marker genes for each sample
    for ( i in sample_names ) {
      temp <- list()
      # try up to three times to run enrichR annotation (fails sometimes)
      attempt <- 1
      while( length(temp) == 0 && attempt <= 3 ) {
        attempt <- attempt + 1
        try(
          temp <- markers_by_sample %>%
            filter(sample == i) %>%
            select("gene") %>%
            t() %>%
            as.vector() %>%
            enrichr(databases = enrichr_dbs)
        )
      }
      # filter results
      for ( j in names(temp) ) {
        length <- temp[[j]] %>% filter(Adjusted.P.value <= adj_p_cutoff) %>% nrow()
        # if there are more than max_terms entries with an adjusted p-value of 1 or less...
        if ( length > max_terms ) {
          temp[[j]] <- temp[[j]] %>%
            top_n(-max_terms, Adjusted.P.value)
        # if there is at least 1 entry with an adjusted p-value of 1 or less...
        } else if ( length > 0 ) {
          temp[[j]] <- temp[[j]] %>%
            filter(Adjusted.P.value <= adj_p_cutoff)
        # remove the curent database
        } else {
          temp[j] <- NULL
        }
      }
      if ( is.null(seurat_object@misc$marker_genes$by_sample_annotation) ) {
        seurat_object@misc$marker_genes$by_sample_annotation <- list()
      }
      seurat_object@misc$marker_genes$by_sample_annotation[[i]] <- temp
    }
  }
  # check if marker genes by cluster are available
  if ( !is.null(seurat_object@misc$marker_genes$by_cluster) ) {
    # if cluster column is already a factor, take the levels from there
    if ( is.factor(seurat_object@meta.data[column_cluster]) ) {
      cluster_names <- levels(seurat_object@meta.data[column_cluster])
    } else {
      cluster_names <- unique(seurat_object@meta.data[column_cluster])
    }
    markers_by_cluster <- seurat_object@misc$marker_genes$by_cluster
    # annotate marker genes for each cluster
    for ( i in cluster_names ) {
      temp <- list()
      # try up to three times to run enrichR annotation (fails sometimes)
      attempt <- 1
      while( length(temp) == 0 && attempt <= 3 ) {
        attempt <- attempt + 1
        try(
          temp <- markers_by_cluster %>%
            filter(cluster == i) %>%
            select("gene") %>%
            t() %>%
            as.vector() %>%
            enrichr(databases = enrichr_dbs)
        )
      }
      # filter results
      for ( j in names(temp) ) {
        length <- temp[[j]] %>% filter(Adjusted.P.value <= adj_p_cutoff) %>% nrow()
        # if there are more than max_terms entries with an adjusted p-value of 1 or less...
        if ( length > max_terms ) {
          temp[[j]] <- temp[[j]] %>%
            top_n(-max_terms, Adjusted.P.value)
        # if there is at least 1 entry with an adjusted p-value of 1 or less...
        } else if ( length > 0 ) {
          temp[[j]] <- temp[[j]] %>%
            filter(Adjusted.P.value <= adj_p_cutoff)
        # remove the curent database
        } else {
          temp[j] <- NULL
        }
      }
      if ( is.null(seurat_object@misc$marker_genes$by_cluster_annotation) ) {
        seurat_object@misc$marker_genes$by_cluster_annotation <- list()
      }
      seurat_object@misc$marker_genes$by_cluster_annotation[[i]] <- temp
    }
  }
  return(seurat_object)
}







