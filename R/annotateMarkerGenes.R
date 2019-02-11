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
  require("enrichR")
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
  temp_seurat <- object
  # check if marker genes are present
  if ( is.null(temp_seurat@misc$marker_genes) ) {
    stop("Please run 'getMarkerGenes()' first.", call. = FALSE)
  }
  # check if marker genes by sample are available
  if ( !is.null(temp_seurat@misc$marker_genes$by_sample) ) {
    # if sample column is already a factor, take the levels from there
    if ( is.factor(temp_seurat@meta.data[[column_sample]]) ) {
      sample_names <- as.character(levels(temp_seurat@meta.data[[column_sample]]))
    } else {
      sample_names <- unique(temp_seurat@meta.data[[column_sample]])
    }
    markers_by_sample <- temp_seurat@misc$marker_genes$by_sample
    # pb = txtProgressBar(min = 0, max = length(sample_names), initial = 0, style = 3) 
    # annotate marker genes for each sample
    for ( i in 1:length(sample_names) ) {
      message(paste0("Get annotation for sample '", sample_names[i], "'"))
      temp <- list()
      # try up to three times to run enrichR annotation (fails sometimes)
      attempt <- 1
      while( length(temp) == 0 && attempt <= 3 ) {
        attempt <- attempt + 1
        try(
          temp <- markers_by_sample %>%
            filter(sample == sample_names[i]) %>%
            select("gene") %>%
            t() %>%
            as.vector() %>%
            enrichR::enrichr(databases = enrichr_dbs)
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
      if ( is.null(temp_seurat@misc$marker_genes$by_sample_annotation) ) {
        temp_seurat@misc$marker_genes$by_sample_annotation <- list()
      }
      temp_seurat@misc$marker_genes$by_sample_annotation[[sample_names[i]]] <- temp
      # setTxtProgressBar(pb, i)
    }
    message("\n")
  }
  # check if marker genes by cluster are available
  if ( !is.null(temp_seurat@misc$marker_genes$by_cluster) ) {
    # if cluster column is already a factor, take the levels from there
    if ( is.factor(temp_seurat@meta.data[[column_cluster]]) ) {
      cluster_names <- as.character(levels(temp_seurat@meta.data[[column_cluster]]))
    } else {
      cluster_names <- unique(temp_seurat@meta.data[[column_cluster]])
    }
    markers_by_cluster <- temp_seurat@misc$marker_genes$by_cluster
    # pb = txtProgressBar(min = 0, max = length(cluster_names), initial = 0, style = 3) 
    # annotate marker genes for each cluster
    for ( i in 1:length(cluster_names) ) {
      message(paste0("Get annotation for cluster ", cluster_names[i]))
      temp <- list()
      # try up to three times to run enrichR annotation (fails sometimes)
      attempt <- 1
      while( length(temp) == 0 && attempt <= 3 ) {
        attempt <- attempt + 1
        try(
          temp <- markers_by_cluster %>%
            filter(cluster == cluster_names[i]) %>%
            select("gene") %>%
            t() %>%
            as.vector() %>%
            enrichR::enrichr(databases = enrichr_dbs)
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
      if ( is.null(temp_seurat@misc$marker_genes$by_cluster_annotation) ) {
        temp_seurat@misc$marker_genes$by_cluster_annotation <- list()
      }
      temp_seurat@misc$marker_genes$by_cluster_annotation[[cluster_names[i]]] <- temp
      # setTxtProgressBar(pb, i)
    }
    message("\n")
  }
  return(temp_seurat)
}







