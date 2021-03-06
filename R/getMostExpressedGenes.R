#' Get most expressed genes for every sample and cluster in Seurat object.
#' @title Get most expressed genes for every sample and cluster in Seurat
#' object.
#' @description This function calculates the most expressed genes for every
#' sample and cluster of the Seurat object.
#' @param object Seurat object.
#' @param column_sample Column in object@meta.data that contains information
#' about sample; defaults to 'sample'.
#' @param column_cluster Column in object@meta.data that contains information
#' about cluster; defaults to 'cluster'.
#' @keywords seurat cerebro
#' @export
#' @import dplyr
#' @examples
#' seurat <- getMostExpressedGenes(
#'   object = seurat,
#'   column_sample = 'sample',
#'   column_cluster = 'cluster'
#' )

getMostExpressedGenes <- function(
  object,
  column_sample = 'sample',
  column_cluster = 'cluster'
) {
  ##--------------------------------------------------------------------------##
  ## create backup of Seurat object (probably not necessary)
  ##--------------------------------------------------------------------------##
  temp_seurat <- object
  ##--------------------------------------------------------------------------##
  ## create slot for results in Seurat object if not already existing
  ##--------------------------------------------------------------------------##
  if ( is.null(temp_seurat@misc$most_expressed_genes) ) {
    temp_seurat@misc$most_expressed_genes <- list()
  }
  ##--------------------------------------------------------------------------##
  ## samples
  ## - check if column_sample was provided and exists in meta data
  ## - get sample names
  ## - check if more than 1 sample is available
  ## - retrieve most expressed genes in parallel with future_lapply
  ## - combine results in a single data frame
  ## - sort by sample (probably not necessary)
  ## - store results in Seurat object
  ##--------------------------------------------------------------------------##
  if ( !is.null(column_sample) && column_sample %in% names(temp_seurat@meta.data) ) {
    #
    if ( is.factor(temp_seurat@meta.data[[column_sample]]) ) {
      sample_names <- as.character(levels(temp_seurat@meta.data[[column_sample]]))
    } else {
      sample_names <- unique(temp_seurat@meta.data[[column_sample]])
      if ( any(is.na(sample_names)) ) {
        number_of_NA_in_sample_column <- which(is.na(sample_names)) %>% length()
        sample_names <- na.omit(sample_names)
        message(
          paste0(
            '[', format(Sys.time(), '%H:%M:%S'), '] Found ', number_of_NA_in_sample_column,
            ' cell(s) without cluster information (NA). These cells will be ignored during this analysis.'
          )
        )
      }
    }
    if ( length(sample_names) > 1 ) {
      most_expressed_genes_by_sample <- data.frame(
        'sample' = character(),
        'gene' = character(),
        'pct' = double(),
        stringsAsFactors = FALSE
      )
      message(paste0('[', format(Sys.time(), '%H:%M:%S'), '] Get most expressed genes by sample...'))
      if ( temp_seurat@version < 3 ) {
        results <- pbapply::pblapply(sample_names, function(x) {
          temp_table <- temp_seurat@raw.data %>%
            as.data.frame(stringsAsFactors = FALSE) %>%
            dplyr::select(which(temp_seurat@meta.data[[column_sample]] == x)) %>%
            mutate(
              sample = x,
              gene = rownames(.),
              rowSums = rowSums(.),
              pct = rowSums / sum(.[1:(ncol(.))]) * 100
            ) %>%
            dplyr::select(c('sample','gene','pct')) %>%
            arrange(-pct) %>%
            head(100)
        })
      } else {
        results <- pbapply::pblapply(sample_names, function(x) {
          temp_table <- temp_seurat@assays$RNA@counts %>%
            as.data.frame(stringsAsFactors = FALSE) %>%
            dplyr::select(which(temp_seurat@meta.data[[column_sample]] == x)) %>%
            mutate(
              sample = x,
              gene = rownames(.),
              rowSums = rowSums(.),
              pct = rowSums / sum(.[1:(ncol(.))]) * 100
            ) %>%
            dplyr::select(c('sample','gene','pct')) %>%
            arrange(-pct) %>%
            head(100)
        })
      }
      most_expressed_genes_by_sample <- do.call(rbind, results) %>%
        mutate(sample = factor(sample, levels = sample_names))
      temp_seurat@misc$most_expressed_genes$by_sample <- most_expressed_genes_by_sample
    }
  }
  #
  ##--------------------------------------------------------------------------##
  ## cluster
  ## - check if column_cluster was provided and exists in meta data
  ## - get cluster names
  ## - check if more than 1 cluster is available
  ## - retrieve most expressed genes in parallel with future_lapply
  ## - combine results in a single data frame
  ## - sort by cluster (probably not necessary)
  ## - store results in Seurat object
  ##--------------------------------------------------------------------------##
  if ( !is.null(column_cluster) && column_cluster %in% names(temp_seurat@meta.data) ) {
    if ( is.factor(temp_seurat@meta.data[[column_cluster]]) ) {
      cluster_names <- as.character(levels(temp_seurat@meta.data[[column_cluster]]))
    } else {
      cluster_names <- sort(unique(temp_seurat@meta.data[[column_cluster]]))
      if ( any(is.na(cluster_names)) ) {
        number_of_NA_in_cluster_column <- which(is.na(cluster_names)) %>% length()
        cluster_names <- na.omit(cluster_names)
        message(
          paste0(
            '[', format(Sys.time(), '%H:%M:%S'), '] Found ', number_of_NA_in_cluster_column,
            ' cell(s) without cluster information (NA). These cells will be ignored during this analysis.'
          )
        )
      }
    }
    if ( length(cluster_names) > 1 ) {
      most_expressed_genes_by_cluster <- data.frame(
          'cluster' = character(),
          'gene' = character(),
          'expr' = double(),
          stringsAsFactors = FALSE
        )
      message(paste0('[', format(Sys.time(), '%H:%M:%S'), '] Get most expressed genes by cluster...'))
      if ( temp_seurat@version < 3 ) {
        results <- pbapply::pblapply(cluster_names, function(x) {
          temp_table <- temp_seurat@raw.data %>%
            as.data.frame(stringsAsFactors = FALSE) %>%
            dplyr::select(which(temp_seurat@meta.data[[column_cluster]] == x)) %>%
            mutate(
              cluster = x,
              gene = rownames(.),
              rowSums = rowSums(.),
              pct = rowSums / sum(.[1:(ncol(.))]) * 100
            ) %>%
            dplyr::select(c('cluster','gene','pct')) %>%
            arrange(-pct) %>%
            head(100)
        })
      } else {
        results <- pbapply::pblapply(cluster_names, function(x) {
          temp_table <- temp_seurat@assays$RNA@counts %>%
            as.data.frame(stringsAsFactors = FALSE) %>%
            dplyr::select(which(temp_seurat@meta.data[[column_cluster]] == x)) %>%
            mutate(
              cluster = x,
              gene = rownames(.),
              rowSums = rowSums(.),
              pct = rowSums / sum(.[1:(ncol(.))]) * 100
            ) %>%
            dplyr::select(c('cluster','gene','pct')) %>%
            arrange(-pct) %>%
            head(100)
        })
      }
      most_expressed_genes_by_cluster <- do.call(rbind, results) %>%
        mutate(cluster = factor(cluster, levels = cluster_names))
      temp_seurat@misc$most_expressed_genes$by_cluster <- most_expressed_genes_by_cluster
    }
  }
  ##--------------------------------------------------------------------------##
  ## return Seurat object
  ##--------------------------------------------------------------------------##
  return(temp_seurat)
}

