#' Get most expressed genes for every sample and cluster in Seurat object.
#'
#' This function calculates the most expressed genes for every sample and cluster of the Seurat object.
#' @param object Seurat object.
#' @param column_sample Column in object@meta.data that contains information about sample; defaults to "sample".
#' @param column_cluster Column in object@meta.data that contains information about cluster; defaults to "cluster".
#' @keywords seurat cerebro
#' @export
#' @examples
#' getMostExpressedGenes(object = seurat)

getMostExpressedGenes <- function(
  object,
  column_sample = "sample",
  column_cluster = "cluster"
) {
  require("dplyr")
  #
  seurat_object <- object
  #
  if ( !is.null(column_sample) && column_sample %in% names(seurat_object@meta.data) ) {
    # if sample column is already a factor, take the levels from there
    if ( is.factor(seurat_object@meta.data[[column_sample]]) ) {
      sample_names <- as.character(levels(seurat_object@meta.data[[column_sample]]))
    } else {
      sample_names <- unique(seurat_object@meta.data[[column_sample]])
    }
    if ( length(sample_names) > 1 ) {
      most_expressed_genes_by_sample <- data.frame(
        "sample" = character(),
        "gene" = character(),
        "pct" = double(),
        stringsAsFactors = FALSE
      )
      message("Get most expressed genes by sample...")
      pb = txtProgressBar(min = 0, max = length(sample_names), initial = 0, style = 3) 
      for ( i in 1:length(sample_names) ) {
        temp_table <- seurat_object@raw.data %>%
          as.data.frame(stringsAsFactors = FALSE) %>%
          select(which(seurat_object@meta.data[[column_sample]] == sample_names[i])) %>%
          mutate(
            sample = sample_names[i],
            gene = rownames(.),
            rowSums = rowSums(.),
            pct = rowSums / sum(.[1:(ncol(.))]) * 100
          ) %>%
          select(c("sample","gene","pct")) %>%
          arrange(-pct) %>%
          head(100)
        most_expressed_genes_by_sample <- rbind(
            most_expressed_genes_by_sample, temp_table
          )
        setTxtProgressBar(pb, i)
      }
      message("\n")
      most_expressed_genes_by_sample <- most_expressed_genes_by_sample %>%
        mutate(sample = factor(sample, levels = sample_names))
    }
  }
  #
  if ( !is.null(column_cluster) && column_cluster %in% names(seurat_object@meta.data) ) {
    if ( is.factor(seurat_object@meta.data[[column_cluster]]) ) {
      cluster_names <- as.character(levels(seurat_object@meta.data[[column_cluster]]))
    } else {
      cluster_names <- unique(seurat_object@meta.data[[column_cluster]])
    }
    if ( length(cluster_names) > 1 ) {
      most_expressed_genes_by_cluster <- data.frame(
          "cluster" = character(),
          "gene" = character(),
          "expr" = double(),
          stringsAsFactors = FALSE
        )
      message("Get most expressed genes by cluster...")
      pb = txtProgressBar(min = 0, max = length(cluster_names), initial = 0, style = 3) 
      for ( i in 1:length(cluster_names) ) {
        temp_table <- seurat_object@raw.data %>%
          as.data.frame(stringsAsFactors = FALSE) %>%
          select(which(seurat_object@meta.data$cluster == cluster_names[i])) %>%
          mutate(
            cluster = cluster_names[i],
            gene = rownames(.),
            rowSums = rowSums(.),
            pct = rowSums / sum(.[1:(ncol(.))]) * 100
          ) %>%
          select(c("cluster","gene","pct")) %>%
          arrange(-pct) %>%
          head(100)
        most_expressed_genes_by_cluster <- rbind(
            most_expressed_genes_by_cluster, temp_table
          )
        setTxtProgressBar(pb, i)
      }
      message("\n")
      most_expressed_genes_by_cluster <- most_expressed_genes_by_cluster %>%
        mutate(cluster = factor(cluster, levels = cluster_names))
    }
  }
  if ( is.null(seurat_object@misc$most_expressed_genes) ) {
    seurat_object@misc$most_expressed_genes <- list()
  }
  seurat_object@misc$most_expressed_genes$by_sample <- most_expressed_genes_by_sample
  seurat_object@misc$most_expressed_genes$by_cluster <- most_expressed_genes_by_cluster
  return(seurat_object)
}

