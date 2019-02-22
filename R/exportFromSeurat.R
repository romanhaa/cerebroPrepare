#' Export Seurat object to Cerebro.
#' @title Export Seurat object to Cerebro.
#' @description This function allows to export a Seurat object to visualize in Cerebro.
#' @param object Seurat object.
#' @param file Where to save the output.
#' @param experiment_name Experiment name.
#' @param organism Organism, e.g. hg (human), mm (mouse), etc.
#' @param column_sample Column in object@meta.data that contains information about sample; defaults to "sample".
#' @param column_cluster Column in object@meta.data that contains information about cluster; defaults to "cluster".
#' @param column_nUMI Column in object@meta.data that contains information about number of transcripts per cell; defaults to "nUMI".
#' @param column_nGene Column in object@meta.data that contains information about number of expressed genes per cell; defaults to "nGene".
#' @param column_cell_cycle_seurat Optional column in object@meta.data that contains information about cell cycle phase based on Regev method (default of Seurat); defaults to NULL.
#' @param column_cell_cycle_cyclone Optional column in object@meta.data that contains information about cell cycle phase based on Cyclone method; defaults to NULL.
#' @param add_all_meta_data If set to TRUE, all further meta data columns will be extracted as well.
#' @return Returns object to be saved and loaded in Cerebro.
#' @keywords seurat cerebro
#' @export
#' @import dplyr
#' @import Seurat
#' @import tidyr
#' @examples
#' exportFromSeurat(object = seurat, file = "PDX_patient_A.cerebro" experiment_name = "PDX_patient_A", organism = "hg")
exportFromSeurat <- function(
  object,
  file,
  experiment_name,
  organism,
  column_sample = "sample",
  column_cluster = "cluster",
  column_nUMI = "nUMI",
  column_nGene = "nGene",
  column_cell_cycle_seurat = NULL,
  column_cell_cycle_cyclone = NULL,
  add_all_meta_data = TRUE
) {
  ##--------------------------------------------------------------------------##
  ## check provided parameters
  ##--------------------------------------------------------------------------##
  if ( (column_sample %in% names(object@meta.data) == FALSE ) ) {
    stop(
      "Column specified in 'column_sample' not found in meta data.",
      call. = FALSE
    )
  }
  if ( (column_cluster %in% names(object@meta.data) == FALSE ) ) {
    stop(
      "Column specified in 'column_cluster' not found in meta data.",
      call. = FALSE
    )
  }
  if ( (column_nUMI %in% names(object@meta.data) == FALSE ) ) {
    stop(
      "Column with number of transcripts per cell ('nUMI') not found in meta data.",
      call. = FALSE
    )
  }
  if ( (column_nGene %in% names(object@meta.data) == FALSE ) ) {
    stop(
      "Column with number of expressed genes per cell ('nGene') not found in meta data.",
      call. = FALSE
    )
  }
  ##--------------------------------------------------------------------------##
  ## colors
  ##--------------------------------------------------------------------------##
  # Dutch palette from flatuicolors.com
  colors_dutch <- c(
      "#FFC312","#C4E538","#12CBC4","#FDA7DF","#ED4C67",
      "#F79F1F","#A3CB38","#1289A7","#D980FA","#B53471",
      "#EE5A24","#009432","#0652DD","#9980FA","#833471",
      "#EA2027","#006266","#1B1464","#5758BB","#6F1E51"
    )
  # Spanish palette from flatuicolors.com
  colors_spanish <- c(
      "#40407a","#706fd3","#f7f1e3","#34ace0","#33d9b2",
      "#2c2c54","#474787","#aaa69d","#227093","#218c74",
      "#ff5252","#ff793f","#d1ccc0","#ffb142","#ffda79",
      "#b33939","#cd6133","#84817a","#cc8e35","#ccae62"
    )
  colors <- c(colors_dutch, colors_spanish)
  cell_cycle_colorset <- setNames(
    c("#45aaf2","#f1c40f","#e74c3c", "#7f8c8d"),
    c("G1","S","G2M","-")
  )
  ##--------------------------------------------------------------------------##
  ## initialize export object
  ##--------------------------------------------------------------------------##
  export <- list(
    experiment = list(
      experiment_name = experiment_name,
      organism = organism
    )
  )
  ##--------------------------------------------------------------------------##
  ## collect some more data if present
  ##--------------------------------------------------------------------------##
  #
  export$experiment$date_of_analysis <- object@misc$experiment$date_of_analysis
  #
  export$parameters <- ifelse(
    is.null(object@misc$parameters), list(), object@misc$parameters)
  #
  export$gene_lists <- ifelse(
    is.null(object@misc$gene_lists), list(), object@misc$gene_lists)
  #
  export$technical_info <- ifelse(
    is.null(object@misc$technical_info), list(), object@misc$technical_info)
  ##--------------------------------------------------------------------------##
  ## samples
  ##--------------------------------------------------------------------------##
  if ( is.factor(object@meta.data[[column_sample]]) ) {
    sample_names <- levels(object@meta.data[[column_sample]])
  } else {
    sample_names <- unique(object@meta.data[[column_sample]])
  }
  export$samples <- list(
    colors = setNames(colors[ 1:length(sample_names) ], sample_names),
    overview = data.frame("sample" = sample_names)
  )
  ##--------------------------------------------------------------------------##
  ## clusters
  ##--------------------------------------------------------------------------##
  if ( is.factor(object@meta.data[[column_cluster]]) ) {
    cluster_names <- levels(object@meta.data[[column_cluster]])
  } else {
    cluster_names <- sort(unique(object@meta.data[[column_cluster]]))
  }
  export$clusters <- list(
    colors = setNames(colors[ 1:length(cluster_names) ], cluster_names),
    overview = data.frame("cluster" = cluster_names)
  )
  ##--------------------------------------------------------------------------##
  ## meta data
  ##--------------------------------------------------------------------------##
  meta_data_columns <- names(object@meta.data)
  export$cells <- data.frame(
    "sample" = factor(object@meta.data[[column_sample]], levels = c(sample_names)),
    "cluster" = factor(object@meta.data[[column_cluster]], levels = c(cluster_names)),
    "nUMI" = object@meta.data[column_nUMI],
    "nGene" = object@meta.data[column_nGene]
  )
  ##--------------------------------------------------------------------------##
  ## samples by cluster
  ##--------------------------------------------------------------------------##
  export$samples$by_cluster <- export$cells %>%
    group_by(sample, cluster) %>%
    summarize(count=n()) %>%
    spread(cluster, count, fill = 0) %>%
    ungroup() %>%
    mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
    dplyr::select(c("sample", "total_cell_count", everything())) %>%
    arrange(factor(sample, levels = sample_names))
  ##--------------------------------------------------------------------------##
  ## clusters by sample
  ##--------------------------------------------------------------------------##
  export$clusters$by_samples <- export$cells %>%
    group_by(cluster, sample) %>%
    summarize(count=n()) %>%
    spread(sample, count, fill = 0) %>%
    ungroup() %>%
    mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
    dplyr::select(c("cluster", "total_cell_count", everything())) %>%
    arrange(factor(cluster, levels = cluster_names))
  meta_data_columns <- meta_data_columns[-which(meta_data_columns == column_sample)]
  meta_data_columns <- meta_data_columns[-which(meta_data_columns == column_cluster)]
  meta_data_columns <- meta_data_columns[-which(meta_data_columns == column_nUMI)]
  meta_data_columns <- meta_data_columns[-which(meta_data_columns == column_nGene)]
  ##--------------------------------------------------------------------------##
  ## cell cycle Seurat (if present)
  ##--------------------------------------------------------------------------##
  if ( !is.null(column_cell_cycle_seurat) && column_cell_cycle_seurat %in% names(object@meta.data) ) {
    export$cells$cell_cycle_seurat <- object@meta.data[[column_cell_cycle_seurat]]
    # by sample
    export$samples$by_cell_cycle_seurat <- export$cells %>%
      group_by(sample, cell_cycle_seurat) %>%
      summarize(count=n()) %>%
      spread(cell_cycle_seurat, count, fill = 0) %>%
      ungroup() %>%
      mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
      dplyr::select(c("sample", "total_cell_count", everything())) %>%
      arrange(factor(sample, levels = sample_names))
    # by cluster
    export$clusters$by_cell_cycle_seurat <- export$cells %>%
      group_by(cluster, cell_cycle_seurat) %>%
      summarize(count=n()) %>%
      spread(cell_cycle_seurat, count, fill = 0) %>%
      ungroup() %>%
      mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
      dplyr::select(c("cluster", "total_cell_count", everything())) %>%
      arrange(factor(cluster, levels = cluster_names))
    meta_data_columns <- meta_data_columns[-which(meta_data_columns == column_cell_cycle_seurat)]
  }
  ##--------------------------------------------------------------------------##
  ## cell cycle Cyclone (if present)
  ##--------------------------------------------------------------------------##
  if ( !is.null(column_cell_cycle_cyclone) && column_cell_cycle_cyclone %in% names(object@meta.data) ) {
    export$cells$cell_cycle_cyclone <- object@meta.data[[column_cell_cycle_cyclone]]
    # by sample
    export$samples$by_cell_cycle_cyclone <- export$cells %>%
      group_by(sample, cell_cycle_cyclone) %>%
      summarize(count=n()) %>%
      spread(cell_cycle_cyclone, count, fill = 0) %>%
      ungroup() %>%
      mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
      dplyr::select(c("sample", "total_cell_count", everything())) %>%
      arrange(factor(sample, levels = sample_names))
    # by cluster
    export$clusters$by_cell_cycle_cyclone <- export$cells %>%
      group_by(cluster, cell_cycle_cyclone) %>%
      summarize(count=n()) %>%
      spread(cell_cycle_cyclone, count, fill = 0) %>%
      ungroup() %>%
      mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
      dplyr::select(c("cluster", "total_cell_count", everything())) %>%
      arrange(factor(cluster, levels = cluster_names))
    meta_data_columns <- meta_data_columns[-which(meta_data_columns == column_cell_cycle_cyclone)]
  }
  ##--------------------------------------------------------------------------##
  ## cell barcode
  ##--------------------------------------------------------------------------##
  if ( !is.null(rownames(as.data.frame(object@meta.data))) ) {
    export$cells$cell_barcode <- rownames(as.data.frame(object@meta.data))
  }
  ##--------------------------------------------------------------------------##
  ## add all other meta data if specified
  ##--------------------------------------------------------------------------##
  if ( add_all_meta_data ) {
    export$cells <- cbind(export$cells, object@meta.data[meta_data_columns])
  }
  ##--------------------------------------------------------------------------##
  ## most expressed genes
  ##--------------------------------------------------------------------------##
  if ( !is.null(object@misc$most_expressed_genes) ) {
    export$most_expressed_genes <- object@misc$most_expressed_genes
  }
  ##--------------------------------------------------------------------------##
  ## marker genes
  ##--------------------------------------------------------------------------##
  if ( !is.null(object@misc$marker_genes) ) {
    export$marker_genes <- object@misc$marker_genes
  }
  ##--------------------------------------------------------------------------##
  ## dimensional reductions
  ##--------------------------------------------------------------------------##
  export$projections <- list()
  projections_available <- names(object@dr)
  projections_available_non_pca <- projections_available[grep(projections_available, pattern = "pca", invert = TRUE)]
  if ( length(projections_available) == 0 ) {
    stop("Warning: No dimensional reductions available.", call. = FALSE)
  } else if ( length(projections_available) == 1 && projections_available_non_pca == 1 ) {
    export$projections[[projections_available]] <- as.data.frame(object@dr[[projections_available]]@cell.embeddings)[,1:2]
    warning("Warning: Only PCA as dimensional reduction found, will export first and second principal components. Consider using tSNE and/or UMAP instead.")
  } else if ( length(projections_available_non_pca) > 0 ) {
    for ( i in projections_available_non_pca ) {
      export$projections[[i]] <- as.data.frame(object@dr[[i]]@cell.embeddings)
    }
  }
  ##--------------------------------------------------------------------------##
  ## cluster tree
  ##--------------------------------------------------------------------------##
  if ( !is.null(object@cluster.tree) ) {
    export$clusters$tree <- object@cluster.tree[[1]]
  }
  ##--------------------------------------------------------------------------##
  ## log-normalized expression
  ##--------------------------------------------------------------------------##
  export$expression <- object@data
  ##--------------------------------------------------------------------------##
  ## save export object to disk
  ##--------------------------------------------------------------------------##
  if ( !file.exists(dirname(file)) ) {
    dir.create(dirname(file), showWarnings = FALSE)
  }
  saveRDS(export, file)
}