#' Extract trajectory from Monocle and add to Seurat object.
#' @title Extract trajectory from Monocle and add to Seurat object.
#' @description This function takes a Monocle object, extracts a trajectory that was calculated, and stores it in the specified Seurat object. Trajectory info (state, pseudotime, projection and tree) will be stored in seurat@misc$trajectory under the specified name.
#' @param monocle Monocle object to extract trajectory from.
#' @param seurat Seurat object to store trajectory in.
#' @param trajectory_name Name of trajectory.
#' @param column_state Name of meta data column that holds info about the state of a cell; defaults to "State".
#' @param column_pseudotime Name of meta data column that holds info about the pseudotime of a cell; defaults to "Pseudotime".
#' @return Returns Seurat object with added trajectory. Trajectory info (state, pseudotime, projection and tree) will be stored in seurat@misc$trajectory under the specified name.
#' @keywords seurat monocle trajectory cerebro
#' @export
#' @examples
#' extractMonocleTrajectory(monocle = monocle, seurat = seurat, name = "trajectory_1")

extractMonocleTrajectory <- function(
    monocle,
    seurat,
    trajectory_name,
    column_state = "State",
    column_pseudotime = "Pseudotime"
  ) {

  ##--------------------------------------------------------------------------##
  ## Check if...
  ## - provided Monocle and Seurat objects are of correct type
  ## - required data is present in Monocle object
  ## - number of cells is equal on Monocle and Seurat objects
  ##--------------------------------------------------------------------------##
  if ( !is(monocle, "CellDataSet") ) {
    stop('Object "monocle" is not of type "CellDataSet".', call. = FALSE)
  }
  if ( !is(seurat, "Seurat") ) {
    stop('Object "seurat" is not of type "Seurat".', call. = FALSE)
  }
  if ( (column_state %in% colnames(monocle@phenoData@data)) == FALSE ) {
    stop(paste0("Specified column for state info ('", column_state, "') could not be found in meta data."), call. = FALSE)
  }
  if ( (column_pseudotime %in% colnames(monocle@phenoData@data)) == FALSE ) {
    stop(paste0("Specified column for pseudotime info ('", column_pseudotime, "') could not be found in meta data."), call. = FALSE)
  }
  if ( length(monocle@minSpanningTree) == 0 ) {
    stop("monocle@minSpanningTree appears to be empty but is required.", call. = FALSE)
  }
  if ( length(monocle@reducedDimK) == 0 ) {
    stop("monocle@reducedDimK appears to be empty but is required.", call. = FALSE)
  }
  if ( length(monocle@reducedDimS) == 0 ) {
    stop("monocle@reducedDimS appears to be empty but is required.", call. = FALSE)
  }
  if ( nrow(monocle@phenoData@data) != nrow(seurat@meta.data) ) {
    stop(paste0("Number of cells in monocle object (", nrow(monocle@phenoData@data), ") must be equal to number of cells in Seurat object (", nrow(seurat@meta.data), ")."), call. = FALSE)
  }
  if ( !is.null(seurat@misc$trajectory[[trajectory_name]]) ) {
    stop(paste0("Trajectory with specified name ('", trajectory_name, "') already exists in seurat@misc$trajectory. Please choose a different name or manually remove data from that slot."), call. = FALSE)
  }

  ##--------------------------------------------------------------------------##
  ## Prepare reducedDimK.
  ##--------------------------------------------------------------------------##
  reduced_dim_K <- base::t(monocle@reducedDimK) %>%
    base::as.data.frame() %>%
    dplyr::rename(dim_1 = 1, dim_2 = 2) %>%
    dplyr::mutate(sample_name = base::rownames(.))

  ##--------------------------------------------------------------------------##
  ## Transform trajectory into plottable edges.
  ##--------------------------------------------------------------------------##
  edges <- monocle@minSpanningTree %>%
    igraph::as_data_frame() %>%
    dplyr::rename(source = "from", target = "to") %>%
    dplyr::left_join(
      reduced_dim_K %>% dplyr::rename(source = "sample_name", source_dim_1 = "dim_1", source_dim_2 = "dim_2"),
      by = "source"
    ) %>%
    dplyr::left_join(
      reduced_dim_K %>% dplyr::rename(target = "sample_name", target_dim_1 = "dim_1", target_dim_2 = "dim_2"),
      by = "target"
    ) %>%
    dplyr::select(source, target, weight, source_dim_1, source_dim_2, target_dim_1, target_dim_2)

  ##--------------------------------------------------------------------------##
  ## Extract state and pseudotime info and cell position in projection.
  ##--------------------------------------------------------------------------##
  trajectory_info <- base::data.frame(
    pseudotime = monocle@phenoData@data[[column_pseudotime]],
    state = monocle@phenoData@data[[column_state]],
    row.names = base::rownames(monocle@phenoData@data)
  )

  trajectory_info <- base::t(monocle@reducedDimS) %>%
    base::as.data.frame() %>%
    dplyr::rename(DR_1 = 1, DR_2 = 2) %>%
    dplyr::mutate(cell = base::rownames(.)) %>%
    dplyr::left_join(trajectory_info %>% dplyr::mutate(cell = base::rownames(.)), by = "cell")

  base::rownames(trajectory_info) <- trajectory_info$cell
  trajectory_info <- trajectory_info %>% dplyr::select(-cell)

  trajectory_info[match(rownames(seurat@meta.data), rownames(trajectory_info)),]

  ##--------------------------------------------------------------------------##
  ## Add trajectory info to Seurat object.
  ##--------------------------------------------------------------------------##
  if ( is.null(seurat@misc$trajectory) ) {
    seurat@misc$trajectory <- list()
  }
  seurat@misc$trajectory[[trajectory_name]] <- list(
    meta = trajectory_info,
    edges = edges
  )

  return(seurat)
}


