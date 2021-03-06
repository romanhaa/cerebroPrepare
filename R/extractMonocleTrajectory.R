#' Extract trajectory from Monocle and add to Seurat object.
#' @title Extract trajectory from Monocle and add to Seurat object.
#' @description This function takes a Monocle object, extracts a trajectory that
#' was calculated, and stores it in the specified Seurat object. Trajectory info
#' (state, pseudotime, projection and tree) will be stored in
#' seurat@misc$trajectory under the specified name.
#' @param monocle Monocle object to extract trajectory from.
#' @param seurat Seurat object to transfer trajectory to.
#' @param trajectory_name Name of trajectory.
#' @param column_state Name of meta data column that holds info about the state
#' of a cell; defaults to 'State'.
#' @param column_pseudotime Name of meta data column that holds info about the
#' pseudotime of a cell; defaults to 'Pseudotime'.
#' @return Returns Seurat object with added trajectory. Trajectory info (state,
#' pseudotime, projection and tree) will be stored in seurat@misc$trajectory
#' under the specified name.
#' @keywords seurat monocle trajectory cerebro
#' @export
#' @examples
#' seurat <- extractMonocleTrajectory(
#'   monocle = monocle,
#'   seurat = seurat,
#'   name = 'trajectory_1'
#'   column_state = 'State',
#'   column_pseudotime = 'Pseudotime'
#' )

extractMonocleTrajectory <- function(
    monocle,
    seurat,
    trajectory_name,
    column_state = 'State',
    column_pseudotime = 'Pseudotime'
  ) {

  ##--------------------------------------------------------------------------##
  ## Check if...
  ## - provided Monocle and Seurat objects are of correct type
  ## - required data is present in Monocle object
  ## - number of cells is equal on Monocle and Seurat objects
  ##--------------------------------------------------------------------------##
  if ( !is(monocle, 'CellDataSet') ) {
    stop("The provided object for 'monocle' is not of type 'CellDataSet'.", call. = FALSE)
  }
  if ( !(class(seurat) %in% c('seurat','Seurat')) ) {
    stop("The provided object for 'seurat' doesn't seem to be a Seurat object.", call. = FALSE)
  }
  if ( (column_state %in% colnames(monocle@phenoData@data)) == FALSE ) {
    stop(paste0("Specified column for state info ('", column_state, "') could not be found in meta data."), call. = FALSE)
  }
  if ( (column_pseudotime %in% colnames(monocle@phenoData@data)) == FALSE ) {
    stop(paste0("Specified column for pseudotime info ('", column_pseudotime, "') could not be found in meta data."), call. = FALSE)
  }
  if ( length(monocle@minSpanningTree) == 0 ) {
    stop('monocle@minSpanningTree appears to be empty but is required.', call. = FALSE)
  }
  if ( length(monocle@reducedDimK) == 0 ) {
    stop('monocle@reducedDimK appears to be empty but is required.', call. = FALSE)
  }
  if ( length(monocle@reducedDimS) == 0 ) {
    stop('monocle@reducedDimS appears to be empty but is required.', call. = FALSE)
  }
  if ( !is.null(seurat@misc$trajectory[[trajectory_name]]) ) {
    stop(paste0("Trajectory with specified name ('", trajectory_name, "') already exists in seurat@misc$trajectory. Please choose a different name or manually remove data from that slot."), call. = FALSE)
  }
  if ( nrow(monocle@phenoData@data) > nrow(seurat@meta.data) ) {
    stop(paste0('Number of cells in monocle object (', nrow(monocle@phenoData@data), ') cannot be larger than number of cells in Seurat object (', nrow(seurat@meta.data), ').'), call. = FALSE)
  }
  if ( length(which(rownames(monocle@phenoData@data) %in% rownames(seurat@meta.data))) != nrow(monocle@phenoData@data) ) {
    stop('Some cells provided in the Monocle object are not present in the Seurat object. This is not supported. Please re-calculate the trajectory with all or a subset of the cells present in the Seurat object.', call. = FALSE)
  }
  if ( nrow(monocle@reducedDimK) > 2 ) {
    stop('Currently only two-dimensional reductions are supported.', call. = FALSE)
  }
  if ( nrow(monocle@phenoData@data) < nrow(seurat@meta.data) ) {
    warning(paste0('There are ', nrow(seurat@meta.data) - nrow(monocle@phenoData@data), ' cells present in the Seurat object but not in the Monocle object. Cells without trajectory information will not be visible in Cerebro.'), call. = FALSE)
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
    dplyr::rename(source = 'from', target = 'to') %>%
    dplyr::left_join(
      .,
      reduced_dim_K %>% dplyr::rename(source = 'sample_name', source_dim_1 = 'dim_1', source_dim_2 = 'dim_2'),
      by = 'source'
    ) %>%
    dplyr::left_join(
      .,
      reduced_dim_K %>% dplyr::rename(target = 'sample_name', target_dim_1 = 'dim_1', target_dim_2 = 'dim_2'),
      by = 'target'
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
    dplyr::left_join(
      .,
      trajectory_info %>% dplyr::mutate(cell = base::rownames(.)),
      by = 'cell'
    ) %>%
    dplyr::left_join(
      seurat@meta.data %>% dplyr::mutate(cell = base::rownames(seurat@meta.data)),
      .,
      by = 'cell'
    ) %>%
    dplyr::select(DR_1,DR_2,pseudotime,state,cell)

  base::rownames(trajectory_info) <- trajectory_info$cell
  trajectory_info <- trajectory_info %>% dplyr::select(-cell)

  ##--------------------------------------------------------------------------##
  ## Add trajectory info to Seurat object.
  ##--------------------------------------------------------------------------##
  if ( is.null(seurat@misc$trajectory) ) {
    seurat@misc$trajectory <- list()
  }
  seurat@misc$trajectory$monocle2[[trajectory_name]] <- list(
    meta = trajectory_info,
    edges = edges
  )

  return(seurat)
}


