library("monocle")
library("dplyr")

monocle <- readRDS("~/Desktop/monocle.rds")

extractTrajectory <- function(object) {
  # required data:
  # - monocle object with
  #   - "State" and "Pseudotime" in phenoData slot
  #   - calculated trajectory
  # required packages:
  # - monocle
  # - Matrix
  # - dplyr
  # - tibble
  # - igraph

  # CHECK IF REQUIRED DATA IS AVAILABLE
  # CHECK IF CELL NUMBER IS EQUAL TO OTHER META DATA IN CEREBRO FILE

  monocle <- object

  min_spanning_tree = monocle::minSpanningTree(monocle)
  reduced_dim_K <- Matrix::t(monocle::reducedDimK(monocle)) %>%
    as.data.frame() %>%
    dplyr::rename(dim_1 = 1, dim_2 = 2) %>%
    dplyr::mutate(sample_name = rownames(.))

  edges <- min_spanning_tree %>%
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
    select(source, target, weight, source_dim_1, source_dim_2, target_dim_1, target_dim_2)

  pheno_data <- data.frame(
    Pseudotime = monocle@phenoData@data$Pseudotime,
    State = monocle@phenoData@data$State,
    row.names = rownames(monocle@phenoData@data)
  )

  pheno_data <- t(monocle::reducedDimS(monocle)) %>%
    as.data.frame() %>%
    dplyr::rename(DR_1 = 1, DR_2 = 2) %>%
    mutate(cell = rownames(.)) %>%
    dplyr::left_join(pheno_data %>% mutate(cell = rownames(.)), by = "cell")

  rownames(pheno_data) <- pheno_data$cell
  pheno_data <- pheno_data %>% select(-cell)

  return(
    list(
      pheno_data = pheno_data,
      edges = edges
    )
  )
}

t <- extractTrajectory(monocle)

p <- ggplot() +
  geom_point(data = t$pheno_data, aes(DR_1, DR_2, color = State)) +
  geom_segment(
    data = t$edges,
    aes(source_dim_1, source_dim_2, xend = target_dim_1, yend = target_dim_2),
    size = 0.75, linetype = "solid", na.rm = TRUE
  ) +
  theme_bw()
ggsave("~/Desktop/by_cell_cycle.pdf", p, height = 6, width = 8)
ggsave("~/Desktop/by_cell_cycle.png", p, height = 6, width = 8)


