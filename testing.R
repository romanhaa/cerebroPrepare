##----------------------------------------------------------------------------##
##
##----------------------------------------------------------------------------##
source("/hpcnfs/scratch/PGP/rhillje/cerebroPrepare/R/getMostExpressedGenes.R")
source("/hpcnfs/scratch/PGP/rhillje/cerebroPrepare/R/getMarkerGenes.R")
source("/hpcnfs/scratch/PGP/rhillje/cerebroPrepare/R/annotateMarkerGenes.R")
source("/hpcnfs/scratch/PGP/rhillje/cerebroPrepare/R/exportFromSeurat.R")

##----------------------------------------------------------------------------##
## Output from Cerebro as input.
## - experiment name already present
## - info about organism already present
## - column names for sample, cluster, nUMI, and nGene correctly named
## - previous cell cycle results not carried over by default, should be fine like this
## - most expressed genes already saved
## - marker genes already identified
## - pathway enrichment already done
##----------------------------------------------------------------------------##
seurat <- readRDS("/hpcnfs/scratch/PGP/rhillje/cerebroPrepare_test_data/from_cerebro/seurat.rds")
t <- exportFromSeurat(
  seurat,
  experiment_name = "from_cerebro",
  organism = "mm",
  column_cell_cycle_regev = "cell_cycle_Regev",
  column_cell_cycle_cyclone = "cell_cycle_Cyclone"
)
saveRDS(t, "/hpcnfs/scratch/PGP/rhillje/cerebroPrepare_test_data/from_cerebro/cerebro_from_cerebro.rds")

##----------------------------------------------------------------------------##
## Naked Seurat object as input.
## - experiment name and organism info not present
## - column names for sample and cluster are different
## - most expressed genes not present
## - marker genes not present
## - enriched pathways not present
## - cell cycle results not present at all
##----------------------------------------------------------------------------##
seurat <- readRDS("/hpcnfs/scratch/PGP/rhillje/cerebroPrepare_test_data/naked/seurat.rds")
seurat <- getMostExpressedGenes(
  seurat,
  column_sample = "sampleID",
  column_cluster = "clusterID"
)
seurat <- getMarkerGenes(
  seurat,
  column_sample = "sampleID",
  column_cluster = "clusterID",
  organism = "mm"
)
seurat <- annotateMarkerGenes(
  seurat,
  column_sample = "sampleID",
  column_cluster = "clusterID"
)
t <- exportFromSeurat(
  seurat,
  experiment_name = "naked_seurat",
  organism = "mm",
  column_sample = "sampleID",
  column_cluster = "clusterID"
)
saveRDS(t, "/hpcnfs/scratch/PGP/rhillje/cerebroPrepare_test_data/naked/cerebro_from_naked.rds")

# alternative
# seurat <- readRDS("/hpcnfs/scratch/PGP/rhillje/cerebroPrepare_test_data/naked/seurat.rds")
# seurat@meta.data$sample <- seurat@meta.data$sampleID
# seurat@meta.data$cluster <- seurat@meta.data$clusterID
# seurat <- getMostExpressedGenes(seurat)
# seurat <- getMarkerGenes(seurat, organism = "mm")
# seurat <- annotateMarkerGenes(seurat)
# t <- exportFromSeurat(seurat, experiment_name = "naked_seurat", organism = "mm")
# saveRDS(t, "/hpcnfs/scratch/PGP/rhillje/cerebroPrepare_test_data/naked/cerebro_from_naked.rds")

##----------------------------------------------------------------------------##
## Output from scRNAseq pipeline as input.
## - experiment name and organism stored in different place
## - column with sample information named differently
## - cell cycle info stored under different column names
## - most expressed genes and marker genes available but stored in different places
##   - convert are call again?
##     - problem: all '.' are replaced with '_', also in column names
##     - conversion might be laborious
## - enriched pathways may be present, not always
##----------------------------------------------------------------------------##
seurat <- readRDS("/hpcnfs/scratch/PGP/rhillje/cerebroPrepare_test_data/from_scRNAseq/seurat.rds")
seurat <- getMostExpressedGenes(
  seurat,
  column_sample = "sampleID"
)
seurat <- getMarkerGenes(
  seurat,
  column_sample = "sampleID",
  organism = "mm"
)
seurat <- annotateMarkerGenes(
  seurat,
  column_sample = "sampleID"
)
t <- exportFromSeurat(
  seurat,
  experiment_name = "from_scRNAseq",
  organism = "mm",
  column_sample = "sampleID",
  column_cell_cycle_regev = "cell.cycle.Regev",
  column_cell_cycle_cyclone = "cell.cycle.Cyclone"
)
saveRDS(t, "/hpcnfs/scratch/PGP/rhillje/cerebroPrepare_test_data/from_scRNAseq/cerebro_from_scRNAseq.rds")












