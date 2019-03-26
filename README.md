# cerebroPrepare

R package with some helper function that prepare single-cell RNA-seq data stored in a Seurat object for visualization in Cerebro.

## Installation

```
library("devtools")
install_github("romanhaa/cerebroPrepare")
```

## How to use

**Required meta data:**

* Experiment name.
* Organism, e.g. 'hg' (human) or 'mm' (mouse).
* Sample.
* Cluster.
* Number of transcripts (usually created by Seurat by default -> nUMI).
* Number of expressed genes (usually created by Seurat by default -> nGene).

**Note:** It is recommended to save sample information in a column called `sample` and cluster information in a column called `cluster`. Otherwise, the respective column names need to specified below.

Prepare data:

```
library("cerebroPrepare")
cerebro <- exportFromSeurat(seurat, "my_experiment.crb")
```

Launch Cerebro and load the RDS file you just exported from R.

To take full advantage of Cerebro, it is recommended to also run the commands below before exporting the data.

```
seurat <- addPercentMtRibo(seurat)
seurat <- getMostExpressedGenes(seurat)
seurat <- getMarkerGenes(seurat)
seurat <- getPathwayEnrichment(seurat)
```

## Credit

* Pathway enrichment in marker genes is done through the enrichR API (<https://github.com/wjawaid/enrichR>). I took the `enrichr` function and modified it to run in parallel (`future_lapply`) and not print status messages.

## Testing

* 3 cases for input:
  * "Naked" Seurat object.
    * Not pre-processed at all, with different column names for sample and column, an extra meta column, and t-SNE and UMAP.
  * Seurat object from our "old" scRNA-seq pipeline.
  * Seurat object from Cerebro pipeline.
    * Should be exactly like the output from the pipeline.
* Other variables to test:
  * Export with and without generating optional data (most expressed genes, marker genes, annotation).

* 4 combinations:
  * object Seurat v2.x + package Seurat v2.x.
  * object Seurat v2.x + package Seurat v3.
  * object Seurat v3 + package Seurat v2.x.
  * object Seurat v3 + package Seurat v3.
