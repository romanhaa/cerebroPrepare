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
* Number of transcripts (usually created by Seurat by default).
* Number of expressed genes (usually created by Seurat by default).

**Note:** It is recommended to save sample information in a column called `sample` and cluster information in a column called `cluster`. Otherwise, the respective column names need to specified below.

Prepare data:

```
library("cerebroPrepare")
cerebro <- exportFromSeurat(seurat)
saveRDS(cerebro, "cerebro_data.rds")
```

Launch Cerebro and load the RDS file you just exported from R.

To take full advantage of Cerebro, it is recommended to also run the commands below before exporting the data.

```
seurat <- getMostExpressedGenes(seurat)
seurat <- getMarkerGenes(seurat)
seurat <- annotateMarkerGenes(seurat)
```

## To do

* Check if functions work.
* Check if it works with "normal" Seurat object and one from Cerebro pipeline.
* Check how Cerebro behaves with different input data.