![LICENSE.md](https://img.shields.io/github/license/romanhaa/cerebroPrepare)
![https://twitter.com/fakechek1](https://img.shields.io/twitter/follow/fakechek1?label=Twitter&style=social)

# cerebroPrepare

R package with some helper function that prepare single-cell RNA-seq data stored in a Seurat object for visualization in Cerebro.
Both Seurat v2 and Seurat v3 are supported.

## Installation

```r
library('BiocManager')
BiocManager::install('romanhaa/cerebroPrepare')
```

## How to use

**Required meta data:**

* Experiment name.
* Organism, e.g. 'hg' (human) or 'mm' (mouse).
* Sample.
* Cluster.
* Number of transcripts (usually created by Seurat by default; `nUMI` / `nCounts_RNA` in Seurat v2 and v3).
* Number of expressed genes (usually created by Seurat by default; `nGene` / `nFeature_RNA` in Seurat v2 and v3).

**Note:** It is recommended to save sample information in a column called `sample` and cluster information in a column called `cluster`. Otherwise, the respective column names need to specified below.

Prepare data:

```r
library('cerebroPrepare')
cerebro <- exportFromSeurat(seurat, 'my_experiment.crb')
```

Launch Cerebro and load the RDS file you just exported from R.

To take full advantage of Cerebro, it is recommended to also run the commands below before exporting the data.

```r
seurat <- addPercentMtRibo(seurat)
seurat <- getMostExpressedGenes(seurat)
seurat <- getMarkerGenes(seurat)
seurat <- getEnrichedPathways(seurat)
```

## Further details

For further details, please refer to the official [Cerebro](https://github.com/romanhaa/Cerebro) repository.

## Credit

* Pathway enrichment in marker genes is done through the enrichR API (<https://github.com/wjawaid/enrichR>). I took the `enrichr` function and modified it to run in parallel (`future_lapply`) and not print status messages.
