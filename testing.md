## Test cases

* 3 cases for input:
  * "Naked" Seurat object.
    * Not pre-processed at all, with different column names for sample and column, an extra meta column, and t-SNE and UMAP.
  * Seurat object from our "old" scRNA-seq pipeline.
  * Seurat object from Cerebro pipeline.
    * Should be exactly like the output from the pipeline.
* Other variables to test:
  * Export with and without generating optional data (most expressed genes, marker genes, annotation).

* In principle there are at least 4 combinations to be tested:
  * object Seurat v2.x + package Seurat v2.x.
  * object Seurat v2.x + package Seurat v3.
  * object Seurat v3 + package Seurat v2.x.
  * object Seurat v3 + package Seurat v3.
* For simplicity, I will only do v2 + v2 and v3 + v3.
