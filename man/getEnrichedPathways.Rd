% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getEnrichedPathways.R
\name{getEnrichedPathways}
\alias{getEnrichedPathways}
\title{Get enriched pathways based on marker genes from EnrichR.}
\usage{
getEnrichedPathways(object, column_sample = "sample",
  column_cluster = "cluster",
  databases = c("GO_Biological_Process_2018",
  "GO_Cellular_Component_2018", "GO_Molecular_Function_2018", "KEGG_2016",
  "WikiPathways_2016", "Reactome_2016", "Panther_2016", "Human_Gene_Atlas",
  "Mouse_Gene_Atlas"), adj_p_cutoff = 0.05, max_terms = 100,
  URL_API = "http://amp.pharm.mssm.edu/Enrichr/enrich")
}
\arguments{
\item{object}{Seurat object.}

\item{column_sample}{Column in object@meta.data that contains information
about sample; defaults to 'sample'.}

\item{column_cluster}{Column in object@meta.data that contains information
about cluster; defaults to 'cluster'.}

\item{databases}{Which databases to query. Use enrichR::listEnrichrDbs() to
check what databases are available.}

\item{adj_p_cutoff}{Cut-off for adjusted p-value of enriched pathways;
defaults to 0.05,}

\item{max_terms}{Save only first n entries of each database; defaults to 100.}

\item{URL_API}{URL to send requests to (Enrichr API). Allows to overwrite
default URL with an alternative taken from the Enrichr website in case the
original is out-of-service; defaults to
'http://amp.pharm.mssm.edu/Enrichr/enrich'.}
}
\description{
This function uses the enrichR API to look for enriched pathways
in marker gene sets of samples and clusters.
}
\details{
Get enriched pathways based on marker genes from EnrichR.
}
\examples{
seurat <- getEnrichedPathways(
  object = seurat,
  column_sample = 'sample',
  column_cluster = 'cluster',
  databases = c('GO_Biological_Process_2018','GO_Cellular_Component_2018'),
  adj_p_cutoff = 0.01,
  max_terms = 100,
  URL_API = 'http://amp.pharm.mssm.edu/Enrichr/enrich'
)
}
\keyword{cerebro}
\keyword{seurat}
