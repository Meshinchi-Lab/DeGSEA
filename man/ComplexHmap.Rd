% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Heatmaps_Function.r
\name{ComplexHmap}
\alias{ComplexHmap}
\title{Create ComplexHeatmap plot}
\usage{
ComplexHmap(
  mat,
  name = "z-scores",
  scale = TRUE,
  threshold = FALSE,
  hmap_anno_obj,
  hmap_anno_obj_genes = NULL,
  space.type = "sRGB",
  color_palette = NULL,
  split = NULL,
  cluster.method = "ward.D2",
  dge_dendrograms.res = NULL,
  samp_dend_order = NULL
)
}
\arguments{
\item{mat}{the normalized, log2 (usually) transformed counts}

\item{name}{is the title}

\item{scale}{boolean whether to scale by row}

\item{threshold}{whether to make all z-scores in a certain range.}

\item{hmap_anno_obj}{from HeatmapAnnotation() function}

\item{hmap_anno_obj_genes}{from HeatmapAnnotation() function}

\item{space.type}{value for color palette, like sRGB or LAB}

\item{color_palette}{colors vector from circlize::colorRamp2}

\item{split}{data.frame with the groups identifying which rows to split.}

\item{cluster.method}{method for dendrogram.}

\item{dge_dendrograms.res}{the object from DeGSEA::dge_dendrograms()}

\item{samp_dend_order}{numeric vector or character vector of column names from the matrix (mat) or the dge_dengrograms.res$TMMCPM matix, in the desired order.}
}
\value{
complexHeatmap object
}
\description{
Create ComplexHeatmap plot
}
\examples{
ex <- c('TBD')
}
