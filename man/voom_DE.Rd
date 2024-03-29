% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DifferentialExpressionPipeline_01.07.19.r
\name{voom_DE}
\alias{voom_DE}
\title{Function for Limma Voom differential expression}
\usage{
voom_DE(
  expnData,
  pheno,
  ref,
  percent,
  logCPM = FALSE,
  trend = FALSE,
  normalization = FALSE,
  GOI = NULL,
  eBayesRobust = FALSE,
  lmMethod = "ls"
)
}
\arguments{
\item{expnData}{a matrix or data frame with the raw counts. Patient IDs as colnames, genes as rownames}

\item{pheno}{a character vector with patient IDs as names, and the status for each in each group(eg pos,neg)}

\item{ref}{a character vector of the reference level for DE. for example ref="No".}

\item{percent}{percent is the fraction (0-1 numberic) of AML samples to include when setting an expression threshold. eg 5\% of AMLs, percent=0.05}

\item{logCPM}{boolean. use logCPM values or not.}

\item{trend}{trend is for using limma trend method with log2 CPMs}

\item{normalization}{method of normalization such as quantile if necessary. should be either FALSE or "qt" so far}

\item{GOI}{a character vector of genes (or numeric vector of row indices) of interest to subset at the end. keeps BH adjuted p-values more accurate.}

\item{eBayesRobust}{boolean}

\item{lmMethod}{least squares (ls)}
}
\value{
list
}
\description{
Function for Limma Voom differential expression
}
\examples{
ex <- c('TBD')

}
