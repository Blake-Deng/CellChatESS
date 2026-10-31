\name{aggregateNetESS}
\alias{aggregateNetESS}
\title{Accelerated CellChat Network Aggregation}
\usage{
aggregateNetESS(object, thresh = 0.05)
}
\arguments{
\item{object}{A CellChat object with \code{object@net$prob} and \code{object@net$pval}.}
\item{thresh}{P-value threshold.}
}
\value{A CellChat object with \code{object@net$count} and \code{object@net$weight}.}
\description{
Rcpp implementation of the default \code{CellChat::aggregateNet()} branch.
}
