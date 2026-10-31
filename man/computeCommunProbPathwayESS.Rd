\name{computeCommunProbPathwayESS}
\alias{computeCommunProbPathwayESS}
\title{Accelerated CellChat Pathway Aggregation}
\usage{
computeCommunProbPathwayESS(object = NULL, net = NULL, pairLR.use = NULL,
  thresh = 0.05)
}
\arguments{
\item{object}{A CellChat object.}
\item{net}{Optional CellChat network list.}
\item{pairLR.use}{Optional ligand-receptor table.}
\item{thresh}{P-value threshold.}
}
\value{A CellChat object, or a list when \code{object} is \code{NULL}.}
\description{
Rcpp helper for aggregating ligand-receptor communication probabilities by pathway. Entries at the threshold are retained. Output dimensions follow the pinned CellChat reference, including dropping the pathway dimension when a single pathway remains.
}
