\name{computeCommunProbESS}
\alias{computeCommunProbESS}
\title{Efficient Signal Streaming for CellChat Communication Probability}
\usage{
computeCommunProbESS(object, type = "triMean", raw.use = TRUE,
  population.size = FALSE, nboot = 100, seed.use = 1L, Kh = 0.5,
  n = 1, distance.use = NULL, algorithm = c("ess", "dense"))
}
\arguments{
\item{object}{A preprocessed single-dataset CellChat object.}
\item{type}{Mean type. Only \code{"triMean"} is supported.}
\item{raw.use}{Use \code{object@data.signaling}; otherwise use projected data.}
\item{population.size}{Only \code{FALSE} is supported.}
\item{nboot}{Number of bootstrap permutations.}
\item{seed.use}{Random seed.}
\item{Kh}{Hill function parameter.}
\item{n}{Hill function exponent.}
\item{distance.use}{Only non-spatial workflows are supported.}
\item{algorithm}{\code{"ess"} runs Efficient Signal Streaming; \code{"dense"} runs the dense Rcpp comparator.}
}
\value{A CellChat object containing communication probability and p-value tensors.}
\description{
Computes triMean-based CellChat communication probabilities using ESS or the dense Rcpp comparator while preserving the CellChat object interface.
}
