\name{computeCommunProbESS}
\alias{computeCommunProbESS}
\title{Efficient Signal Streaming for CellChat Communication Probability}
\usage{
computeCommunProbESS(object, type = "triMean", raw.use = TRUE,
  population.size = FALSE, nboot = 100, seed.use = 1L, Kh = 0.5,
  n = 1, distance.use = NULL, threads = 1L,
  algorithm = c("ess", "dense"))
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
\item{threads}{Positive integer controlling native threads for
\code{algorithm = "ess"}. One uses the serial fast path; values greater than
one use OpenMP when available. The effective count is capped by \code{nboot}.}
\item{algorithm}{\code{"ess"} runs Efficient Signal Streaming; \code{"dense"} runs the dense Rcpp comparator.}
}
\details{
The signaling expression matrix must be a dense numeric matrix or a numeric sparse \pkg{Matrix} representation used by CellChat, and must contain only finite, non-negative values. Inputs containing negative values, \code{NA}, \code{NaN}, or infinite values are rejected before inference. This restriction is required for expression-support pruning to preserve zero contributions. Sparse inputs are checked through stored values without densification; referenced rows are subsequently projected and converted to dense native input. Cell identities must match expression columns, required complex subunits must be present, and gene names must be unique. Direct genes use their group summaries; named complexes retain the geometric mean, including single-subunit complexes. Only the declared non-spatial, single-dataset RNA subset is supported.
}
\value{A CellChat object containing communication probability and p-value tensors.}
\description{
Computes triMean-based CellChat communication probabilities using ESS or the dense Rcpp comparator while preserving the CellChat object interface.
}
