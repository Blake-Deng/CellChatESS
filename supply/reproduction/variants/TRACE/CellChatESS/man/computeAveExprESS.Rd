\name{computeAveExprESS}
\alias{computeAveExprESS}
\title{Accelerated Group Average Expression}
\usage{
computeAveExprESS(object, features = NULL, group.by = NULL,
  type = "triMean", slot.name = c("data.signaling", "data"),
  data.use = NULL)
}
\arguments{
\item{object}{A CellChat object.}
\item{features}{Optional features to summarize.}
\item{group.by}{Optional metadata column.}
\item{type}{Mean type. Only \code{"triMean"} is currently supported.}
\item{slot.name}{CellChat data slot.}
\item{data.use}{Optional expression matrix.}
}
\value{A feature-by-group average expression matrix.}
\description{
Rcpp implementation of triMean group average expression used by CellChat.
}
