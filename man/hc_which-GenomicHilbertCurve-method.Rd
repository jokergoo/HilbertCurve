\name{hc_which-GenomicHilbertCurve-method}
\alias{hc_which,GenomicHilbertCurve-method}
\title{
Query regions
}
\description{
Query regions
}
\usage{
\S4method{hc_which}{GenomicHilbertCurve}(object, ix, iy)
}
\arguments{

  \item{object}{a \code{\link{GenomicHilbertCurve-class}} object}
  \item{ix}{A single position on x-axis.}
  \item{iy}{A single position on y-axis.}

}
\details{
Values of \code{ix} and \code{iy} should be integers and take values in [1, 2^level].
}
\value{
A data frame with three columns \code{chr}, \code{start} and \code{end}. The value corresponds to the genomic ranges.
}
\examples{
# There is no example
NULL

}
