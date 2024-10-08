\name{hc_which-HilbertCurve-method}
\alias{hc_which,HilbertCurve-method}
\title{
Query regions
}
\description{
Query regions
}
\usage{
\S4method{hc_which}{HilbertCurve}(object, ix, iy)
}
\arguments{

  \item{object}{A \code{\link{HilbertCurve-class}} object.}
  \item{ix}{A single position on x-axis.}
  \item{iy}{A single position on y-axis.}

}
\details{
Values of \code{ix} and \code{iy} should be integers and take values in [1, 2^level].
}
\value{
A data frame with two columns \code{start} and \code{end}. The value corresponds to the range in data.
}
\examples{
# There is no example
NULL

}
