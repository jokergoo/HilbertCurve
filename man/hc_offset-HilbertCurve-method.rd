\name{hc_offset-HilbertCurve-method}
\alias{hc_offset,HilbertCurve-method}
\alias{hc_offset}
\title{
Adjust positions
}
\description{
Adjust positions
}
\usage{
\S4method{hc_offset}{HilbertCurve}(object, x)
}
\arguments{

  \item{object}{A \code{\link{HilbertCurve-class}} object.}
  \item{x}{positions.}

}
\details{
Since internally positions are transformed to positive integers, if input positions
are specified as negative values when initializing the Hilbert curve, a shift will be recorded
internally and positions are transformed to positive value automatically.
}
\value{
A positive numeric value.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
hc = HilbertCurve(-100, 100)
hc_offset(hc, c(-100, -50, 0, 50, 100))
}
