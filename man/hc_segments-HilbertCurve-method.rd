\name{hc_segments-HilbertCurve-method}
\alias{hc_segments,HilbertCurve-method}
\title{
Add line segments to Hilbert curve
}
\description{
Add line segments to Hilbert curve
}
\usage{
\S4method{hc_segments}{HilbertCurve}(object, ir = NULL, x1 = NULL, x2 = NULL,
    gp = gpar(lty = 1, lwd = 1, col = 1))
}
\arguments{

  \item{object}{A \code{\link{HilbertCurve-class}} object.}
  \item{ir}{an \code{\link[IRanges:IRanges-constructor]{IRanges}} object which specifies the input intervals.}
  \item{x1}{if start positions are not integers, they can be set by \code{x1}.}
  \item{x2}{if end positions are not integers, they can be set by \code{x2}.}
  \item{gp}{graphic parameters for lines. It should be specified by \code{\link[grid]{gpar}}. Note you cannot set \code{linejoin} and \code{lineend}.}

}
\value{
A data frame which contains coordinates (in the 2D space) of segments.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
hc = HilbertCurve(1, 100, level = 4, reference = TRUE)

x = sort(sample(100, 20))
s = x[1:10*2 - 1]
e = x[1:10*2]
require(IRanges)
ir = IRanges(s, e)

hc_segments(hc, ir)
}
