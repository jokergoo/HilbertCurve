\name{hc_segments-HilbertCurve-method}
\alias{hc_segments,HilbertCurve-method}
\title{
Add line segments to Hilbert curve
}
\description{
Add line segments to Hilbert curve
}
\usage{
\S4method{hc_segments}{HilbertCurve}(object, ir, x1 = NULL, x2 = NULL,
    gp = gpar(lty = 1, lwd = 1, col = 1))
}
\arguments{

  \item{object}{A \code{\link{HilbertCurve-class}} object.}
  \item{ir}{a \code{\link[IRanges]{IRanges}} object.}
  \item{x1}{if positions are not integers, they can be set by \code{x1} and \code{x2}.}
  \item{x2}{if positions are not integers, they can be set by \code{x1} and \code{x2}.}
  \item{gp}{graphical parameters for lines. It should be specified by \code{\link[grid]{gpar}}.}

}
\value{
A data frame which contains coordinates for segments.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
hc = HilbertCurve(1, 100, level = 4, reference = TRUE)
x = sort(sample(100, 20))
 s = x[1:10*2 - 1]
 e = x[1:10*2]
 ir = IRanges(s, e)
hc_segments(hc, ir)
}
