\name{hc_segments-GenomicHilbertCurve-method}
\alias{hc_segments,GenomicHilbertCurve-method}
\title{
Add line segments to Hilbert curve
}
\description{
Add line segments to Hilbert curve
}
\usage{
\S4method{hc_segments}{GenomicHilbertCurve}(object, gr, gp = gpar(lty = 1, lwd = 1, col = 1))
}
\arguments{

  \item{object}{a \code{\link{GenomicHilbertCurve-class}} object}
  \item{gr}{a \code{\link[GenomicRanges]{GRanges}} object}
  \item{gp}{pass to \code{\link{hc_segments,HilbertCurve-method}}}

}
\details{
It is basically a wrapper of \code{\link{hc_segments,HilbertCurve-method}}.
}
\value{
refer to \code{\link{hc_segments,HilbertCurve-method}}
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
require(circlize)
 bed = generateRandomBed(nr = 100)
 gr = GRanges(seqnames = bed[[1]], ranges = IRanges(bed[[2]], bed[[3]]))
 hc = GenomicHilbertCurve()
 hc_segments(hc, gr, gp = gpar(col = rand_color(length(gr))))
}
