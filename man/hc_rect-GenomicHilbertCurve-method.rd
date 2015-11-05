\name{hc_rect-GenomicHilbertCurve-method}
\alias{hc_rect,GenomicHilbertCurve-method}
\title{
Add rectangles on Hilbert curve
}
\description{
Add rectangles on Hilbert curve
}
\usage{
\S4method{hc_rect}{GenomicHilbertCurve}(object, gr, gp = gpar(fill = "red", col = "red"),
    mean_mode = c("w0", "absolute", "weighted"))
}
\arguments{

  \item{object}{a \code{\link{GenomicHilbertCurve-class}} object}
  \item{gr}{a \code{\link[GenomicRanges]{GRanges}} object}
  \item{gp}{pass to \code{\link{hc_rect,HilbertCurve-method}}}
  \item{mean_mode}{pass to \code{\link{hc_rect,HilbertCurve-method}}}

}
\details{
It is basically a wrapper of \code{\link{hc_rect,HilbertCurve-method}}.
}
\value{
refer to \code{\link{hc_rect,HilbertCurve-method}}
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
require(circlize)
 bed = generateRandomBed(nr = 100)
 gr = GRanges(seqnames = bed[[1]], ranges = IRanges(bed[[2]], bed[[3]]))
 hc = GenomicHilbertCurve()
 hc_rect(hc, gr, gp = gpar(fill = rand_color(length(gr))))
}
