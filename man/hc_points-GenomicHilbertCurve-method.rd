\name{hc_points-GenomicHilbertCurve-method}
\alias{hc_points,GenomicHilbertCurve-method}
\title{
Add points to the Hilbert curve
}
\description{
Add points to the Hilbert curve
}
\usage{
\S4method{hc_points}{GenomicHilbertCurve}(object, gr,
    np = max(c(2, 10 - hc_level(object))), size = unit(1, "char"),
    pch = 1, gp = gpar(), mean_mode = c("w0", "absolute", "weighted"),
    shape = "circle")
}
\arguments{

  \item{object}{a \code{\link{GenomicHilbertCurve-class}} object}
  \item{gr}{a \code{\link[GenomicRanges]{GRanges}} object}
  \item{np}{pass to \code{\link{hc_points,HilbertCurve-method}}}
  \item{size}{pass to \code{\link{hc_points,HilbertCurve-method}}}
  \item{pch}{pass to \code{\link{hc_points,HilbertCurve-method}}}
  \item{gp}{pass to \code{\link{hc_points,HilbertCurve-method}}}
  \item{mean_mode}{pass to \code{\link{hc_points,HilbertCurve-method}}}
  \item{shape}{pass to \code{\link{hc_points,HilbertCurve-method}}}

}
\details{
It is basically a wrapper of \code{\link{hc_points,HilbertCurve-method}}.
}
\value{
refer to \code{\link{hc_points,HilbertCurve-method}}
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
require(circlize)
bed = generateRandomBed(nr = 100)
gr = GRanges(seqnames = bed[[1]], ranges = IRanges(bed[[2]], bed[[3]]))
hc = GenomicHilbertCurve()
hc_points(hc, gr, gp = gpar(fill = rand_color(length(gr))))
}
