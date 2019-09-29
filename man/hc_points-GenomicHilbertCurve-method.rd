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
  \item{gr}{a \code{\link[GenomicRanges:GRanges-class]{GRanges}} object which contains the genomic regions to be mapped to the curve}
  \item{np}{pass to \code{\link{hc_points,HilbertCurve-method}}}
  \item{size}{size of points when \code{np <= 1}, pass to \code{\link{hc_points,HilbertCurve-method}}}
  \item{pch}{shape of the points when \code{np <= 1}, pass to \code{\link{hc_points,HilbertCurve-method}}}
  \item{gp}{graphic parameters of the points when \code{np <= 1}, pass to \code{\link{hc_points,HilbertCurve-method}}}
  \item{mean_mode}{pass to \code{\link{hc_points,HilbertCurve-method}}}
  \item{shape}{shape of the points when \code{np >= 2}, pass to \code{\link{hc_points,HilbertCurve-method}}}

}
\details{
It is basically a wrapper of \code{\link{hc_points,HilbertCurve-method}}.
}
\value{
Refer to \code{\link{hc_points,HilbertCurve-method}}
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
require(circlize)
require(GenomicRanges)
bed = generateRandomBed(nr = 100)
gr = GRanges(seqnames = bed[[1]], ranges = IRanges(bed[[2]], bed[[3]]))
hc = GenomicHilbertCurve()
hc_points(hc, gr, gp = gpar(fill = rand_color(length(gr))))
}
