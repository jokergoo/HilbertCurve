\name{hc_polygon-GenomicHilbertCurve-method}
\alias{hc_polygon,GenomicHilbertCurve-method}
\title{
Add text to Hilbert curve
}
\description{
Add text to Hilbert curve
}
\usage{
\S4method{hc_polygon}{GenomicHilbertCurve}(object, gr, gp = gpar(),
    end_type = c("average", "expanding", "shrinking"))
}
\arguments{

  \item{object}{a \code{\link{GenomicHilbertCurve-class}} object}
  \item{gr}{a \code{\link[GenomicRanges:GRanges-class]{GRanges}} object which contains the genomic regions to be mapped to the curve}
  \item{gp}{pass to \code{\link{hc_polygon,HilbertCurve-method}}}
  \item{end_type}{pass to \code{\link{hc_polygon,HilbertCurve-method}}}

}
\details{
It is basically a wrapper of \code{\link{hc_polygon,HilbertCurve-method}}.
}
\value{
Refer to \code{\link{hc_polygon,HilbertCurve-method}}
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
require(circlize)
require(GenomicRanges)
bed = generateRandomBed(nr = 20)
gr = GRanges(seqnames = bed[[1]], ranges = IRanges(bed[[2]], bed[[3]]))
hc = GenomicHilbertCurve()
hc_polygon(hc, gr)
}
