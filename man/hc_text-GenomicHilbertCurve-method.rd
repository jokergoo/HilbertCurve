\name{hc_text-GenomicHilbertCurve-method}
\alias{hc_text,GenomicHilbertCurve-method}
\title{
Add text to Hilbert curve
}
\description{
Add text to Hilbert curve
}
\usage{
\S4method{hc_text}{GenomicHilbertCurve}(object, gr, labels, gp = gpar(), ...)
}
\arguments{

  \item{object}{a \code{\link{GenomicHilbertCurve-class}} object}
  \item{gr}{a \code{\link[GenomicRanges]{GRanges}} object which contains the genomic regions to be mapped to the curve}
  \item{labels}{pass to \code{\link{hc_text,HilbertCurve-method}}}
  \item{gp}{pass to \code{\link{hc_text,HilbertCurve-method}}}
  \item{...}{pass to \code{\link{hc_text,HilbertCurve-method}}}

}
\details{
It is basically a wrapper of \code{\link{hc_text,HilbertCurve-method}}.
}
\value{
Refer to \code{\link{hc_text,HilbertCurve-method}}
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
require(circlize)
bed = generateRandomBed(nr = 20)
gr = GRanges(seqnames = bed[[1]], ranges = IRanges(bed[[2]], bed[[3]]))
hc = GenomicHilbertCurve()
hc_text(hc, gr, labels = sample(letters, nrow(bed), replace = TRUE))
}
