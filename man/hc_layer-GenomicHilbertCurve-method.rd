\name{hc_layer-GenomicHilbertCurve-method}
\alias{hc_layer,GenomicHilbertCurve-method}
\title{
Add a new layer to the Hilbert curve
}
\description{
Add a new layer to the Hilbert curve
}
\usage{
\S4method{hc_layer}{GenomicHilbertCurve}(object, gr, col = "red", border = NA,
    mean_mode = c("w0", "absolute", "weighted"), grid_line = 0,
    grid_line_col = "black", overlay = default_overlay)
}
\arguments{

  \item{object}{a \code{\link{GenomicHilbertCurve-class}} object}
  \item{gr}{a \code{\link[GenomicRanges]{GRanges}} object which contains the genomic regions to be mapped to the curve}
  \item{col}{pass to \code{\link{hc_layer,HilbertCurve-method}}}
  \item{border}{colors for the borders of every regions. Set it to \code{NA} if borders are suppressed.}
  \item{mean_mode}{pass to \code{\link{hc_layer,HilbertCurve-method}}}
  \item{grid_line}{pass to \code{\link{hc_layer,HilbertCurve-method}}}
  \item{grid_line_col}{pass to \code{\link{hc_layer,HilbertCurve-method}}}
  \item{overlay}{pass to \code{\link{hc_layer,HilbertCurve-method}}}

}
\details{
It is basically a wrapper of \code{\link{hc_layer,HilbertCurve-method}}.
}
\value{
Refer to \code{\link{hc_layer,HilbertCurve-method}}
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
require(circlize)
bed = generateRandomBed()
gr = GRanges(seqnames = bed[[1]], ranges = IRanges(bed[[2]], bed[[3]]))
hc = GenomicHilbertCurve(mode = "pixel")
hc_layer(hc, gr, col = rand_color(length(gr)))
}
