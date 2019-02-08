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
  \item{gr}{a \code{\link[GenomicRanges:GRanges-class]{GRanges}} object which contains the genomic regions to be mapped to the curve}
  \item{col}{a scalar or a vector of colors which correspond to regions in \code{gr}, pass to \code{\link{hc_layer,HilbertCurve-method}}}
  \item{border}{a scalar or a vector of colors which correspond to the borders of regions. Set it to \code{NA} if borders are suppressed.}
  \item{mean_mode}{Under 'pixel' mode, each pixel represents a small window. This argument provides methods to summarize value for the small window if the input genomic regions can not completely overlap with the window,  pass to \code{\link{hc_layer,HilbertCurve-method}}}
  \item{grid_line}{whether add grid lines to show blocks of the Hilber curve, pass to \code{\link{hc_layer,HilbertCurve-method}}}
  \item{grid_line_col}{color for the grid lines, pass to \code{\link{hc_layer,HilbertCurve-method}}}
  \item{overlay}{a self-defined function which defines how to overlay new layer to the plot, pass to \code{\link{hc_layer,HilbertCurve-method}}}

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
require(GenomicRanges)
bed = generateRandomBed()
gr = GRanges(seqnames = bed[[1]], ranges = IRanges(bed[[2]], bed[[3]]))
hc = GenomicHilbertCurve(mode = "pixel", level = 9)
hc_layer(hc, gr, col = rand_color(length(gr)))
}
