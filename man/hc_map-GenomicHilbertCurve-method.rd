\name{hc_map-GenomicHilbertCurve-method}
\alias{hc_map,GenomicHilbertCurve-method}
\alias{hc_map}
\title{
Draw a map which represents positions of different genomic categories on the curve
}
\description{
Draw a map which represents positions of different genomic categories on the curve
}
\usage{
\S4method{hc_map}{GenomicHilbertCurve}(object, level = 7, fill = NULL,
    labels = names(object@background), labels_gp = gpar(),
    add = FALSE, ...)
}
\arguments{

  \item{object}{a \code{\link{GenomicHilbertCurve-class}} object}
  \item{level}{Since a map does not need to have high resolution, a value of around 6 would be enough.  If \code{add} is set to \code{TRUE}, \code{level} will be enforced to have the same level in the current Hilbert curve.}
  \item{fill}{colors for different genomic categories (current you cannot adjust the style of the borders)}
  \item{labels}{label for each genomic category. By default, if the category is unique for different chromosome, the label will be the chromosome name, or else they will be the combination of chromosme names and positions. It is ignored if the curve is under 'pixel' mode.}
  \item{labels_gp}{graphic settings for labels}
  \item{add}{whether add the map to the current curve or draw it in a new graphic device. Notice if \code{add} is set to \code{TRUE}, you should set \code{fill} with transparency so that it will not hide your original plot.}
  \item{...}{pass to \code{\link{GenomicHilbertCurve}}.}

}
\details{
When multiple genomic categories are drawn into one single Hilbert curve, a map which shows the positions
of different categories on the curve is necessary to correspond to the graphics on the curve.
}
\value{
A \code{\link{GenomicHilbertCurve-class}} object
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
# add it in the same graphic device
hc_map(hc, fill = rand_color(24, transparency = 0.5), add = TRUE)

# or open a new graphic device
hc_map(hc, fill = rand_color(24))
}
