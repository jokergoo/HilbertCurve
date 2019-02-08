\name{hc_map-GenomicHilbertCurve-method}
\alias{hc_map,GenomicHilbertCurve-method}
\alias{hc_map}
\title{
Draw a map which represents positions of different chromosomes on the curve
}
\description{
Draw a map which represents positions of different chromosomes on the curve
}
\usage{
\S4method{hc_map}{GenomicHilbertCurve}(object, level = 7,
    fill = rand_color(length(background), transparency = 0.5), border = NA,
    labels = names(object@background), show_labels = TRUE, labels_gp = gpar(),
    add = FALSE, ...)
}
\arguments{

  \item{object}{a \code{\link{GenomicHilbertCurve-class}} object}
  \item{level}{Since a map does not need to have high resolution, a value of around 7 would be enough.  If \code{add} is set to \code{TRUE}, \code{level} will be enforced to have the same level in the current Hilbert curve.}
  \item{fill}{colors for different chromosomes, or more generally, for different 'seqnames'.}
  \item{border}{colors for the borders of chromosomes. Set it to \code{NA} if borders are suppressed.}
  \item{labels}{label for each chromosome, or more generally, for different 'sequences'}
  \item{show_labels}{whether show text labels}
  \item{labels_gp}{graphic settings for labels}
  \item{add}{whether add the map to the current curve or draw it in a new graphic device. Notice if \code{add} is set to \code{TRUE}, you should set \code{fill} with transparency so that it will not hide your original plot.}
  \item{...}{pass to \code{\link{GenomicHilbertCurve}}. It is only used if you want the map to be plotted in a new graphic device.}

}
\details{
When multiple genomic categories (e.g. chromosomes) are drawn into one single Hilbert curve, a map which shows the positions
of categories on the curve is necessary to distinguish different genomic categories.

Under "pixel" mode, if the map is directly added to the Hilbert curve, no chromosome name is drawn. The chromosome names
are only drawn if the map is plotted in a new graphic device or added to the Hilbert curve under "normal" mode.

Just be careful if you directly overlay the map to the curve that the color of the map does not affect the original
plot too much.
}
\value{
A \code{\link{GenomicHilbertCurve-class}} object
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
# add it in the same graphic device
hc_map(hc, fill = rand_color(24, transparency = 0.5), add = TRUE)

# add the map only with borders
hc = GenomicHilbertCurve()
hc_points(hc, gr, gp = gpar(fill = rand_color(length(gr))))
hc_map(hc, fill = NA, border = "grey", add = TRUE)

# or open a new graphic device
hc_map(hc, fill = rand_color(24))
}
