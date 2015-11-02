\name{hc_map-GenomicHilbertCurve-method}
\alias{hc_map,GenomicHilbertCurve-method}
\alias{hc_map}
\title{
Draw a map which represents positions of different genomic categories
}
\description{
Draw a map which represents positions of different genomic categories
}
\usage{
\S4method{hc_map}{GenomicHilbertCurve}(object, level = 7, fill = NULL,
    labels = names(object@background), labels_gp = gpar(),
    add = FALSE, ...)
}
\arguments{

  \item{object}{a \code{\link{GenomicHilbertCurve-class}} object}
  \item{level}{Since a map does not need to have too high resolution, a value of around 6 would be enough. It is ignored if \code{add} is set to \code{TRUE}.}
  \item{fill}{colors for differnet genomic categories (current you cannot adjust the style of the borders)}
  \item{labels}{labels of each genomic categories. By default, if the categories is unique for different chromosome,the labels will be the chromosome name or they will be the combination of chromosme names and positions.It is ignored if \code{add} is set to \code{TRUE} and the curve is under 'pixel' mode.}
  \item{labels_gp}{graphic settings for labels}
  \item{add}{whether add the map to the current curve or draw it in a new curve. But notice if \code{add} is set to \code{TRUE},you should set \code{fill} with transparency.}

}
\details{
When multiple genomic categories are draw into one single Hilbert curve, a map which shows the position
of different categories on the curve is necessary to correspond to the graphics on the curve.
}
\examples{
# There is no example
NULL

}
