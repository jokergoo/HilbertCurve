\name{HilbertCurve}
\alias{HilbertCurve}
\title{
Initialize a Hilbert curve
}
\description{
Initialize a Hilbert curve
}
\usage{
HilbertCurve(s, e, level = 4, mode = c("normal", "pixel"),
    reference = FALSE, reference_gp = gpar(lty = 3, col = "#999999"),
    arrow = TRUE, zoom = NULL, newpage = TRUE,
    background = "white", title = NULL, title_gp = gpar(fontsize = 16),
    start_from = c("bottomleft", "topleft", "bottomright", "topright"),
    legend = list())
}
\arguments{

  \item{s}{position that will be mapped to the start of the Hilbert curve.}
  \item{e}{position that will be mapped to the end of the Hilbert curve.}
  \item{level}{level of the Hilbert curve. There will by \code{4^level - 1} segments in the curve.}
  \item{mode}{"normal" mode is used for low \code{level} value and "pixel" mode is always used for high \code{level} value, so the "normal" mode is always for low-resolution visualization while "pixel" mode is used for high-resolution visualization. See 'details' for explanation.}
  \item{reference}{whether add reference line on the plot. Only works under 'normal' mode. The reference line is only used for letting users know how the curve folds.}
  \item{reference_gp}{graphic settings for the reference line. It should be specified by \code{\link[grid]{gpar}}.}
  \item{arrow}{whether add arrows on the reference line. Only works under 'normal' mode.}
  \item{zoom}{Internally, position are stored as integer values. To better map the data to the Hilbert curve,  the original positions are zoomed according to the range and the level of Hilbert curve. E.g. if  the curve visualizes data ranging from 1 to 2 but level of the curve is set to 4, the positions will be zoomed by ~x2000 so that values like 1.5, 1.555 can be mapped to the curve with more accuracy. You don't need to care the zooming thing,  proper zooming factor is calculated automatically.}
  \item{newpage}{whether call \code{\link[grid]{grid.newpage}} to draw on a new graphic device.}
  \item{background}{background color, currently please don't modify this value.}
  \item{title}{title of the plot.}
  \item{title_gp}{graphic parameters for the title. It should be specified by \code{\link[grid]{gpar}}.}
  \item{start_from}{which corner on the plot should the curve starts?}
  \item{legend}{a \code{\link[grid]{grob}} object or a list of \code{\link[grid]{grob}} objects. You can construct a \code{\link[ComplexHeatmap]{ColorMapping-class}} object and generate a legend, see in the Example section.}

}
\details{
This funciton initializes a Hilbert curve with level \code{level} which corresponds 
to the range between \code{s} and \code{e}.

Under 'normal' mode, there is a visible Hilbert curve which plays like a folded axis and
different low-level graphics can be added on according to the coordinate. 
It only works nice if the level of the Hilbert curve is small (say less than 6).

When the level is high (e.g. > 10), the whole 2D space will be almost completely filled by the curve and
it is impossible to add or visualize e.g. points on the curve. In this case, the 'pixel'
mode visualizes each tiny 'segment' as a pixel and maps values to colors. So the Hilbert
curve with level 11 will generate a PNG figure with 2048x2048 resolution. This is extremely
useful for visualize genomic data. E.g. If we make a Hilbert curve for human chromosome 1 with
level 11, then each pixel can represent 60bp (\code{249250621/2048/2048}) which is of very high resolution.

Under 'pixel' mode, if the current device is an interactive deivce, every time a new layer is added, 
the image will be add to the interactive device as a rastered image. But still you can use \code{\link{hc_png,HilbertCurve-method}}
to export the plot as PNG file.

Notice, \code{s} and \code{e} are not necessarily to be integers, it can be any values (e.g. numeric or even negative values).
}
\value{
A \code{\link{HilbertCurve-class}} object.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
HilbertCurve(1, 100, reference = TRUE)
HilbertCurve(1, 100, level = 5)
HilbertCurve(1, 100, title = "title")
HilbertCurve(1, 100, start_from = "topleft")

require(ComplexHeatmap)
cm = ColorMapping(colors = c("red", "blue"), levels = c("a", "b"))
legend = color_mapping_legend(cm, plot = FALSE, title = "foo")
hc = HilbertCurve(1, 100, title = "title", legend = legend)
hc_segments(hc, x1 = 20, x2 = 40)
}
