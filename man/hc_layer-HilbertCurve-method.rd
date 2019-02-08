\name{hc_layer-HilbertCurve-method}
\alias{hc_layer,HilbertCurve-method}
\title{
Add a new layer to the Hilbert curve
}
\description{
Add a new layer to the Hilbert curve
}
\usage{
\S4method{hc_layer}{HilbertCurve}(object, ir, x1 = NULL, x2 = x1, col = "red", border = NA,
    mean_mode = c("w0", "absolute", "weighted"), grid_line = 0,
    grid_line_col = "black", overlay = default_overlay)
}
\arguments{

  \item{object}{A \code{\link{HilbertCurve-class}} object.}
  \item{ir}{an \code{\link[IRanges:IRanges-constructor]{IRanges}} object which specifies the input intervals.}
  \item{x1}{if start positions are not integers, they can be set by \code{x1}.}
  \item{x2}{if end positions are not integers, they can be set by \code{x2}.}
  \item{col}{a scalar or a vector of colors which correspond to intervals in \code{ir} (or \code{x1} and \code{x2}).}
  \item{border}{a scalar or a vector of colors for the borders of intervals. Set it to \code{NA} if borders are suppressed.}
  \item{mean_mode}{Under 'pixel' mode, each pixel represents a small window. This argument provides methods to summarize value for the small window if the input intervals can not completely overlap with the window.  See explanation in \code{\link{hc_points,HilbertCurve-method}}.}
  \item{grid_line}{whether add grid lines to show blocks of the Hilber curve.  It should be an integer number and there will be \code{2^(grid_line-1)-1} horizontal and vertical grid lines.}
  \item{grid_line_col}{color for the grid lines}
  \item{overlay}{a self-defined function which defines how to overlay new layer to the plot. By default it is \code{\link{default_overlay}}. Let's assume the red channel for the layers which are already in the plot is \code{r0}, the red channel for the new layer is \code{r} and the alpha channel is \code{alpha}, the overlayed color is calculated as \code{r*alpha + r0*(1-alpha)}. This self-defined function should accept 7 arguments which are: vectors of r, g, b channels which correspond to the layers that are already in the plot, and r, g, b, alpha channels which corresponds to the new layer. All the  values passed into are between 0 to 1. The returned value for this function should be a list which contains r, g, b channels which correspond to the overlayed colors. Note that these 7 arguments only correspond to the pixels which are covered by the new layer.}

}
\details{
This function only works under 'pixel' mode.

Under "pixel" mode, color is the only graphic representation of values in the input intervals. 
To make a more precise and robust color mapping, users may consider \code{\link[circlize]{colorRamp2}} to create
a color mapping function.

If you want to add more than one layers to the curve, remember to set colors with transparency.

\code{overlay} argument is useful for changing color themes for the overlapped areas, please refer to the vignette
to see examples of how to swith color themes in easy ways.
}
\value{
No value is returned.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
hc = HilbertCurve(1, 100, level = 9, mode = "pixel")

x = sort(sample(100, 20))
s = x[1:10*2 - 1]
e = x[1:10*2]
require(IRanges)
ir = IRanges(s, e)

hc_layer(hc, ir)

hc = HilbertCurve(1, 100, level = 9, mode = "pixel")
hc_layer(hc, ir, grid_line = 3)

hc = HilbertCurve(1, 100, level = 9, mode = "pixel")
hc_layer(hc, ir, border = "black")
}
