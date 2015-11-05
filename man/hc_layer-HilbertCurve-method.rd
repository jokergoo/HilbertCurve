\name{hc_layer-HilbertCurve-method}
\alias{hc_layer,HilbertCurve-method}
\title{
Add a new layer to the Hilbert curve
}
\description{
Add a new layer to the Hilbert curve
}
\usage{
\S4method{hc_layer}{HilbertCurve}(object, ir, x1 = NULL, x2 = x1, col = "red",
    mean_mode = c("w0", "absolute", "weighted"), grid_line = 0,
    grid_line_col = "black", overlay = default_overlay)
}
\arguments{

  \item{object}{A \code{\link{HilbertCurve-class}} object.}
  \item{ir}{a \code{\link[IRanges]{IRanges}} object.}
  \item{x1}{if positions are not integers, they can be set by \code{x1} and \code{x2}.}
  \item{x2}{if positions are not integers, they can be set by \code{x1} and \code{x2}.}
  \item{col}{colors corresponding to intervals in \code{ir} (or combinations of \code{x1} and \code{x2}).}
  \item{mean_mode}{when a segment in the curve overlaps with intervals in \code{ir}, how to calculate  the mean values for this segment. See explanation in \code{\link{hc_points}}.}
  \item{grid_line}{whether add grid lines to show blocks of the Hilber curve.  It should be an integer number and there will be \code{2^(grid_line-1)-1} grid lines horizontal and vertical.}
  \item{grid_line_col}{color for the grid lines}
  \item{overlay}{a function which calculates the overlayed colors. Let's assume the r channel for the layers which are already added is \code{r0}, the r channel for the new layer is \code{r} and the alpha channel is \code{alpha}, the overlayed color is \code{r*alpha + r0*(1-alpha)}. This self-defined function should accept 7 arguments which are: vectors of r, g, b channels which correspond to the layers that are already added, and r, g, b, alpha channels which corresponds to the new layer. All the  values passed into are between 0 to 1. The returned value for this function should be a list which contains r, g, b channels which correspond to the overlayed colors. Note that these 7 arguments only correspond to the pixels which are covered by the new layer.}

}
\details{
If you want to add more than one layers to the curve, remember to set colors with transparency.

This function only works under 'pixel' mode.
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
ir = IRanges(s, e)
hc_layer(hc, ir)
hc = HilbertCurve(1, 100, level = 9, mode = "pixel")
hc_layer(hc, ir, grid_line = 3)
}
