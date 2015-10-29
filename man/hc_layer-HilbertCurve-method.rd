\name{hc_layer-HilbertCurve-method}
\alias{hc_layer,HilbertCurve-method}
\alias{hc_layer}
\title{
Add a new layer to the Hilbert curve
}
\description{
Add a new layer to the Hilbert curve
}
\usage{
\S4method{hc_layer}{HilbertCurve}(object, ir, x1 = NULL, x2 = NULL, col = "red",
    mean_mode = c("w0", "absolute", "weighted"), grid_line = 0)
}
\arguments{

  \item{object}{A \code{\link{HilbertCurve-class}} object.}
  \item{ir}{a \code{\link[IRanges]{IRanges}} object.}
  \item{x1}{if positions are not integers, they can be set by \code{x1} and \code{x2}.}
  \item{x2}{if positions are not integers, they can be set by \code{x1} and \code{x2}.}
  \item{col}{colors corresponding to intervals in \code{ir} (or combinations of \code{x1} and \code{x2}).}
  \item{mean_mode}{when a segment in the curve overlaps with intervals in \code{ir}, how to calculate the mean values for this segment. See explanation in \code{\link{hc_points}}.}
  \item{grid_line}{whether add grid lines to show blocks of the Hilber curve. It should be an integer number and there will be \code{2^(grid_line-1)-1}grid lines horizontal and vertical.}

}
\details{
If you want to add more than one layers to the curve, remember to set colors transparent.

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
