\name{hc_layer-HilbertCurve-method}
\alias{hc_layer,HilbertCurve-method}
\alias{hc_layer}
\title{
Add a new layer on the Hilbert curve

}
\description{
Add a new layer on the Hilbert curve

}
\usage{
\S4method{hc_layer}{HilbertCurve}(object, ir, x1 = NULL, x2 = NULL, col = "red",
    mean_mode = c("w0", "absolute", "weighted"))}
\arguments{

  \item{object}{A \code{\link{HilbertCurve-class}} object}
  \item{ir}{a \code{\link[IRanges]{IRanges}} object}
  \item{x1}{if positions are not integer, the position can be set by \code{x1} and \code{x2}}
  \item{x2}{if positions are not integer, the position can be set by \code{x1} and \code{x2}}
  \item{col}{colors corresponding to each interval in \code{ir}}
  \item{mean_mode}{when a window is not perfectkt matched to one region in \code{ir}, how to calculate the mean values in this window. See explanation in \code{\link{hc_points}}.}
}
\details{
If you want to add more than one layers on the plot, remember to set transparent colors.

}
\value{
No value is returned

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
hc_save(hc, file = "test.png")

}
