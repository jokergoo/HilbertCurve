\name{hc_rect-HilbertCurve-method}
\alias{hc_rect,HilbertCurve-method}
\alias{hc_rect}
\title{
Add rectangles on Hilbert curve

}
\description{
Add rectangles on Hilbert curve

}
\usage{
\S4method{hc_rect}{HilbertCurve}(object, ir, x1 = NULL, x2 = NULL, gp = gpar(fill = "red", col = "red"),
    mean_mode = c("w0", "absolute", "weighted"))}
\arguments{

  \item{object}{A \code{\link{HilbertCurve-class}} object}
  \item{ir}{a \code{\link[IRanges]{IRanges}} object}
  \item{x1}{if positions are not integer, the position can be set by \code{x1} and \code{x2}}
  \item{x2}{if positions are not integer, the position can be set by \code{x1} and \code{x2}}
  \item{gp}{graphical parameters for rectangles}
  \item{mean_mode}{when a window is not perfectkt matched to one region in \code{ir}, how to calculate the mean values in this window. See explanation in \code{\link{hc_points}}.}
}
\value{
A data frame which contains coordinates for rectangles.

}
\author{
Zuguang Gu <z.gu@dkfz.de>

}
\examples{
hc = HilbertCurve(1, 100, level = 4, reference = TRUE)

x = sort(sample(100, 20))
s = x[1:10*2 - 1]
e = x[1:10*2]
ir = IRanges(s, e)
hc_rect(hc, ir)

}
