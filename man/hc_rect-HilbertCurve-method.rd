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

  \item{object}{A \code{\link{HilbertCurve-class}} object.}
  \item{ir}{a \code{\link[IRanges]{IRanges}} object.}
  \item{x1}{if positions are not integers, they can be set by \code{x1} and \code{x2}.}
  \item{x2}{if positions are not integers, they can be set by \code{x1} and \code{x2}.}
  \item{gp}{graphical parameters for rectangles. It should be specified by \code{\link[grid]{gpar}}.}
  \item{mean_mode}{when a segment in the curve overlaps with intervals in \code{ir}, how to calculate the mean values for this segment. See explanation in \code{\link{hc_points}}.}
}
\details{
You cannot set the width or height of the rectangles. Rectangles are always located
at the turning points and have width or height equal to the length of the segments.

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
