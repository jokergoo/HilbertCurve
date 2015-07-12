\name{hc_points-HilbertCurve-method}
\alias{hc_points,HilbertCurve-method}
\alias{hc_points}
\title{
Add points to the Hilbert curve

}
\description{
Add points to the Hilbert curve

}
\usage{
\S4method{hc_points}{HilbertCurve}(object, ir, x1 = NULL, x2 = NULL, np = max(c(2, 10 - hc_level(object))),
    size = unit(1, "char"), pch = 1, gp = gpar(), mean_mode = c("w0", "absolute", "weighted"),
    shape = c("circle", "square", "triangle", "hexagon", "star"))}
\arguments{

  \item{object}{A \code{\link{HilbertCurve-class}} object}
  \item{ir}{a \code{\link[IRanges]{IRanges}} object}
  \item{x1}{if positions are not integer, the position can be set by \code{x1} and \code{x2}}
  \item{x2}{if positions are not integer, the position can be set by \code{x1} and \code{x2}}
  \item{np}{number of points (a circle or a square, ...) that are put in a segment}
  \item{size}{size of the points}
  \item{pch}{shape of points, used for points if \code{np >= 2}}
  \item{gp}{graphical parameters for points}
  \item{mean_mode}{when a window is not perfectkt matched to one region in \code{ir}, how to calculate the mean values in this window. See 'Details' section for a detailed explanation.}
  \item{shape}{shape of points, used for points if \code{np <= 1}}
}
\details{
If \code{np} is set to a value less than 2 or \code{NULL}, points will be added at every middle point in \code{ir}.
If \code{np} is set to a value larger or equal to 2, a list of e.g. circles are put at every segment in \code{ir},
so, longer segments will have more circles on it.

Following illustrates different settings for \code{mean_mode}:

  \preformatted{
       4      5      2     values in ir
    ++++++   +++   +++++   gr
      ================     window (16bp)

    absolute: (4 + 5 + 2)/3
    weighted: (4*4 + 5*3 + 2*3)/(4 + 3 + 3)
    w0:       (4*4 + 5*3 + 2*3)/16
  }

}
\value{
A data frame which contains coordinates for points.

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

hc_points(hc, ir)

hc = HilbertCurve(1, 100, level = 4, reference = TRUE)
hc_points(hc, x1 = c(1.5, 50.5), x2 = c(10.5, 60.5))

require(circlize)
value = runif(length(ir))
col_fun = colorRamp2(range(value), c("white", "red"))
hc = HilbertCurve(1, 100, level = 4, reference = TRUE)
hc_points(hc, ir, np = 3, shape = "star", gp = gpar(fill = col_fun(value)))

hc = HilbertCurve(1, 100, level = 4, reference = TRUE)
hc_points(hc, ir, np = 0)

hc = HilbertCurve(1, 100, level = 4, reference = TRUE)
hc_points(hc, np = 0, x1 = c(1.5, 50.5), x2 = c(10.5, 60.5))
hc_points(hc, np = 0, x1 = 70.5, gp = gpar(col = "red"))

}
