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
\S4method{hc_points}{HilbertCurve}(object, ir, np = max(c(2, 10 - hc_level(object))),
    size = unit(1, "char"), pch = 1, gp = gpar(), mean_mode = c("w0", "absolute", "weighted"),
    shape = c("circle", "square", "triangle", "hexagon", "star"))}
\arguments{

  \item{object}{A \code{\link{HilbertCurve-class}} object}
  \item{ir}{a \code{\link[IRanges]{IRanges}} object}
  \item{np}{number of points (a circle or a square, ...) that are put in a segment}
  \item{size}{size of the points}
  \item{pch}{shape of points, used for points if \code{np >= 2}}
  \item{gp}{graphical parameters for points}
  \item{mean_mode}{-mean_mode}
  \item{shape}{shape of points, used for points if \code{np <= 1}}
}
\details{
If \code{np} is set to a value less than 2 or \code{NULL}, points will be added at every middle point in \code{ir}.
If \code{np} is set to a value larger or equal to 2, a list of e.g. circles are put at every segment in \code{ir},
so, longer segments will have more circles on it.

}
\value{
A data frame which contains coordinates for points.

}
\author{
Zuguang Gu <z.gu@dkfz.de>

}
