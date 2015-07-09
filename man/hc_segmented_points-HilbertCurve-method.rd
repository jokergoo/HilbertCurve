\name{hc_segmented_points-HilbertCurve-method}
\alias{hc_segmented_points,HilbertCurve-method}
\alias{hc_segmented_points}
\title{
Add points to the Hilbert curve

}
\description{
Add points to the Hilbert curve

}
\usage{
\S4method{hc_segmented_points}{HilbertCurve}(object, ir, gp = gpar(), np = max(c(2, 10 - hc_level(object))),
    mean_mode = c("w0", "absolute", "weighted"), shape = c("circle", "square", "triangle", "hexagon", "star"))}
\arguments{

  \item{object}{A \code{\link{HilbertCurve-class}} object}
  \item{ir}{a \code{\link[IRanges]{IRanges}} object}
  \item{np}{number of points (a circle or a square, ...) that are put in a segment}
  \item{gp}{graphical parameters for points}
  \item{mean_mode}{-mean_mode}
  \item{shape}{shape of points, used for points if \code{np <= 1}}
}
\details{
A list of e.g. circles are put at every segment in \code{ir},
so, longer segments will have more circles on it.

This function is used internally.

}
\value{
A data frame which contains coordinates for points.

}
\author{
Zuguang Gu <z.gu@dkfz.de>

}
