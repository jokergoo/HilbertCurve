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
\S4method{hc_segmented_points}{HilbertCurve}(object, ir, x1 = NULL, x2 = NULL, gp = gpar(), np = max(c(2, 10 - hc_level(object))),
    mean_mode = c("w0", "absolute", "weighted"),
    shape = c("circle", "square", "triangle", "hexagon", "star"))}
\arguments{

  \item{object}{A \code{\link{HilbertCurve-class}} object}
  \item{ir}{a \code{\link[IRanges]{IRanges}} object}
  \item{x1}{if positions are not integer, the position can be set by \code{x1} and \code{x2}}
  \item{x2}{if positions are not integer, the position can be set by \code{x1} and \code{x2}}
  \item{np}{number of points (a circle or a square, ...) that are put in a segment}
  \item{gp}{graphical parameters for points}
  \item{mean_mode}{when a window is not perfectkt matched to one region in \code{ir}, how to calculate the mean values in this window. See explanation in \code{\link{hc_points}}.}
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
\examples{
# see hc_points
NULL

}
