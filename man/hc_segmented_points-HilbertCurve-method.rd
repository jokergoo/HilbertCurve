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
\S4method{hc_segmented_points}{HilbertCurve}(object, ir = NULL, x1 = NULL, x2 = NULL, gp = gpar(),
    np = max(c(2, 10 - hc_level(object))),
    mean_mode = c("w0", "absolute", "weighted", "max_freq"),
    shape = "circle")
}
\arguments{

  \item{object}{A \code{\link{HilbertCurve-class}} object.}
  \item{ir}{an \code{\link[IRanges:IRanges-constructor]{IRanges}} object which specifies the input intervals.}
  \item{x1}{if start positions are not integers, they can be set by \code{x1}.}
  \item{x2}{if end positions are not integers, they can be set by \code{x2}.}
  \item{np}{number of points (a circle or a square, ...) that are put in a segment.}
  \item{gp}{graphic parameters for points. It should be specified by \code{\link[grid]{gpar}}. The size of the points can be set here because the size of points are determined by \code{np} argument.}
  \item{mean_mode}{when a segment in the curve overlaps with intervals in \code{ir}, how to calculate  the mean values for this segment. See explanation in \code{\link{hc_points}}.}
  \item{shape}{shape of points. Possible values are "circle", "square", "triangle", "hexagon", "star".}

}
\details{
Every segment on the curve is split by \code{np} points.

This function is used internally, please use \code{\link{hc_points,HilbertCurve-method}} directly.
}
\value{
A data frame which contains coordinates (in the 2D space) of points.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# see documentation of hc_points
NULL
}
