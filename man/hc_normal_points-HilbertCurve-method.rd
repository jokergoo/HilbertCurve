\name{hc_normal_points-HilbertCurve-method}
\alias{hc_normal_points,HilbertCurve-method}
\alias{hc_normal_points}
\title{
Add points to the Hilbert curve
}
\description{
Add points to the Hilbert curve
}
\usage{
\S4method{hc_normal_points}{HilbertCurve}(object, ir = NULL, x1 = NULL, x2 = x1, gp = gpar(),
    pch = 1, size = unit(1, "char"))
}
\arguments{

  \item{object}{A \code{\link{HilbertCurve-class}} object.}
  \item{ir}{an \code{\link[IRanges:IRanges-constructor]{IRanges}} object which specifies the input intervals.}
  \item{x1}{if start positions are not integers, they can be set by \code{x1}.}
  \item{x2}{if end positions are not integers, they can be set by \code{x2}.}
  \item{size}{size of the points. It should be a \code{\link[grid]{unit}} object, pass to \code{\link[grid]{grid.points}}.}
  \item{pch}{shape of points, pass to \code{\link[grid]{grid.points}}.}
  \item{gp}{graphic parameters for points. It should be specified by \code{\link[grid]{gpar}}.}

}
\details{
Points are added at the middle of the intervals in \code{ir} (or \code{x1} and \code{x2}),
so there is only one point for each interval.

This function is used internally. Please use \code{\link{hc_points,HilbertCurve-method}} instead.
}
\value{
A data frame which contains coordinates (in the 2D space) of points.
}
\seealso{
\code{\link{hc_points,HilbertCurve-method}}
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# see documentation of hc_points
NULL
}
