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
\S4method{hc_normal_points}{HilbertCurve}(object, ir, x1 = NULL, x2 = NULL, gp = gpar(),
    pch = 1, size = unit(1, "char"))
}
\arguments{

  \item{object}{A \code{\link{HilbertCurve-class}} object.}
  \item{ir}{a \code{\link[IRanges]{IRanges}} object.}
  \item{x1}{if positions are not integers, they can be set by \code{x1} and \code{x2}.}
  \item{x2}{if positions are not integers, they can be set by \code{x1} and \code{x2}.}
  \item{size}{size of the points. It should be a \code{\link[grid]{unit}} object.}
  \item{pch}{shape of points.}
  \item{gp}{graphical parameters for points. It should be specified by \code{\link[grid]{gpar}}.}

}
\details{
Points are added at the middle of the intervals in \code{ir}.

This function is used internally.
}
\value{
A data frame which contains coordinates for points.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# see documentation of hc_points
NULL
}
