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
\S4method{hc_normal_points}{HilbertCurve}(object, ir, gp = gpar(), pch = 1, size = unit(1, "char"))}
\arguments{

  \item{object}{A \code{\link{HilbertCurve-class}} object}
  \item{ir}{a \code{\link[IRanges]{IRanges}} object}
  \item{size}{size of the points}
  \item{pch}{shape of points, used for points if \code{np >= 2}}
  \item{gp}{graphical parameters for points}
}
\details{
Points are added at every middle point in \code{ir}.

This function is used internally.

}
\value{
A data frame which contains coordinates for points.

}
\author{
Zuguang Gu <z.gu@dkfz.de>

}
