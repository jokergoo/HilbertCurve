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
\S4method{hc_segmented_points}{HilbertCurve}(object, ir, x1 = NULL, x2 = NULL, gp = gpar(),
    np = max(c(2, 10 - hc_level(object))),
    mean_mode = c("w0", "absolute", "weighted"),
    shape = "circle")
}
\arguments{

  \item{object}{A \code{\link{HilbertCurve-class}} object.}
  \item{ir}{a \code{\link[IRanges]{IRanges}} object.}
  \item{x1}{if positions are not integers, they can be set by \code{x1} and \code{x2}.}
  \item{x2}{if positions are not integers, they can be set by \code{x1} and \code{x2}.}
  \item{np}{number of points (a circle or a square, ...) that are put in a segment.}
  \item{gp}{graphical parameters for points. It should be specified by \code{\link[grid]{gpar}}.}
  \item{mean_mode}{when a segment in the curve overlaps with intervals in \code{ir}, how to calculate  the mean values for this segment. See explanation in \code{\link{hc_points}}.}
  \item{shape}{shape of points. Possible values are "circle", "square", "triangle", "hexagon", "star".}

}
\details{
Every segment that overlaps to \code{ir} will be segmented into \code{np} parts
 and a circle (or star, ...) is put on every 'small segments'.

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
