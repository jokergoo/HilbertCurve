\name{hc_segments-HilbertCurve-method}
\alias{hc_segments,HilbertCurve-method}
\alias{hc_segments}
\title{
Add line segments to Hilbert curve

}
\description{
Add line segments to Hilbert curve

}
\usage{
\S4method{hc_segments}{HilbertCurve}(object, ir, gp = gpar())}
\arguments{

  \item{object}{A \code{\link{HilbertCurve-class}} object}
  \item{ir}{a \code{\link[IRanges]{IRanges}} object}
  \item{gp}{graphical parameters for rectangles}
}
\value{
A data frame which contains coordinates for segments.

}
\author{
Zuguang Gu <z.gu@dkfz.de>

}
