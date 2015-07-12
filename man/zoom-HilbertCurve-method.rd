\name{zoom-HilbertCurve-method}
\alias{zoom,HilbertCurve-method}
\alias{zoom}
\title{
Zoom original positions

}
\description{
Zoom original positions

}
\usage{
\S4method{zoom}{HilbertCurve}(object, x)}
\arguments{

  \item{object}{A \code{\link{HilbertCurve-class}} object}
  \item{x}{positions}
}
\details{
Internally, position are stored as integer values. To increase the resolution
of the data that maps to the Hilbert curve, the original position would be zoom
according to the range of the position and the level of Hilbert curve. E.g. if 
the curve visualizes data ranging from 1 to 100 but level of the curve is set to 8,
the positions will be zoomed by ~x1000 so that values link 1.1, 1.111 can be mapped
to the curve with more accuracy.

The function is used internally.

}
\value{
Zoomed positions

}
\author{
Zuguang Gu <z.gu@dkfz.de>

}
\examples{
NULL

}
