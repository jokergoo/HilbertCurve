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
\S4method{zoom}{HilbertCurve}(object, x)
}
\arguments{

  \item{object}{A \code{\link{HilbertCurve-class}} object.}
  \item{x}{positions.}

}
\details{
Internally, position are stored as integer values. To increase the resolution
 of the data that maps to the Hilbert curve, the original position would be zoom
 according to the range of the position and the level of Hilbert curve. E.g. if 
 the curve visualizes data ranging from 1 to 2 but level of the curve is set to 4,
 the positions will be zoomed by ~x2000 so that values link 1.5, 1.555 can be mapped
 to the curve with more accuracy.

The function is used internally.
}
\value{
A numeric vector which is zoomed positions.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
hc = HilbertCurve(1, 2)
 zoom(hc, 1.5)
}
