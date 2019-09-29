\name{zoom-HilbertCurve-method}
\alias{zoom,HilbertCurve-method}
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
  \item{x}{original positions.}

}
\details{
Internally, position are stored as integer values. To better map the data to the Hilbert curve, 
the original positions are zoomed
according to the range and the level of Hilbert curve. E.g. if 
the curve visualizes data ranging from 1 to 2 but level of the curve is set to 4,
the positions will be zoomed by ~x2000 so that values like 1.5, 1.555 can be mapped
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
\alias{zoom}
