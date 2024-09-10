\name{unzoom-HilbertCurve-method}
\alias{unzoom,HilbertCurve-method}
\title{
Transform zoomed positions to their original values
}
\description{
Transform zoomed positions to their original values
}
\usage{
\S4method{unzoom}{HilbertCurve}(object, x)
}
\arguments{

  \item{object}{A \code{\link{HilbertCurve-class}} object.}
  \item{x}{positions.}

}
\details{
This is a reverse function of \code{\link{zoom,HilbertCurve-method}}.

The function is used internally.
}
\value{
A numeric vector of original positions.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
hc = HilbertCurve(1, 2)
z = zoom(hc, 1.5)
unzoom(hc, z)
}
\alias{unzoom}
