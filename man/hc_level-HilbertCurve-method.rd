\name{hc_level-HilbertCurve-method}
\alias{hc_level,HilbertCurve-method}
\alias{hc_level}
\title{
Level of the Hilbert curve

}
\description{
Level of the Hilbert curve

}
\usage{
\S4method{hc_level}{HilbertCurve}(object)}
\arguments{

  \item{object}{A \code{\link{HilbertCurve-class}} object}
}
\value{
The level of the Hilbert curve

}
\author{
Zuguang Gu <z.gu@dkfz.de>

}
\examples{
hc = HilbertCurve(1, 100)
hc_level(hc)

hc = HilbertCurve(1, 100, level = 5)
hc_level(hc)

}
