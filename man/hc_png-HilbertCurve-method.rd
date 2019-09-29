\name{hc_png-HilbertCurve-method}
\alias{hc_png,HilbertCurve-method}
\title{
Save Hilbert curve as a PNG figure
}
\description{
Save Hilbert curve as a PNG figure
}
\usage{
\S4method{hc_png}{HilbertCurve}(object, file = "HilbertCurve.png")
}
\arguments{

  \item{object}{A \code{\link{HilbertCurve-class}} object.}
  \item{file}{file name. If the suffix of the file name is not \code{.png},  it will be added automatically no matter you like it or not.}

}
\details{
A PNG figure with resolution of \code{2^level x 2^level} is generated.

Only the body of the Hilbert curve will be written to PNG file.

This function only works under 'pixel' mode.
}
\value{
No value is returned.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
hc = HilbertCurve(1, 100, level = 9, mode = "pixel")

x = sort(sample(100, 20))
s = x[1:10*2 - 1]
e = x[1:10*2]
require(IRanges)
ir = IRanges(s, e)

hc_layer(hc, ir)
hc_png(hc, file = "test.png")
}
\alias{hc_png}
