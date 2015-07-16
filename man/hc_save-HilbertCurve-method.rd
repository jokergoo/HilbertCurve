\name{hc_save-HilbertCurve-method}
\alias{hc_save,HilbertCurve-method}
\alias{hc_save}
\title{
Save Hilbert curve as a PNG figure

}
\description{
Save Hilbert curve as a PNG figure

}
\usage{
\S4method{hc_save}{HilbertCurve}(object, file = "Rplot.png", grid = 0)}
\arguments{

  \item{object}{A \code{\link{HilbertCurve-class}} object.}
  \item{file}{file name. If the suffix of the file name is not \code{.png}, it will be added automatically no matter you like it or not.}
  \item{grid}{whether add grid lines to show blocks of the Hilber curve. It should be an integer number and should not exceed the level of the curve.}
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
ir = IRanges(s, e)

hc_layer(hc, ir)
hc_save(hc, file = "test2.png", grid = 2)
hc_save(hc, file = "test3.png", grid = 3)
hc_save(hc, file = "test4.png", grid = 4)

}
