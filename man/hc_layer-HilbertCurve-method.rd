\name{hc_layer-HilbertCurve-method}
\alias{hc_layer,HilbertCurve-method}
\alias{hc_layer}
\title{
Add a new layer on the Hilbert curve

}
\description{
Add a new layer on the Hilbert curve

}
\usage{
\S4method{hc_layer}{HilbertCurve}(object, ir, col = "red", mean_mode = c("w0", "absolute", "weighted"))}
\arguments{

  \item{object}{A \code{\link{HilbertCurve-class}} object}
  \item{ir}{a \code{\link[IRanges]{IRanges}} object}
  \item{col}{colors corresponding to each interval in \code{ir}}
  \item{mean_mode}{-mean_mode}
}
\details{
If you want to add more than one layers on the plot, remember to set transparent colors.

}
\value{
No value is returned

}
\author{
Zuguang Gu <z.gu@dkfz.de>

}
