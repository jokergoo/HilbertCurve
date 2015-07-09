\name{hc_rect-HilbertCurve-method}
\alias{hc_rect,HilbertCurve-method}
\alias{hc_rect}
\title{
Add rectangles on Hilbert curve

}
\description{
Add rectangles on Hilbert curve

}
\usage{
\S4method{hc_rect}{HilbertCurve}(object, ir, gp = gpar(fill = "red", col = "red"),
    mean_mode = c("w0", "absolute", "weighted"))}
\arguments{

  \item{object}{A \code{\link{HilbertCurve-class}} object}
  \item{ir}{a \code{\link[IRanges]{IRanges}} object}
  \item{gp}{graphical parameters for rectangles}
  \item{mean_mode}{-mean_mode}
}
\value{
A data frame which contains coordinates for rectangles.

}
\author{
Zuguang Gu <z.gu@dkfz.de>

}
