\name{hc_rect-HilbertCurve-method}
\alias{hc_rect,HilbertCurve-method}
\title{
Add rectangles on Hilbert curve
}
\description{
Add rectangles on Hilbert curve
}
\usage{
\S4method{hc_rect}{HilbertCurve}(object, ir = NULL, x1 = NULL, x2 = NULL,
    gp = gpar(fill = "red"),
    mean_mode = c("w0", "absolute", "weighted"))
}
\arguments{

  \item{object}{A \code{\link{HilbertCurve-class}} object.}
  \item{ir}{an \code{\link[IRanges:IRanges-constructor]{IRanges}} object which specifies the input intervals.}
  \item{x1}{if start positions are not integers, they can be set by \code{x1}.}
  \item{x2}{if end positions are not integers, they can be set by \code{x2}.}
  \item{gp}{graphic parameters for rectangles. It should be specified by \code{\link[grid]{gpar}}. Note you cannot set \code{linejoin} and \code{lineend}.}
  \item{mean_mode}{when a segment in the curve can not be overlapped with intervals in \code{ir}, how to calculate  the mean values for this segment. See explanation in \code{\link{hc_points,HilbertCurve-method}}.}

}
\details{
Rectangles are put if a segment in the Hilbert curve overlaps with the input intervals. 
You cannot set the width or height of the rectangles. It is always fixed (actually it is a square).

It can be thought as the low-resolution version of \code{\link{hc_layer,HilbertCurve-method}}.
}
\value{
A data frame which contains coordinates (in the 2D space) of rectangles.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
hc = HilbertCurve(1, 100, level = 4, reference = TRUE)

x = sort(sample(100, 20))
s = x[1:10*2 - 1]
e = x[1:10*2]
require(IRanges)
ir = IRanges(s, e)
hc_rect(hc, ir)
}
