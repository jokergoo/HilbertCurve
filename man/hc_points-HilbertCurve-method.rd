\name{hc_points-HilbertCurve-method}
\alias{hc_points,HilbertCurve-method}
\title{
Add points to the Hilbert curve
}
\description{
Add points to the Hilbert curve
}
\usage{
\S4method{hc_points}{HilbertCurve}(object, ir, x1 = NULL, x2 = x1,
    np = max(c(2, 10 - hc_level(object))), size = unit(1, "char"),
    pch = 1, gp = gpar(), mean_mode = c("w0", "absolute", "weighted"),
    shape = "circle")
}
\arguments{

  \item{object}{A \code{\link{HilbertCurve-class}} object.}
  \item{ir}{a \code{\link[IRanges]{IRanges}} object.}
  \item{x1}{if start positions are not integers, they can be set by \code{x1}.}
  \item{x2}{if end positions are not integers, they can be set by \code{x2}.}
  \item{np}{number of points (a circle or a square, ...) that are put in a segment. \code{np} controls the mode of how to add the points to the curve. See 'Details' section.}
  \item{size}{size of the points. It should be a \code{\link[grid]{unit}} object. Only works if \code{np <= 1}}
  \item{pch}{shape of points, used for points if \code{np <= 1}.}
  \item{gp}{graphic parameters for points. It should be specified by \code{\link[grid]{gpar}}.}
  \item{mean_mode}{when \code{np >= 2}, each segment on the curve is split into \code{np} windows and each window actually represents an interval in the axis. When overlapping input intervals to the windows on the curve and when the window can not completely cover the input intervals, some averaging method should be applied to get a more accurate estimation for the value in the window. Here the HilbertCurve package provides three modes: "w0", "weighted" and "absolute" which calculate the mean value in the window with respect to different scenarios. See 'Details' section and the vignette for more informative explanation.}
  \item{shape}{shape of points, used for points if \code{np >= 2}. Possible values are "circle", "square", "triangle", "hexagon", "star".}

}
\details{
If \code{np} is set to 1 or \code{NULL}, points will be added in the middle for each interval in \code{ir} (or \code{x1}, \code{x2}).

If \code{np} is set to a value larger or equal to 2, every segment on the curve will be split by \code{np} points (e.g. circles).
In this case, each point actually represent a window on the curve and when the window is not fully covered by
the input intervals, there are three different metrics to average the values in the window.

Following illustrates different settings for \code{mean_mode}:

  \preformatted{
       100    80     60    values in ir
    ++++++   +++   +++++   ir
      ================     window (width = 16)
        4     3     3      overlap

    absolute: (100 + 80 + 60)/3
    weighted: (100*4 + 80*3 + 60*3)/(4 + 3 + 3)
    w0:       (100*4 + 80*3 + 60*3 + 0*6)/16  }

So which mode to use depends on specific scenario. If the background is not of interest, \code{\link{absolute}} and \code{\link{weighted}}
modes may be proper and if the value also needs to be averaged with background, \code{\link{w0}} is the proper choice.

Graphic parameters is always represented as numeric values (e.g. colors can be converted into numeric RGB values) 
and they will be averaged according to above rules.

Internally, it will depatch to \code{\link{hc_normal_points,HilbertCurve-method}} or \code{\link{hc_segmented_points,HilbertCurve-method}}
depending on the value of \code{np}.
}
\value{
A data frame which contains coordinates for points.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
hc = HilbertCurve(1, 100, level = 4, reference = TRUE)

x = sort(sample(100, 20))
s = x[1:10*2 - 1]
e = x[1:10*2]
ir = IRanges(s, e)

hc_points(hc, ir)

hc = HilbertCurve(1, 100, level = 4, reference = TRUE)
hc_points(hc, x1 = c(1.5, 50.5), x2 = c(10.5, 60.5))

require(circlize)
value = runif(length(ir))
col_fun = colorRamp2(range(value), c("white", "red"))
hc = HilbertCurve(1, 100, level = 4, reference = TRUE)
hc_points(hc, ir, np = 3, shape = "star", gp = gpar(fill = col_fun(value)))

hc = HilbertCurve(1, 100, level = 4, reference = TRUE)
hc_points(hc, ir, np = 0)

hc = HilbertCurve(1, 100, level = 4, reference = TRUE)
hc_points(hc, np = 0, x1 = c(1.5, 50.5), x2 = c(10.5, 60.5))
hc_points(hc, np = 0, x1 = 70.5, gp = gpar(col = "red"))
}
