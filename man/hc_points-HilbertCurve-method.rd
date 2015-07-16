\name{hc_points-HilbertCurve-method}
\alias{hc_points,HilbertCurve-method}
\alias{hc_points}
\title{
Add points to the Hilbert curve

}
\description{
Add points to the Hilbert curve

}
\usage{
\S4method{hc_points}{HilbertCurve}(object, ir, x1 = NULL, x2 = NULL, np = max(c(2, 10 - hc_level(object))),
    size = unit(1, "char"), pch = 1, gp = gpar(), mean_mode = c("w0", "absolute", "weighted"),
    shape = c("circle", "square", "triangle", "hexagon", "star"))}
\arguments{

  \item{object}{A \code{\link{HilbertCurve-class}} object.}
  \item{ir}{a \code{\link[IRanges]{IRanges}} object.}
  \item{x1}{if positions are not integers, they can be set by \code{x1} and \code{x2}.}
  \item{x2}{if positions are not integers, they can be set by \code{x1} and \code{x2}.}
  \item{np}{number of points (a circle or a square, ...) that are put in a segment. \code{np} controlsthe mode of how to add the points to the curve. See 'details' section.}
  \item{size}{size of the points. It should be a \code{\link[grid]{unit}} object. Only works if \code{np < 2}}
  \item{pch}{shape of points, used for points if \code{np < 2}.}
  \item{gp}{graphical parameters for points. It should be specified by \code{\link[grid]{gpar}}.}
  \item{mean_mode}{when a segment in the curve overlaps with intervals in \code{ir}, how to calculate the mean values for this segment (such as the RGB colors). See 'Details' section for a detailed explanation.}
  \item{shape}{shape of points, used for points if \code{np >= 2}.}
}
\details{
If \code{np} is set to a value less than 2 or \code{NULL}, points will be added at the middle points in \code{ir} (or \code{x1}, \code{x2}).
If \code{np} is set to a value larger or equal to 2, every segment that overlaps to \code{ir} will be segmented into \code{np} parts
and a circle (or star, ...) is put on every 'small segments'.

Following illustrates different settings for \code{mean_mode}:

  \preformatted{
       100    80     60    values in ir (e.g. red compoment for colors)
    ++++++   +++   +++++   ir
      ================     window (width = 16)
        4     3     3      overlap

    absolute: (100 + 80 + 60)/3
    weighted: (100*4 + 80*3 + 60*3)/(4 + 3 + 3)
    w0:       (100*4 + 80*3 + 60*3)/16
  }

So use of the mode depends on specific scenario. For example, if \code{ir} corresponds to positions of genes,
then the mode of \code{w0} is perhaps a good choise. If \code{ir} corresponds to positions of CpG sites which is
has width of 1 and most of the time is sparse in genomic windows, then \code{absolute} is a correct choice.

Graphical parameters can be set as a vector and they will be averaged according to above rules.

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
