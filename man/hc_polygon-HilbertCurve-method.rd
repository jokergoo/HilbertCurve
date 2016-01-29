\name{hc_polygon-HilbertCurve-method}
\alias{hc_polygon,HilbertCurve-method}
\title{
Add areas to Hilbert curve
}
\description{
Add areas to Hilbert curve
}
\usage{
\S4method{hc_polygon}{HilbertCurve}(object, ir, x1 = NULL, x2 = NULL,
    gp = gpar(),
    end_type = c("average", "expanding", "shrinking"))
}
\arguments{

  \item{object}{A \code{\link{HilbertCurve-class}} object.}
  \item{ir}{a \code{\link[IRanges]{IRanges}} object which specifies the input intervals.}
  \item{x1}{if start positions are not integers, they can be set by \code{x1}.}
  \item{x2}{if end positions are not integers, they can be set by \code{x2}.}
  \item{gp}{graphic parameters. It should be specified by \code{\link[grid]{gpar}}.}
  \item{end_type}{since two ends of a continuous interval do not necessarily completely overlap with  the hilbert curve segments, this argument controls how to determine the ends of the interval which will be presented on the curve.  \code{average}: if the end covers more than half of the segment, the whole segment is included and if the end covers less than half of hte segment, the segment is removed; \code{expanding}: segments are included as long as they are overlapped; \code{shrinking}: segments are removed if they are not completely covered.}

}
\details{
Drawing polygons are quite similar as drawing rectangles in the context of Hilbert curve.
The major difference is: 1) Colors for the ends of the interval can change if they are not completely
covered by the segments, and 2) polygons can have borders.

Basically polygons are used to mark areas in the Hilbert curve.
}
\value{
No value is returned.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
ir = IRanges(10, 40)

hc = HilbertCurve(1, 100, level = 4, reference = TRUE)
hc_segments(hc, ir)
hc_polygon(hc, ir, gp = gpar(fill = "#FF000080", col = 1))
}
