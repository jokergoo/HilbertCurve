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
  \item{ir}{a \code{\link[IRanges]{IRanges}} object.}
  \item{x1}{if start positions are not integers, they can be set by \code{x1}.}
  \item{x2}{if end positions are not integers, they can be set by \code{x2}.}
  \item{gp}{graphic parameters for lines. It should be specified by \code{\link[grid]{gpar}}.}
  \item{end_type}{since two ends of a continuous curve do not necessary completely overlap with  the hilbert curve segments, this argument determines how to include the ends of the curve. \code{average}: if the end covers more than half of the segment, the whole segment is included and if the end covers less than half of hte segment, the segment is removed; \code{expanding}: segments are included as long as they overlap; \code{shrinking}: segments are removed if they are not completely covered.}

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
