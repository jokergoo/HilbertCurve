\name{HilbertCurve-class}
\docType{class}
\alias{HilbertCurve-class}
\title{
The HilbertCurve class
}
\description{
The HilbertCurve class
}
\details{
Hilbert curve (\url{https://en.wikipedia.org/wiki/Hilbert_curve} ) is a type of space-filling curves
that folds one-dimensional axis into a two-dimensional space, but still keeps the locality.
It has advantages to visualize data with long axis with high resolution.

This package aims to provide an easy and flexible way to visualize data through Hilbert curve.
The implementation and example figures are based on following sources:

\itemize{
  \item \url{https://mkweb.bcgsc.ca/hilbert/}
  \item \url{https://corte.si/posts/code/hilbert/portrait/index.html}
  \item \url{https://bioconductor.org/packages/devel/bioc/html/HilbertVis.html}
}
}
\section{Methods}{
The \code{\link{HilbertCurve-class}} provides following methods:

\itemize{
  \item \code{\link{HilbertCurve}}: constructor method;
  \item \code{\link{hc_points,HilbertCurve-method}}: add points;
  \item \code{\link{hc_segments,HilbertCurve-method}}: add lines;
  \item \code{\link{hc_rect,HilbertCurve-method}}: add rectangles;
  \item \code{\link{hc_polygon,HilbertCurve-method}}: add polygons;
  \item \code{\link{hc_text,HilbertCurve-method}}: add text;
  \item \code{\link{hc_layer,HilbertCurve-method}}: add layers, works under "pixel" mode;
  \item \code{\link{hc_png,HilbertCurve-method}}: save plot as PNG format, works under "pixel" mode.
}}
\seealso{
The \code{\link{GenomicHilbertCurve-class}} inherits \code{\link{HilbertCurve-class}} and is designed specifically for handling genomic data.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
NULL
}
