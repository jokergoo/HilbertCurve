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
that fold one dimensional axis into a two dimensional space, but with still keeping the locality.
It has advantages to visualize data with long axis in following two aspects:
1. greatly improve resolution for the visualization;
2. easy to visualization clusters because generally they will also be close in the Hilbert curve. 

This package aims to provide a easy and flexible way to visualize data through Hilbert curve.
The implementation and example figures are based on following sources:

\itemize{
  \item \url{http://mkweb.bcgsc.ca/hilbert/}
  \item \url{http://corte.si/posts/code/hilbert/portrait/index.html}
  \item \url{http://bioconductor.org/packages/devel/bioc/html/HilbertVis.html}
}

}
\section{Methods}{
The \code{\link{HilbertCurve-class}} provides following methods:

\itemize{
  \item \code{\link{HilbertCurve}}: constructor method
  \item \code{\link{hc_points,HilbertCurve-method}}: add points
  \item \code{\link{hc_segments,HilbertCurve-method}}: add lines
  \item \code{\link{hc_rect,HilbertCurve-method}}: add rectangles
  \item \code{\link{hc_layer,HilbertCurve-method}}: add layers
  \item \code{\link{hc_save,HilbertCurve-method}}: save plot as png format
}

}
\author{
Zuguang Gu <z.gu@dkfz.de>

}
\examples{
NULL

}
