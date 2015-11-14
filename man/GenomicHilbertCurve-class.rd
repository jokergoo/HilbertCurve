\name{GenomicHilbertCurve-class}
\docType{class}
\alias{GenomicHilbertCurve-class}
\title{
The GenomicHilbertCurve class
}
\description{
The GenomicHilbertCurve class
}
\details{
The \code{\link{GenomicHilbertCurve-class}} is inherited from the \code{\link{HilbertCurve-class}}. Basically
the structure of this class is almost the same as the \code{\link{HilbertCurve-class}} but with several
additional slots added to facilitate visualizing genomic data.
}
\section{Methods}{
The \code{\link{GenomicHilbertCurve-class}} provides following methods:

\itemize{
  \item \code{\link{GenomicHilbertCurve}}: constructor method;
  \item \code{\link{hc_points,GenomicHilbertCurve-method}}: add points;
  \item \code{\link{hc_segments,GenomicHilbertCurve-method}}: add lines;
  \item \code{\link{hc_rect,GenomicHilbertCurve-method}}: add rectangles;
  \item \code{\link{hc_text,GenomicHilbertCurve-method}}: add text;
  \item \code{\link{hc_layer,GenomicHilbertCurve-method}}: add layers undel "pixel" mode;
  \item \code{\link{hc_map,GenomicHilbertCurve-method}}: show the map of different categories on the curve.
}}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
NULL
}
