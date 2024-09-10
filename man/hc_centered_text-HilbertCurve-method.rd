\name{hc_centered_text-HilbertCurve-method}
\alias{hc_centered_text,HilbertCurve-method}
\title{
Add text to the center of the block
}
\description{
Add text to the center of the block
}
\usage{
\S4method{hc_centered_text}{HilbertCurve}(object, ir = NULL, labels, x1 = NULL, x2 = NULL, gp = gpar(), ...)
}
\arguments{

  \item{object}{A \code{\link{HilbertCurve-class}} object.}
  \item{ir}{an \code{\link[IRanges:IRanges-constructor]{IRanges}} object that contains positions which correspond to text. The middle points of the intervals will be the positions of the text.}
  \item{labels}{text corresponding to intervals in \code{ir}.}
  \item{x1}{if start positions are not integers, they can be set by \code{x1}.}
  \item{x2}{if end positions are not integers, they can be set by \code{x2}.}
  \item{gp}{graphic parameters for text. It should be specified by \code{\link[grid]{gpar}}.}
  \item{...}{pass to \code{\link[grid]{grid.text}}. E.g. you can set text justification by \code{just} here.}

}
\details{
If the interval is long enough that it represents as a block in the 2D space, the corresponding
label is put approximately at center (or at least inside) of the block.

Please use \code{\link{hc_text,HilbertCurve-method}} directly.
}
\value{
\code{NULL}
}
\seealso{
It is used in \code{\link{hc_map,GenomicHilbertCurve-method}} to put chromosome names in the center
of chromosomes.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
hc = HilbertCurve(1, 10)
hc_rect(hc, x1 = c(1, 3, 7), x2 = c(3, 7, 10), gp = gpar(fill = 2:5))
hc_centered_text(hc, x1 = 1, x2 = 3, labels = "A")
hc_centered_text(hc, x1 = 3, x2 = 7, labels = "B")
hc_centered_text(hc, x1 = 7, x2 = 10, labels = "C")
}
\alias{hc_centered_text}
