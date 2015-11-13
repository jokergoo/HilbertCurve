\name{hc_centered_text-HilbertCurve-method}
\alias{hc_centered_text,HilbertCurve-method}
\alias{hc_centered_text}
\title{
Add text to the center of the block
}
\description{
Add text to the center of the block
}
\usage{
\S4method{hc_centered_text}{HilbertCurve}(object, ir, labels, x1 = NULL, x2 = NULL, gp = gpar(), ...)
}
\arguments{

  \item{object}{A \code{\link{HilbertCurve-class}} object.}
  \item{ir}{a \code{\link[IRanges]{IRanges}} object that contains positions of text. Basically, the middle point of the interval will be the position of the text.}
  \item{labels}{text corresponding to intervals in \code{ir}.}
  \item{x1}{if start positions are not integers, they can be set by \code{x1}.}
  \item{x2}{if end positions are not integers, they can be set by \code{x2}.}
  \item{gp}{graphical parameters for text. It should be specified by \code{\link[grid]{gpar}}.}
  \item{...}{pass to \code{\link[grid]{grid.text}}. E.g. you can set text justification by \code{just} here.}

}
\details{
If the interval is long enough that it represents as a block in the curve, the corresponding
label is put approximately at center of the block.

It is quite experimental and only used internally.
}
\value{
\code{NULL}
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
