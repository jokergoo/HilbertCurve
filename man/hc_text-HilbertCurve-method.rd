\name{hc_text-HilbertCurve-method}
\alias{hc_text,HilbertCurve-method}
\alias{hc_text}
\title{
Add text to Hilbert curve

}
\description{
Add text to Hilbert curve

}
\usage{
\S4method{hc_text}{HilbertCurve}(object, ir, labels, x1 = NULL, x2 = NULL, gp = gpar(), ...)}
\arguments{

  \item{object}{A \code{\link{HilbertCurve-class}} object}
  \item{ir}{a \code{\link[IRanges]{IRanges}} object that contains positions of text. If interval has width larger than 1,the middle point of the interval will be the position of the text.}
  \item{labels}{text corresponding the \code{ir}.}
  \item{x1}{if positions are not integer, the position can be set by \code{x1} and \code{x2}}
  \item{x2}{if positions are not integer, the position can be set by \code{x1} and \code{x2}}
  \item{gp}{graphical parameters for text}
  \item{...}{pass to \code{\link[grid]{grid.text}}. You can set \code{just} for text here}
}
\value{
A data frame which contains coordinates for text.

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

labels = sample(letters, length(ir), replace = TRUE)
hc_text(hc, ir, labels = labels)

}
