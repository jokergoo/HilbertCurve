\name{hc_text-HilbertCurve-method}
\alias{hc_text,HilbertCurve-method}
\title{
Add text to Hilbert curve
}
\description{
Add text to Hilbert curve
}
\usage{
\S4method{hc_text}{HilbertCurve}(object, ir, labels, x1 = NULL, x2 = x1, gp = gpar(), ...)
}
\arguments{

  \item{object}{A \code{\link{HilbertCurve-class}} object.}
  \item{ir}{an \code{\link[IRanges:IRanges-constructor]{IRanges}} object that contains positions which correspond to text. The middle point of the interval will be the position of the text.}
  \item{labels}{text corresponding to intervals in \code{ir}.}
  \item{x1}{if start positions are not integers, they can be set by \code{x1}.}
  \item{x2}{if end positions are not integers, they can be set by \code{x2}.}
  \item{gp}{graphic parameters for text. It should be specified by \code{\link[grid]{gpar}}.}
  \item{...}{pass to \code{\link[grid]{grid.text}}. E.g. you can set text justification by \code{just} here.}

}
\details{
The text is added correspoding to the middle of each interval in \code{ir}.
}
\value{
A data frame which contains coordinates (in the 2D space) of text.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
hc = HilbertCurve(1, 100, level = 4, reference = TRUE)

x = sort(sample(100, 20))
s = x[1:10*2 - 1]
e = x[1:10*2]
require(IRanges)
ir = IRanges(s, e)

labels = sample(letters, length(ir), replace = TRUE)
hc_text(hc, ir, labels = labels)
}
