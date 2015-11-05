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
  \item{x1}{if positions are not integers, they can be set by \code{x1} and \code{x2}.}
  \item{x2}{if positions are not integers, they can be set by \code{x1} and \code{x2}.}
  \item{gp}{graphical parameters for text. It should be specified by \code{\link[grid]{gpar}}.}
  \item{...}{pass to \code{\link[grid]{grid.text}}. E.g. you can set text justification by \code{just} here.}

}
\details{
Internally used.
}
\value{
\code{NULL}
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
NULL
}
