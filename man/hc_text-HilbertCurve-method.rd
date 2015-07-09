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
\S4method{hc_text}{HilbertCurve}(object, ir, text, gp = gpar(), ...)}
\arguments{

  \item{object}{A \code{\link{HilbertCurve-class}} object}
  \item{ir}{a \code{\link[IRanges]{IRanges}} object that contains positions of text. If interval has width larger than 1,the middle point of the interval will be the position of the text.}
  \item{text}{text corresponding the \code{ir}.}
  \item{gp}{graphical parameters for text}
  \item{...}{pass to \code{\link[grid]{grid.text}}. You can set \code{just} for text here}
}
\value{
A data frame which contains coordinates for text.

}
\author{
Zuguang Gu <z.gu@dkfz.de>

}
