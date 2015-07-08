\name{HilbertCurve}
\alias{HilbertCurve}
\title{
Initialize a Hilbert curve

}
\description{
Initialize a Hilbert curve

}
\usage{
HilbertCurve(s, e, level = 4, mode = c("normal", "pixel"),
    reference = FALSE, zoom = NULL, newpage = TRUE, background = "white")}
\arguments{

  \item{s}{start of Hilbert curve, should be an integer}
  \item{e}{end of Hilbert curve, should be an integer}
  \item{level}{order of Hilbert curve. There will by \code{4^level} segments in the Hilbert curve.}
  \item{mode}{make it like a normal R plot or write the plot directly into png file.}
  \item{reference}{add reference information on the plot}
  \item{zoom}{zooming of data ranges}
  \item{newpage}{whether call \code{\link[grid]{newpage}}`}
  \item{background}{background color}
}
