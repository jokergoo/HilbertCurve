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
    reference = FALSE, arrow = TRUE, zoom = NULL, newpage = TRUE,
    background = "white", title = NULL, title_gp = gpar(fontsize = 16),
    legend = list())}
\arguments{

  \item{s}{start of the Hilbert curve, should be an integer}
  \item{e}{end of the Hilbert curve, should be an integer}
  \item{level}{level of the Hilbert curve. There will by \code{4^level} segments in the Hilbert curve.}
  \item{mode}{make it like a normal R plot or write the plot directly into png file.}
  \item{reference}{add reference line on the plot}
  \item{arrow}{whther draw arrows on the reference line}
  \item{zoom}{zooming of the position ranges}
  \item{newpage}{whether call \code{\link[grid]{grid.newpage}}`}
  \item{background}{background color}
  \item{title}{title of the plot}
  \item{title_gp}{graphical parameters for title}
  \item{legend}{a \code{\link[grid]{grob}} object or a list of \code{\link[grid]{grob}} objects.}
}
\value{
A \code{\link{HilbertCurve-class}} object.

}
\author{
Zuguang Gu <z.gu@dkfz.de>

}
