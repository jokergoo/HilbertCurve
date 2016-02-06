\name{default_overlay}
\alias{default_overlay}
\title{
Default color overlay for adding new layers
}
\description{
Default color overlay for adding new layers
}
\usage{
default_overlay(r0, g0, b0, r, g, b, alpha = 1)
}
\arguments{

  \item{r0}{red channel for the layers that are already in the plot.}
  \item{g0}{green channel for the layers that are already in the plot.}
  \item{b0}{blue channel for the layers that are already in the plot.}
  \item{r}{red channel for the new layer}
  \item{g}{green channel for the new layer}
  \item{b}{blue channel for the new layer}
  \item{alpha}{alpha channel for the new layer}

}
\details{
The default overlay is (take red channel for example) \code{r*alpha + r0*(1-alpha)}.
}
\value{
A list which contains overlayed RGB colors.
}
\seealso{
Color overlay function is always used in \code{\link{hc_layer,HilbertCurve-method}} or \code{\link{hc_layer,GenomicHilbertCurve-method}}.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# red (1, 0, 0) overlay to the grey (0.5, 0.5, 0.5) with 0.5 transparency
default_overlay(1, 0, 0, 0.5, 0.5, 0.5, 0.5)
}
