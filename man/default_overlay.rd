\name{default_overlay}
\alias{default_overlay}
\title{
Default overlay for adding new layers
}
\description{
Default overlay for adding new layers
}
\usage{
default_overlay(r0, g0, b0, r, g, b, alpha = rep(1, length(r0)))
}
\arguments{

  \item{r0}{r channel for the layers that are already added.}
  \item{g0}{g channel for the layers that are already added.}
  \item{b0}{b channel for the layers that are already added.}
  \item{r}{r channel for the new layer}
  \item{g}{g channel for the new layer}
  \item{b}{b channel for the new layer}
  \item{alpha}{alpha channel for the new layer}

}
\details{
The default overlay is (take r channel for example) \code{r*alpha + r0*(1-alpha)}.
}
\value{
A list contains overlayed RGB channel.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
NULL
}
