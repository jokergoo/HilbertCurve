\name{GenomicHilbertCurve}
\alias{GenomicHilbertCurve}
\title{
Initialize a Hilbert curve specifically for genomic data
}
\description{
Initialize a Hilbert curve specifically for genomic data
}
\usage{
GenomicHilbertCurve(chr = paste0("chr", c(1:22, "X", "Y")), species = "hg19",
    background = NULL, ...)
}
\arguments{

  \item{chr}{a vector of chromosome names. Note it should have 'chr' prefix. This argument will be ignored when \code{background} is set.}
  \item{species}{abbreviation of species, e.g. 'hg19' or 'mm10'}
  \item{background}{the background can be privided as a 'GenomicRanges::GRanges' object.  It is not very well supported. It assumes that all regions are subset of background regions.}
  \item{...}{pass to \code{\link{HilbertCurve}}}

}
\details{
If the background region contains more than one categories (e.g. more than one chromosomes), 
 they are concatenated on a same Hilbert curve.
}
\value{
A \code{\link{GenomicHilbertCurve-class}} object
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
require(circlize)
 bed = generateRandomBed()
 gr = GRanges(seqnames = bed[[1]], ranges = IRanges(bed[[2]], bed[[3]]))
 hc = GenomicHilbertCurve()
 hc_points(hc, gr)
hc = GenomicHilbertCurve(chr = c("chr1", "chr2"))
 hc_points(hc, gr)
background = GRanges(seqnames = "chr1", ranges = IRanges(1, 10000000))
 hc = GenomicHilbertCurve(background = background)
 hc_points(hc, gr)
}
