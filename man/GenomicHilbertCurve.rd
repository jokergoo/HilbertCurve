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
  \item{species}{abbreviation of species, e.g. 'hg19' or 'mm10'. \code{\link[circlize]{read.chromInfo}} is used to retrieve the chromosome information.}
  \item{background}{the background can be provided as a \code{\link[GenomicRanges]{GRanges}} object. Chromosomes should be unique across rows. Or more generally, the 'seqnames' should be different.}
  \item{...}{common arguments in \code{\link{HilbertCurve}} can be used here.}

}
\details{
Multiple chromosomes can be visualized in a same Hilbert curve. All chromosomes
are concatenated on after the other based on the order which is specified.

Since chromosomes will have irregular shapes on the curve, under 'pixel' mode, 
users can set \code{border} option in \code{\link{hc_map,GenomicHilbertCurve-method}} to highlight 
borders of chromosomes to identify their locations on the curve.
}
\value{
A \code{\link{GenomicHilbertCurve-class}} object
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
require(circlize)
require(GenomicRanges)
bed = generateRandomBed()
gr = GRanges(seqnames = bed[[1]], ranges = IRanges(bed[[2]], bed[[3]]))
hc = GenomicHilbertCurve()
hc_points(hc, gr)

hc = GenomicHilbertCurve(chr = c("chr1", "chr2"))
hc_points(hc, gr)

bg = GRanges(seqnames = c("chr1", "chr2"), 
    ranges = IRanges(c(1,10000000), c(10000000,20000000)))
hc = GenomicHilbertCurve(background = bg, level = 6)
hc_points(hc, gr, gp = gpar(fill = rand_color(length(gr))))
hc_map(hc, fill = NA, border = "grey", add = TRUE)
}
