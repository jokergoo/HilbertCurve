
library(IRanges)
library(HilbertVis)

hc = HilbertCurve(1, 1600, level = 8)

x = sort(sample(1600, 120))
s = x[1:60*2 - 1]
e = x[1:60*2]
i = sample(60, 10)
ir = IRanges(s[i], e[i])

#hc_points(ir)
hc_rect(hc, ir)



hc = HilbertCurve(1, 1600, level = 8, mode = "pixel")
hc_layer(hc, ir)
hc_save(hc)


library(TxDb.Hsapiens.UCSC.hg19.lincRNAsTranscripts)
g = genes(TxDb.Hsapiens.UCSC.hg19.lincRNAsTranscripts)
g = g[seqnames(g) == "chr1"]
ir = ranges(g)

hc = HilbertCurve(1, 249250621, level = 5, reference = TRUE)
w = width(ir)
w[w > 50000] = 50000
size = w/20000
hc_points(hc, ir, np = 1, pch = 16, size = unit(size, "mm"))
