
# == title
# The HilbertCurve class
#
# == details
# The `GenomicHilbertCurve-class` inherits the `HilbertCurve-class`. Basically
# the structure of this class is same as the `HilbertCurve-class` but with several
# additional slots to faciliated for visualizing genomic data.
#
# == Methods
# The `GenomicHilbertCurve-class` provides following methods:
#
# - `GenomicHilbertCurve`: constructor method;
# - `hc_points,GenomicHilbertCurve-method`: add points;
# - `hc_segments,GenomicHilbertCurve-method`: add lines;
# - `hc_rect,GenomicHilbertCurve-method`: add rectangles;
# - `hc_text,GenomicHilbertCurve-method`: add text;
# - `hc_layer,GenomicHilbertCurve-method`: add layers;
# - `hc_map,GenomicHilbertCurve-method`: show the map of chromosomes.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# NULL
#
GenomicHilbertCurve = setClass("GenomicHilbertCurve",
	slots = c(getClass("Heatmap")@slots,
		      list(background = "GRanges"))
	contains = "HilbertCurve"
)

# == title
# Initialize a Hilbert curve specifically for genomic data
#
# == param
# -chr a vector of chromosome names. Note it should have 'chr' prefix
# -species species
# -background the background can be privided as a 'GenomicRanges::GRanges' object.
# -... pass to `HilbertCurve`
#
# == details
# If more than one chromosomes are selected, they are concatenated on a same Hilbert curve.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# require(circlize)
# bed = generateRandomBed()
# gr = GRanges(seqnames = bed[[1]], ranges = IRanges(bed[[2]], bed[[3]]))
# hc = GenomicHilbertCurve()
# hc_points(hc, gr)
#
# hc = GenomicHilbertCurve(chr = c("chr1", "chr2'"))
# hc_points(hc, gr)
GenomicHilbertCurve = function(chr = paste0("chr", c(1:22, "X", "Y")), species = "hg19", 
	background = NULL, ...) {

	if(is.null(background)) {
		chr = unique(chr)
		chromInfo = read.chromInfo(species = species)
		chr.len = chromInfo$chr.len[chr]
		background = GRanges(seqnames = chr, 
			                 ranges = IRanges(rep(1, length(chr)),
				                              chr.len))
	} else {
		if(any(duplicated(as.vector(seqnames(background))))) {
			stop("Currently, seqnames in `background` should be unique.")
		}
	}
	
    hc = HilbertCurve(1, sum(chr.len), ...)

    hc2 = new("GenomicHilbertCurve")
    for(sn in slotnames(hc)) {
    	slot(hc2, sn) = slot(hc, sn)
    }
    hc@background = background
    return(hc)
}

# == 
setMethod(f = "hc_points",
	signature = "GenomicHilbertCurve",
	definition = function(object, gr, 
	np = max(c(2, 10 - hc_level(object))), size = unit(1, "char"), 
	pch = 1, gp = gpar(), mean_mode = c("w0", "absolute", "weighted"),
	shape = c("circle", "square", "triangle", "hexagon", "star")) {

	mtch = findOverlaps(gr, object@background)
	gr = gr[mtch[, 1]]
	if(length(size) > 1) {
		size = size[mtch[, 1]]
	}
	if(length(pch) > 1) {
		pch = pch[mtch[, 1]] 
	}
	gp = recycle_gp(gp, nrow(gr))
	gp = subset_gp(gp, mtch[, 1])

	df = merge_into_one_chr(gr, object@background)

	callNextMethod(object, x1 = df[,1], x2 = df[,2], size = size, pch = pch, gp = gp, 
		mean_mode = mean_mode, shape = shape)
})


setMethod(f = "hc_rect",
	signature = "GenomicHilbertCurve",
	definition = function(object, gr, gp = gpar(fill = "red", col = "red"), 
	mean_mode = c("w0", "absolute", "weighted")) {

	mtch = findOverlaps(gr, object@background)
	gr = gr[mtch[, 1]]
	gp = recycle_gp(gp, nrow(gr))
	gp = subset_gp(gp, mtch[, 1])

	df = merge_into_one_chr(gr, object@background)

	callNextMethod(object, x1 = df[,1], x2 = df[,2], gp = gp, mean_mode = mean_mode)
})

setMethod(f = "hc_segments",
	signature = "GenomicHilbertCurve",
	definition = function(object, gr, gp = gpar(lty = 1, lwd = 1, col = 1)) {

	mtch = findOverlaps(gr, object@background)
	gr = gr[mtch[, 1]]
	gp = recycle_gp(gp, nrow(gr))
	gp = subset_gp(gp, mtch[, 1])

	df = merge_into_one_chr(gr, object@background)

	callNextMethod(object, x1 = df[,1], x2 = df[,2], gp = gp)
})

setMethod(f = "hc_text",
	signature = "GenomicHilbertCurve",
	definition = function(object, gr, labels, gp = gpar(), ...) {

	mtch = findOverlaps(gr, object@background)
	gr = gr[mtch[, 1]]
	gp = recycle_gp(gp, nrow(gr))
	gp = subset_gp(gp, mtch[, 1])

	df = merge_into_one_chr(gr, object@background)

	callNextMethod(object, x1 = df[,1], x2 = df[,2], labels, gp = gp)
})

setMethod(f = "hc_layer",
	signature = "GenomicHilbertCurve",
	definition = function(object, gr, col = "red", 
	mean_mode = c("w0", "absolute", "weighted"), grid_line = 0) {

	mtch = findOverlaps(gr, object@background)
	gr = gr[mtch[, 1]]
	col = col[mtch[,1 ]]

	df = merge_into_one_chr(gr, object@background)

	callNextMethod(object, x1 = df[,1], x2 = df[,2], col = col, mean_mode = mean_mode, 
		grid_line = grid_line)
})

setMethod(f = "hc_map",
	signature = "GenomicHilbertCurve",
	definition = function(object, level = 6, fill = NULL, labels_gp = gpar(fontsize = 16)) {
		background = object@background_names
		if(is.null(fill)) {
			fill = rand_color(length(background))
		}
		hc = GenomicHilbertCurve(background = background, level = level)
		hc_rect(hc, background, gp = gpar(fill = fill, col = NA))
		hc_text(hc, as.vector(seqnames(background)), gp = labels_gp)
		return(object)
})

merge_into_one_chr = function(gr, background) {
	background_names = as.vector(seqnames(background))
	background_length = width(background)

	term = numeric()
	for(i in seq_along(background_names)) {
		term[i] = cumsum(background_length[seq_len(i-1)])
	}
	term[1] = 0
	names(term) = background_names

	x1 = start(gr) + term[as.vector(seqnames(gr))]
	x2 = end(gr) + term[as.vector(seqnames(gr))]
	return(data.frame(x1, x2))
}


recycle_gp = function(gp, n = 1) {
	g = lapply(gp, function(x) {
		if(length(x) == 1 && n > 1) {
			rep(x, n)
		} else {
			x
		}
	})
	class(g) = "gpar"
	return(g)
}

subset_gp = function(gp, i = 1) {
	g = lapply(gp, function(x) {
		x[i]
	})
	class(g) = "gpar"
	return(g)
}
