
# == title
# The GenomicHilbertCurve class
#
# == details
# The `GenomicHilbertCurve-class` is inherited from the `HilbertCurve-class`. Basically
# the structure of this class is almost the same as the `HilbertCurve-class` but with several
# additional slots added to facilitate visualizing genomic data.
#
# == Methods
# The `GenomicHilbertCurve-class` provides following methods:
#
# - `GenomicHilbertCurve`: constructor method;
# - `hc_points,GenomicHilbertCurve-method`: add points;
# - `hc_segments,GenomicHilbertCurve-method`: add lines;
# - `hc_rect,GenomicHilbertCurve-method`: add rectangles;
# - `hc_text,GenomicHilbertCurve-method`: add text;
# - `hc_layer,GenomicHilbertCurve-method`: add layers undel "pixel" mode;
# - `hc_map,GenomicHilbertCurve-method`: show the map of different categories on the curve.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# NULL
#
GenomicHilbertCurve = setClass("GenomicHilbertCurve",
	slots = c(getClass("HilbertCurve")@slots,
		      list(background = "GRanges")),
	contains = "HilbertCurve"
)

# == title
# Initialize a Hilbert curve specifically for genomic data
#
# == param
# -chr a vector of chromosome names. Note it should have 'chr' prefix. This argument will be ignored
#      when ``background`` is set.
# -species abbreviation of species, e.g. 'hg19' or 'mm10'
# -background the background can be provided as a 'GenomicRanges::GRanges' object. 
#         It is not very well supported. It assumes that all regions which will be specified in 
#         `hc_layer,GenomicHilbertCurve-method` are subsets of background regions.
# -... pass to `HilbertCurve`
#
# == details
# If the background region contains more than one categories (e.g. more than one chromosomes), 
# they are concatenated onto a same Hilbert curve.
#
# == value
# A `GenomicHilbertCurve-class` object
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
# hc = GenomicHilbertCurve(chr = c("chr1", "chr2"))
# hc_points(hc, gr)
#
# background = GRanges(seqnames = "chr1", ranges = IRanges(1, 10000000))
# hc = GenomicHilbertCurve(background = background)
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
		names(background) = chr
	} else {
		if(length(background) > 50) {
			warning("Too many background regions, it will be hard to interpret the plot.")
		}
		if(!any(duplicated(as.vector(seqnames(background))))) {
			names(background) = seqnames(background)
		} else {
			names(background) = paste0(seqnames(background), ":", start(background), "-", end(background), sep = "")
		}
	}

	chr.len = as.numeric(width(background))
    hc = HilbertCurve(1, sum(chr.len), ...)

    hc2 = new("GenomicHilbertCurve")
    for(sn in slotNames(hc)) {
    	slot(hc2, sn) = slot(hc, sn)
    }
    hc2@background = background
    return(hc2)
}

# == title
# Add points to the Hilbert curve
#
# == param
# -object a `GenomicHilbertCurve-class` object
# -gr a `GenomicRanges::GRanges` object which contains the genomic regions to be mapped to the curve
# -np pass to `hc_points,HilbertCurve-method`
# -size pass to `hc_points,HilbertCurve-method`
# -pch pass to `hc_points,HilbertCurve-method`
# -gp pass to `hc_points,HilbertCurve-method`
# -mean_mode pass to `hc_points,HilbertCurve-method`
# -shape pass to `hc_points,HilbertCurve-method`
#
# == details
# It is basically a wrapper of `hc_points,HilbertCurve-method`.
#
# == value
# Refer to `hc_points,HilbertCurve-method`
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# require(circlize)
# bed = generateRandomBed(nr = 100)
# gr = GRanges(seqnames = bed[[1]], ranges = IRanges(bed[[2]], bed[[3]]))
# hc = GenomicHilbertCurve()
# hc_points(hc, gr, gp = gpar(fill = rand_color(length(gr))))
setMethod(f = "hc_points",
	signature = "GenomicHilbertCurve",
	definition = function(object, gr, 
	np = max(c(2, 10 - hc_level(object))), size = unit(1, "char"), 
	pch = 1, gp = gpar(), mean_mode = c("w0", "absolute", "weighted"),
	shape = "circle") {

	if(is.data.frame(gr)) {
		gr = GRanges(seqnames = gr[[1]], ranges = IRanges(gr[[2]], gr[[3]]))
	}
	if(!inherits(gr, "GRanges")) {
		stop("`gr` should be a `GRanges` object or a data frame.")
	}

	mtch = as.matrix(findOverlaps(gr, object@background))
	gr = gr[mtch[, 1]]
	if(length(size) > 1) {
		size = size[mtch[, 1]]
	}
	if(length(pch) > 1) {
		pch = pch[mtch[, 1]] 
	}
	gp = recycle_gp(gp, length(gr))
	gp = subset_gp(gp, mtch[, 1])

	df = merge_into_one_chr(gr, object@background)
	callNextMethod(object, x1 = df[,1], x2 = df[,2], np = np, size = size, pch = pch, gp = gp, 
		mean_mode = mean_mode, shape = shape)
})

# == title
# Add rectangles on Hilbert curve
#
# == param
# -object a `GenomicHilbertCurve-class` object
# -gr a `GenomicRanges::GRanges` object which contains the genomic regions to be mapped to the curve
# -gp pass to `hc_rect,HilbertCurve-method`
# -mean_mode pass to `hc_rect,HilbertCurve-method`
#
# == details
# It is basically a wrapper of `hc_rect,HilbertCurve-method`.
#
# == value
# Refer to `hc_rect,HilbertCurve-method`
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# require(circlize)
# bed = generateRandomBed(nr = 100)
# gr = GRanges(seqnames = bed[[1]], ranges = IRanges(bed[[2]], bed[[3]]))
# hc = GenomicHilbertCurve()
# hc_rect(hc, gr, gp = gpar(fill = rand_color(length(gr))))
setMethod(f = "hc_rect",
	signature = "GenomicHilbertCurve",
	definition = function(object, gr, gp = gpar(fill = "red", col = "red"), 
	mean_mode = c("w0", "absolute", "weighted")) {

	if(is.data.frame(gr)) {
		gr = GRanges(seqnames = gr[[1]], ranges = IRanges(gr[[2]], gr[[3]]))
	}
	if(!inherits(gr, "GRanges")) {
		stop("`gr` should be a `GRanges` object or a data frame.")
	}

	mtch = as.matrix(findOverlaps(gr, object@background))
	gr = gr[mtch[, 1]]
	gp = recycle_gp(gp, length(gr))
	gp = subset_gp(gp, mtch[, 1])
	df = merge_into_one_chr(gr, object@background)

	callNextMethod(object, x1 = df[,1], x2 = df[,2], gp = gp, mean_mode = mean_mode)
})

# == title
# Add line segments to Hilbert curve
#
# == param
# -object a `GenomicHilbertCurve-class` object
# -gr a `GenomicRanges::GRanges` object which contains the genomic regions to be mapped to the curve
# -gp pass to `hc_segments,HilbertCurve-method`
#
# == details
# It is basically a wrapper of `hc_segments,HilbertCurve-method`.
#
# == value
# Refer to `hc_segments,HilbertCurve-method`
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# require(circlize)
# bed = generateRandomBed(nr = 100)
# gr = GRanges(seqnames = bed[[1]], ranges = IRanges(bed[[2]], bed[[3]]))
# hc = GenomicHilbertCurve()
# hc_segments(hc, gr, gp = gpar(col = rand_color(length(gr))))
#
setMethod(f = "hc_segments",
	signature = "GenomicHilbertCurve",
	definition = function(object, gr, gp = gpar(lty = 1, lwd = 1, col = 1)) {

	if(is.data.frame(gr)) {
		gr = GRanges(seqnames = gr[[1]], ranges = IRanges(gr[[2]], gr[[3]]))
	}
	if(!inherits(gr, "GRanges")) {
		stop("`gr` should be a `GRanges` object or a data frame.")
	}

	mtch = as.matrix(findOverlaps(gr, object@background))
	gr = gr[mtch[, 1]]
	gp = recycle_gp(gp, length(gr))
	gp = subset_gp(gp, mtch[, 1])

	df = merge_into_one_chr(gr, object@background)

	callNextMethod(object, x1 = df[,1], x2 = df[,2], gp = gp)
})

# == title
# Add text to Hilbert curve
#
# == param
# -object a `GenomicHilbertCurve-class` object
# -gr a `GenomicRanges::GRanges` object which contains the genomic regions to be mapped to the curve
# -labels pass to `hc_text,HilbertCurve-method`
# -gp pass to `hc_text,HilbertCurve-method`
# -... pass to `hc_text,HilbertCurve-method`
#
# == details
# It is basically a wrapper of `hc_text,HilbertCurve-method`.
#
# == value
# Refer to `hc_text,HilbertCurve-method`
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# require(circlize)
# bed = generateRandomBed(nr = 20)
# gr = GRanges(seqnames = bed[[1]], ranges = IRanges(bed[[2]], bed[[3]]))
# hc = GenomicHilbertCurve()
# hc_text(hc, gr, labels = letters[1:20])
#
setMethod(f = "hc_text",
	signature = "GenomicHilbertCurve",
	definition = function(object, gr, labels, gp = gpar(), ...) {

	if(is.data.frame(gr)) {
		gr = GRanges(seqnames = gr[[1]], ranges = IRanges(gr[[2]], gr[[3]]))
	}
	if(!inherits(gr, "GRanges")) {
		stop("`gr` should be a `GRanges` object or a data frame.")
	}

	mtch = as.matrix(findOverlaps(gr, object@background))
	gr = gr[mtch[, 1]]
	labels = labels[mtch[, 1]]
	gp = recycle_gp(gp, length(gr))
	gp = subset_gp(gp, mtch[, 1])

	df = merge_into_one_chr(gr, object@background)

	callNextMethod(object, x1 = df[,1], x2 = df[,2], labels = labels, gp = gp)
})

# == title
# Add text to Hilbert curve
#
# == param
# -object a `GenomicHilbertCurve-class` object
# -gr a `GenomicRanges::GRanges` object which contains the genomic regions to be mapped to the curve
# -gp pass to `hc_text,HilbertCurve-method`
# -end_type pass to `hc_text,HilbertCurve-method`
# -... pass to `hc_text,HilbertCurve-method`
#
# == details
# It is basically a wrapper of `hc_polygon,HilbertCurve-method`.
#
# == value
# Refer to `hc_polygon,HilbertCurve-method`
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# require(circlize)
# bed = generateRandomBed(nr = 20)
# gr = GRanges(seqnames = bed[[1]], ranges = IRanges(bed[[2]], bed[[3]]))
# hc = GenomicHilbertCurve()
# hc_polygon(hc, gr)
#
setMethod(f = "hc_polygon",
	signature = "GenomicHilbertCurve",
	definition = function(object, gr, gp = gpar(), 
		end_type = c("average", "expanding", "shrinking"), ...) {

	if(is.data.frame(gr)) {
		gr = GRanges(seqnames = gr[[1]], ranges = IRanges(gr[[2]], gr[[3]]))
	}
	if(!inherits(gr, "GRanges")) {
		stop("`gr` should be a `GRanges` object or a data frame.")
	}

	mtch = as.matrix(findOverlaps(gr, object@background))
	gr = gr[mtch[, 1]]
	gp = recycle_gp(gp, length(gr))
	gp = subset_gp(gp, mtch[, 1])
	df = merge_into_one_chr(gr, object@background)

	callNextMethod(object, x1 = df[,1], x2 = df[,2], gp = gp, end_type = end_type)
})

# == title
# Add a new layer to the Hilbert curve
#
# == param
# -object a `GenomicHilbertCurve-class` object
# -gr a `GenomicRanges::GRanges` object which contains the genomic regions to be mapped to the curve
# -col pass to `hc_layer,HilbertCurve-method`
# -border colors for the borders of every regions. Set it to ``NA`` if borders are suppressed.
# -mean_mode pass to `hc_layer,HilbertCurve-method`
# -grid_line pass to `hc_layer,HilbertCurve-method`
# -grid_line_col pass to `hc_layer,HilbertCurve-method`
# -overlay pass to `hc_layer,HilbertCurve-method`
#
# == details
# It is basically a wrapper of `hc_layer,HilbertCurve-method`.
#
# == value
# Refer to `hc_layer,HilbertCurve-method`
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# require(circlize)
# bed = generateRandomBed()
# gr = GRanges(seqnames = bed[[1]], ranges = IRanges(bed[[2]], bed[[3]]))
# hc = GenomicHilbertCurve(mode = "pixel")
# hc_layer(hc, gr, col = rand_color(length(gr)))
#
setMethod(f = "hc_layer",
	signature = "GenomicHilbertCurve",
	definition = function(object, gr, col = "red", border = NA,
	mean_mode = c("w0", "absolute", "weighted"), grid_line = 0,
	grid_line_col = "black", overlay = default_overlay) {

	if(is.data.frame(gr)) {
		gr = GRanges(seqnames = gr[[1]], ranges = IRanges(gr[[2]], gr[[3]]))
	}
	if(!inherits(gr, "GRanges")) {
		stop("`gr` should be a `GRanges` object or a data frame.")
	}

	mtch = as.matrix(findOverlaps(gr, object@background))
	if(length(col) == 1) col = rep(col, length(gr))
	gr = gr[mtch[, 1]]
	col = col[mtch[,1 ]]

	df = merge_into_one_chr(gr, object@background)

	callNextMethod(object, x1 = df[,1], x2 = df[,2], col = col, border = border, mean_mode = mean_mode, 
		grid_line = grid_line, grid_line_col = grid_line_col, overlay = overlay)
})

# == title
# Draw a map which represents positions of different genomic categories on the curve
#
# == param
# -object a `GenomicHilbertCurve-class` object
# -level Since a map does not need to have high resolution, a value of around 6 would be enough. 
#        If ``add`` is set to ``TRUE``, ``level`` will be enforced to have the same level in the current Hilbert curve.
# -fill colors for different genomic categories (current you cannot adjust the style of the borders)
# -border colors for the borders of every regions. Set it to ``NA`` if borders are suppressed.
# -labels label for each genomic category. By default, if the category is unique for different chromosome,
#          the label will be the chromosome name, or else they will be the combination of chromosme names and positions.
#          It is ignored if the curve is under 'pixel' mode.
# -labels_gp graphic settings for labels
# -add whether add the map to the current curve or draw it in a new graphic device. Notice if ``add`` is set to ``TRUE``,
#      you should set ``fill`` with transparency so that it will not hide your original plot.
# -... pass to `GenomicHilbertCurve`.
#
# == details
# When multiple genomic categories are drawn into one single Hilbert curve, a map which shows the positions
# of different categories on the curve is necessary to correspond to the graphics on the curve.
#
# == value
# A `GenomicHilbertCurve-class` object
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# require(circlize)
# bed = generateRandomBed(nr = 100)
# gr = GRanges(seqnames = bed[[1]], ranges = IRanges(bed[[2]], bed[[3]]))
# hc = GenomicHilbertCurve()
# hc_points(hc, gr, gp = gpar(fill = rand_color(length(gr))))
# # add it in the same graphic device
# hc_map(hc, fill = rand_color(24, transparency = 0.5), add = TRUE)
#
# # or open a new graphic device
# hc_map(hc, fill = rand_color(24))
setMethod(f = "hc_map",
	signature = "GenomicHilbertCurve",
	definition = function(object, level = 7, fill = NULL, border = NA,
	labels = names(object@background), labels_gp = gpar(),
	add = FALSE, ...) {

	background = object@background
	df = merge_into_one_chr(background, object@background)

	if(is.null(fill)) {
		fill = rand_color(length(background), transparency = 0.5)
	}
	if(add) {
		if(object@MODE == "pixel") {
			hc_layer(object, background, col = fill, border = border)
		} else {
			hc_polygon(object, background, gp = gpar(fill = fill, col = border))
			hc_centered_text(object, x1 = df[, 1], x2 = df[, 2], labels = labels, gp = labels_gp)
		}
	} else {
		hc = GenomicHilbertCurve(background = background, level = level, start_from = object@start_from, ...)
		hc_polygon(hc, background, gp = gpar(fill = fill, col = border))
		hc_centered_text(hc, x1 = df[, 1], x2 = df[, 2], labels = labels, gp = labels_gp)
	}
	return(invisible(object))
})


merge_into_one_chr = function(gr, background) {
	background_names = as.vector(seqnames(background))
	background_length = as.numeric(width(background))

	mtch = findOverlaps(gr, background, select = "arbitrary")
	gr$category = seqnames(background)[mtch]

	term = cumsum(background_length)
	term = c(0, term[-length(term)])
	names(term) = background_names

	x1 = start(gr) + term[as.vector(gr$category)]
	x2 = end(gr) + term[as.vector(gr$category)]
	return(data.frame(x1, x2))
}
