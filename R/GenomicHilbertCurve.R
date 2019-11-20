
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
# - `hc_polygon,GenomicHilbertCurve-method`: add poygons;
# - `hc_text,GenomicHilbertCurve-method`: add text;
# - `hc_layer,GenomicHilbertCurve-method`: add layers undel "pixel" mode;
# - `hc_map,GenomicHilbertCurve-method`: show the map of different categories on the curve. Works both for "normal" and "pixel" mode
#
# The usage of above functions are almost same as those functions for the `HilbertCurve-class`
# except that the second argument which specifies the intervals should be a `GenomicRanges::GRanges` object.
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
# -species abbreviation of species, e.g. 'hg19' or 'mm10'. `circlize::read.chromInfo` is used to retrieve
#          the chromosome information.
# -background the background can be provided as a `GenomicRanges::GRanges` object. Chromosomes should be unique
#          across rows. Or more generally, the 'seqnames' should be different.
# -... common arguments in `HilbertCurve` can be used here.
#
# == details
# Multiple chromosomes can be visualized in a same Hilbert curve. All chromosomes
# are concatenated on after the other based on the order which is specified.
#
# Since chromosomes will have irregular shapes on the curve, under 'pixel' mode, 
# users can set ``border`` option in `hc_map,GenomicHilbertCurve-method` to highlight 
# borders of chromosomes to identify their locations on the curve.
#
# == value
# A `GenomicHilbertCurve-class` object
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# require(circlize)
# require(GenomicRanges)
# bed = generateRandomBed()
# gr = GRanges(seqnames = bed[[1]], ranges = IRanges(bed[[2]], bed[[3]]))
# hc = GenomicHilbertCurve()
# hc_points(hc, gr)
#
# hc = GenomicHilbertCurve(chr = c("chr1", "chr2"))
# hc_points(hc, gr)
#
# bg = GRanges(seqnames = c("chr1", "chr2"), 
#     ranges = IRanges(c(1,10000000), c(10000000,20000000)))
# hc = GenomicHilbertCurve(background = bg, level = 6)
# hc_points(hc, gr, gp = gpar(fill = rand_color(length(gr))))
# hc_map(hc, fill = NA, border = "grey", add = TRUE)
#
GenomicHilbertCurve = function(chr = paste0("chr", c(1:22, "X", "Y")), species = "hg19", 
	background = NULL, ...) {
	
	if(missing(chr)) chr = usable_chromosomes(species)

	if(is.null(background)) {
		chr = unique(chr)
		chromInfo = read.chromInfo(species = species)
		chr.len = chromInfo$chr.len[chr]
		background = GRanges(seqnames = chr, 
			                 ranges = IRanges(rep(1, length(chr)),
				                              chr.len))
		names(background) = chr
	} else {
		background = reduce(background)
		if(!any(duplicated(as.vector(seqnames(background))))) {
			if(is.null(names(background))) {
				names(background) = seqnames(background)
			}
		} else {
			stop("Chromosomes cannot be duplicated in background regions.")
		}
	}

	len = as.numeric(width(background))
    hc = HilbertCurve(1, sum(len), ...)

    hc2 = new("GenomicHilbertCurve")
    for(sn in slotNames(hc)) {
    	slot(hc2, sn) = slot(hc, sn)
    }
    hc2@background = background
    return(hc2)
}


usable_chromosomes = function(species) {
	if(is.null(species)) return(NULL)

	switch(gsub("\\d+$", "", species),
		"hg" = paste0("chr", c(1:22, "X", "Y")),
		"mm" = paste0("chr", c(1:19, "X", "Y")),
		"rn" = paste0("chr", c(1:20, "X", "Y")),
		"dm" = paste0("chr", c("2L", "2R", "3L", "3R", "4", "X")),
		"ce" = paste0("chr", c("I", "II", "III", "IV", "V", "X")),
		"sacCer" = paste0("chr", c("I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X", "XI", "XII", "XIII", "XIV", "XV")),
		NULL
	)
}

# == title
# Add points to the Hilbert curve
#
# == param
# -object a `GenomicHilbertCurve-class` object
# -gr a `GenomicRanges::GRanges` object which contains the genomic regions to be mapped to the curve
# -np pass to `hc_points,HilbertCurve-method`
# -size size of points when ``np <= 1``, pass to `hc_points,HilbertCurve-method`
# -pch shape of the points when ``np <= 1``, pass to `hc_points,HilbertCurve-method`
# -gp graphic parameters of the points when ``np <= 1``, pass to `hc_points,HilbertCurve-method`
# -mean_mode pass to `hc_points,HilbertCurve-method`
# -shape shape of the points when ``np >= 2``, pass to `hc_points,HilbertCurve-method`
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
# require(GenomicRanges)
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
# require(GenomicRanges)
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
# require(GenomicRanges)
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
# -centered_by how to define the "center" of the interval represented in Hilbert curve. Pass to `hc_text,HilbertCurve-method`.
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
# require(GenomicRanges)
# bed = generateRandomBed(nr = 20)
# gr = GRanges(seqnames = bed[[1]], ranges = IRanges(bed[[2]], bed[[3]]))
# hc = GenomicHilbertCurve()
# hc_text(hc, gr, labels = sample(letters, nrow(bed), replace = TRUE))
#
setMethod(f = "hc_text",
	signature = "GenomicHilbertCurve",
	definition = function(object, gr, labels, gp = gpar(), 
	centered_by = c("interval", "polygon"), ...) {

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

	callNextMethod(object, x1 = df[,1], x2 = df[,2], labels = labels, gp = gp, centered_by = centered_by, ...)
})

# == title
# Add text to Hilbert curve
#
# == param
# -object a `GenomicHilbertCurve-class` object
# -gr a `GenomicRanges::GRanges` object which contains the genomic regions to be mapped to the curve
# -gp pass to `hc_polygon,HilbertCurve-method`
# -end_type pass to `hc_polygon,HilbertCurve-method`
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
# require(GenomicRanges)
# bed = generateRandomBed(nr = 20)
# gr = GRanges(seqnames = bed[[1]], ranges = IRanges(bed[[2]], bed[[3]]))
# hc = GenomicHilbertCurve()
# hc_polygon(hc, gr)
#
setMethod(f = "hc_polygon",
	signature = "GenomicHilbertCurve",
	definition = function(object, gr, gp = gpar(), 
	end_type = c("average", "expanding", "shrinking")) {

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

	end_type = match.arg(end_type)[1]
	callNextMethod(object, x1 = df[,1], x2 = df[,2], gp = gp, end_type = end_type)
})

# == title
# Add a new layer to the Hilbert curve
#
# == param
# -object a `GenomicHilbertCurve-class` object
# -gr a `GenomicRanges::GRanges` object which contains the genomic regions to be mapped to the curve
# -col a scalar or a vector of colors which correspond to regions in ``gr``, pass to `hc_layer,HilbertCurve-method`
# -border a scalar or a vector of colors which correspond to the borders of regions. Set it to ``NA`` if borders are suppressed.
# -mean_mode Under 'pixel' mode, each pixel represents a small window. This argument provides methods
#            to summarize value for the small window if the input genomic regions can not completely overlap with the window, 
#            pass to `hc_layer,HilbertCurve-method`
# -grid_line whether add grid lines to show blocks of the Hilber curve, pass to `hc_layer,HilbertCurve-method`
# -grid_line_col color for the grid lines, pass to `hc_layer,HilbertCurve-method`
# -overlay a self-defined function which defines how to overlay new layer to the plot, pass to `hc_layer,HilbertCurve-method`
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
# require(GenomicRanges)
# bed = generateRandomBed()
# gr = GRanges(seqnames = bed[[1]], ranges = IRanges(bed[[2]], bed[[3]]))
# hc = GenomicHilbertCurve(mode = "pixel", level = 9)
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
# Draw a map which represents positions of different chromosomes on the curve
#
# == param
# -object a `GenomicHilbertCurve-class` object
# -level Since a map does not need to have high resolution, a value of around 7 would be enough. 
#        If ``add`` is set to ``TRUE``, ``level`` will be enforced to have the same level in the current Hilbert curve.
# -fill colors for different chromosomes, or more generally, for different 'seqnames'.
# -border colors for the borders of chromosomes. Set it to ``NA`` if borders are suppressed.
# -labels label for each chromosome, or more generally, for different 'sequences'
# -show_labels whether show text labels
# -labels_gp graphic settings for labels
# -add whether add the map to the current curve or draw it in a new graphic device. Notice if ``add`` is set to ``TRUE``,
#      you should set ``fill`` with transparency so that it will not hide your original plot.
# -... pass to `GenomicHilbertCurve`. It is only used if you want the map to be plotted in a new graphic device.
#
# == details
# When multiple genomic categories (e.g. chromosomes) are drawn into one single Hilbert curve, a map which shows the positions
# of categories on the curve is necessary to distinguish different genomic categories.
#
# Under "pixel" mode, if the map is directly added to the Hilbert curve, no chromosome name is drawn. The chromosome names
# are only drawn if the map is plotted in a new graphic device or added to the Hilbert curve under "normal" mode.
#
# Just be careful if you directly overlay the map to the curve that the color of the map does not affect the original
# plot too much.
#
# == value
# A `GenomicHilbertCurve-class` object
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# require(circlize)
# require(GenomicRanges)
# bed = generateRandomBed(nr = 100)
# gr = GRanges(seqnames = bed[[1]], ranges = IRanges(bed[[2]], bed[[3]]))
# hc = GenomicHilbertCurve()
# hc_points(hc, gr, gp = gpar(fill = rand_color(length(gr))))
# # add it in the same graphic device
# hc_map(hc, fill = rand_color(24, transparency = 0.5), add = TRUE)
# 
# # add the map only with borders
# hc = GenomicHilbertCurve()
# hc_points(hc, gr, gp = gpar(fill = rand_color(length(gr))))
# hc_map(hc, fill = NA, border = "grey", add = TRUE)
#
# # or open a new graphic device
# hc_map(hc, fill = rand_color(24))
setMethod(f = "hc_map",
	signature = "GenomicHilbertCurve",
	definition = function(object, level = 7, 
	fill = rand_color(length(background), transparency = 0.5), border = NA,
	labels = names(object@background), show_labels = TRUE, labels_gp = gpar(),
	add = FALSE, ...) {

	background = object@background
	df = merge_into_one_chr(background, object@background)

	if(add) {
		if(object@MODE == "pixel") {
			oi = .ENV$I_PLOT
			seekViewport(paste0("hilbert_curve_", .ENV$I_PLOT))

			hc2 = GenomicHilbertCurve(mode = "normal", chr = unique(as.vector(seqnames(background))), 
				level = level, newpage = FALSE, zoom = object@ZOOM, start_from = object@start_from, first_seg = object@first_seg,
				padding = unit(0, "mm"))
			hc_map(hc2, add = TRUE, labels = labels, fill = fill, border = border, show_labels = show_labels, 
				labels_gp = labels_gp)
			seekViewport(name = paste0("hilbert_curve_", oi, "_global"))
			upViewport()
			
		} else {
			hc_polygon(object, background, gp = gpar(fill = fill, col = border))
			if(show_labels) {
				hc_centered_text(object, x1 = df[, 1], x2 = df[, 2], labels = labels, gp = labels_gp)
			}
		}
	} else {
		hc = GenomicHilbertCurve(background = background, level = level, start_from = object@start_from, ...)
		hc_polygon(hc, background, gp = gpar(fill = fill, col = border))
		if(show_labels) {
			hc_centered_text(hc, x1 = df[, 1], x2 = df[, 2], labels = labels, gp = labels_gp)
		}
	}
	return(invisible(object))
})

# chromosomes are unique across rows
# the only thing that needs to notice is the start site in background may not be zero
merge_into_one_chr = function(gr, background) {
	background_names = names(background)
	background_length = as.numeric(width(background))

	term = cumsum(background_length)
	term = c(0, term[-length(term)])
	names(term) = background_names

	mtch = as.matrix(findOverlaps(gr, background))

	x1 = rep(-1, length(gr))
	x2 = rep(-1, length(gr))

	category = names(background)[mtch[, 2]]
	offset = start(background)[mtch[, 2]] - 1

	l = start(gr)[mtch[,1]] < start(background)[mtch[,2]] & end(gr)[mtch[,1]] >= start(background)[mtch[,2]]
	x1[mtch[, 1]] = ifelse(l, start(background)[mtch[,2]], start(gr)[mtch[, 1]]) - offset + term[as.vector(category)]
	l = start(gr)[mtch[,1]] <= end(background)[mtch[,2]] & end(gr)[mtch[,1]] > end(background)[mtch[,2]]
	x2[mtch[, 1]] = ifelse(l, end(background)[mtch[,2]], end(gr)[mtch[, 1]]) - offset + term[as.vector(category)]
	
	l = start(gr)[mtch[,1]] > end(background)[mtch[,2]]
	x1[mtch[, 1][l]] = -1
	x2[mtch[, 1][l]] = -1
	l = end(gr)[mtch[,1]] < start(background)[mtch[,2]]
	x1[mtch[, 1][l]] = -1
	x2[mtch[, 1][l]] = -1

	return(data.frame(x1, x2))
}
