
.ENV = new.env()

.ENVI_PLOT = 0

get_plot_index = function() {
	.ENV$I_PLOT
}

increase_plot_index = function() {
	.ENV$I_PLOT = .ENV$I_PLOT + 1
}


# == title
# The HilbertCurve class
#
# == 
HilbertCurve = setClass("HilbertCurve",
	slots = list(
		BINS = "IRanges",
		POS = "data.frame",
		RANGE = "numeric",
		ZOOM = "numeric",
		MODE = "character",
		RGB = "environment",
		LEVEL = "integer"
))

# == title
# Zoom original positions
#
# == param
# -object A `HilbertCurve-class` object
# -x positions
#
# == details
# The function is used internally
#
# == value
# Zoomed positions
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
setMethod(f = "zoom",
	signature = "HilbertCurve",
	definition = function(object, x) {
	x * object@ZOOM
})

# == title
# Transform zoomed positions to their original values
#
# == param
# -object A `HilbertCurve-class` object
# -x positions
#
# == details
# The function is used internally
#
# == value
# Original positioins
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
setMethod(f = "unzoom",
	signature = "HilbertCurve",
	definition = function(object, x) {
	x / object@ZOOM
})

# == title
# Initialize a Hilbert curve
#
# == param
# -s start of the Hilbert curve, should be an integer
# -e end of the Hilbert curve, should be an integer
# -level order of the Hilbert curve. There will by ``4^level`` segments in the Hilbert curve.
# -mode make it like a normal R plot or write the plot directly into png file.
# -reference add reference information on the plot
# -zoom zooming of the position ranges
# -newpage whether call `grid::grid.newpage``
# -background background color
#
# == value
# A `HilbertCurve-class` object.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
HilbertCurve = function(s, e, level = 4, mode = c("normal", "pixel"),
	reference = FALSE, zoom = NULL, newpage = TRUE, background = "white") {

	hc = new("HilbertCurve")
	level = as.integer(level)

	if(is.null(zoom)) zoom = ceiling(10*(4^level-1)/(e - s))
	hc@ZOOM = zoom

	pos = hilbertCurve(level)
	n = 4^level  # number of points

	breaks = round(seq(zoom(hc, s), zoom(hc, e), length = n)) 

	# n - 1 rows
	bins = IRanges(start = breaks[-length(breaks)], end = breaks[-1]) # number of segments

	# x and y correspond to the end point of each interval in ir
	hc@BINS = bins

	# position for two end points for every segment (e.g. interval in bins)
	hc@POS = data.frame(x1 = pos$x[-n], y1 = pos$y[-n], x2 = pos$x[-1], y2 = pos$y[-1])
	hc@RANGE = c(breaks[1], breaks[length(breaks)])
	hc@LEVEL = level

	mode = match.arg(mode)[1]
	hc@MODE = mode

	if(mode == "normal") {
		if(newpage) grid.newpage()

		increase_plot_index()
		pushViewport(viewport(name = paste0("hilbert_curve_", get_plot_index()), xscale = c(-0.5, 2^level - 0.5), yscale = c(-0.5, sqrt(n)-0.5)))
		
		h = round(convertHeight(unit(1, "npc"), "mm", valueOnly = TRUE))
		w = round(convertWidth(unit(1, "npc"), "mm", valueOnly = TRUE))
		if(h != w) {
			warning(sprintf("Hilbert curve has height %imm and width %imm, but it should be put in a squared viewport.\n", h, w))
		}

		if(reference) {
			grid.segments(hc@POS$x1, hc@POS$y1, hc@POS$x2, hc@POS$y2, default.units = "native", gp = gpar(lty = 3, col = "#CCCCCC"))
			
			# grid.points(hc@POS$x1[1], hc@POS$y1[1], default.units = "native", gp = gpar(col = "#CCCCCC", cex = 0.5))
			# grid.points(hc@POS$x2, hc@POS$y2, default.units = "native", gp = gpar(col = "#CCCCCC", cex = 0.5))

			# grid.text(round(unzoom(hc, start(bins)[1])), hc@POS$x1[1], hc@POS$y1[1], default.units = "native", gp = gpar(col = "#999999", cex = 0.5))
			# grid.text(round(unzoom(hc, end(bins))), hc@POS$x2, hc@POS$y2, default.units = "native", gp = gpar(col = "#999999", cex = 0.5))
		
			grid_arrows(hc@POS$x1, hc@POS$y1, (hc@POS$x1+hc@POS$x2)/2, (hc@POS$y1+hc@POS$y2)/2, only.head = TRUE, arrow_gp = gpar(fill = "#CCCCCC", col = NA))
		}

		upViewport()
	} else {
		background = background[1]
		background = col2rgb(background) / 255
		red = matrix(background[1], nrow = 2^level, ncol = 2^level)
		green = matrix(background[2], nrow = 2^level, ncol = 2^level)
		blue = matrix(background[3], nrow = 2^level, ncol = 2^level)
		hc@RGB$red = red
		hc@RGB$green = green
		hc@RGB$blue = blue
	}

	return(invisible(hc))
}

# == title
# Print the HilbertCurve object
# 
# == param
# -object A `HilbertCurve-class` object
#
# == value
# No value is returned
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
setMethod(f = "show",
	signature = "HilbertCurve",
	definition = function(object) {
	cat("A Hilbert curve with order", object@LEVEL, "\n")
	cat("mode:", object@MODE, "\n")
})

# == title
# Level of the Hilbert curve
#
# == param
# -object A `HilbertCurve-class` object
#
# == value
# The level of the Hilbert curve
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
setMethod(f = "hc_level",
	signature = "HilbertCurve",
	definition = function(object) {
	object@LEVEL
})

# == title
# Add points to the Hilbert curve
#
# == param
# -object A `HilbertCurve-class` object
# -ir a `IRanges::IRanges` object
# -np number of points (a circle or a square, ...) that are put in a segment
# -size size of the points
# -pch shape of points, used for points if ``np >= 2``
# -gp graphical parameters for points
# -mean_mode
# -shape shape of points, used for points if ``np <= 1``
#
# == details
# If ``np`` is set to a value less than 2 or ``NULL``, points will be added at every middle point in ``ir``.
# If ``np`` is set to a value larger or equal to 2, a list of e.g. circles are put at every segment in ``ir``,
# so, longer segments will have more circles on it.
#
# == value
# A data frame which contains coordinates for points.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
setMethod(f = "hc_points",
	signature = "HilbertCurve",
	definition = function(object, ir, np = max(c(2, 10 - hc_level(object))), 
	size = unit(1, "char"), pch = 1, gp = gpar(), mean_mode = c("w0", "absolute", "weighted"),
	shape = c("circle", "square", "triangle", "hexagon", "star")) {

	if(object@MODE == "pixel") {
		stop("`hc_points()` can only be used under 'normal' mode.")
	}

	if(np >= 2) {
		hc_segmented_points(object, ir, gp = gp, np = np, mean_mode = mean_mode, shape = shape)
	} else {
		hc_normal_points(object, ir, gp = gp, pch = pch, size = size)
	}

})

# == title
# Add points to the Hilbert curve
#
# == param
# -object A `HilbertCurve-class` object
# -ir a `IRanges::IRanges` object
# -size size of the points
# -pch shape of points, used for points if ``np >= 2``
# -gp graphical parameters for points
#
# == details
# Points are added at every middle point in ``ir``.
#
# This function is used internally.
#
# == value
# A data frame which contains coordinates for points.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
setMethod(f = "hc_normal_points",
	signature = "HilbertCurve",
	definition = function(object, ir, gp = gpar(), pch = 1, size = unit(1, "char")) {

	ir = IRanges(start = zoom(object, (start(ir) + end(ir))/2),
	             end = zoom(object, (start(ir) + end(ir))/2))

	seekViewport(name = paste0("hilbert_curve_", get_plot_index()))

	od = order(ir)
	ir = ir[od]

	mtch = as.matrix(findOverlaps(object@BINS, ir))
	l = duplicated(mtch[,2])
	mtch = mtch[!l, , drop = FALSE]

	r1 = object@BINS[mtch[,1]]
	pos = object@POS[mtch[,1], ]
	r = pintersect(r1, ir[mtch[,2]])

	sr1 = start(r1)
	sr = start(r)
	er1 = end(r1)
	er = end(r)

	x1 = pos$x2 - (pos$x2 - pos$x1)*(er1 - sr)/(er1 - sr1)
	y1 = pos$y2 - (pos$y2 - pos$y1)*(er1 - sr)/(er1 - sr1)

	grid.points(x1, y1, default.units = "native", gp = gp, pch = pch, size = size)	
	
	df = data.frame(x = x1, y = y1, stringsAsFactors = FALSE)
	return(invisible(df))
})

# == title
# Add points to the Hilbert curve
#
# == param
# -object A `HilbertCurve-class` object
# -ir a `IRanges::IRanges` object
# -np number of points (a circle or a square, ...) that are put in a segment
# -gp graphical parameters for points
# -mean_mode
# -shape shape of points, used for points if ``np <= 1``
#
# == details
# A list of e.g. circles are put at every segment in ``ir``,
# so, longer segments will have more circles on it.
#
# This function is used internally.
#
# == value
# A data frame which contains coordinates for points.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
setMethod(f = "hc_segmented_points",
	signature = "HilbertCurve",
	definition = function(object, ir, gp = gpar(), np = max(c(2, 10 - hc_level(object))),
	mean_mode = c("w0", "absolute", "weighted"), shape = c("circle", "square", "triangle", "hexagon", "star")) {

	ir = IRanges(start = zoom(object, start(ir)),
	             end = zoom(object, end(ir)))

	col = normalize_gp("col", gp$col, length(ir))
	fill = normalize_gp("fill", gp$fill, length(ir))

	seekViewport(name = paste0("hilbert_curve_", get_plot_index()))

	mtch = as.matrix(findOverlaps(object@BINS, ir))
	r = pintersect(object@BINS[mtch[,1]], ir[mtch[,2]])
	l = width(r) > 1
	r = r[l]
	mtch = mtch[l, , drop = FALSE]

	r1 = object@BINS[mtch[,1]]
	pos = object@POS[mtch[,1], ]
	col = col[mtch[,2]]
	fill = fill[mtch[,2]]

	sr1 = start(r1)
	sr = start(r)
	er1 = end(r1)
	er = end(r)

	x1 = pos$x2 - (pos$x2 - pos$x1)*(er1 - sr)/(er1 - sr1)
	y1 = pos$y2 - (pos$y2 - pos$y1)*(er1 - sr)/(er1 - sr1)
	x2 = pos$x2 - (pos$x2 - pos$x1)*(er1 - er)/(er1 - sr1)
	y2 = pos$y2 - (pos$y2 - pos$y1)*(er1 - er)/(er1 - sr1)

	#grid.segments(x1, y1, x2, y2, default.units = "native", gp = gpar(lwd = 2, col = "red"))
	shape = match.arg(shape)[1]
	xx = numeric()
	yy = numeric()
	for(i in seq_along(x1)) {
		if(pos$x1[i] == pos$x2[i]) {
			x = rep(pos$x1[i], np)
		} else if(pos$x1[i] > pos$x2[i]) {
			x = seq(pos$x1[i], pos$x2[i], by = -1/(np-1))
		} else {
			x = seq(pos$x1[i], pos$x2[i], by = 1/(np-1))
		}
		if(pos$y1[i] == pos$y2[i]) {
			y = rep(pos$y1[i], np)
		} else if(pos$y1[i] > pos$y2[i]) {
			y = seq(pos$y1[i], pos$y2[i], by = -1/(np-1))
		} else {
			y = seq(pos$y1[i], pos$y2[i], by = 1/(np-1))
		}

		l = x >= min(c(x1[i], x2[i])) & x <= max(c(x1[i], x2[i])) &
		    y >= min(c(y1[i], y2[i])) & y <= max(c(y1[i], y2[i]))
		x = x[l]
		y = y[l]

		if(length(x)) {
			r = 1/(np-1)/2
			if(shape == "circle") {
				grid.circle(x, y, r = 1/(np-1)/2, default.units = "native", gp = gpar(fill = fill[i], col = col[i]))
			} else if(shape == "square") {
				grid.rect(x, y, width = 2*r, height = 2*r, default.units = "native", gp = gpar(fill = fill[i], col = col[i]))
			} else if(shape == "triangle") {
				grid.polygon(c(x, x-2*r*tan(pi/6), x+2*r*tan(pi/6)), 
					c(y+r, y-r, y-r), id = rep(seq_along(x), 3), default.units = "native", gp = gpar(fill = fill[i], col = col[i]))
			} else if(shape == "hexagon") {
				grid.polygon(rep(cos(0:5 * pi/3)*r, length(x)) + rep(x, each = 6),
					         rep(sin(0:5 * pi/3)*r, length(y)) + rep(y, each = 6),
					         id = rep(seq_along(x), each = 6),
					         default.units = "native", gp = gpar(fill = fill[i], col = col[i]))
			} else if(shape == "star") {
				r2 = cos(pi/10)*tan(pi/20)/(tan(pi/20) + tan(pi/10)) * r
				grid.polygon(rep(sin(0:9*pi/5 + pi*0)*rep(c(r,r2), 5), length(x)) + rep(x, each = 10),
					         rep(cos(0:9*pi/5 + pi*0)*rep(c(r,r2), 5), length(y)) + rep(y, each = 10),
					         id = rep(seq_along(x), each = 10),
					         default.units = "native", gp = gpar(fill = fill[i], col = col[i]))
			}

			xx = c(xx, x)
			yy = c(yy, y)
		}
	}

	df = data.frame(x = xx, y = yy, r = rep(1/(np-1)/2, length(xx)))
	return(invisible(df))
})

average_in_window = function(window, ir, mtch, v, mean_mode, empty_v = 0) {
 	intersect = pintersect(window[mtch[,1]], ir[mtch[,2]])
 	v = v[mtch[,2]]

	if(mean_mode == "w0") {
		p = (width(intersect)-1)/(width(window[mtch[,1]])-1)
		x = tapply(p*v, mtch[, 1], sum)
		sum_p = tapply(p, mtch[, 1], sum)
		x2 = (1-sum_p)*empty_v
		x = x + x2
	} else if(mean_mode == "absolute") {
		x = tapply(v, mtch[, 1], mean)
		
	} else {
		w = width(intersect)-1
		x = tapply(w*v, mtch[, 1], sum) / tapply(w, mtch[, 1], sum)
	}

	x = structure(as.vector(x), names = dimnames(x)[[1]])
	
	return(x)
}

# == title
# Add rectangles on Hilbert curve
#
# == param
# -object A `HilbertCurve-class` object
# -ir a `IRanges::IRanges` object
# -gp graphical parameters for rectangles
# -mean_mode
#
# == value
# A data frame which contains coordinates for rectangles.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
setMethod(f = "hc_rect",
	signature = "HilbertCurve",
	definition = function(object, ir, gp = gpar(fill = "red", col = "red"), 
	mean_mode = c("w0", "absolute", "weighted")) {

	if(object@MODE == "pixel") {
		stop("`hc_rect()` can only be used under 'normal' mode.")
	}

	ir = IRanges(start = zoom(object, start(ir)),
	             end = zoom(object, end(ir)))

	seekViewport(name = paste0("hilbert_curve_", get_plot_index()))

	s = start(object@BINS)
	e = end(object@BINS)
	mid = round((s + e)/2)
	window = IRanges(start = c(s[1], mid), end = c(mid, e[length(e)]))
		
	mtch = as.matrix(findOverlaps(window, ir))
	
	# colors correspond to mean value for each window
	mean_mode = match.arg(mean_mode)[1]

	fill = normalize_gp("fill", gp$fill, length(ir))

	rgb = col2rgb(fill, alpha = TRUE)

	r = average_in_window(window, ir, mtch, rgb[1, ], mean_mode, 255)
	g = average_in_window(window, ir, mtch, rgb[2, ], mean_mode, 255)
	b = average_in_window(window, ir, mtch, rgb[3, ], mean_mode, 255)
	alpha = rep(max(rgb[4, ]), length(r))

	index = as.numeric(names(r))
	fill2 = rgb(red = r, green = g, blue = b, alpha = alpha, maxColorValue = 255)

	pos = object@POS
	n = 4^object@LEVEL

	if(n %in% index) {
		ind = setdiff(index, n)
		grid.rect(pos$x1[ind], pos$y1[ind], width = 1, height = 1, default.units = "native", gp = gpar(fill = fill2[-length(fill2)], col = NA, lineend = "butt", linejoin = "mitre"))
		df = data.frame(x = pos$x1[ind], y = pos$y1[ind])

		grid.rect(pos$x2[n], pos$y2[n], width = 1, height = 1, default.units = "native", gp = gpar(fill = fill2[length(fill2)], col = NA, lineend = "butt", linejoin = "mitre"))
		df = data.frame(x = c(df$x, pos$x2[n]), y = c(df$y, pos$y2[n]))
	} else {
		grid.rect(pos$x1[index], pos$y1[index], width = 1, height = 1, default.units = "native", gp = gpar(fill = fill2, col = NA, lineend = "butt", linejoin = "mitre"))
		df = data.frame(x = pos$x1[index], y = pos$y1[index])
	}

	return(invisible(df))

})

normalize_gp = function(name = NULL, value = NULL, length = NULL) {
	if(is.null(value)) value = get.gpar(name)[[name]]
	if(length(value) == 1) value = rep(value, length)
	return(value)
}

# == title
# Add line segments to Hilbert curve
#
# == param
# -object A `HilbertCurve-class` object
# -ir a `IRanges::IRanges` object
# -gp graphical parameters for rectangles
#
# == value
# A data frame which contains coordinates for segments.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
setMethod(f = "hc_segments",
	signature = "HilbertCurve",
	definition = function(object, ir, gp = gpar()) {

	if(object@MODE == "pixel") {
		stop("`hc_segments()` can only be used under 'normal' mode.")
	}

	ir = IRanges(start = zoom(object, start(ir)),
	             end = zoom(object, end(ir)))

	seekViewport(name = paste0("hilbert_curve_", get_plot_index()))

	mtch = as.matrix(findOverlaps(object@BINS, ir))
	r = pintersect(object@BINS[mtch[,1]], ir[mtch[,2]])
	l = width(r) > 1
	r = r[l]
	mtch = mtch[l, , drop = FALSE]

	r1 = object@BINS[mtch[,1]]
	pos = object@POS[mtch[,1], ]

	sr1 = start(r1)
	sr = start(r)
	er1 = end(r1)
	er = end(r)

	x1 = pos$x2 - (pos$x2 - pos$x1)*(er1 - sr)/(er1 - sr1)
	y1 = pos$y2 - (pos$y2 - pos$y1)*(er1 - sr)/(er1 - sr1)
	x2 = pos$x2 - (pos$x2 - pos$x1)*(er1 - er)/(er1 - sr1)
	y2 = pos$y2 - (pos$y2 - pos$y1)*(er1 - er)/(er1 - sr1)

	gp$lineend = "butt"
	grid.segments(x1, y1, x2, y2, default.units = "native", gp = gp)

	df = data.frame(x1 = x1, y1 = y1, x2 = x2, y2 = y2)
	return(invisible(df))
	
})

# == title
# Add text to Hilbert curve
#
# == param
# -object A `HilbertCurve-class` object
# -ir a `IRanges::IRanges` object that contains positions of text. If interval has width larger than 1,
#     the middle point of the interval will be the position of the text.
# -text text corresponding the ``ir``.
# -gp graphical parameters for text
# -... pass to `grid::grid.text`. You can set ``just`` for text here
#
# == value
# A data frame which contains coordinates for text.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
setMethod(f = "hc_text",
	signature = "HilbertCurve",
	definition = function(object, ir, text, gp = gpar(), ...) {

	if(object@MODE == "pixel") {
		stop("`hc_text()` can only be used under 'normal' mode.")
	}

	ir = IRanges(start = zoom(object, (start(ir) + end(ir))/2),
	             end = zoom(object, (start(ir) + end(ir))/2))

	seekViewport(name = paste0("hilbert_curve_", get_plot_index()))

	od = order(ir)
	ir = ir[od]
	text = text[od]

	mtch = as.matrix(findOverlaps(object@BINS, ir))
	l = duplicated(mtch[,2])
	mtch = mtch[!l, , drop = FALSE]

	r1 = object@BINS[mtch[,1]]
	pos = object@POS[mtch[,1], ]
	text = text[mtch[,2]]
	r = pintersect(r1, ir[mtch[,2]])

	sr1 = start(r1)
	sr = start(r)
	er1 = end(r1)
	er = end(r)

	x1 = pos$x2 - (pos$x2 - pos$x1)*(er1 - sr)/(er1 - sr1)
	y1 = pos$y2 - (pos$y2 - pos$y1)*(er1 - sr)/(er1 - sr1)

	grid.text(text, x1, y1, default.units = "native", gp = gp, ...)	
	
	df = data.frame(x = x1, y = y1, text = text, stringsAsFactors = FALSE)
	return(invisible(df))
})

# add arrow at (x2, y2) and x1 -> x2
grid_arrows = function(x1, y1, x2, y2, length = unit(2, "mm"), angle = 15, only.head = FALSE, 
	arrow_gp = gpar(), lines_gp = gpar()) {

	length = convertUnit(length*2, "native", valueOnly = TRUE) - convertUnit(length, "native", valueOnly = TRUE)

	X1 = length*tan(angle/180*pi)
	X2 = rep(0, length = length(X1))
	X3 = -length*tan(angle/180*pi)

	Y1 = rep(-length/2, length = length(X1))
	Y2 = rep(length/2, length = length(X1))
	Y3 = rep(-length/2, length = length(X1))

	
	theta = ifelse(x2 >= x1, pi/2 + atan((y2-y1)/(x2-x1)) + pi, pi/2 + atan((y2-y1)/(x2-x1)))
	mat = matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), 2, 2)
	
	a1 = X1*cos(theta) - Y1*sin(theta) + x2
	a2 = X2*cos(theta) - Y2*sin(theta) + x2
	a3 = X3*cos(theta) - Y3*sin(theta) + x2
	b1 = X1*sin(theta) + Y1*cos(theta) + y2
	b2 = X2*sin(theta) + Y2*cos(theta) + y2
	b3 = X3*sin(theta) + Y3*cos(theta) + y2

	grid.polygon(c(a1, a2, a3), c(b1, b2, b3), default.units = "native", id = rep(seq_along(a1), times = 3), gp = arrow_gp)
	if(!only.head) grid.segments(x1, y1, x2, y2, default.units = "native", gp = lines_gp)
	#grid.points(x2, y2, gp = gpar(col = "red"))
}

# == title
# Add a new layer on the Hilbert curve
#
# == param
# -object A `HilbertCurve-class` object
# -ir a `IRanges::IRanges` object
# -col colors corresponding to each interval in ``ir``
# -mean_mode
#
# == details
# If you want to add more than one layers on the plot, remember to set transparent colors.
#
# == value
# No value is returned
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
setMethod(f = "hc_layer",
	signature = "HilbertCurve",
	definition = function(object, ir, col = "red", mean_mode = c("w0", "absolute", "weighted")) {

	if(object@MODE == "normal") {
		stop("`hc_layer()` can only be used under 'pixel' mode.")
	}

	ir = IRanges(start = zoom(object, start(ir)),
	             end = zoom(object, end(ir)))

	s = start(object@BINS)
	e = end(object@BINS)
	mid = round((s + e)/2)
	window = IRanges(start = c(s[1], mid), end = c(mid, e[length(e)]))
		
	mtch = as.matrix(findOverlaps(window, ir))
	
	# colors correspond to mean value for each window
	mean_mode = match.arg(mean_mode)[1]

	if(length(col) == 1) col = rep(col, length(ir))

	rgb = col2rgb(col, alpha = TRUE)

	r = average_in_window(window, ir, mtch, rgb[1, ], mean_mode, 255)/255
	g = average_in_window(window, ir, mtch, rgb[2, ], mean_mode, 255)/255
	b = average_in_window(window, ir, mtch, rgb[3, ], mean_mode, 255)/255
	alpha = rep(max(rgb[4, ]), length(r))/255
	r[r > 1] = 1
	g[g > 1] = 1
	b[b > 1] = 1
	alpha[alpha > 1] = 1

	col_index = c(object@POS$x1, object@POS$x2[4^object@LEVEL-1]) + 1
	row_index = c(object@POS$y1, object@POS$y2[4^object@LEVEL-1]) + 1
	col_index = col_index[as.numeric(names(r))]
	row_index = row_index[as.numeric(names(r))]

	ind = row_index + (col_index - 1)*nrow(object@RGB$red)
	object@RGB$red[ind] = r*alpha + object@RGB$red[ind] * (1-alpha)
	object@RGB$green[ind] = g*alpha + object@RGB$green[ind] * (1-alpha)
	object@RGB$blue[ind] = b*alpha + object@RGB$blue[ind] * (1-alpha)

	return(invisible(NULL))
})

# == title
# Save Hilbert curve as PNG figure
#
# == param
# -object A `HilbertCurve-class` object
# -file file name
# -grid add grid lines
#
# == value
# No value is returned
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
setMethod(f = "hc_save",
	signature = "HilbertCurve",
	definition = function(object, file = "Rplot.png", grid = 0) {

	red = object@RGB$red
	green = object@RGB$green
	blue = object@RGB$blue
	red = red[seq(nrow(red), 1), , drop = FALSE]
	green = green[seq(nrow(green), 1), , drop = FALSE]
	blue = blue[seq(nrow(blue), 1), , drop = FALSE]
	
	if(grid > 1) {
		n = ncol(red)
		for(i in 1:(grid-1)) {
			red[round(n/grid*i), ] = 0
			green[round(n/grid*i), ] = 0
			blue[round(n/grid*i), ] = 0
			
			red[, round(n/grid*i)] = 0
			green[, round(n/grid*i)] = 0
			blue[, round(n/grid*i)] = 0
		}
	}

	img = array(dim = c(2^object@LEVEL, 2^object@LEVEL, 3))
	img[, , 1] = red
	img[, , 2] = green
	img[, , 3] = blue

	if(!grepl("\\.png", file, ignore.case = TRUE)) {
		file = paste0(file, ".png")
	}

	writePNG(img, target = file)
})








