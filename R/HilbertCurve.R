
.ENV = new.env()

.ENV$I_PLOT = 0

get_plot_index = function() {
	.ENV$I_PLOT
}

increase_plot_index = function() {
	.ENV$I_PLOT = .ENV$I_PLOT + 1
}


# == title
# The HilbertCurve class
#
# == details
# Hilbert curve (https://en.wikipedia.org/wiki/Hilbert_curve ) is a type of space-filling curves
# that folds one-dimensional axis into a two-dimensional space, but still keeps the locality.
# It has advantages to visualize data with long axis with high resolution.
#
# This package aims to provide an easy and flexible way to visualize data through Hilbert curve.
# The implementation and example figures are based on following sources:
#
# - http://mkweb.bcgsc.ca/hilbert/
# - http://corte.si/posts/code/hilbert/portrait/index.html
# - http://bioconductor.org/packages/devel/bioc/html/HilbertVis.html
#
# == Methods
# The `HilbertCurve-class` provides following methods:
#
# - `HilbertCurve`: constructor method;
# - `hc_points,HilbertCurve-method`: add points;
# - `hc_segments,HilbertCurve-method`: add lines;
# - `hc_rect,HilbertCurve-method`: add rectangles;
# - `hc_polygon,HilbertCurve-method`: add polygons;
# - `hc_text,HilbertCurve-method`: add text;
# - `hc_layer,HilbertCurve-method`: add layers, works under "pixel" mode;
# - `hc_png,HilbertCurve-method`: save plot as PNG format, works under "pixel" mode.
#
# == seealso
# The `GenomicHilbertCurve-class` inherits `HilbertCurve-class` and is designed specifically for handling genomic data.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# NULL
#
HilbertCurve = setClass("HilbertCurve",
	slots = list(
		BINS = "IRanges",
		POS = "data.frame",
		RANGE = "numeric",
		ZOOM = "numeric",
		MODE = "character",
		RGB = "environment",
		LEVEL = "integer",
		OFFSET = "numeric",
		start_from = "character"
))

# == title
# Zoom original positions
#
# == param
# -object A `HilbertCurve-class` object.
# -x original positions.
#
# == details
# Internally, position are stored as integer values. To better map the data to the Hilbert curve, 
# the original positions are zoomed
# according to the range and the level of Hilbert curve. E.g. if 
# the curve visualizes data ranging from 1 to 2 but level of the curve is set to 4,
# the positions will be zoomed by ~x2000 so that values like 1.5, 1.555 can be mapped
# to the curve with more accuracy.
#
# The function is used internally.
#
# == value
# A numeric vector which is zoomed positions.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# hc = HilbertCurve(1, 2)
# zoom(hc, 1.5)
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
# -object A `HilbertCurve-class` object.
# -x positions.
#
# == details
# This is a reverse function of `zoom,HilbertCurve-method`.
#
# The function is used internally.
#
# == value
# A numeric vector of original positions
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# hc = HilbertCurve(1, 2)
# z = zoom(hc, 1.5)
# unzoom(hc, z)
#
setMethod(f = "unzoom",
	signature = "HilbertCurve",
	definition = function(object, x) {
	x / object@ZOOM
})

# == title
# Adjust positions
#
# == param
# -object A `HilbertCurve-class` object.
# -x positions.
#
# == details
# Since internally positions are transformed to positive integers, if input positions
# are specified as negative values when initialize the Hilbert curve, a shift will be recorded
# internally and positions are transformed to positive value automatically.
#
# == value
# A positive numeric value
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# hc = HilbertCurve(-100, 100)
# hc_offset(hc, c(-100, -50, 0, 50, 100))
#
setMethod(f = "hc_offset",
	signature = "HilbertCurve",
	definition = function(object, x) {
	x - object@OFFSET
})

# == title
# Initialize a Hilbert curve
#
# == param
# -s position that will be mapped as the start of the Hilbert curve.
# -e position that will be mapped as the end of the Hilbert curve.
# -level iteration level of the Hilbert curve. There will by ``4^level - 1`` segments in the curve.
# -mode "normal" mode is used for low ``level`` value and "pixel" mode is always used for high ``level`` value,
#       so the "normal" mode is always for low-resolution visualization while "pixel" mode is used for high-resolution visualization.
#       See 'details' for explanation.
# -reference whether add reference line on the plot. Only works under 'normal' mode. The reference line
#            is only used for illustrating how the curve folds.
# -reference_gp graphic settings for the reference line. It should be specified by `grid::gpar`.
# -arrow whether add arrows on the reference line. Only works under 'normal' mode.
# -zoom Internally, position are stored as integer values. To better map the data to the Hilbert curve, 
#       the original positions are zoomed
#       according to the range and the level of Hilbert curve. E.g. if 
#       the curve visualizes data ranging from 1 to 2 but level of the curve is set to 4,
#       the positions will be zoomed by ~x2000 so that values like 1.5, 1.555 can be mapped
#       to the curve with more accuracy. You don't need to care the zooming thing, 
#       proper zooming factor is calculated automatically.
# -newpage whether call `grid::grid.newpage` to draw on a new graphic device.
# -background_col background color.
# -title title of the plot.
# -title_gp graphic parameters for the title. It should be specified by `grid::gpar`.
# -start_from which corner on the plot should the curve starts?
# -first_seg the orientation of the first segment
# -legend a `grid::grob` object or a list of `grid::grob` objects. You can construct a `ComplexHeatmap::ColorMapping-class`
#         object and generate a legend, see the Example section.
#
# == details
# This funciton initializes a Hilbert curve with level ``level`` which corresponds 
# to the range between ``s`` and ``e``.
#
# Under 'normal' mode, there is a visible Hilbert curve which plays like a folded axis and
# different low-level graphics can be added afterwards according to the coordinates. 
# It works nice if the level of the Hilbert curve is small (say less than 6). 
#
# When the level is high (e.g. > 10), the whole 2D space will be almost completely filled by the curve and
# it is impossible to add or visualize e.g. points on the curve. In this case, the 'pixel'
# mode visualizes each tiny 'segment' as a pixel and maps values to colors. Internally, the whole plot
# is represented as an RGB matrix and every time a new layer is added to the plot, the RGB matrix
# will be updated according to the color overlay. When all the layers are added, normally a PNG figure is generated
# directly from the RGB matrix. So the Hilbert
# curve with level 11 will generate a PNG figure with 2048x2048 resolution. This is extremely
# useful for visualize genomic data. E.g. If we make a Hilbert curve for human chromosome 1 with
# level 11, then each pixel can represent 60bp (``249250621/2048/2048``) which is of very high resolution.
#
# Under 'pixel' mode, if the current device is an interactive deivce, every time a new layer is added, 
# the image will be add to the interactive device as a rastered image. But still you can use `hc_png,HilbertCurve-method`
# to export the plot as PNG file.
#
# To make it short and clear, under "normal" mode, you use following low-level graphic functions:
# 
# - `hc_points,HilbertCurve-method`
# - `hc_segments,HilbertCurve-method`
# - `hc_rect,HilbertCurve-method`
# - `hc_polygon,HilbertCurve-method`
# - `hc_text,HilbertCurve-method`
#
# And under "pixel" mode, you can use following functions:
#
# - `hc_layer,HilbertCurve-method`
# - `hc_png,HilbertCurve-method`
#
# Notice, ``s`` and ``e`` are not necessarily to be integers, it can be any values (e.g. numeric or even negative values).
#
# == value
# A `HilbertCurve-class` object.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# HilbertCurve(1, 100, reference = TRUE)
# HilbertCurve(1, 100, level = 5, reference = TRUE)
# HilbertCurve(1, 100, title = "title", reference = TRUE)
# HilbertCurve(1, 100, start_from = "topleft", reference = TRUE)
#
# # plot with one legend
# require(ComplexHeatmap)
# cm = ColorMapping(colors = c("red", "blue"), levels = c("a", "b"))
# legend = color_mapping_legend(cm, plot = FALSE, title = "foo")
# hc = HilbertCurve(1, 100, title = "title", legend = legend)
# hc_segments(hc, x1 = 20, x2 = 40)
#
# # plot with more than one legend
# require(circlize)
# cm1 = ColorMapping(colors = c("red", "blue"), levels = c("a", "b"))
# legend1 = color_mapping_legend(cm1, plot = FALSE, title = "foo")
# cm2 = ColorMapping(col_fun = colorRamp2(c(-1, 0, 1), c("green", "white", "red")))
# legend2 = color_mapping_legend(cm2, plot = FALSE, title = "bar")
# hc = HilbertCurve(1, 100, title = "title", legend = list(legend1, legend2))
# hc_segments(hc, x1 = 20, x2 = 40)
#
HilbertCurve = function(s, e, level = 4, mode = c("normal", "pixel"),
	reference = FALSE, reference_gp = gpar(lty = 3, col = "#999999"), 
	arrow = TRUE, zoom = NULL, newpage = TRUE, 
	background_col = "transparent", title = NULL, title_gp = gpar(fontsize = 16), 
	start_from = c("bottomleft", "topleft", "bottomright", "topright"),
	first_seg = c("horizontal", "vertical"), legend = list()) {

	hc = new("HilbertCurve")
	level = as.integer(level)

	if(s > e) {
		stop("`s` must be smaller than `e`.")
	}
	if(s < 0) {
		hc@OFFSET = s
	} else {
		hc@OFFSET = 0
	}
	s = hc_offset(hc, s)
	e = hc_offset(hc, e)

	if(is.null(zoom)) zoom = 10*(4^level-1)/(e - s)
	hc@ZOOM = zoom

	pos = hilbertCurve(level)
	start_from = match.arg(start_from)[1]
	first_seg = match.arg(first_seg)[1]
	
	if(start_from == "bottomleft") {
		if(first_seg == "horizontal") {
			# default, no transformation
		} else if(first_seg == "vertical") {
			pos = rotate(pos, degree = 90)
			pos = flip(pos, direction = "horizontal")
		}
	} else if(start_from == "topleft") {
		if(first_seg == "horizontal") {
			pos = rotate(pos, degree = 180)
			pos = flip(pos, direction = "horizontal")
		} else {
			pos = rotate(pos, degree = 270)
		}
	} else if(start_from == "bottomright") {
		if(first_seg == "horizontal") {
			pos = flip(pos, direction = "horizontal")
		} else {
			pos = rotate(pos, degree = 90)
		}
	} else if(start_from == "topright") {
		if(first_seg == "horizontal") {
			pos = rotate(pos, degree = 180)
		} else {
			pos = rotate(pos, degree = 270)
			pos = flip(pos, direction = "horizontal")
		}
	}

	hc@start_from = start_from
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

	
	if(newpage) grid.newpage()
	increase_plot_index()

	# create a 2x2 layout
	if(length(title) == 0) {
		title_height = unit(0, "mm")
	} else {
		title_height = grobHeight(textGrob(title, gp = title_gp)) + unit(5, "mm")
	}

	if(length(legend) == 0) {
		legend_width = unit(0, "mm")
	} else {
		if(inherits(legend, "grob")) legend = list(legend)
		legend_width = max(do.call("unit.c", lapply(legend, grobWidth))) + unit(4, "mm")
	}

	pushViewport(viewport(layout = grid.layout(nrow = 2, ncol = 2, widths = unit.c(unit(1, "npc") - legend_width, legend_width), 
			                                                       heights = unit.c(title_height, unit(1, "npc") - title_height)),
	                     name = paste0("hilbert_curve_", get_plot_index(), "global")))

	title_gp = validate_gpar(title_gp, default = list(fontsize = 16))
	if(length(title) != 0) {
		pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
		grid.text(title, gp = title_gp)
		upViewport()
	}
	if(length(legend) != 0) {
		gap = unit(2, "mm")
		pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2))
		# draw the list of legend
		legend_height = sum(do.call("unit.c", lapply(legend, grobHeight))) + gap*(length(legend)-1)
		y = unit(0.5, "npc") + legend_height*0.5 
		for(i in seq_along(legend)) {
			pushViewport(viewport(x = unit(2, "mm"), y = y, height = grobHeight(legend[[i]]), width = grobWidth(legend[[i]]), just = c("left", "top")))
			grid.draw(legend[[i]])
			upViewport()
			y = y - gap - grobHeight(legend[[i]])
		}

		upViewport()
	}

	size = unit(1, "snpc")
	pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))

	pushViewport(viewport(name = paste0("hilbert_curve_", get_plot_index()), xscale = c(-0.5, 2^level - 0.5), yscale = c(-0.5, sqrt(n)-0.5), width = size, height = size))
	grid.rect(gp = gpar(fill = background_col, col = NA))

	reference_gp = validate_gpar(reference_gp, default = list(lty = 3, col = "#999999"))
	if(mode == "normal") {
		if(reference) {
			grid.segments(hc@POS$x1, hc@POS$y1, hc@POS$x2, hc@POS$y2, default.units = "native", gp = reference_gp)
			
			# grid.points(hc@POS$x1[1], hc@POS$y1[1], default.units = "native", gp = gpar(col = "#CCCCCC", cex = 0.5))
			# grid.points(hc@POS$x2, hc@POS$y2, default.units = "native", gp = gpar(col = "#CCCCCC", cex = 0.5))

			# grid.text(round(unzoom(hc, start(bins)[1])), hc@POS$x1[1], hc@POS$y1[1], default.units = "native", gp = gpar(col = "#999999", cex = 0.5))
			# grid.text(round(unzoom(hc, end(bins))), hc@POS$x2, hc@POS$y2, default.units = "native", gp = gpar(col = "#999999", cex = 0.5))
		
			if(arrow) grid_arrows(hc@POS$x1, hc@POS$y1, (hc@POS$x1+hc@POS$x2)/2, (hc@POS$y1+hc@POS$y2)/2, only.head = TRUE, arrow_gp = gpar(fill = reference_gp$col, col = NA))
		}
		upViewport(3)
	} else {
		background_col = background_col[1]
		background_col = col2rgb(background_col) / 255
		red = matrix(background_col[1], nrow = 2^level, ncol = 2^level)
		green = matrix(background_col[2], nrow = 2^level, ncol = 2^level)
		blue = matrix(background_col[3], nrow = 2^level, ncol = 2^level)
		hc@RGB$red = red
		hc@RGB$green = green
		hc@RGB$blue = blue

		add_raster(hc@RGB)  # already jump to top vp

	}

	# message(strwrap("If your regions are mostly very small, consider to set 'mean_mode' to 'absolute' when using `hc_points()`, `hc_rect()` and `hc_layer()`."))

	return(invisible(hc))
}

rotate = function(pos, degree = 0) {
	theta = degree/180 * pi
	x = pos[, 1]
	y = pos[, 2]

	center = c(mean(x), mean(y))
	x = x - center[1]
	y = y - center[2]
	new_x = x*cos(theta) - y*sin(theta)
	new_y = x*sin(theta) + y*cos(theta)
	new_x = round(new_x + center[1])
	new_y = round(new_y + center[2])

	return(data.frame(x = new_x, y = new_y))
}

flip = function(pos, direction = c("horizontal", "vertical")) {
	x = pos[, 1]
	y = pos[, 2]

	center = c(mean(x), mean(y))
	new_x = x - center[1]
	new_y = y - center[2]
	if("horizontal" %in% direction) {
		new_x = -new_x
	}
	if("vertical" %in% direction) {
		new_y = -new_y
	}
	new_x = round(new_x + center[1])
	new_y = round(new_y + center[2])

	return(data.frame(x = new_x, y = new_y))
}

# == title
# Print the HilbertCurve object
# 
# == param
# -object A `HilbertCurve-class` object.
#
# == value
# No value is returned.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# HilbertCurve(1, 100)
#
setMethod(f = "show",
	signature = "HilbertCurve",
	definition = function(object) {
	cat("A Hilbert curve with level", object@LEVEL, "\n")
	cat("mode:", object@MODE, "\n")
})

# == title
# Level of the Hilbert curve
#
# == param
# -object A `HilbertCurve-class` object.
#
# == value
# The level of the Hilbert curve.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# hc = HilbertCurve(1, 100)
# hc_level(hc)
#
# hc = HilbertCurve(1, 100, level = 5)
# hc_level(hc)
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
# -object A `HilbertCurve-class` object.
# -ir an `IRanges::IRanges` object which specifies the input intervals.
# -x1 if start positions are not integers, they can be set by ``x1``.
# -x2 if end positions are not integers, they can be set by ``x2``.
# -np number of points (a circle or a square, ...) that are put in a segment. ``np`` controls
#     the mode of how to add the points to the curve. See 'Details' section.
# -size size of the points. It should be a `grid::unit` object. Only works if ``np <= 1``
# -pch shape of points, used for points if ``np <= 1``.
# -gp graphic parameters for points. It should be specified by `grid::gpar`.
# -mean_mode when ``np >= 2``, each segment on the curve is split into ``np`` windows and each window
#           actually represents an small interval in the axis. When overlapping input intervals to the windows
#           on the curve and when the window can not completely cover the input intervals, some averaging method should
#           be applied to get a more accurate estimation for the value in the window. Here the HilbertCurve package
#           provides three modes: "w0", "weighted" and "absolute" which calculate the mean value in the window
#           with respect to different scenarios. See 'Details' section and the vignette for more informative explanation.
# -shape shape of points, used for points if ``np >= 2``. Possible values are "circle", "square", "triangle", "hexagon", "star".
#
# == details
# If ``np`` is set to 1 or ``NULL``, points will be added in the middle for each interval in ``ir`` (or ``x1``, ``x2``).
#
# If ``np`` is set to a value larger or equal to 2, every segment on the curve is split by ``np`` points (e.g. circles).
# In this case, each point actually represent a window on the curve and when the window is not fully covered by
# the input intervals, there are three different metrics to average the values in the window.
#
# Following illustrates different settings for ``mean_mode``:
#
#        100    80     60    values in ir
#     ++++++   +++   +++++   ir
#       ================     window (width = 16)
#         4     3     3      overlap
#
#     absolute: (100 + 80 + 60)/3
#     weighted: (100*4 + 80*3 + 60*3)/(4 + 3 + 3)
#     w0:       (100*4 + 80*3 + 60*3 + 0*6)/16
#
# So which mode to use depends on specific scenario. If the background is not of interest, ``absolute`` and ``weighted``
# modes may be proper and if the value also needs to be averaged with background, ``w0`` is the proper choice. Section "Averaging models"
# in the vignette gives a more detailed explanation for this argument.
#
# If ``np >= 2``, the value of ``np`` also controls the size of points.
# 
# Graphic parameters is always represented as numeric values (e.g. colors can be converted into numeric RGB values) 
# and they will be averaged according to above rules.
#
# Internally, it will depatch to `hc_normal_points,HilbertCurve-method` or `hc_segmented_points,HilbertCurve-method`
# depending on the value of ``np``.
#
# == value
# A data frame which contains coordinates (in the 2D space) of points.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# hc = HilbertCurve(1, 100, level = 4, reference = TRUE)
#
# x = sort(sample(100, 20))
# s = x[1:10*2 - 1]
# e = x[1:10*2]
# require(IRanges)
# ir = IRanges(s, e)
#
# hc_points(hc, ir)
#
# hc = HilbertCurve(1, 100, level = 4, reference = TRUE)
# hc_points(hc, x1 = c(1.5, 50.5), x2 = c(10.5, 60.5))
#
# require(circlize)
# value = runif(length(ir))
# col_fun = colorRamp2(range(value), c("white", "red"))
# hc = HilbertCurve(1, 100, level = 4, reference = TRUE)
# hc_points(hc, ir, np = 3, shape = "star", gp = gpar(fill = col_fun(value)))
# 
# hc = HilbertCurve(1, 100, level = 4, reference = TRUE)
# hc_points(hc, ir, np = 0)
#
# hc = HilbertCurve(1, 100, level = 4, reference = TRUE)
# hc_points(hc, np = 0, x1 = c(1.5, 50.5), x2 = c(10.5, 60.5))
# hc_points(hc, np = 0, x1 = 70.5, gp = gpar(col = "red"))
#
setMethod(f = "hc_points",
	signature = "HilbertCurve",
	definition = function(object, ir, x1 = NULL, x2 = x1, 
	np = max(c(2, 10 - hc_level(object))), size = unit(1, "char"), 
	pch = 1, gp = gpar(), mean_mode = c("w0", "absolute", "weighted"),
	shape = "circle") {

	if(object@MODE == "pixel") {
		stop("`hc_points()` can only be used under 'normal' mode.")
	}

	if(is.null(np)) np = 1

	if(np >= 2) {
		if(missing(ir)) {
			hc_segmented_points(object, x1 = x1, x2 = x2, gp = gp, np = np, mean_mode = mean_mode, shape = shape)
		} else {
			hc_segmented_points(object, ir, gp = gp, np = np, mean_mode = mean_mode, shape = shape)
		}
	} else {
		if(missing(ir)) {
			hc_normal_points(object, x1 = x1, x2 = x2, gp = gp, pch = pch, size = size)
		} else {
			hc_normal_points(object, ir, gp = gp, pch = pch, size = size)
		}
	}

})

# == title
# Add points to the Hilbert curve
#
# == param
# -object A `HilbertCurve-class` object.
# -ir an `IRanges::IRanges` object which specifies the input intervals.
# -x1 if start positions are not integers, they can be set by ``x1``.
# -x2 if end positions are not integers, they can be set by ``x2``.
# -size size of the points. It should be a `grid::unit` object, pass to `grid::grid.points`.
# -pch shape of points, pass to `grid::grid.points`.
# -gp graphic parameters for points. It should be specified by `grid::gpar`.
#
# == details
# Points are added at the middle of the intervals in ``ir`` (or ``x1`` and ``x2``),
# so there is only one point for each interval.
#
# This function is used internally. Please use `hc_points,HilbertCurve-method` instead.
#
# == value
# A data frame which contains coordinates (in the 2D space) of points.
#
# == seealso
# `hc_points,HilbertCurve-method`
# 
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# # see documentation of hc_points
# NULL
#
setMethod(f = "hc_normal_points",
	signature = "HilbertCurve",
	definition = function(object, ir, x1 = NULL, x2 = x1, gp = gpar(), 
	pch = 1, size = unit(1, "char")) {

	if(object@MODE == "pixel") {
        stop("`hc_points()` can only be used under 'normal' mode.")
    }

    if(!missing(ir)) {
    	if(!inherits(ir, "IRanges")) {
    		stop("It seems you want to specify positions by one or two numeric vectors, specify them with argument names. (`x1 = ..., x2 = ...`)")
    	}
    }

	if(missing(ir) || is.null(ir)) {
		if(is.null(x1) ) {
			stop("You should either specify `ir`, or `x1` or `x1` and `x2`.")
		} else if(!is.null(x1) && !is.null(x2)) {
			x1 = hc_offset(object, x1)
			x2 = hc_offset(object, x2)
			ir = IRanges(start = round(zoom(object, (x1+x2)/2)),
				         end = round(zoom(object, (x1+x2)/2)))
		}
	} else {
		ir = IRanges(hc_offset(object, start(ir)),
			         hc_offset(object, end(ir)))
		ir = IRanges(start = zoom(object, (start(ir) + end(ir))/2),
		             end = zoom(object, (start(ir) + end(ir))/2))
	}

	od = order(ir)
	ir = ir[od]

	mtch = as.matrix(findOverlaps(object@BINS, ir))
	if(nrow(mtch) == 0) {
		return(invisible(NULL))
	}
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

	seekViewport(name = paste0("hilbert_curve_", get_plot_index()))

	grid.points(x1, y1, default.units = "native", gp = gp, pch = pch, size = size)	

	seekViewport(name = paste0("hilbert_curve_", get_plot_index(), "global"))
	upViewport()
	
	df = data.frame(x = x1, y = y1, stringsAsFactors = FALSE)
	return(invisible(df))
})

# == title
# Add points to the Hilbert curve
#
# == param
# -object A `HilbertCurve-class` object.
# -ir an `IRanges::IRanges` object which specifies the input intervals.
# -x1 if start positions are not integers, they can be set by ``x1``.
# -x2 if end positions are not integers, they can be set by ``x2``.
# -np number of points (a circle or a square, ...) that are put in a segment.
# -gp graphic parameters for points. It should be specified by `grid::gpar`.
#     The size of the points can be set here because the size of points are determined by ``np`` argument.
# -mean_mode when a segment in the curve overlaps with intervals in ``ir``, how to calculate 
#       the mean values for this segment. See explanation in `hc_points`.
# -shape shape of points. Possible values are "circle", "square", "triangle", "hexagon", "star".
#
# == details
# Every segment on the curve is split by ``np`` points.
#
# This function is used internally, please use `hc_points,HilbertCurve-method` directly.
#
# == value
# A data frame which contains coordinates (in the 2D space) of points.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# # see documentation of hc_points
# NULL
#
setMethod(f = "hc_segmented_points",
    signature = "HilbertCurve",
    definition = function(object, ir, x1 = NULL, x2 = NULL, gp = gpar(), 
    np = max(c(2, 10 - hc_level(object))),
    mean_mode = c("w0", "absolute", "weighted"), 
    shape = "circle") {

    if(object@MODE == "pixel") {
        stop("`hc_points()` can only be used under 'normal' mode.")
    }

    if(!missing(ir)) {
    	if(!inherits(ir, "IRanges")) {
    		stop("It seems you want to specify positions by one or two numeric vectors, specify them with argument names. (`x1 = ..., x2 = ...`)")
    	}
    }

    if(missing(ir) || is.null(ir)) {
        if(is.null(x1) || is.null(x2)) {
            stop("You should either specify `ir`, or `x1` and `x2`.")
        } else {
            x1 = hc_offset(object, x1)
            x2 = hc_offset(object, x2)
            ir = IRanges(start = round(zoom(object, x1)),
                         end = round(zoom(object, x2)))
        }
    } else {
        ir = IRanges(hc_offset(object, start(ir)),
                     hc_offset(object, end(ir)))
        ir = IRanges(start = zoom(object, start(ir)),
                     end = zoom(object, end(ir)))
    }
 
    if(length(shape) == 1) {
        shape = rep(shape, length(ir))
    }

    bin_width = width(object@BINS[1])
    bin_radius = bin_width/(2*(np-1))
    breaks = seq(object@RANGE[1] - bin_radius, object@RANGE[2] + bin_radius, length = (4^object@LEVEL)*(np-1)+1)
    radius = 1/(2*(np-1))
    s = breaks[-length(breaks)]
    e = breaks[-1]
    s[1] = s[1] + bin_radius
    e[length(e)] = e[length(e)] - bin_radius
    window = IRanges(start = round(s), end = round(e))

    mid = apply(object@POS, 1, function(v) {
        x = seq(v[1], v[3], length = np)
        y = seq(v[2], v[4], length = np)
        data.frame(x = x[-1], y = y[-1])
    })
    mid = do.call("rbind", mid)
    mid = rbind(c(object@POS[1, 1], object@POS[1, 2]), mid)

    mtch = as.matrix(findOverlaps(window, ir))
    if(nrow(mtch) == 0) {
    	return(invisible(NULL))
    }
    
    # colors correspond to mean value for each window
    mean_mode = match.arg(mean_mode)[1]

    gp = validate_gpar(gp, default = list(col = "black", fill = "transparent"))
    if(length(gp$col) == 1) gp$col = rep(gp$col, length(ir))
    if(length(gp$fill) == 1) gp$fill = rep(gp$fill, length(ir))
    fill = gp$fill
    rgb_fill = col2rgb(fill, alpha = TRUE)
    col = gp$col
    rgb_col = col2rgb(col, alpha = TRUE)

    rgb_mat = average_in_window(window, ir, mtch, list(rgb_fill[1, ], rgb_fill[2, ], rgb_fill[3, ], rgb_col[1, ], rgb_col[2, ], rgb_col[3, ]), mean_mode, 255)

    alpha_fill = rep(max(rgb_fill[4, ]), nrow(rgb_mat))
    fill2 = rgb(red = rgb_mat[,1], green = rgb_mat[,2], blue = rgb_mat[,3], alpha = alpha_fill, maxColorValue = 255)
	alpha_col = rep(max(rgb_col[4, ]), nrow(rgb_mat))
    col2 = rgb(red = rgb_mat[,4], green = rgb_mat[,5], blue = rgb_mat[,6], alpha = alpha_col, maxColorValue = 255)

    index = as.numeric(rownames(rgb_mat))
    shape = tapply(shape[mtch[,2]], mtch[, 1], function(x) x[1])
    mid = mid[unique(mtch[, 1]), , drop = FALSE]

    seekViewport(name = paste0("hilbert_curve_", get_plot_index()))
    for(s in unique(shape)) {
        l = shape == s
        x = mid[l, 1]
        y = mid[l, 2]
        if(s == "circle") {
            grid.circle(x, y, r = radius, default.units = "native", gp = gpar(fill = fill2[l], col = col2[l]))
        } else if(s == "square") {
            grid.rect(x, y, width = 2*radius, height = 2*radius, default.units = "native", gp = gpar(fill = fill2[l], col = col2[l]))
        } else if(s == "triangle") {
            grid.polygon(c(x, x-2*radius*tan(pi/6), x+2*radius*tan(pi/6)), 
                c(y+radius, y-radius, y-radius), id = rep(seq_along(x), 3), default.units = "native", gp = gpar(fill = fill2[l], col = col2[l]))
        } else if(s == "hexagon") {
            grid.polygon(rep(cos(0:5 * pi/3)*radius, length(x)) + rep(x, each = 6),
                         rep(sin(0:5 * pi/3)*radius, length(y)) + rep(y, each = 6),
                         id = rep(seq_along(x), each = 6),
                         default.units = "native", gp = gpar(fill = fill2[l], col = col2[l]))
        } else if(s == "star") {
            r2 = cos(pi/10)*tan(pi/20)/(tan(pi/20) + tan(pi/10)) * radius
            grid.polygon(rep(sin(0:9*pi/5 + pi*0)*rep(c(radius,r2), 5), length(x)) + rep(x, each = 10),
                         rep(cos(0:9*pi/5 + pi*0)*rep(c(radius,r2), 5), length(y)) + rep(y, each = 10),
                         id = rep(seq_along(x), each = 10),
                         default.units = "native", gp = gpar(fill = fill2[l], col = col2[l]))
        }
    }

    seekViewport(name = paste0("hilbert_curve_", get_plot_index(), "global"))
    upViewport()

    df = cbind(mid, radius = rep(radius, nrow(mid)))
    return(invisible(df))
})

# after zooming, ir may overlap
# it returns a vector having same length as window
# window = IRanges(1, 10)
# ir = IRanges(c(1, 4, 7), c(2, 5, 8))
# mtch = as.matrix(findOverlaps(window, ir))
average_in_window = function(window, ir, mtch, v, mean_mode, empty_v = 0) {

	if(!is.list(v)) {
		v = list(v)
	}

	if(length(empty_v) == 1) {
		empty_v = rep(empty_v, length(v))
	}

	v = do.call("cbind", v)

 	intersect = pintersect(window[mtch[,1]], ir[mtch[,2]])
 	v = v[mtch[,2], , drop = FALSE]
 	n = nrow(v)

 	ind_list = as.list(tapply(seq_len(n), mtch[, 1], function(x) x))

	if(mean_mode == "w0") {
		w = width(intersect)

		ir2 = reduce(ir, min.gapwidth = 0)
		mtch2 = as.matrix(findOverlaps(window, ir2))
		intersect2 = pintersect(window[mtch2[, 1]], ir2[mtch2[, 2]])

		width_intersect = tapply(width(intersect2), mtch2[, 1], sum)
		ind = unique(mtch2[, 1])
		width_setdiff = width(window[ind]) - width_intersect

		w2 = width(window[ind])
		u = matrix(nrow = length(ind_list), ncol = ncol(v))
		rownames(u) = names(ind_list)
		for(i in seq_along(ind_list)) {
			ind = ind_list[[i]]
			x = colSums(v[ind, , drop = FALSE]*w[ind])/sum(w[ind])
			u[i, ] = (x*width_intersect[i] + empty_v*width_setdiff[i])/w2[i]
		}

	} else if(mean_mode == "absolute") {
		u = matrix(nrow = length(ind_list), ncol = ncol(v))
		rownames(u) = names(ind_list)
		for(i in seq_along(ind_list)) {
			u[i, ] = colMeans(v[ind_list[[i]], , drop = FALSE])
		}
		
	} else {
		w = width(intersect)
		u = matrix(nrow = length(ind_list), ncol = ncol(v))
		rownames(u) = names(ind_list)
		for(i in seq_along(ind_list)) {
			ind = ind_list[[i]]
			u[i, ] = colSums(v[ind, , drop = FALSE]*w[ind])/sum(w[ind])
		}
	}

	return(u)
}

# == title
# Add rectangles on Hilbert curve
#
# == param
# -object A `HilbertCurve-class` object.
# -ir an `IRanges::IRanges` object which specifies the input intervals.
# -x1 if start positions are not integers, they can be set by ``x1``.
# -x2 if end positions are not integers, they can be set by ``x2``.
# -gp graphic parameters for rectangles. It should be specified by `grid::gpar`. Note you cannot set ``linejoin`` and ``lineend``.
# -mean_mode when a segment in the curve can not be overlapped with intervals in ``ir``, how to calculate 
#       the mean values for this segment. See explanation in `hc_points,HilbertCurve-method`.
#
# == details
# Rectangles are put if a segment in the Hilbert curve overlaps with the input intervals. 
# You cannot set the width or height of the rectangles. It is always fixed (actually it is a square).
#
# It can be thought as the low-resolution version of `hc_layer,HilbertCurve-method`.
#
# == value
# A data frame which contains coordinates (in the 2D space) of rectangles.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# hc = HilbertCurve(1, 100, level = 4, reference = TRUE)
#
# x = sort(sample(100, 20))
# s = x[1:10*2 - 1]
# e = x[1:10*2]
# require(IRanges)
# ir = IRanges(s, e)
# hc_rect(hc, ir)
#
setMethod(f = "hc_rect",
	signature = "HilbertCurve",
	definition = function(object, ir, x1 = NULL, x2 = NULL, 
	gp = gpar(fill = "red"), 
	mean_mode = c("w0", "absolute", "weighted")) {

	if(object@MODE == "pixel") {
		stop("`hc_rect()` can only be used under 'normal' mode.")
	}

	if(!missing(ir)) {
    	if(!inherits(ir, "IRanges")) {
    		stop("It seems you want to specify positions by one or two numeric vectors, specify them with argument names. (`x1 = ..., x2 = ...`)")
    	}
    }

	if(missing(ir) || is.null(ir)) {
		if(is.null(x1) || is.null(x2)) {
			stop("You should either specify `ir`, or `x1` and `x2`.")
		} else {
			x1 = hc_offset(object, x1)
			x2 = hc_offset(object, x2)
			ir = IRanges(start = round(zoom(object, x1)),
				         end = round(zoom(object, x2)))
		}
	} else {
		ir = IRanges(hc_offset(object, start(ir)),
			         hc_offset(object, end(ir)))
		ir = IRanges(start = zoom(object, start(ir)),
	                 end = zoom(object, end(ir)))
	}

	s = start(object@BINS)
	e = end(object@BINS)
	mid = round((s + e)/2)
	window = IRanges(start = c(s[1], mid), end = c(mid, e[length(e)]))
		
	mtch = as.matrix(findOverlaps(window, ir))
	if(nrow(mtch) == 0) {
		return(invisible(NULL))
	}
	
	# colors correspond to mean value for each window
	mean_mode = match.arg(mean_mode)[1]

	gp = validate_gpar(gp, default = list(fill = "red"))
	if(length(gp$fill) == 1) gp$fill = rep(gp$fill, length(ir))
	fill = gp$fill

	rgb = col2rgb(fill, alpha = TRUE)
	rgb_mat = average_in_window(window, ir, mtch, list(rgb[1, ], rgb[2, ], rgb[3, ]), mean_mode, 255)
	r = rgb_mat[, 1]
	g = rgb_mat[, 2]
	b = rgb_mat[, 3]
	alpha = rep(max(rgb[4, ]), length(r))

	index = as.numeric(rownames(rgb_mat))
	fill2 = rgb(red = r, green = g, blue = b, alpha = alpha, maxColorValue = 255)

	pos = object@POS
	n = 4^object@LEVEL

	seekViewport(name = paste0("hilbert_curve_", get_plot_index()))

	if(n %in% index) {
		ind = setdiff(index, n)
		grid.rect(pos$x1[ind], pos$y1[ind], width = 1, height = 1, default.units = "native", gp = gpar(fill = fill2[-length(fill2)], col = NA, lineend = "butt", linejoin = "mitre"))
		df = data.frame(x = pos$x1[ind], y = pos$y1[ind])

		grid.rect(pos$x2[n-1], pos$y2[n-1], width = 1, height = 1, default.units = "native", gp = gpar(fill = fill2[length(fill2)], col = NA, lineend = "butt", linejoin = "mitre"))
		df = data.frame(x = c(df$x, pos$x2[n]), y = c(df$y, pos$y2[n]))
	} else {
		grid.rect(pos$x1[index], pos$y1[index], width = 1, height = 1, default.units = "native", gp = gpar(fill = fill2, col = NA, lineend = "butt", linejoin = "mitre"))
		df = data.frame(x = pos$x1[index], y = pos$y1[index])
	}

	seekViewport(name = paste0("hilbert_curve_", get_plot_index(), "global"))
	upViewport()

	return(invisible(df))

})


# == title
# Add line segments to Hilbert curve
#
# == param
# -object A `HilbertCurve-class` object.
# -ir an `IRanges::IRanges` object which specifies the input intervals.
# -x1 if start positions are not integers, they can be set by ``x1``.
# -x2 if end positions are not integers, they can be set by ``x2``.
# -gp graphic parameters for lines. It should be specified by `grid::gpar`. Note you cannot set ``linejoin`` and ``lineend``.
#
# == value
# A data frame which contains coordinates (in the 2D space) of segments.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# hc = HilbertCurve(1, 100, level = 4, reference = TRUE)
#
# x = sort(sample(100, 20))
# s = x[1:10*2 - 1]
# e = x[1:10*2]
# require(IRanges)
# ir = IRanges(s, e)
#
# hc_segments(hc, ir)
#
setMethod(f = "hc_segments",
	signature = "HilbertCurve",
	definition = function(object, ir, x1 = NULL, x2 = NULL, 
	gp = gpar(lty = 1, lwd = 1, col = 1)) {

	if(object@MODE == "pixel") {
		stop("`hc_segments()` can only be used under 'normal' mode.")
	}

	if(!missing(ir)) {
    	if(!inherits(ir, "IRanges")) {
    		stop("It seems you want to specify positions by one or two numeric vectors, specify them with argument names. (`x1 = ..., x2 = ...`)")
    	}
    }

	if(missing(ir) || is.null(ir)) {
		if(is.null(x1) || is.null(x2)) {
			stop("You should either specify `ir`, or `x1` and `x2`.")
		} else {
			x1 = hc_offset(object, x1)
			x2 = hc_offset(object, x2)
			ir = IRanges(start = round(zoom(object, x1)),
				         end = round(zoom(object, x2)))
		}
	} else {
		ir = IRanges(hc_offset(object, start(ir)),
			         hc_offset(object, end(ir)))
		ir = IRanges(start = zoom(object, start(ir)),
	                 end = zoom(object, end(ir)))
	}

	gp = validate_gpar(gp, default = list(lty = 1, lwd = 1, col = 1))
	if(length(gp$lty) == 1) gp$lty = rep(gp$lty, length(ir))
	if(length(gp$lwd) == 1) gp$lwd = rep(gp$lwd, length(ir))
	if(length(gp$col) == 1) gp$col = rep(gp$col, length(ir))

	mtch = as.matrix(findOverlaps(object@BINS, ir))
	if(nrow(mtch) == 0) {
		return(invisible(NULL))
	}
	r = pintersect(object@BINS[mtch[,1]], ir[mtch[,2]])
	# l = width(r) > 1
	# r = r[l]
	# mtch = mtch[l, , drop = FALSE]

	seekViewport(name = paste0("hilbert_curve_", get_plot_index()))
	df = tapply(seq_len(nrow(mtch)), mtch[, 2], function(ind) {
		i1 = mtch[ind, 1]  ## index for BINS
		i2 = mtch[ind[1], 2]  ## index for ir

		r = r[ind]  # intersected part
		r1 = object@BINS[i1]  # BIN part
		pos = object@POS[i1, ,drop = FALSE]  # pos on  the 2D space

		sr1 = start(r1)
		sr = start(r)
		er1 = end(r1)
		er = end(r)

		#er[er == er1] = er[er == er1] + 1

		x1 = pos$x2 - (pos$x2 - pos$x1)*(er1 - sr)/(er1 - sr1)
		y1 = pos$y2 - (pos$y2 - pos$y1)*(er1 - sr)/(er1 - sr1)
		x2 = pos$x2 - (pos$x2 - pos$x1)*(er1 - er)/(er1 - sr1)
		y2 = pos$y2 - (pos$y2 - pos$y1)*(er1 - er)/(er1 - sr1)

		df = data.frame(x = c(x1, x2[length(x2)]), y = c(y1, y2[length(y2)]))
		grid.lines(df$x, df$y, default.units = "native", 
			gp = gpar(lwd = gp$lwd[i2], lty = gp$lty[i2], col = gp$col[i2], lineend ="butt", linejoin = "mitre"))
		
		df$which = rep(i2, nrow(df))

		return(df)
	})

	# gp$lty = gp$lty[mtch[, 2]]
	# gp$lwd = gp$lwd[mtch[, 2]]
	# gp$col = gp$col[mtch[, 2]]

	# r1 = object@BINS[mtch[,1]]
	# pos = object@POS[mtch[,1], ]

	# sr1 = start(r1)
	# sr = start(r)
	# er1 = end(r1)
	# er = end(r)

	# er[er == er1] = er[er == er1] + 1

	# x1 = pos$x2 - (pos$x2 - pos$x1)*(er1 - sr)/(er1 - sr1)
	# y1 = pos$y2 - (pos$y2 - pos$y1)*(er1 - sr)/(er1 - sr1)
	# x2 = pos$x2 - (pos$x2 - pos$x1)*(er1 - er)/(er1 - sr1)
	# y2 = pos$y2 - (pos$y2 - pos$y1)*(er1 - er)/(er1 - sr1)

	# gp$lineend = "butt"
	# grid.segments(x1, y1, x2, y2, default.units = "native", gp = gp)

	# df = data.frame(x1 = x1, y1 = y1, x2 = x2, y2 = y2)

	seekViewport(name = paste0("hilbert_curve_", get_plot_index(), "global"))
	upViewport()

	return(invisible(df))
	
})

# == title
# Add text to Hilbert curve
#
# == param
# -object A `HilbertCurve-class` object.
# -ir an `IRanges::IRanges` object that contains positions which correspond to text.
#     The middle point of the interval will be the position of the text.
# -labels text corresponding to intervals in ``ir``.
# -x1 if start positions are not integers, they can be set by ``x1``.
# -x2 if end positions are not integers, they can be set by ``x2``.
# -gp graphic parameters for text. It should be specified by `grid::gpar`.
# -... pass to `grid::grid.text`. E.g. you can set text justification by ``just`` here.
#
# == details
# The text is added correspoding to the middle of each interval in ``ir``.
#
# == value
# A data frame which contains coordinates (in the 2D space) of text.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# hc = HilbertCurve(1, 100, level = 4, reference = TRUE)
#
# x = sort(sample(100, 20))
# s = x[1:10*2 - 1]
# e = x[1:10*2]
# require(IRanges)
# ir = IRanges(s, e)
#
# labels = sample(letters, length(ir), replace = TRUE)
# hc_text(hc, ir, labels = labels)
#
setMethod(f = "hc_text",
	signature = "HilbertCurve",
	definition = function(object, ir, labels, x1 = NULL, x2 = x1, gp = gpar(), ...) {

	if(object@MODE == "pixel") {
		stop("`hc_text()` can only be used under 'normal' mode.")
	}

	if(!missing(ir)) {
    	if(!inherits(ir, "IRanges")) {
    		stop("It seems you want to specify positions by one or two numeric vectors, specify them with argument names. (`x1 = ..., x2 = ...`)")
    	}
    }

	if(missing(ir) || is.null(ir)) {
		if(is.null(x1) ) {
			stop("You should either specify `ir`, or `x1` or `x1` and `x2`.")
		} else if(!is.null(x1) && !is.null(x2)) {
			x1 = hc_offset(object, x1)
			x2 = hc_offset(object, x2)
			ir = IRanges(start = round(zoom(object, (x1+x2)/2)),
				         end = round(zoom(object, (x1+x2)/2)))
		}
	} else {
		ir = IRanges(hc_offset(object, start(ir)),
			         hc_offset(object, end(ir)))
		ir = IRanges(start = zoom(object, (start(ir) + end(ir))/2),
		             end = zoom(object, (start(ir) + end(ir))/2))
	}

	od = order(ir)
	ir = ir[od]
	labels = labels[od]

	mtch = as.matrix(findOverlaps(object@BINS, ir))
	if(nrow(mtch) == 0) {
		return(invisible(NULL))
	}

	l = duplicated(mtch[,2])
	mtch = mtch[!l, , drop = FALSE]

	r1 = object@BINS[mtch[,1]]
	pos = object@POS[mtch[,1], ]
	labels = labels[mtch[,2]]
	r = pintersect(r1, ir[mtch[,2]])

	sr1 = start(r1)
	sr = start(r)
	er1 = end(r1)
	er = end(r)

	x1 = pos$x2 - (pos$x2 - pos$x1)*(er1 - sr)/(er1 - sr1)
	y1 = pos$y2 - (pos$y2 - pos$y1)*(er1 - sr)/(er1 - sr1)

	seekViewport(name = paste0("hilbert_curve_", get_plot_index()))

	grid.text(labels, x1, y1, default.units = "native", gp = gp, ...)	
	
	df = data.frame(x = x1, y = y1, labels = labels, stringsAsFactors = FALSE)

	seekViewport(name = paste0("hilbert_curve_", get_plot_index(), "global"))
	upViewport()
	
	return(invisible(df))
})

# == title
# Add text to the center of the block
#
# == param
# -object A `HilbertCurve-class` object.
# -ir an `IRanges::IRanges` object that contains positions which correspond to text.
#     The middle points of the intervals will be the positions of the text.
# -labels text corresponding to intervals in ``ir``.
# -x1 if start positions are not integers, they can be set by ``x1``.
# -x2 if end positions are not integers, they can be set by ``x2``.
# -gp graphic parameters for text. It should be specified by `grid::gpar`.
# -... pass to `grid::grid.text`. E.g. you can set text justification by ``just`` here.
#
# == details
# If the interval is long enough that it represents as a block in the 2D space, the corresponding
# label is put approximately at center (or at least inside) of the block.
#
# It is quite experimental and only used internally.
#
# == value
# ``NULL``
#
# == seealso
# It is basically used in `hc_map,GenomicHilbertCurve-method` to put chromosome names in the center
# of chromosomes.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# hc = HilbertCurve(1, 10)
# hc_rect(hc, x1 = c(1, 3, 7), x2 = c(3, 7, 10), gp = gpar(fill = 2:5))
# hc_centered_text(hc, x1 = 1, x2 = 3, labels = "A")
# hc_centered_text(hc, x1 = 3, x2 = 7, labels = "B")
# hc_centered_text(hc, x1 = 7, x2 = 10, labels = "C")
setMethod(f = "hc_centered_text",
	signature = "HilbertCurve",
	definition = function(object, ir, labels, x1 = NULL, x2 = NULL, gp = gpar(), ...) {

	if(object@MODE == "pixel") {
		stop("`hc_text()` can only be used under 'normal' mode.")
	}

	if(!missing(ir)) {
    	if(!inherits(ir, "IRanges")) {
    		stop("It seems you want to specify positions by one or two numeric vectors, specify them with argument names. (`x1 = ..., x2 = ...`)")
    	}
    }

	if(missing(ir) || is.null(ir)) {
		if(is.null(x1) || is.null(x2)) {
			stop("You should either specify `ir`, or `x1` and `x2`.")
		} else {
			x1 = hc_offset(object, x1)
			x2 = hc_offset(object, x2)
			ir = IRanges(start = round(zoom(object, x1)),
				         end = round(zoom(object, x2)))
		}
	} else {
		ir = IRanges(hc_offset(object, start(ir)),
			         hc_offset(object, end(ir)))
		ir = IRanges(start = zoom(object, start(ir)),
	                 end = zoom(object, end(ir)))
	}

	mtch = as.matrix(findOverlaps(object@BINS, ir))
	if(nrow(mtch) == 0) {
		return(invisible(NULL))
	}

	seekViewport(name = paste0("hilbert_curve_", get_plot_index()))

	lapply(seq_along(unique(mtch[, 2])), function(i) {
		ind = mtch[mtch[, 2] == i, 1]
		pos = object@POS[ind, ]
		pos = cbind(pos[, 1:2], pos[nrow(pos), 3:4])
		h1 = hist(pos[, 1], plot = FALSE)
		i1 = which(h1$counts >= quantile(h1$counts, 0.6))
		ind_ir = reduce(IRanges(i1, i1))
		selected = ind_ir[which.max(width(ind_ir))[1]]
		i1 = start(selected)
		i2 = end(selected)
		x = mean(h1$mids[i1:i2])

		h1 = hist(pos[, 2], plot = FALSE)
		i1 = which(h1$counts >= quantile(h1$counts, 0.6))
		ind_ir = reduce(IRanges(i1, i1))
		selected = ind_ir[which.max(width(ind_ir))[1]]
		i1 = start(selected)
		i2 = end(selected)
		y = mean(h1$mids[i1:i2])
		
		grid.text(labels[i], x, y, default.units = "native", gp = gp, ...)
	})

	seekViewport(name = paste0("hilbert_curve_", get_plot_index(), "global"))
	upViewport()
	
	return(invisible(NULL))
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
# Add a new layer to the Hilbert curve
#
# == param
# -object A `HilbertCurve-class` object.
# -ir an `IRanges::IRanges` object which specifies the input intervals.
# -x1 if start positions are not integers, they can be set by ``x1``.
# -x2 if end positions are not integers, they can be set by ``x2``.
# -col a scalar or a vector of colors which correspond to intervals in ``ir`` (or ``x1`` and ``x2``).
# -border a scalar or a vector of colors for the borders of intervals. Set it to ``NA`` if borders are suppressed.
# -mean_mode Under 'pixel' mode, each pixel represents a small window. This argument provides methods
#            to summarize value for the small window if the input intervals can not completely overlap with the window. 
#            See explanation in `hc_points,HilbertCurve-method`.
# -grid_line whether add grid lines to show blocks of the Hilber curve. 
#        It should be an integer number and there will be ``2^(grid_line-1)-1``
#        horizontal and vertical grid lines.
# -grid_line_col color for the grid lines
# -overlay a self-defined function which defines how to overlay new layer to the plot. By default it is `default_overlay`. Let's assume the red channel for the layers
#          which are already in the plot is ``r0``, the red channel for the new layer is ``r`` and the alpha channel
#          is ``alpha``, the overlayed color is calculated as ``r*alpha + r0*(1-alpha)``. This self-defined function
#          should accept 7 arguments which are: vectors of r, g, b channels which correspond to the layers
#          that are already in the plot, and r, g, b, alpha channels which corresponds to the new layer. All the 
#          values passed into are between 0 to 1. The returned value for this function should be a list which contains
#          r, g, b channels which correspond to the overlayed colors. Note that these 7 arguments
#          only correspond to the pixels which are covered by the new layer.
#
# == details
# This function only works under 'pixel' mode.
#
# Under "pixel" mode, color is the only graphic representation of values in the input intervals. 
# To make a more precise and robust color mapping, users may consider `circlize::colorRamp2` to create
# a color mapping function.
#
# If you want to add more than one layers to the curve, remember to set colors with transparency.
#
# ``overlay`` argument is useful for changing color themes for the overlapped areas, please refer to the vignette
# to see examples of how to swith color themes in easy ways.
#
# == value
# No value is returned.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# hc = HilbertCurve(1, 100, level = 9, mode = "pixel")
#
# x = sort(sample(100, 20))
# s = x[1:10*2 - 1]
# e = x[1:10*2]
# require(IRanges)
# ir = IRanges(s, e)
#
# hc_layer(hc, ir)
#
# hc = HilbertCurve(1, 100, level = 9, mode = "pixel")
# hc_layer(hc, ir, grid_line = 3)
#
# hc = HilbertCurve(1, 100, level = 9, mode = "pixel")
# hc_layer(hc, ir, border = "black")
#
setMethod(f = "hc_layer",
	signature = "HilbertCurve",
	definition = function(object, ir, x1 = NULL, x2 = x1, col = "red", border = NA,
	mean_mode = c("w0", "absolute", "weighted"), grid_line = 0,
	grid_line_col = "black", overlay = default_overlay) {

	if(object@MODE == "normal") {
		stop("`hc_layer()` can only be used under 'pixel' mode.")
	}

	if(!missing(ir)) {
    	if(!inherits(ir, "IRanges")) {
    		stop("It seems you want to specify positions by one or two numeric vectors, specify them as `x1 = ..., x2 = ...`")
    	}
    }

	if(missing(ir) || is.null(ir)) {
		if(is.null(x1) || is.null(x2)) {
			stop("You should either specify `ir`, or `x1` and `x2`.")
		} else {
			x1 = hc_offset(object, x1)
			x2 = hc_offset(object, x2)
			ir = IRanges(start = round(zoom(object, x1)),
				         end = round(zoom(object, x2)))
		}
	} else {
		ir = IRanges(hc_offset(object, start(ir)),
			         hc_offset(object, end(ir)))
		ir = IRanges(start = zoom(object, start(ir)),
	                 end = zoom(object, end(ir)))
	}

	s = start(object@BINS)
	e = end(object@BINS)
	mid = round((s + e)/2)
	window = IRanges(start = c(s[1], mid), end = c(mid, e[length(e)]))
		
	mtch = as.matrix(findOverlaps(window, ir))
	if(nrow(mtch) == 0) {
		return(invisible(NULL))
	}
	
	# colors correspond to mean value for each window
	mean_mode = match.arg(mean_mode)[1]

	if(length(col) == 1) col = rep(col, length(ir))

	rgb = col2rgb(col, alpha = TRUE)

	rgb_mat = average_in_window(window, ir, mtch, list(rgb[1, ], rgb[2, ], rgb[3, ], rgb[4, ]), mean_mode, c(255, 255, 255, 0))
	r = rgb_mat[, 1]/255
	g = rgb_mat[, 2]/255
	b = rgb_mat[, 3]/255
	alpha = rgb_mat[, 4]/255
	r[r > 1] = 1
	g[g > 1] = 1
	b[b > 1] = 1
	alpha[alpha > 1] = 1

	col_index = c(object@POS$x1, object@POS$x2[4^object@LEVEL-1]) + 1
	row_index = c(object@POS$y1, object@POS$y2[4^object@LEVEL-1]) + 1
	col_index = col_index[as.numeric(rownames(rgb_mat))]
	row_index = row_index[as.numeric(rownames(rgb_mat))]

	ind = row_index + (col_index - 1)*nrow(object@RGB$red)
	# object@RGB$red[ind] = r*alpha + object@RGB$red[ind] * (1-alpha)
	# object@RGB$green[ind] = g*alpha + object@RGB$green[ind] * (1-alpha)
	# object@RGB$blue[ind] = b*alpha + object@RGB$blue[ind] * (1-alpha)

	lt = overlay(object@RGB$red[ind],
		         object@RGB$green[ind],
		         object@RGB$blue[ind],
		         r, g, b, alpha)
	object@RGB$red[ind] = lt[[1]]
	object@RGB$green[ind] = lt[[2]]
	object@RGB$blue[ind] = lt[[3]]

	if(length(border) == 1) {
		if(is.logical(border) && !is.na(border)) {
			if(border) {
				border = "black"
			} else {
				border = NA
			}
		}
	}
	
	if(any(!is.na(border))) {
		mtch = as.matrix(findOverlaps(object@BINS, ir))
		rp = pintersect(object@BINS[mtch[,1]], ir[mtch[,2]])
		if(length(border) == 1) border = rep(border, length(ir))

		mean_bin_width = mean(width(object@BINS))
		tapply(seq_len(nrow(mtch)), mtch[, 2], function(ind) {
			i1 = mtch[ind, 1]  ## index for BINS
			i2 = mtch[ind[1], 2]  ## index for ir
			n = length(ind)
			
			rp = rp[ind]  # intersected part
			r1 = object@BINS[i1]  # BIN part
			pos = object@POS[i1, ,drop = FALSE]  # pos on  the 2D space

			removed_index = NULL
			if(width(rp)[1] < mean_bin_width/2) {
				removed_index = c(removed_index, 1)
			} else {
				start(rp[1]) = start(r1[1])
			}
			if(width(rp)[n] < mean_bin_width/2) {
				removed_index = c(removed_index, n)
			} else {
				end(rp[n]) = end(r1[n])
			}
			if(length(removed_index)) {
				removed_index = unique(removed_index)
				rp = rp[-removed_index]
				r1 = r1[-removed_index]
				pos = pos[-removed_index, , drop = FALSE]
			}

			if(length(rp) == 0) {
				return(NULL)
			}

			df = data.frame(x = pos$x1, y = pos$y1)
			is_border = which_border(df$x, df$y, hc_level(object)^2)

			is_border = is_border$left_border | is_border$right_border | is_border$top_border | is_border$bottom_border
			border_df = df[is_border, , drop = FALSE]

			ind = border_df$y + 1 + (border_df$x + 1 - 1)*nrow(object@RGB$red)

			border_rgb = col2rgb(border[i2], alpha = TRUE)[, 1]/255
			
			lt = overlay(object@RGB$red[ind],
				         object@RGB$green[ind],
				         object@RGB$blue[ind],
				         border_rgb[1], border_rgb[2], border_rgb[3], border_rgb[4])
			object@RGB$red[ind] = lt[[1]]
			object@RGB$green[ind] = lt[[2]]
			object@RGB$blue[ind] = lt[[3]]
			
		})
	}

	grid_line = as.integer(grid_line)
	if(grid_line > object@LEVEL) {
		stop("`grid_line` can not exceed the level of the curve.")
	}

	if(grid_line > 1) {
		grid = 2^(grid_line-1)
		grid_line_col = col2rgb(grid_line_col)/255
		n = ncol(object@RGB$red)
		for(i in 1:(grid-1)) {
			object@RGB$red[round(n/grid*i), ] = grid_line_col[1]
			object@RGB$green[round(n/grid*i), ] = grid_line_col[2]
			object@RGB$blue[round(n/grid*i), ] = grid_line_col[3]
			
			object@RGB$red[, round(n/grid*i)] = grid_line_col[1]
			object@RGB$green[, round(n/grid*i)] = grid_line_col[2]
			object@RGB$blue[, round(n/grid*i)] = grid_line_col[3]
		}
	}

	add_raster(object@RGB)

	return(invisible(NULL))
})

# == title
# Save Hilbert curve as a PNG figure
#
# == param
# -object A `HilbertCurve-class` object.
# -file file name. If the suffix of the file name is not ``.png``, 
#        it will be added automatically no matter you like it or not.
#
# == details
# A PNG figure with resolution of ``2^level x 2^level`` is generated.
#
# Only the body of the Hilbert curve will be written to PNG file.
#
# This function only works under 'pixel' mode.
#
# == value
# No value is returned.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# hc = HilbertCurve(1, 100, level = 9, mode = "pixel")
#
# x = sort(sample(100, 20))
# s = x[1:10*2 - 1]
# e = x[1:10*2]
# require(IRanges)
# ir = IRanges(s, e)
#
# hc_layer(hc, ir)
# hc_png(hc, file = "test.png")
#
setMethod(f = "hc_png",
	signature = "HilbertCurve",
	definition = function(object, file = "HilbertCurve.png") {

	if(object@MODE == "normal") {
		stop("`hc_png()` can only be used under 'pixel' mode.")
	}

	red = object@RGB$red
	green = object@RGB$green
	blue = object@RGB$blue
	red = red[seq(nrow(red), 1), , drop = FALSE]
	green = green[seq(nrow(green), 1), , drop = FALSE]
	blue = blue[seq(nrow(blue), 1), , drop = FALSE]

	img = array(dim = c(2^object@LEVEL, 2^object@LEVEL, 3))
	img[, , 1] = red
	img[, , 2] = green
	img[, , 3] = blue

	if(!grepl("\\.png", file, ignore.case = TRUE)) {
		file = paste0(file, ".png")
	}

	writePNG(img, target = file)
})

add_raster = function(lt_rgb) {
	#if(dev.interactive()) {
		red = lt_rgb$red
		green = lt_rgb$green
		blue = lt_rgb$blue
		red = red[seq(nrow(red), 1), , drop = FALSE]
		green = green[seq(nrow(green), 1), , drop = FALSE]
		blue = blue[seq(nrow(blue), 1), , drop = FALSE]

		img = array(dim = c(dim(red), 3))
		img[, , 1] = red
		img[, , 2] = green
		img[, , 3] = blue

		seekViewport(paste0("hilbert_curve_", get_plot_index()))
		grid.raster(img, x = unit(0.5, "npc"), y = unit(0.5, "npc"), width = unit(1, "npc"), height = unit(1, "npc"))
		seekViewport(name = paste0("hilbert_curve_", get_plot_index(), "global"))
		upViewport()
	#}
}

# == title
# Default color overlay for adding new layers
#
# == param
# -r0  red channel for the layers that are already in the plot.
# -g0 green channel for the layers that are already in the plot.
# -b0 blue channel for the layers that are already in the plot.
# -r red channel for the new layer
# -g green channel for the new layer
# -b blue channel for the new layer
# -alpha alpha channel for the new layer
#
# == details
# The default overlay is (take red channel for example) ``r*alpha + r0*(1-alpha)``.
#
# == value
# A list which contains overlayed RGB colors.
#
# == seealso
# Color overlay function is always used in `hc_layer,HilbertCurve-method` or `hc_layer,GenomicHilbertCurve-method`.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# # red (1, 0, 0) overlay to the grey (0.5, 0.5, 0.5) with 0.5 transparency
# default_overlay(1, 0, 0, 0.5, 0.5, 0.5, 0.5)
default_overlay = function(r0, g0, b0, r, g, b, alpha = 1) {
	list(r = r*alpha + r0*(1-alpha),
		 g = g*alpha + g0*(1-alpha),
		 b = b*alpha + b0*(1-alpha))
}
