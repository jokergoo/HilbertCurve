
# given a list of points (with x and y coordinate),
# test whether the points are one the border of this groups of points
# x: 0..(n-1), y: 0..(n-1)
which_border = function(x, y, n) {

	min_x_pos = min(x)
	max_x_pos = max(x)
	min_y_pos = min(y)
	max_y_pos = max(y)

	mat = matrix(NA, nrow = max_x_pos - min_x_pos + 1 + 2,
		             ncol = max_y_pos - min_y_pos + 1 + 2)

	x = x - min_x_pos + 1 + 1 # start from 2
	y = y - min_y_pos + 1 + 1
	
	nr = nrow(mat)
	ind = x + (y - 1) * nr
	mat[ind] = 1

	left_x = mat[x - 1 + (y - 1)*nr]
	left_y = mat[x     + (y - 1)*nr]
	right_x = mat[x + 1 + (y - 1)*nr]
	right_y = mat[x     + (y - 1)*nr]
	bottom_x = mat[x + (y - 1)*nr]
	bottom_y = mat[x + (y - 1 - 1)*nr]
	top_x = mat[x + (y - 1)*nr]
	top_y = mat[x + (y - 1 + 1)*nr]

	return(list(left_border = is.na(left_x) | is.na(left_y),
		        right_border = is.na(right_x) | is.na(right_y),
		        bottom_border = is.na(bottom_x) | is.na(bottom_y),
		        top_border = is.na(top_x) | is.na(top_y)))
}

# == title
# Add polygons to Hilbert curve
#
# == param
# -object a `HilbertCurve-class` object.
# -ir an `IRanges::IRanges` object which specifies the input intervals.
# -x1 if start positions are not integers, they can be set by ``x1``.
# -x2 if end positions are not integers, they can be set by ``x2``.
# -gp graphic parameters. It should be specified by `grid::gpar`.
# -end_type since two ends of a continuous interval do not necessarily completely overlap with 
#      the Hilbert curve segments, this argument controls how to determine the
#      ends of the interval which will be presented on the curve. 
#      ``average``: if the end covers more than half of the segment,
#      the whole segment is included and if the end covers less than half of the segment,
#      the segment is removed; ``expanding``: segments are included as long as they are overlapped;
#      ``shrinking``: segments are removed if they are not completely covered.
#
# == detail
# Drawing polygons are quite visually similar as drawing rectangles.
# The major differences are: 1) for rectangles, colors for the ends of the interval can change if they are not completely
# covered by the Hilbert curve segments, and 2) polygons can have borders.
#
# Normally polygons are used to mark areas in the Hilbert curve.
#
# == value
# No value is returned.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# require(IRanges)
# ir = IRanges(10, 40)
#
# hc = HilbertCurve(0, 100, level = 4, reference = TRUE)
# hc_segments(hc, ir)
# hc_text(hc, x1 = 10:40, labels = 10:40)
# hc_polygon(hc, ir, gp = gpar(fill = "#FF000080", col = 1))
#
setMethod(f = "hc_polygon",
	signature = "HilbertCurve",
	definition = function(object, ir = NULL, x1 = NULL, x2 = NULL, 
	gp = gpar(), end_type = c("expanding", "average", "shrinking")) {

	end_type = match.arg(end_type)[1]

	if(object@MODE == "pixel") {
		oi = .ENV$I_PLOT
		seekViewport(paste0("hilbert_curve_", .ENV$I_PLOT))

		hc2 = HilbertCurve(s = object@data_range[1], e = object@data_range[2], mode = "normal", 
			level = min(object@LEVEL, 9), newpage = FALSE, zoom = object@ZOOM, legend = FALSE,
			start_from = object@start_from, first_seg = object@first_seg, padding = unit(0, "mm"), legend = FALSE)
		hc_polygon(hc2, ir = ir, x1 = x1, x2 = x2, gp = gp, end_type = end_type)
		seekViewport(name = paste0("hilbert_curve_", oi, "_global"))
		upViewport()
		return(invisible(NULL))
	}
	seekViewport(name = paste0("hilbert_curve_", get_plot_index()))
	polygons = get_polygons(object, ir = ir, x1 = x1, x2 = x2, end_type = end_type)

	gp = validate_gpar(gp, default = list(lty = 1, lwd = 1, col = 1, fill = "transparent"))

	if(length(gp$lty) == 1) gp$lty = rep(gp$lty, length(polygons))
	if(length(gp$lwd) == 1) gp$lwd = rep(gp$lwd, length(polygons))
	if(length(gp$col) == 1) gp$col = rep(gp$col, length(polygons))
	if(length(gp$fill) == 1) gp$fill = rep(gp$fill, length(polygons))

	for(i in seq_along(polygons)) {
		if(nrow(polygons[[i]]) >= 2) {
			grid.polygon(polygons[[i]][, 1], polygons[[i]][, 2], default.units = "native", gp = gpar(fill = gp$fill[i], 
				col = gp$col[i], lty = gp$lty[i], lwd = gp$lwd[i], lineend = "butt", linejoin = "mitre"))
		}
	}
	seekViewport(name = paste0("hilbert_curve_", get_plot_index(), "_global"))
	upViewport()

	return(invisible(NULL))
	
})


get_polygons = function(object, ir = NULL, x1 = NULL, x2 = NULL, end_type = c("average", "expanding", "shrinking")) {
	
	if(!is.null(ir)) {
    	if(!inherits(ir, "IRanges")) {
    		stop("It seems you want to specify positions by one or two numeric vectors, specify them with argument names. (`x1 = ..., x2 = ...`)")
    	}
    }
    if(is.null(ir)) {
		if(is.null(x1) || is.null(x2)) {
			stop("You should either specify `ir`, or `x1` and `x2`.")
		} else {
			x1 = hc_offset(object, x1)
			x2 = hc_offset(object, x2)
			start = round(zoom(object, x1))
			end = round(zoom(object, x2))
			if(length(start) > 1) {
				l = start[2:length(start)] == end[1:(length(end)-1)] & x2[2:length(start)] > x1[1:(length(end)-1)]
				start[which(l) + 1] = start[which(l) + 1] + 1
			}
		}
	} else {
		ir = IRanges(hc_offset(object, start(ir)),
			         hc_offset(object, end(ir)))
		start = zoom(object, start(ir))
	    end = zoom(object, end(ir))
	    if(length(start) > 1) {
			l = start[2:length(start)] == end[1:(length(end)-1)] & start(ir)[2:length(start)] > end(ir)[1:(length(end)-1)]
			start[which(l) + 1] = start[which(l) + 1] + 1
		}
	}
	ir = IRanges(start = start, end = end)

	if(length(ir) > 200) {
		cat("It is not suggested to use `hc_polygon()` on more than 200 intervals.\n")
	}

	mtch = as.matrix(findOverlaps(object@BINS, ir))
	if(nrow(mtch) == 0) {
		return(invisible(NULL))
	}
	r = pintersect(object@BINS[mtch[,1]], ir[mtch[,2]])

	mean_bin_width = mean(width(object@BINS))
	end_type = match.arg(end_type)[1]

	tapply(seq_len(nrow(mtch)), mtch[, 2], function(ind) {
		i1 = mtch[ind, 1]  ## index for BINS
		i2 = mtch[ind[1], 2]  ## index for ir
		n = length(ind)

		r = r[ind]  # intersected part
		r1 = object@BINS[i1]  # BIN part
		pos = object@POS[i1, ,drop = FALSE]  # pos on  the 2D space

		if(end_type == "expanding") {
			start(r[1]) = start(r1[1])
			end(r[n]) = end(r1[n])
		} else if(end_type == "shrinking") {
			removed_index = NULL
			if(width(r)[1] < mean_bin_width) {
				removed_index = c(removed_index, 1)
			}
			if(width(r)[n] < mean_bin_width) {
				removed_index = c(removed_index, n)
			}
			if(length(removed_index)) {
				removed_index = unique(removed_index)
				r = r[-removed_index]
				r1 = r1[-removed_index]
				pos = pos[-removed_index, , drop = FALSE]
			}
		} else if(end_type == "average") {
			removed_index = NULL
			if(width(r)[1] < mean_bin_width/2) {
				removed_index = c(removed_index, 1)
			} else {
				start(r[1]) = start(r1[1])
			}
			if(width(r)[n] < mean_bin_width/2) {
				removed_index = c(removed_index, n)
			} else {
				end(r[n]) = end(r1[n])
			}
			if(length(removed_index)) {
				removed_index = unique(removed_index)
				r = r[-removed_index]
				r1 = r1[-removed_index]
				pos = pos[-removed_index, , drop = FALSE]
			}
		}

		if(length(r) == 0) {
			return(matrix(nrow = 0, ncol = 2))
		}

		if(max(i1) == length(object@BINS)) {
			df = data.frame(x = c(pos$x1, pos$x2[nrow(pos)]), y = c(pos$y1, pos$y2[nrow(pos)]))
		} else {
			df = data.frame(x = c(pos$x1), y = c(pos$y1))
		}

		border = which_border(df$x, df$y, hc_level(object)^2)
		seg_mat = matrix(nrow = sum(border$left_border) + sum(border$right_border) + sum(border$top_border) + sum(border$bottom_border), ncol = 4)
		colnames(seg_mat) = c("x0", "y0", "x1", "y1")

		l = border$left_border
		k = 0
		if(sum(l)) {
			seg_mat[k + 1:sum(l), ] = cbind(x0 = df$x[l] - 0.5, y0 = df$y[l] - 0.5, x1 = df$x[l] - 0.5, y1 = df$y[l] + 0.5)
			k = k + sum(l) 
		}
		l = border$right_border
		if(sum(l)) {
			seg_mat[k + 1:sum(l), ] = cbind(x0 = df$x[l] + 0.5, y0 = df$y[l] - 0.5, x1 = df$x[l] + 0.5, y1 = df$y[l] + 0.5)
			k = k + sum(l) 
		}
		l = border$top_border
		if(sum(l)) {
			seg_mat[k + 1:sum(l), ] = cbind(x0 = df$x[l] - 0.5, y0 = df$y[l] + 0.5, x1 = df$x[l] + 0.5, y1 = df$y[l] + 0.5)
			k = k + sum(l) 
		}
		l = border$bottom_border
		if(sum(l)) {
			seg_mat[k + 1:sum(l), ] = cbind(x0 = df$x[l] - 0.5, y0 = df$y[l] - 0.5, x1 = df$x[l] + 0.5, y1 = df$y[l] - 0.5)
			k = k + sum(l) 
		}
		seg_mat = seg_mat[!is.na(seg_mat[, 1]),  , drop = FALSE]
		return(reorder_polygon_segments(seg_mat))
	})
}

reorder_polygon_segments = function(pos) {

	pos2 = matrix(nrow = nrow(pos), ncol = ncol(pos))
	colnames(pos2) = colnames(pos)

	if(nrow(pos) == 0) {
		return(pos2)
	}

	pos2[1, ] = pos[1, ,drop = FALSE]
	pos[1, 1] = -Inf
	pos[1, 2] = -Inf
	pos[1, 3] = -Inf
	pos[1, 4] = -Inf

	if(nrow(pos) == 1) {
		return(pos2)
	}

	k = 1
	for(i in seq_len(nrow(pos))) {
		x = pos2[k, "x1"]
		y = pos2[k, "y1"]
		k = k + 1

		ind = which(pos[, "x0"] == x & pos[, "y0"] == y)
		if(length(ind)) {
			ind = ind[1]
			pos2[k, ] = pos[ind, ]
			# qqcat("[@{x}, @{y}] -> [@{pos2[k, 3]}, @{pos2[k, 4]}]\n")
			pos[ind, 1] = -Inf
			pos[ind, 2] = -Inf
			pos[ind, 3] = -Inf
			pos[ind, 4] = -Inf
		} else {
			ind = which(pos[, "x1"] == x & pos[, "y1"] == y)
			if(length(ind)) {
				ind = ind[1]
				pos2[k, ] = cbind(x0 = pos[ind, 3], y0 = pos[ind, 4], x1 = pos[ind, 1], y1 = pos[ind, 2])
				# qqcat("[@{x}, @{y}] -> [@{pos2[k, 3]}, @{pos2[k, 4]}]\n")
				pos[ind, 1] = -Inf
				pos[ind, 2] = -Inf
				pos[ind, 3] = -Inf
				pos[ind, 4] = -Inf
			} else {
				stop("polygon is not connected.")
			}
		}
		if(k == nrow(pos2)) {
			break
		}
	}
	return(pos2)
}



