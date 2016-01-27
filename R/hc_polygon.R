
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
# Add areas to Hilbert curve
#
# == param
# -object A `HilbertCurve-class` object.
# -ir a `IRanges::IRanges` object.
# -x1 if start positions are not integers, they can be set by ``x1``.
# -x2 if end positions are not integers, they can be set by ``x2``.
# -gp graphic parameters for lines. It should be specified by `grid::gpar`.
# -end_type since two ends of a continuous curve do not necessary completely overlap with 
#      the hilbert curve segments, this argument determines how to include the
#      ends of the curve. ``average``: if the end covers more than half of the segment,
#      the whole segment is included and if the end covers less than half of hte segment,
#      the segment is removed; ``expanding``: segments are included as long as they overlap;
#      ``shrinking``: segments are removed if they are not completely covered.
#
# == value
# No value is returned.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
#
# ir = IRanges(10, 40)
#
# hc = HilbertCurve(1, 100, level = 4, reference = TRUE)
# hc_segments(hc, ir)
# hc_polygon(hc, ir, gp = gpar(fill = "#FF000080", col = 1))
#
setMethod(f = "hc_polygon",
	signature = "HilbertCurve",
	definition = function(object, ir, x1 = NULL, x2 = NULL, 
	gp = gpar(),
	end_type = c("average", "expanding", "shrinking")) {

	if(object@MODE == "pixel") {
		stop("`hc_polygon()` can only be used under 'normal' mode.")
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

	gp = validate_gpar(gp, default = list(lty = 1, lwd = 1, col = 1, fill = "transparent"))

	if(length(gp$lty) == 1) gp$lty = rep(gp$lty, length(ir))
	if(length(gp$lwd) == 1) gp$lwd = rep(gp$lwd, length(ir))
	if(length(gp$col) == 1) gp$col = rep(gp$col, length(ir))
	if(length(gp$fill) == 1) gp$fill = rep(gp$fill, length(ir))

	mtch = as.matrix(findOverlaps(object@BINS, ir))
	if(nrow(mtch) == 0) {
		return(invisible(NULL))
	}
	r = pintersect(object@BINS[mtch[,1]], ir[mtch[,2]])

	mean_bin_width = mean(width(object@BINS))
	end_type = match.arg(end_type)[1]

	seekViewport(name = paste0("hilbert_curve_", get_plot_index()))
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
			return(NULL)
		}

		sr1 = start(r1)
		sr = start(r)
		er1 = end(r1)
		er = end(r)

		x1 = pos$x2 - (pos$x2 - pos$x1)*(er1 - sr)/(er1 - sr1)
		y1 = pos$y2 - (pos$y2 - pos$y1)*(er1 - sr)/(er1 - sr1)
		x2 = pos$x2 - (pos$x2 - pos$x1)*(er1 - er)/(er1 - sr1)
		y2 = pos$y2 - (pos$y2 - pos$y1)*(er1 - er)/(er1 - sr1)

		if(max(i1) == length(object@BINS)) {
			df = data.frame(x = c(pos$x1, pos$x2[nrow(pos)]), y = c(pos$y1, pos$y2[nrow(pos)]))
		} else {
			df = data.frame(x = c(pos$x1), y = c(pos$y1))
		}
		grid.rect(df$x, df$y, width = 1, height = 1, default.units = "native", gp = gpar(fill = gp$fill[i2], 
			col = NA, lineend = "butt", linejoin = "mitre"))

		border = which_border(df$x, df$y, hc_level(object)^2)
		l = border$left_border
		if(sum(l)) grid.segments(df$x[l] - 0.5, df$y[l] - 0.5, df$x[l] - 0.5, df$y[l] + 0.5, default.units = "native",
			gp = gpar(lwd = gp$lwd[i2], lty = gp$lty[i2], col = gp$col[i2], lineend ="butt", linejoin = "mitre"))
		l = border$right_border
		if(sum(l)) grid.segments(df$x[l] + 0.5, df$y[l] - 0.5, df$x[l] + 0.5, df$y[l] + 0.5, default.units = "native",
			gp = gpar(lwd = gp$lwd[i2], lty = gp$lty[i2], col = gp$col[i2], lineend ="butt", linejoin = "mitre"))
		l = border$top_border
		if(sum(l)) grid.segments(df$x[l] - 0.5, df$y[l] + 0.5, df$x[l] + 0.5, df$y[l] + 0.5, default.units = "native",
			gp = gpar(lwd = gp$lwd[i2], lty = gp$lty[i2], col = gp$col[i2], lineend ="butt", linejoin = "mitre"))
		l = border$bottom_border
		if(sum(l)) grid.segments(df$x[l] - 0.5, df$y[l] - 0.5, df$x[l] + 0.5, df$y[l] - 0.5, default.units = "native",
			gp = gpar(lwd = gp$lwd[i2], lty = gp$lty[i2], col = gp$col[i2], lineend ="butt", linejoin = "mitre"))
	})

	seekViewport(name = paste0("hilbert_curve_", get_plot_index(), "global"))
	upViewport()

	return(invisible(NULL))
	
})
