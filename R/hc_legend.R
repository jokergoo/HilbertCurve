
hc_legend = function(start_from, first_segment) {
	level = 2
	n = 4^level
	pos = hilbert_curve_cpp(level)
	pos = data.frame(x = pos[[1]], y = pos[[2]])
	
	if(start_from == "bottomleft") {
		if(first_segment == "horizontal") {
			# default, no transformation
		} else if(first_segment == "vertical") {
			pos = rotate(pos, degree = 90)
			pos = flip(pos, direction = "horizontal")
		}
	} else if(start_from == "topleft") {
		if(first_segment == "horizontal") {
			pos = rotate(pos, degree = 180)
			pos = flip(pos, direction = "horizontal")
		} else {
			pos = rotate(pos, degree = 270)
		}
	} else if(start_from == "bottomright") {
		if(first_segment == "horizontal") {
			pos = flip(pos, direction = "horizontal")
		} else {
			pos = rotate(pos, degree = 90)
		}
	} else if(start_from == "topright") {
		if(first_segment == "horizontal") {
			pos = rotate(pos, degree = 180)
		} else {
			pos = rotate(pos, degree = 270)
			pos = flip(pos, direction = "horizontal")
		}
	}

	pos = data.frame(x1 = pos$x[-n], y1 = pos$y[-n], x2 = pos$x[-1], y2 = pos$y[-1])

	vp = viewport(x = 0, y = 0, just = c("left", "bottom"), xscale = c(0, 2^level - 1), yscale = c(0, sqrt(n)-1), width = unit(0.8, "cm"), height = unit(0.8, "cm"))
	gb_hc = gTree(children = gList(
		        segmentsGrob(pos$x1, pos$y1, pos$x2, pos$y2, default.units = "native"),
		        pointsGrob(pos$x1[1], pos$y1[1], default.units = "native", pch = 16, size = unit(4, "pt")),
				segmentsGrob(pos$x1[n-1], pos$y1[n-1], pos$x2[n-1], pos$y2[n-1], 
					default.units = "native", arrow = arrow(angle = 15, length = unit(6, "pt"), type = "closed"), gp = gpar(fill = "black"))
			), vp = vp)

	title_height = convertHeight(grobHeight(textGrob("Direction", gp = gpar(fontsize = 10, fontface = "bold"))), unitTo = "mm") + unit(2, "mm")
	title_width = convertWidth(grobWidth(textGrob("Direction", gp = gpar(fontsize = 10, fontface = "bold"))), unitTo = "mm")
	gb_title = textGrob("Direction", x = 0, y = 1, just = c("left", "top"), gp = gpar(fontsize = 10, fontface = "bold"), 
		vp = viewport(x = 0, y = unit(1, "npc"), just = c("left", "top"), width = unit(0.8, "cm"), height = title_height))
	
	gb = gTree(children = gList(gb_title, gb_hc), cl = "hc_legend", 
		vp = viewport(width = title_width + unit(4, "pt"), height = unit(0.8, "cm") + title_height))

	attr(gb, "width") = gb$vp$width
	attr(gb, "height") = gb$vp$height
	gb
}

