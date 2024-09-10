validate_gpar = function(gp, default) {
	for(nm in names(default)) {
		if(length(gp[[nm]]) == 0) {
			gp[[nm]] = default[[nm]]
		}
	}
	return(gp)
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

validate_input = function(object, ir, x1, x2) {
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
            ir = IRanges(start = round(zoom(object, x1)),
                         end = round(zoom(object, x2)))
        }
    } else {
        ir = IRanges(hc_offset(object, start(ir)),
                     hc_offset(object, end(ir)))
        ir = IRanges(start = zoom(object, start(ir)),
                     end = zoom(object, end(ir)))
    }
    return(ir)
}

# == title
# Whether the color is white
#
# == param
# -r Red channel.
# -g Green channel.
# -b Blue channel.
# -maxColorValue 1 or 255.
#
is_white = function(r, g, b, maxColorValue = 1) {
	r == maxColorValue & g == maxColorValue & b == maxColorValue
}



is_RStudio_current_dev = function() {
    dv = names(dev.list())
    if(length(dv) < 2) {
        FALSE
    } else {
        n = length(dv)
        if(dv[n-1] == "RStudioGD") {
            TRUE
        } else {
            FALSE
        }
    }
}


stop_wrap = function (...) {
    x = paste0(...)
    x = paste(strwrap(x), collapse = "\n")
    stop(x, call. = FALSE)
}

message_wrap = function (...)  {
    x = paste0(...)
    x = paste(strwrap(x), collapse = "\n")
    message(x)
}

warning_wrap = function (...)  {
    x = paste0(...)
    x = paste(strwrap(x), collapse = "\n")
    warning(x)
}

