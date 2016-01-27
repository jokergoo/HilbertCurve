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
