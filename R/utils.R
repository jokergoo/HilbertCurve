validate_gpar = function(gp, default) {
	for(nm in names(default)) {
		if(length(gp[[nm]]) == 0) {
			gp[[nm]] = default[[nm]]
		}
	}
	return(gp)
}

