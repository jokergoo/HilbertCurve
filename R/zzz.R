.onAttach = function(libname, pkgname) {
    version = packageDescription(pkgname, fields = "Version")

  	msg = paste0("========================================
", pkgname, " version ", version, "
Bioconductor page: http://bioconductor.org/packages/HilbertCurve/
Github page: https://github.com/jokergoo/HilbertCurve
Documentation: http://bioconductor.org/packages/HilbertCurve/

If you use it in published research, please cite:
Gu, Z. HilbertCurve: an R/Bioconductor package for high-resolution 
  visualization of genomic data. Bioinformatics 2016.
========================================
")	

    packageStartupMessage(msg)
}
