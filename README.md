[![Build Status](https://travis-ci.org/jokergoo/HilbertCurve.svg)](https://travis-ci.org/jokergoo/HilbertCurve) 
[![codecov](https://img.shields.io/codecov/c/github/jokergoo/HilbertCurve.svg)](https://codecov.io/github/jokergoo/HilbertCurve) 
[![bioc](http://www.bioconductor.org/shields/downloads/HilbertCurve.svg)](https://bioconductor.org/packages/stats/bioc/HilbertCurve/) 
[![bioc](http://mcube.nju.edu.cn/cgi-bin/zuguanggu/bioc_download.pl?package=HilbertCurve)](https://bioconductor.org/packages/stats/bioc/HilbertCurve/) 
[![bioc](http://www.bioconductor.org/shields/years-in-bioc/HilbertCurve.svg)](http://bioconductor.org/packages/devel/bioc/html/HilbertCurve.html)


## HilbertCurve

[Hilbert curve](https://en.wikipedia.org/wiki/Hilbert_curve) is a type of space-filling curves
that fold one dimensional axis into a two dimensional space, but with still keeping the locality.
It has advantages to visualize data with long axis in following two aspects:

1. greatly improve resolution for the visualization;
2. easy to visualize clusters because generally data points in the cluster will also be close in the Hilbert curve. 

This package aims to provide an easy and flexible way to visualize data through Hilbert curve.
The implementation and example figures are based on following sources:

- http://mkweb.bcgsc.ca/hilbert/
- http://corte.si/posts/code/hilbert/portrait/index.html
- http://bioconductor.org/packages/devel/bioc/html/HilbertVis.html

### Citation

Zuguang Gu, Roland Eils, and Matthias Schlesner, 
[HilbertCurve: an R/Bioconductor package for high-resolution visualization of genomic data.](https://doi.org/10.1093/bioinformatics/btw161)
Bioinformatics 2016

### Install

The package is at [Bioconductor](http://bioconductor.org/packages/devel/bioc/html/HilbertCurve.html) now
and you can install the newest version by:

```r
library(devtools)
install_github("jokergoo/ComplexHeatmap")  # in order to get the newest version of ComplexHeatmap
install_github("jokergoo/HilbertCurve")
```

### Usage

Basically, there are two steps to make a Hilbert curve.

1. Initialize the curve and also map the one-dimensional axis to the curve.
2. add low-level graphics by `hc_points()`, `hc_segments()`, ... by giving the positions of the graphics.

```r
hc = HilbertCurve(1, 100, level = 4)
hc_points(hc, ...)
hc_segments(hc, ...)
hc_rect(hc, ...)
hc_text(hc, ...)
```

There is another 'pixel' mode which provides a high resolution for visualizing genomic data by the Hilbert curve.

```r
hc = HilbertCurve(1, 100000000000, level = 10)
hc_layer(hc, ...) # this can be repeated several times to add multiple layers on the curve
hc_png(hc, ...)
```

### Examples

Rainbow color spectrum:

![1](https://cloud.githubusercontent.com/assets/449218/12678993/f184c4de-c6a1-11e5-8c8c-ed3ed938c487.png)

Chinese dynasty:

![2](https://cloud.githubusercontent.com/assets/449218/12678995/f18981cc-c6a1-11e5-8b66-6222bed67c63.png)

GC percent and genes on chromosome 1:

![3](https://cloud.githubusercontent.com/assets/449218/12678996/f18a6646-c6a1-11e5-9e0b-c99cc7a93f0e.png)

Association between H3K36me3 histone modification and gene bodies:

![4](https://cloud.githubusercontent.com/assets/449218/12678992/f1848320-c6a1-11e5-8225-e6fef169f29b.png)

Methylation on chromosome 1:

![5](https://cloud.githubusercontent.com/assets/449218/12678994/f186827e-c6a1-11e5-884a-b9135f24146e.png)

Copy number alterations in 22 chromosomes:

![6](https://cloud.githubusercontent.com/assets/449218/12678997/f18e405e-c6a1-11e5-9478-3d8fdc4bc834.png)

### License

MIT @ Zuguang Gu
