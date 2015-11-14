[![Build Status](https://travis-ci.org/jokergoo/HilbertCurve.svg)](https://travis-ci.org/jokergoo/HilbertCurve)

# HilbertCurve

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

## Install

The package is at [Bioconductor](http://bioconductor.org/packages/devel/bioc/html/HilbertCurve.html) now
and you can install the newest version by:

```r
library(devtools)
install_github("jokergoo/ComplexHeatmap")  # in order to get the newest version of ComplexHeatmap
install_github("jokergoo/HilbertCurve")
```

## Usage

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

## Examples

Show rainbow color spectrum:

![download](https://cloud.githubusercontent.com/assets/449218/8629460/e8e62eca-275b-11e5-9249-e871194f953f.png)

Genes on chromosome 1:

![download 1](https://cloud.githubusercontent.com/assets/449218/8629458/e8dcaa1c-275b-11e5-95c4-c36a9e1629ce.png)

Methylation on chromosome 1:

![download 2](https://cloud.githubusercontent.com/assets/449218/8629459/e8e1fc2e-275b-11e5-99e0-ed7cfd9aa0ad.png)

GC percent and genes on chromosome 1:

![download 3](https://cloud.githubusercontent.com/assets/449218/8629472/01f441d6-275c-11e5-8760-4d9412c4754d.png)
