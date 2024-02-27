<!DOCTYPE html>
<!-- Generated by pkgdown: do not edit by hand --><html lang="en"><head><meta http-equiv="Content-Type" content="text/html; charset=UTF-8"><meta charset="utf-8"><meta http-equiv="X-UA-Compatible" content="IE=edge"><meta name="viewport" content="width=device-width, initial-scale=1.0"><title>Initialize a Hilbert curve specifically for genomic data — GenomicHilbertCurve • HilbertCurve</title><!-- jquery --><script src="https://cdnjs.cloudflare.com/ajax/libs/jquery/3.4.1/jquery.min.js" integrity="sha256-CSXorXvZcTkaix6Yvo6HppcZGetbYMGWSFlBw8HfCJo=" crossorigin="anonymous"></script><!-- Bootstrap --><link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/twitter-bootstrap/3.4.1/css/bootstrap.min.css" integrity="sha256-bZLfwXAP04zRMK2BjiO8iu9pf4FbLqX6zitd+tIvLhE=" crossorigin="anonymous"><script src="https://cdnjs.cloudflare.com/ajax/libs/twitter-bootstrap/3.4.1/js/bootstrap.min.js" integrity="sha256-nuL8/2cJ5NDSSwnKD8VqreErSWHtnEP9E7AySL+1ev4=" crossorigin="anonymous"></script><!-- bootstrap-toc --><link rel="stylesheet" href="../bootstrap-toc.css"><script src="../bootstrap-toc.js"></script><!-- Font Awesome icons --><link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/5.12.1/css/all.min.css" integrity="sha256-mmgLkCYLUQbXn0B1SRqzHar6dCnv9oZFPEC1g1cwlkk=" crossorigin="anonymous"><link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/5.12.1/css/v4-shims.min.css" integrity="sha256-wZjR52fzng1pJHwx4aV2AO3yyTOXrcDW7jBpJtTwVxw=" crossorigin="anonymous"><!-- clipboard.js --><script src="https://cdnjs.cloudflare.com/ajax/libs/clipboard.js/2.0.6/clipboard.min.js" integrity="sha256-inc5kl9MA1hkeYUt+EC3BhlIgyp/2jDIyBLS6k3UxPI=" crossorigin="anonymous"></script><!-- headroom.js --><script src="https://cdnjs.cloudflare.com/ajax/libs/headroom/0.11.0/headroom.min.js" integrity="sha256-AsUX4SJE1+yuDu5+mAVzJbuYNPHj/WroHuZ8Ir/CkE0=" crossorigin="anonymous"></script><script src="https://cdnjs.cloudflare.com/ajax/libs/headroom/0.11.0/jQuery.headroom.min.js" integrity="sha256-ZX/yNShbjqsohH1k95liqY9Gd8uOiE1S4vZc+9KQ1K4=" crossorigin="anonymous"></script><!-- pkgdown --><link href="../pkgdown.css" rel="stylesheet"><script src="../pkgdown.js"></script><meta property="og:title" content="Initialize a Hilbert curve specifically for genomic data — GenomicHilbertCurve"><meta property="og:description" content="Initialize a Hilbert curve specifically for genomic data"><!-- mathjax --><script src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.5/MathJax.js" integrity="sha256-nvJJv9wWKEm88qvoQl9ekL2J+k/RWIsaSScxxlsrv8k=" crossorigin="anonymous"></script><script src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.5/config/TeX-AMS-MML_HTMLorMML.js" integrity="sha256-84DKXVJXs0/F8OTMzX4UR909+jtl4G7SPypPavF+GfA=" crossorigin="anonymous"></script><!--[if lt IE 9]>
<script src="https://oss.maxcdn.com/html5shiv/3.7.3/html5shiv.min.js"></script>
<script src="https://oss.maxcdn.com/respond/1.4.2/respond.min.js"></script>
<![endif]--></head><body data-spy="scroll" data-target="#toc">
    

    <div class="container template-reference-topic">
      <header><div class="navbar navbar-default navbar-fixed-top" role="navigation">
  <div class="container">
    <div class="navbar-header">
      <button type="button" class="navbar-toggle collapsed" data-toggle="collapse" data-target="#navbar" aria-expanded="false">
        <span class="sr-only">Toggle navigation</span>
        <span class="icon-bar"></span>
        <span class="icon-bar"></span>
        <span class="icon-bar"></span>
      </button>
      <span class="navbar-brand">
        <a class="navbar-link" href="../index.html">HilbertCurve</a>
        <span class="version label label-default" data-toggle="tooltip" data-placement="bottom" title="">1.33.1</span>
      </span>
    </div>

    <div id="navbar" class="navbar-collapse collapse">
      <ul class="nav navbar-nav"><li class="dropdown">
  <a href="#" class="dropdown-toggle" data-toggle="dropdown" role="button" data-bs-toggle="dropdown" aria-expanded="false">
    Articles
     
    <span class="caret"></span>
  </a>
  <ul class="dropdown-menu" role="menu"><li class="dropdown-header">Vignettes</li>
    <li>
      <a href="../articles/hc_general.html">Making 2D Hilbert Curve</a>
    </li>
    <li>
      <a href="../articles/hc_genome.html">GenomicHilbertCurve: specific for genomic data</a>
    </li>
    <li class="divider">
    <li class="dropdown-header">More applications</li>
    <li>
      <a href="../articles/cpan.html">Visualize CPAN modules with Hilbert curve</a>
    </li>
  </ul></li>
<li>
  <a href="../reference/index.html">Reference</a>
</li>
      </ul><ul class="nav navbar-nav navbar-right"><li>
  <a href="https://github.com/jokergoo/HilbertCurve/" class="external-link">
    <span class="fab fa-github fa-lg"></span>
     
  </a>
</li>
      </ul></div><!--/.nav-collapse -->
  </div><!--/.container -->
</div><!--/.navbar -->

      

      </header><div class="row">
  <div class="col-md-9 contents">
    <div class="page-header">
    <h1>Initialize a Hilbert curve specifically for genomic data</h1>
    
    <div class="hidden name"><code>GenomicHilbertCurve.rd</code></div>
    </div>

    <div class="ref-description">
    <p>Initialize a Hilbert curve specifically for genomic data</p>
    </div>

    <div id="ref-usage">
    <div class="sourceCode"><pre class="sourceCode r"><code><span><span class="fu">GenomicHilbertCurve</span><span class="op">(</span>chr <span class="op">=</span> <span class="fu"><a href="https://rdrr.io/r/base/paste.html" class="external-link">paste0</a></span><span class="op">(</span><span class="st">"chr"</span>, <span class="fu"><a href="https://rdrr.io/r/base/c.html" class="external-link">c</a></span><span class="op">(</span><span class="fl">1</span><span class="op">:</span><span class="fl">22</span>, <span class="st">"X"</span>, <span class="st">"Y"</span><span class="op">)</span><span class="op">)</span>, species <span class="op">=</span> <span class="st">"hg19"</span>,</span>
<span>    background <span class="op">=</span> <span class="cn">NULL</span>, <span class="va">...</span><span class="op">)</span></span></code></pre></div>
    </div>

    <div id="arguments">
    <h2>Arguments</h2>
    <dl><dt>chr</dt>
<dd><p>a vector of chromosome names. Note it should have 'chr' prefix. This argument will be ignored when <code>background</code> is set.</p></dd>

  <dt>species</dt>
<dd><p>abbreviation of species, e.g. 'hg19' or 'mm10'. <code><a href="https://rdrr.io/pkg/circlize/man/read.chromInfo.html" class="external-link">read.chromInfo</a></code> is used to retrieve the chromosome information.</p></dd>

  <dt>background</dt>
<dd><p>the background can be provided as a <code><a href="https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html" class="external-link">GRanges</a></code> object. Chromosomes should be unique across rows. Or more generally, the 'seqnames' should be different.</p></dd>

  <dt>...</dt>
<dd><p>common arguments in <code><a href="HilbertCurve.html">HilbertCurve</a></code> can be used here.</p></dd>


</dl></div>
    <div id="details">
    <h2>Details</h2>
    <p>Multiple chromosomes can be visualized in a same Hilbert curve. All chromosomes
are concatenated on after the other based on the order which is specified.</p>
<p>Since chromosomes will have irregular shapes on the curve, under 'pixel' mode, 
users can set <code>border</code> option in <code><a href="hc_map-GenomicHilbertCurve-method.html">hc_map,GenomicHilbertCurve-method</a></code> to highlight 
borders of chromosomes to identify their locations on the curve.</p>
    </div>
    <div id="value">
    <h2>Value</h2>
    

<p>A <code><a href="GenomicHilbertCurve-class.html">GenomicHilbertCurve-class</a></code> object</p>
    </div>
    <div id="author">
    <h2>Author</h2>
    <p>Zuguang Gu &lt;z.gu@dkfz.de&gt;</p>
    </div>

    <div id="ref-examples">
    <h2>Examples</h2>
    <div class="sourceCode"><pre class="sourceCode r"><code><span class="r-in"><span><span class="kw"><a href="https://rdrr.io/r/base/library.html" class="external-link">require</a></span><span class="op">(</span><span class="va"><a href="https://github.com/jokergoo/circlize" class="external-link">circlize</a></span><span class="op">)</span></span></span>
<span class="r-msg co"><span class="r-pr">#&gt;</span> Loading required package: circlize</span>
<span class="r-msg co"><span class="r-pr">#&gt;</span> ========================================</span>
<span class="r-msg co"><span class="r-pr">#&gt;</span> circlize version 0.4.15</span>
<span class="r-msg co"><span class="r-pr">#&gt;</span> CRAN page: https://cran.r-project.org/package=circlize</span>
<span class="r-msg co"><span class="r-pr">#&gt;</span> Github page: https://github.com/jokergoo/circlize</span>
<span class="r-msg co"><span class="r-pr">#&gt;</span> Documentation: https://jokergoo.github.io/circlize_book/book/</span>
<span class="r-msg co"><span class="r-pr">#&gt;</span> </span>
<span class="r-msg co"><span class="r-pr">#&gt;</span> If you use it in published research, please cite:</span>
<span class="r-msg co"><span class="r-pr">#&gt;</span> Gu, Z. circlize implements and enhances circular visualization</span>
<span class="r-msg co"><span class="r-pr">#&gt;</span>   in R. Bioinformatics 2014.</span>
<span class="r-msg co"><span class="r-pr">#&gt;</span> </span>
<span class="r-msg co"><span class="r-pr">#&gt;</span> This message can be suppressed by:</span>
<span class="r-msg co"><span class="r-pr">#&gt;</span>   suppressPackageStartupMessages(library(circlize))</span>
<span class="r-msg co"><span class="r-pr">#&gt;</span> ========================================</span>
<span class="r-in"><span><span class="kw"><a href="https://rdrr.io/r/base/library.html" class="external-link">require</a></span><span class="op">(</span><span class="va"><a href="https://bioconductor.org/packages/GenomicRanges" class="external-link">GenomicRanges</a></span><span class="op">)</span></span></span>
<span class="r-msg co"><span class="r-pr">#&gt;</span> Loading required package: GenomicRanges</span>
<span class="r-msg co"><span class="r-pr">#&gt;</span> Loading required package: stats4</span>
<span class="r-msg co"><span class="r-pr">#&gt;</span> Loading required package: BiocGenerics</span>
<span class="r-msg co"><span class="r-pr">#&gt;</span> </span>
<span class="r-msg co"><span class="r-pr">#&gt;</span> Attaching package: ‘BiocGenerics’</span>
<span class="r-msg co"><span class="r-pr">#&gt;</span> The following objects are masked from ‘package:stats’:</span>
<span class="r-msg co"><span class="r-pr">#&gt;</span> </span>
<span class="r-msg co"><span class="r-pr">#&gt;</span>     IQR, mad, sd, var, xtabs</span>
<span class="r-msg co"><span class="r-pr">#&gt;</span> The following objects are masked from ‘package:base’:</span>
<span class="r-msg co"><span class="r-pr">#&gt;</span> </span>
<span class="r-msg co"><span class="r-pr">#&gt;</span>     Filter, Find, Map, Position, Reduce, anyDuplicated, aperm, append,</span>
<span class="r-msg co"><span class="r-pr">#&gt;</span>     as.data.frame, basename, cbind, colnames, dirname, do.call,</span>
<span class="r-msg co"><span class="r-pr">#&gt;</span>     duplicated, eval, evalq, get, grep, grepl, intersect, is.unsorted,</span>
<span class="r-msg co"><span class="r-pr">#&gt;</span>     lapply, mapply, match, mget, order, paste, pmax, pmax.int, pmin,</span>
<span class="r-msg co"><span class="r-pr">#&gt;</span>     pmin.int, rank, rbind, rownames, sapply, setdiff, sort, table,</span>
<span class="r-msg co"><span class="r-pr">#&gt;</span>     tapply, union, unique, unsplit, which.max, which.min</span>
<span class="r-msg co"><span class="r-pr">#&gt;</span> Loading required package: S4Vectors</span>
<span class="r-wrn co"><span class="r-pr">#&gt;</span> <span class="warning">Warning: </span>package ‘S4Vectors’ was built under R version 4.3.2</span>
<span class="r-msg co"><span class="r-pr">#&gt;</span> </span>
<span class="r-msg co"><span class="r-pr">#&gt;</span> Attaching package: ‘S4Vectors’</span>
<span class="r-msg co"><span class="r-pr">#&gt;</span> The following object is masked from ‘package:utils’:</span>
<span class="r-msg co"><span class="r-pr">#&gt;</span> </span>
<span class="r-msg co"><span class="r-pr">#&gt;</span>     findMatches</span>
<span class="r-msg co"><span class="r-pr">#&gt;</span> The following objects are masked from ‘package:base’:</span>
<span class="r-msg co"><span class="r-pr">#&gt;</span> </span>
<span class="r-msg co"><span class="r-pr">#&gt;</span>     I, expand.grid, unname</span>
<span class="r-msg co"><span class="r-pr">#&gt;</span> Loading required package: IRanges</span>
<span class="r-msg co"><span class="r-pr">#&gt;</span> Loading required package: GenomeInfoDb</span>
<span class="r-in"><span><span class="va">bed</span> <span class="op">=</span> <span class="fu"><a href="https://rdrr.io/pkg/circlize/man/generateRandomBed.html" class="external-link">generateRandomBed</a></span><span class="op">(</span><span class="op">)</span></span></span>
<span class="r-in"><span><span class="va">gr</span> <span class="op">=</span> <span class="fu"><a href="https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html" class="external-link">GRanges</a></span><span class="op">(</span>seqnames <span class="op">=</span> <span class="va">bed</span><span class="op">[[</span><span class="fl">1</span><span class="op">]</span><span class="op">]</span>, ranges <span class="op">=</span> <span class="fu"><a href="https://rdrr.io/pkg/IRanges/man/IRanges-constructor.html" class="external-link">IRanges</a></span><span class="op">(</span><span class="va">bed</span><span class="op">[[</span><span class="fl">2</span><span class="op">]</span><span class="op">]</span>, <span class="va">bed</span><span class="op">[[</span><span class="fl">3</span><span class="op">]</span><span class="op">]</span><span class="op">)</span><span class="op">)</span></span></span>
<span class="r-in"><span><span class="va">hc</span> <span class="op">=</span> <span class="fu">GenomicHilbertCurve</span><span class="op">(</span><span class="op">)</span></span></span>
<span class="r-in"><span><span class="fu"><a href="hc_points-dispatch.html">hc_points</a></span><span class="op">(</span><span class="va">hc</span>, <span class="va">gr</span><span class="op">)</span></span></span>
<span class="r-plt img"><img src="GenomicHilbertCurve-1.png" alt="" width="700" height="433"></span>
<span class="r-in"><span></span></span>
<span class="r-in"><span><span class="va">hc</span> <span class="op">=</span> <span class="fu">GenomicHilbertCurve</span><span class="op">(</span>chr <span class="op">=</span> <span class="fu"><a href="https://rdrr.io/r/base/c.html" class="external-link">c</a></span><span class="op">(</span><span class="st">"chr1"</span>, <span class="st">"chr2"</span><span class="op">)</span><span class="op">)</span></span></span>
<span class="r-in"><span><span class="fu"><a href="hc_points-dispatch.html">hc_points</a></span><span class="op">(</span><span class="va">hc</span>, <span class="va">gr</span><span class="op">)</span></span></span>
<span class="r-plt img"><img src="GenomicHilbertCurve-2.png" alt="" width="700" height="433"></span>
<span class="r-in"><span></span></span>
<span class="r-in"><span><span class="va">bg</span> <span class="op">=</span> <span class="fu"><a href="https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html" class="external-link">GRanges</a></span><span class="op">(</span>seqnames <span class="op">=</span> <span class="fu"><a href="https://rdrr.io/r/base/c.html" class="external-link">c</a></span><span class="op">(</span><span class="st">"chr1"</span>, <span class="st">"chr2"</span><span class="op">)</span>, </span></span>
<span class="r-in"><span>    ranges <span class="op">=</span> <span class="fu"><a href="https://rdrr.io/pkg/IRanges/man/IRanges-constructor.html" class="external-link">IRanges</a></span><span class="op">(</span><span class="fu"><a href="https://rdrr.io/r/base/c.html" class="external-link">c</a></span><span class="op">(</span><span class="fl">1</span>,<span class="fl">10000000</span><span class="op">)</span>, <span class="fu"><a href="https://rdrr.io/r/base/c.html" class="external-link">c</a></span><span class="op">(</span><span class="fl">10000000</span>,<span class="fl">20000000</span><span class="op">)</span><span class="op">)</span><span class="op">)</span></span></span>
<span class="r-in"><span><span class="va">hc</span> <span class="op">=</span> <span class="fu">GenomicHilbertCurve</span><span class="op">(</span>background <span class="op">=</span> <span class="va">bg</span>, level <span class="op">=</span> <span class="fl">6</span><span class="op">)</span></span></span>
<span class="r-in"><span><span class="fu"><a href="hc_points-dispatch.html">hc_points</a></span><span class="op">(</span><span class="va">hc</span>, <span class="va">gr</span>, gp <span class="op">=</span> <span class="fu">gpar</span><span class="op">(</span>fill <span class="op">=</span> <span class="fu"><a href="https://rdrr.io/pkg/circlize/man/rand_color.html" class="external-link">rand_color</a></span><span class="op">(</span><span class="fu"><a href="https://rdrr.io/r/base/length.html" class="external-link">length</a></span><span class="op">(</span><span class="va">gr</span><span class="op">)</span><span class="op">)</span><span class="op">)</span><span class="op">)</span></span></span>
<span class="r-in"><span><span class="fu"><a href="hc_map-GenomicHilbertCurve-method.html">hc_map</a></span><span class="op">(</span><span class="va">hc</span>, fill <span class="op">=</span> <span class="cn">NA</span>, border <span class="op">=</span> <span class="st">"grey"</span>, add <span class="op">=</span> <span class="cn">TRUE</span><span class="op">)</span></span></span>
<span class="r-plt img"><img src="GenomicHilbertCurve-3.png" alt="" width="700" height="433"></span>
</code></pre></div>
    </div>
  </div>
  <div class="col-md-3 hidden-xs hidden-sm" id="pkgdown-sidebar">
    <nav id="toc" data-toggle="toc" class="sticky-top"><h2 data-toc-skip>Contents</h2>
    </nav></div>
</div>


      <footer><div class="copyright">
  <p></p><p>Developed by Zuguang Gu.</p>
</div>

<div class="pkgdown">
  <p></p><p>Site built with <a href="https://pkgdown.r-lib.org/" class="external-link">pkgdown</a> 2.0.7.</p>
</div>

      </footer></div>

  


  <style>nav[data-toggle='toc'] .nav .nav {display: block;}</style></body></html>

