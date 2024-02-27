<!DOCTYPE html>
<!-- Generated by pkgdown: do not edit by hand --><html lang="en"><head><meta http-equiv="Content-Type" content="text/html; charset=UTF-8"><meta charset="utf-8"><meta http-equiv="X-UA-Compatible" content="IE=edge"><meta name="viewport" content="width=device-width, initial-scale=1.0"><title>Add points to the Hilbert curve — hc_normal_points-HilbertCurve-method • HilbertCurve</title><!-- jquery --><script src="https://cdnjs.cloudflare.com/ajax/libs/jquery/3.4.1/jquery.min.js" integrity="sha256-CSXorXvZcTkaix6Yvo6HppcZGetbYMGWSFlBw8HfCJo=" crossorigin="anonymous"></script><!-- Bootstrap --><link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/twitter-bootstrap/3.4.1/css/bootstrap.min.css" integrity="sha256-bZLfwXAP04zRMK2BjiO8iu9pf4FbLqX6zitd+tIvLhE=" crossorigin="anonymous"><script src="https://cdnjs.cloudflare.com/ajax/libs/twitter-bootstrap/3.4.1/js/bootstrap.min.js" integrity="sha256-nuL8/2cJ5NDSSwnKD8VqreErSWHtnEP9E7AySL+1ev4=" crossorigin="anonymous"></script><!-- bootstrap-toc --><link rel="stylesheet" href="../bootstrap-toc.css"><script src="../bootstrap-toc.js"></script><!-- Font Awesome icons --><link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/5.12.1/css/all.min.css" integrity="sha256-mmgLkCYLUQbXn0B1SRqzHar6dCnv9oZFPEC1g1cwlkk=" crossorigin="anonymous"><link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/5.12.1/css/v4-shims.min.css" integrity="sha256-wZjR52fzng1pJHwx4aV2AO3yyTOXrcDW7jBpJtTwVxw=" crossorigin="anonymous"><!-- clipboard.js --><script src="https://cdnjs.cloudflare.com/ajax/libs/clipboard.js/2.0.6/clipboard.min.js" integrity="sha256-inc5kl9MA1hkeYUt+EC3BhlIgyp/2jDIyBLS6k3UxPI=" crossorigin="anonymous"></script><!-- headroom.js --><script src="https://cdnjs.cloudflare.com/ajax/libs/headroom/0.11.0/headroom.min.js" integrity="sha256-AsUX4SJE1+yuDu5+mAVzJbuYNPHj/WroHuZ8Ir/CkE0=" crossorigin="anonymous"></script><script src="https://cdnjs.cloudflare.com/ajax/libs/headroom/0.11.0/jQuery.headroom.min.js" integrity="sha256-ZX/yNShbjqsohH1k95liqY9Gd8uOiE1S4vZc+9KQ1K4=" crossorigin="anonymous"></script><!-- pkgdown --><link href="../pkgdown.css" rel="stylesheet"><script src="../pkgdown.js"></script><meta property="og:title" content="Add points to the Hilbert curve — hc_normal_points-HilbertCurve-method"><meta property="og:description" content="Add points to the Hilbert curve"><!-- mathjax --><script src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.5/MathJax.js" integrity="sha256-nvJJv9wWKEm88qvoQl9ekL2J+k/RWIsaSScxxlsrv8k=" crossorigin="anonymous"></script><script src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.5/config/TeX-AMS-MML_HTMLorMML.js" integrity="sha256-84DKXVJXs0/F8OTMzX4UR909+jtl4G7SPypPavF+GfA=" crossorigin="anonymous"></script><!--[if lt IE 9]>
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
    <h1>Add points to the Hilbert curve</h1>
    
    <div class="hidden name"><code>hc_normal_points-HilbertCurve-method.rd</code></div>
    </div>

    <div class="ref-description">
    <p>Add points to the Hilbert curve</p>
    </div>

    <div id="ref-usage">
    <div class="sourceCode"><pre class="sourceCode r"><code><span><span class="co"># S4 method for HilbertCurve</span></span>
<span><span class="fu">hc_normal_points</span><span class="op">(</span><span class="va">object</span>, ir <span class="op">=</span> <span class="cn">NULL</span>, x1 <span class="op">=</span> <span class="cn">NULL</span>, x2 <span class="op">=</span> <span class="va">x1</span>, gp <span class="op">=</span> <span class="fu"><a href="https://rdrr.io/r/grid/gpar.html" class="external-link">gpar</a></span><span class="op">(</span><span class="op">)</span>,</span>
<span>    pch <span class="op">=</span> <span class="fl">1</span>, size <span class="op">=</span> <span class="fu"><a href="https://rdrr.io/r/grid/unit.html" class="external-link">unit</a></span><span class="op">(</span><span class="fl">1</span>, <span class="st">"char"</span><span class="op">)</span><span class="op">)</span></span></code></pre></div>
    </div>

    <div id="arguments">
    <h2>Arguments</h2>
    <dl><dt>object</dt>
<dd><p>A <code><a href="HilbertCurve-class.html">HilbertCurve-class</a></code> object.</p></dd>

  <dt>ir</dt>
<dd><p>an <code><a href="https://rdrr.io/pkg/IRanges/man/IRanges-constructor.html" class="external-link">IRanges</a></code> object which specifies the input intervals.</p></dd>

  <dt>x1</dt>
<dd><p>if start positions are not integers, they can be set by <code>x1</code>.</p></dd>

  <dt>x2</dt>
<dd><p>if end positions are not integers, they can be set by <code>x2</code>.</p></dd>

  <dt>size</dt>
<dd><p>size of the points. It should be a <code><a href="https://rdrr.io/r/grid/unit.html" class="external-link">unit</a></code> object, pass to <code><a href="https://rdrr.io/r/grid/grid.points.html" class="external-link">grid.points</a></code>.</p></dd>

  <dt>pch</dt>
<dd><p>shape of points, pass to <code><a href="https://rdrr.io/r/grid/grid.points.html" class="external-link">grid.points</a></code>.</p></dd>

  <dt>gp</dt>
<dd><p>graphic parameters for points. It should be specified by <code><a href="https://rdrr.io/r/grid/gpar.html" class="external-link">gpar</a></code>.</p></dd>


</dl></div>
    <div id="details">
    <h2>Details</h2>
    <p>Points are added at the middle of the intervals in <code>ir</code> (or <code>x1</code> and <code>x2</code>),
so there is only one point for each interval.</p>
<p>This function is used internally. Please use <code><a href="hc_points-HilbertCurve-method.html">hc_points,HilbertCurve-method</a></code> instead.</p>
    </div>
    <div id="value">
    <h2>Value</h2>
    

<p>A data frame which contains coordinates (in the 2D space) of points.</p>
    </div>
    <div id="see-also">
    <h2>See also</h2>
    <div class="dont-index"><p><code><a href="hc_points-HilbertCurve-method.html">hc_points,HilbertCurve-method</a></code></p></div>
    </div>
    <div id="author">
    <h2>Author</h2>
    <p>Zuguang Gu &lt;z.gu@dkfz.de&gt;</p>
    </div>

    <div id="ref-examples">
    <h2>Examples</h2>
    <div class="sourceCode"><pre class="sourceCode r"><code><span class="r-in"><span><span class="co"># see documentation of hc_points</span></span></span>
<span class="r-in"><span><span class="cn">NULL</span></span></span>
<span class="r-out co"><span class="r-pr">#&gt;</span> NULL</span>
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

