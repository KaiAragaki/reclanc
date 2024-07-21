
<!-- README.md is generated from README.Rmd. Please edit that file -->

# reclanc

<!-- badges: start -->

[![codecov](https://codecov.io/gh/KaiAragaki/reclanc/graph/badge.svg?token=Z1AquW3Awo)](https://codecov.io/gh/KaiAragaki/reclanc)
<!-- badges: end -->

reclanc is a revival of [ClaNC (Classification of microarrays to nearest
centroids), by Alan R.
Dabney](https://doi.org/10.1093/bioinformatics/bti681).

Since the source has been lost (at least to my knowledge), the code
comes from [here](https://github.com/naikai/sake/blob/master/R/clanc.R)
with heavy modification.

reclanc is a nearest-centroid classifier for expression data. It tends
to be a little more sensitive and accurate than similar models like PAM.

Besides its mere existence, reclanc differs slightly from the original
ClaNC package in a few ways:

- reclanc supports a wide variety of inputs (`data.frame`, `matrix`,
  `formula`, `recipe`, `ExpressionSet`, and `SummarizedExperiment`)
- reclanc plays nicely with [tidymodels](https://www.tidymodels.org/),
  and offloads things like making folds to `rsample` and tuning to
  `tune` (see [this
  vignette](https://kaiaragaki.github.io/reclanc/articles/case-study.html)
  for how to leverage tidymodels with reclanc).
- Provides a prediction method based on correlation, rather than
  distance - useful for predicting classes from data from different
  sequencing modalities

## Installation

You can install the development version of reclanc like so:

``` r
# install.packages("pak")
pak::pak("KaiAragaki/reclanc")
```

## How to use it

``` r
library(reclanc)
lapply(synthetic_expression, head) # dummy data
#> $expression
#>        sample1  sample2  sample3  sample4  sample5  sample6  sample7  sample8
#> gene1 8.097529 7.119188 7.304400 7.554689 7.953206 7.714925 7.512700 8.597547
#> gene2 8.641837 9.400416 8.500865 8.878687 8.318438 8.728683 7.812591 7.638167
#> gene3 3.436236 4.317915 3.435193 3.515755 3.024976 4.762209 5.048956 2.006646
#> gene4 4.368008 5.212750 4.618249 4.201365 3.195294 4.707750 5.126769 6.178658
#> gene5 2.423974 3.563816 4.062362 2.163278 2.021435 2.813873 0.000000 4.652358
#> gene6 5.371205 5.919809 4.366915 4.805534 4.834856 5.622157 3.883531 3.593082
#>        sample9 sample10 sample11 sample12
#> gene1 6.475641 7.648858 8.637526 7.345038
#> gene2 8.110285 7.906104 7.424728 7.927039
#> gene3 2.739211 3.111668 3.161077 4.306611
#> gene4 5.170265 4.259578 5.872855 6.159023
#> gene5 1.532242 3.399823 3.691250 1.932937
#> gene6 4.246205 4.637316 3.575837 2.730452
#> 
#> $classes
#> [1] A A A A A A
#> Levels: A B
```

``` r
centroids <- clanc(
  synthetic_expression$expression,
  classes = synthetic_expression$classes,
  active = 5
)
centroids
#> <clanc> 
#> $centroids
#>    class   gene expression pooled_sd active prior
#> 1      A gene13   8.936462 0.3418472      5   0.5
#> 2      A gene21   7.379940 0.5279636      5   0.5
#> 3      A  gene2   8.744821 0.3147537      5   0.5
#> 4      A gene74   4.028507 0.4940783      5   0.5
#> 5      A gene41   4.328516 0.6317005      5   0.5
#> 6      A gene66   6.124761 0.5883218      5   0.5
#> 7      A gene24   4.307301 0.7214700      5   0.5
#> 8      A gene95   6.288173 0.4462475      5   0.5
#> 9      A gene94   7.777318 0.5375914      5   0.5
#> 10     A gene52   3.743798 0.5173769      5   0.5
#> 11     B gene13   9.938137 0.3418472      5   0.5
#> 12     B  gene2   8.273987 0.3147537      5   0.5
#> 13     B gene21   6.584681 0.5279636      5   0.5
#> 14     B gene41   5.518354 0.6317005      5   0.5
#> 15     B gene74   3.226598 0.4940783      5   0.5
#> 16     B gene24   3.370467 0.7214700      5   0.5
#> 17     B gene66   7.008174 0.5883218      5   0.5
#> 18     B gene94   8.422255 0.5375914      5   0.5
#> 19     B gene95   5.703161 0.4462475      5   0.5
#> 20     B gene52   2.438579 0.5173769      5   0.5
```

For information on basic usage, see
[this](https://kaiaragaki.github.io/reclanc/articles/using-reclanc.html)
vignette. For a case study, as well as how to optimize the `active`
parameter, see
[this](https://kaiaragaki.github.io/reclanc/articles/case-study.html)
vignette.

## How it works

You can find a gentle introduction to how reclanc works
[here](https://kai.rbind.io/posts/projects-reclanc/) and a more in-depth
and statistically rigorous definition of how this algorithm works in the
[original paper](https://doi.org/10.1093/bioinformatics/bti681).

## Citation

Citation for the original ClaNC paper:

Alan R. Dabney, Classification of microarrays to nearest centroids,
Bioinformatics, Volume 21, Issue 22, November 2005, Pages 4148–4154,
<https://doi.org/10.1093/bioinformatics/bti681>
