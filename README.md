
<!-- README.md is generated from README.Rmd. Please edit that file -->

# stmoran

<!-- badges: start -->
<!-- badges: end -->

The goal of stmoran is to make the spatio-temporal variants of global
and local Moran’s proposed by (Gao et al. 2019) (spatio-temporal LISA)
available in R. It borrows code from and works as an extension to Roger
Bivand’s `spdep` package. I put this together because I needed it: it
has not undergone extensive testing and is likely to contain bugs and
errors. Please let me know if you encounter any.

## Installation

You can install the development version of stmoran from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("nbtjaden/stmoran")
```

You will also need a working installation of `spdep`:

``` r
install.packages("spdep")
```

## Usage

The functions `moran_st` and `localmoran_st` are the spatio-temporal
equivalents of `spdep::moran` and `spdep::localmoran`. The main
difference in using them is that `x` is expected to be XXXX instead of
YYYYYY.

``` r
library(stmoran)
## add examples
```

## References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-gao_measuring_2019" class="csl-entry">

Gao, Yong, Jing Cheng, Haohan Meng, and Yu Liu. 2019. “Measuring
Spatio-Temporal Autocorrelation in Time Series Data of Collective Human
Mobility.” *Geo-Spatial Information Science* 22 (3): 166–73.
<https://doi.org/10.1080/10095020.2019.1643609>.

</div>

</div>
