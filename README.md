
<!-- README.md is generated from README.Rmd. Please edit that file -->
codonr
======

The goal of codonr is to ...

Installation
------------

You can install the released version of codonr from [GitHub](https://github.com/santiago1234/codonr) with:

``` r
devtools::install_github("santiago1234/codonr")
```

Example
-------

This is a basic example which shows you how to solve a common problem:

-   Count codon in dna sequences

``` r
library(codonr)

dna_seq <- "ATGACAGAGGGGGGACAGTTAATG"
count_codons(dna_seq)
#> # A tibble: 7 x 2
#>   codon     n
#>   <chr> <int>
#> 1 ACA       1
#> 2 ATG       2
#> 3 CAG       1
#> 4 GAG       1
#> 5 GGA       1
#> 6 GGG       1
#> 7 TTA       1
```
