
<!-- README.md is generated from README.Rmd. Please edit that file -->

# hlapro

<!-- badges: start -->

[![check-standard](https://github.com/lcreteig/hlapro/actions/workflows/check-standard.yaml/badge.svg)](https://github.com/lcreteig/hlapro/actions/workflows/check-standard.yaml)
[![Codecov test
coverage](https://codecov.io/gh/lcreteig/hlapro/branch/main/graph/badge.svg)](https://app.codecov.io/gh/lcreteig/hlapro?branch=main)
[![CRAN
status](https://www.r-pkg.org/badges/version/hlapro)](https://CRAN.R-project.org/package=hlapro)
<!-- badges: end -->

The goal of hlapro is to provide some tooling to work with [Human
Leukocyte Antigen
(HLA)](https://en.wikipedia.org/wiki/Human_leukocyte_antigen) data, in
the context of transplantation research:

- HLA typings
- single bead assays (e.g. Luminex) for detecting HLA-antibodies

The package is currently still in (heavy) development, but ultimately
the goal is to support the following functionality:

- **HLA typings**:
  - Obtaining HLA mismatches between transplant donor and recipient
  - Formatting into [GL Strings](https://glstring.org)
  - Downscaling (reducing) typing resolution
  - Upscaling (imputing) typing resolution, based on haplotype frequency
    tables
- **HLA antibody assays**:
  - Determining which beads in the assay are positive, based on
    automated implementations of different methods
  - Assessing the presence of donor-specific antibodies (DSA)

## Installation

You can install the development version of hlapro from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("lcreteig/hlapro")
```

## Usage

``` r
library(hlapro)
```

### Mismatches

Get mismatched HLAs from donor and recipient typings

``` r
donor_typing <- "A1 A2 B5"
recipient_typing <- "A1 A3 B5 B12"

get_mismatches(donor_typing, recipient_typing)
#> [1] "A2"
```

### Extracting alleles

Extract HLA alleles from a typing string

``` r
typing <- "A1 A2 B7 B8 Cw3 DQ5 DQ8 DR4 DR11 DR52 DR53"
extract_alleles_str(typing)
#>    A_1    A_2    B_1    B_2    C_1    C_2 DPB1_1 DPB1_2 DQA1_1 DQA1_2 DQB1_1 
#>    "1"    "2"    "7"    "8"    "3"     NA     NA     NA     NA     NA    "5" 
#> DQB1_2 DRB1_1 DRB1_2 DRB._1 DRB._2 
#>    "8"    "4"   "11"   "52"   "53"

df <- tidyr::tibble(typing = "A1 A2 B7 B8 Cw3 DQ5 DQ8 DR4 DR11 DR52 DR53")
extract_alleles_df(df, typing, loci = c("A", "B", "C"))
#> Joining with `by = join_by(typing)`
#> Joining with `by = join_by(typing)`
#> # A tibble: 1 × 7
#>   typing                                     A_1   A_2   B_1   B_2   C_1   C_2  
#>   <chr>                                      <chr> <chr> <chr> <chr> <chr> <chr>
#> 1 A1 A2 B7 B8 Cw3 DQ5 DQ8 DR4 DR11 DR52 DR53 1     2     7     8     3     ""
```

### Validating alleles

Check whether alleles are well-formed

``` r
validate_allele(c("A2", "A*01:AABJE", "A*24:02:01:02L", "not-an-HLA"))
#> [1]  TRUE  TRUE  TRUE FALSE
```

### Allele resolution

Determine whether an allele is of low/intermediate/high resolution

``` r
get_resolution(c("A2", "A*01:AABJE", "B*42:08"))
#> [1] "low"          "intermediate" "high"
```

### Downscaling to serological equivalents

Get the serological equivalent of an allele as defined by the [ETRL
HLA](https://etrl.eurotransplant.org/resources/hla-tables/) conversion
tables

``` r
etrl_lookup(c("B*15:79N", "B*15:YETY", "B*15:01:16", "B*15:02", "B*15:85"))
#> # A tibble: 5 × 4
#>   Allele  `ET MatchDeterminantSplit` `ET MatchDeterminantBroad` Public
#>   <chr>   <chr>                      <chr>                      <chr> 
#> 1 <NA>    <NA>                       <NA>                       <NA>  
#> 2 B*15:XX <NA>                       B15                        <NA>  
#> 3 B*15:01 B62                        B15                        Bw6   
#> 4 B*15:02 B75                        B15                        Bw6   
#> 5 B*15:XX <NA>                       B15                        <NA>
```

For alleles in serological notation, conversion to broads and public
epitopes is also supported

``` r
splits <- c("A24", "A25")
get_broad(splits)
#> [1] "A9"  "A10"
```

``` r
b_s <- c("B14", "B63")
get_public(b_s)
#> [1] "Bw6" "Bw4"
```

## Other packages

There’s many other implementations with partly overlapping goals (some
of which hlapro might depend on in the future).

- R:
  - [hlatools](https://github.com/gschofl/hlatools) provides access to
    the [IPD-IMGT/HLA database](https://www.ebi.ac.uk/ipd/imgt/hla/)
  - [hlaR](https://cran.r-project.org/web/packages/hlaR) does cleaning
    and imputation of HLA typings
  - [immunotation](http://bioconductor.org/packages/release/bioc/html/immunotation.html)
    formats HLA typings and can work with haplotype/allele frequencies
- Python:
  - [py-ard](https://github.com/nmdp-bioinformatics/py-ard) reduces HLA
    typing resolution
  - [pyglstring](https://github.com/nmdp-bioinformatics/pyglstring)
    checks whether GL Strings are well formatted
  - [ALLAN](https://github.com/lgragert/hla-who-to-unos) converts
    between serological and molecular HLA notation

## Why `hlapro`?

This package is developed mainly for use in the
[**PRO**CARE](https://doi.org/10.1111/tan.13581) project, which (among
others) aims to **pro**file (HLA-)antibodies and their role in kidney
transplant rejection. hlapro aims to be a one-stop shop for all the
analyses in the second iteration of PROCARE.
