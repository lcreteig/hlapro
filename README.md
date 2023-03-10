
<!-- README.md is generated from README.Rmd. Please edit that file -->

# hlapro

<!-- badges: start -->
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
  - Assessing the present of donor-specific antibodies (DSA)

## Installation

You can install the development version of hlapro from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("lcreteig/hlapro")
```

## Example

Right now, the only (very basic!) functionality included in hlapro is to
get mismatched HLAs from donor and recipient typings:

``` r
library(hlapro)

donor_typing <- "A1 A2 B5"
recipient_typing <- "A1 A3 B5 B12"

get_mismatches(donor_typing, recipient_typing)
#> [1] "A2"
```

## Other packages

There’s many other implementations with partly overlapping goals (some
of which hlapro might depend on in the future):

- [hlatools](https://github.com/gschofl/hlatools) (R) provides access to
  the [IPD-IMGT/HLA database](https://www.ebi.ac.uk/ipd/imgt/hla/)
- [hlaR](https://cran.r-project.org/web/packages/hlaR) (R) does cleaning
  and imputation of HLA typings
- [immunotation](http://bioconductor.org/packages/release/bioc/html/immunotation.html) (R)
  formats HLA typings and can work with haplotype/allele frequencies
- [py-ard](https://github.com/nmdp-bioinformatics/py-ard) (python)
  reduces HLA typing resolution
- [pyglstring](https://github.com/nmdp-bioinformatics/pyglstring) checks
  whether GL Strings are well formatted

## Why `hlapro`?

This package is developed mainly for use in the
[**PRO**CARE](https://doi.org/10.1111/tan.13581) project, which (among
others) aims to **pro**file (HLA-)antibodies and their role in kidney
transplant rejection. hlapro aims to be a one-stop shop for all the
analyses in the 2nd iteration of PROCARE.
