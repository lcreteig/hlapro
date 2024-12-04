
<!-- README.md is generated from README.Rmd. Please edit that file -->

# hlapro <a href="https://lcreteig.github.io/hlapro/"><img src="man/figures/logo.png" align="right" height="139" alt="hlapro website" /></a>

<!-- badges: start -->

[![R-CMD-check](https://github.com/lcreteig/hlapro/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/lcreteig/hlapro/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/lcreteig/hlapro/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/lcreteig/hlapro/actions/workflows/pkgdown.yaml)
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

### Typings

#### Mismatches

Get mismatched HLAs from donor and recipient typings

``` r
donor_typing <- "A1 A2 B5"
recipient_typing <- "A1 A3 B5 B12"

get_mismatches(donor_typing, recipient_typing)
#> [1] "A2"
```

#### Extracting alleles

Extract HLA alleles from a typing string

``` r
typing <- "A1 A2 B7 B8 Cw3 DQ5 DQ8 DR4 DR11 DR52 DR53"
extract_alleles_str(typing)
#>    A_1    A_2    B_1    B_2    C_1    C_2 DPA1_1 DPA1_2 DPB1_1 DPB1_2 DQA1_1 
#>    "1"    "2"    "7"    "8"    "3"     NA     NA     NA     NA     NA     NA 
#> DQA1_2 DQB1_1 DQB1_2 DRB1_1 DRB1_2 DRB._1 DRB._2 
#>     NA    "5"    "8"    "4"   "11"   "52"   "53"

df <- tidyr::tibble(typing = "A1 A2 B7 B8 Cw3 DQ5 DQ8 DR4 DR11 DR52 DR53")
extract_alleles_df(df, typing, loci = c("A", "B", "C"))
#> Joining with `by = join_by(typing)`
#> Joining with `by = join_by(typing)`
#> # A tibble: 1 × 7
#>   typing                                     A_1   A_2   B_1   B_2   C_1   C_2  
#>   <chr>                                      <chr> <chr> <chr> <chr> <chr> <chr>
#> 1 A1 A2 B7 B8 Cw3 DQ5 DQ8 DR4 DR11 DR52 DR53 1     2     7     8     3     ""
```

If there’s more than two alleles for a given locus in the typing, you’ll
receive a warning:

``` r
typing <- "A1 A2 A3 B7 B8 Cw1"
extract_alleles_str(typing)
#> Warning: One or more loci found with more than 2 alleles.
#> ✖ `extract_alleles_str()` will only pick the first two.
#> ℹ Use `hlapro::count_alleles()` to find out more.
#>    A_1    A_2    B_1    B_2    C_1    C_2 DPA1_1 DPA1_2 DPB1_1 DPB1_2 DQA1_1 
#>    "1"    "2"    "7"    "8"    "1"     NA     NA     NA     NA     NA     NA 
#> DQA1_2 DQB1_1 DQB1_2 DRB1_1 DRB1_2 DRB._1 DRB._2 
#>     NA     NA     NA     NA     NA     NA     NA
```

Use `count_alleles()` to easily inspect the number of alleles per locus:

``` r
count_alleles(typing)
#>    A    B    C DPA1 DPB1 DQA1 DQB1 DRB1 DRB. 
#>    3    2    1    0    0    0    0    0    0
```

#### Cleaning allele strings

Correcting some common formatting issues

``` r
allele_vec <- c(
  " A'*1 ", # spurious whitespace, punctuation, no leading zero, missing XX code
  "A10(25)", # redundant typing (includes both split and broad)
  "C*03:01/02" # locus and allele group left off in ambiguity
)

clean_hla(allele_vec)
#> [1] "A*01:XX"         "A25"             "C*03:01/C*03:02"
```

`clean_hla()` can also convert typings in the pre-2010 (v2) nomenclature
to the updated (v3) nomenclature:

``` r
clean_hla(c("Cw*030205", "A*2416", "DPB1*02BYVD", "B*35UMU"))
#> [1] "C*03:02:05"   "A*31:08"      "DPB1*02:FNWG" "B*35:WRE"
```

#### Validating alleles

Check whether alleles are well-formed

``` r
validate_allele(c("A2", "A*01:AABJE", "A*24:02:01:02L", "not-an-HLA"))
#> [1]  TRUE  TRUE  TRUE FALSE
```

#### Allele resolution

Determine whether an allele is of low/intermediate/high resolution

``` r
get_resolution(c("A2", "A*01:AABJE", "B*42:08"))
#> [1] "low"          "intermediate" "high"
```

Use the extended mode to get more information:

``` r
get_resolution(c("A2", "A*24:XX", "A*01:AB", "B*42:08", "A*01:01:01:01"),
  extended = TRUE
)
#> [1] "serology - broad"    "molecular - split"   "intermediate"       
#> [4] "high - second field" "high - fourth field"
```

#### Downscaling to serological equivalents

There’s two ways to do this:

1.  Using the [IPD-IMGT/HLA
    database](https://www.ebi.ac.uk/ipd/imgt/hla/) as implemented by the
    [NMDP](https://network.nmdp.org/services-support/bioinformatics-immunobiology/tools)
    in [py-ard](https://github.com/nmdp-bioinformatics/py-ard):

``` r
ard <- db_initialize(path.expand("~/ipd_db")) # (download &) initialize database
reduce_to_serology(
  ard,
  c("B*15:79N", "B*15:25/B*15:61", "B*15:02", "B*15:01:16")
)
#> [1] ""            "B15/B62/B70" "B15/B70/B75" "B15/B62"
```

2.  Get the serological equivalents of an allele as defined by the [ETRL
    HLA](https://etrl.eurotransplant.org/resources/hla-tables/)
    conversion tables. This is less extensive and precise, but is
    nonetheless relevant for the Eurotransplant region.

``` r
get_serology(c("B*15:79N", "B*15:25/B*15:61", "B*15:02", "B*15:01:16"))
#> [1] NA    "B15" "B75" "B62"
```

The latter also supports lookup of broads or splits specifically:

``` r
alleles <- c("A*01:01", "A*25:76:02")
get_split(alleles)
#> [1] NA    "A25"
```

``` r
splits <- c("A24", "A*25:76:02")
get_broad(splits)
#> [1] "A9"  "A10"
```

And whether an allele has the Bw4 or Bw6 epitope:

``` r
b_s <- c("B14", "B63", "B*40:05", "A*01:01")
get_public(b_s)
#> [1] "Bw6" "Bw4" "Bw6" NA
```

#### “Reducing” to two-field resolution

For many applications it suffices to work with HLA alleles at the
protein-level. `py-ard` can also leverage the IPD-IMGT/HLA database to
accomplish this:

``` r
reduce_to_field2(ard, c("A*01:01:01:01"))
#> [1] "A*01:01"
```

The above case is arguably trivial, but many others are not:

``` r
reduce_to_field2(ard, "A*01:04:01:01N") # shortnull
#> [1] "A*01:04N"
reduce_to_field2(ard, "B*15:AH") # MACs
#> [1] "B*15:01/B*15:07"
reduce_to_field2(ard, "Cw10") # serology
#> [1] "C*03:02/C*03:04/C*03:04Q/C*03:06/C*03:26/C*03:28/C*03:46"
```

#### Upscaling from serological equivalents to 2-field high resolution

Depends on the haplotype frequencies released by the NMDP, which must be
downloaded from [here](https://frequency.nmdp.org) after logging in and
accepting the license.

``` r
upscale_typings(
  filepath = "~/Downloads/A~C~B~DRB3-4-5~DRB1~DQB1.xlsx",
  typing = "A24 A28 B35 B61 DR4 DR11"
) |>
  dplyr::select(unphased_geno, dplyr::starts_with("haplo"))
#> # A tibble: 2 × 5
#>   unphased_geno              haplo_freq_1 haplo_freq_2 haplo_rank_1 haplo_rank_2
#>   <chr>                             <dbl>        <dbl>        <dbl>        <dbl>
#> 1 A*24:02g A*68:01g B*35:03…    0.000585     0.000209           240          665
#> 2 A*24:02g A*68:01g B*35:03…    0.0000849    0.0000597         1565         2058
```

Also able to upscale multiple typings at once, for instance in a
dataframe:

``` r
typing_df <- tidyr::tibble(
  id = c("001", "002"),
  input_typings = c(
    "A24 A28 B35 B61 DR4 DR11",
    "A2 A3 B52 B35 Cw4 DR11 DR52 DQ3"
  )
)
typing_df |>
  dplyr::mutate(geno_df = upscale_typings(
    "~/Downloads/A~C~B~DRB3-4-5~DRB1~DQB1.xlsx",
    input_typings,
    as_list = TRUE
  )) |>
  tidyr::unnest(geno_df)
#> # A tibble: 6 × 13
#>   id    input_typings id_unphased_geno unphased_geno unphased_freq unphased_prob
#>   <chr> <chr>                    <int> <chr>                 <dbl>         <dbl>
#> 1 001   A24 A28 B35 …              364 A*24:02g A*6…  0.000000254          0.121
#> 2 001   A24 A28 B35 …              364 A*24:02g A*6…  0.000000254          0.121
#> 3 002   A2 A3 B52 B3…               17 A*02:01g A*0…  0.0000000542         0.151
#> 4 002   A2 A3 B52 B3…               17 A*02:01g A*0…  0.0000000542         0.151
#> 5 002   A2 A3 B52 B3…               17 A*02:01g A*0…  0.0000000542         0.151
#> 6 002   A2 A3 B52 B3…               17 A*02:01g A*0…  0.0000000542         0.151
#> # ℹ 7 more variables: multiplier <dbl>, phased_freq <dbl>, phased_prob <dbl>,
#> #   haplo_freq_1 <dbl>, haplo_freq_2 <dbl>, haplo_rank_1 <dbl>,
#> #   haplo_rank_2 <dbl>
```

#### Working with multiple allele codes

py-ard can once again be used to lookup and decode [Multiple Allele
Codes (MACs)](https://hml.nmdp.org/MacUI/)

``` r
mac_lookup(ard, "B*08:01/B*08:19N/B*08:109")
#> [1] "B*08:YETY"
```

``` r
mac_decode(ard, "B*08:YETY")
#> [1] "B*08:01/B*08:19N/B*08:109"
```

#### Converting to and from GL Strings

Typing data often comes in a data frame like this:

``` r
typing_df <- tidyr::tibble(
  id = c("001", "002"),
  A_1 = c("A*01:01", "A*02:01"), A_2 = c("A*03:01", "A*29:02"),
  B_1 = c("B*08:01", "B*07:02"), B_2 = c("B*07:02", NA_character_),
  C_1 = c("C*07:01", "C*05:01"), C_2 = c("C*07:02", NA_character_)
)
typing_df
#> # A tibble: 2 × 7
#>   id    A_1     A_2     B_1     B_2     C_1     C_2    
#>   <chr> <chr>   <chr>   <chr>   <chr>   <chr>   <chr>  
#> 1 001   A*01:01 A*03:01 B*08:01 B*07:02 C*07:01 C*07:02
#> 2 002   A*02:01 A*29:02 B*07:02 <NA>    C*05:01 <NA>
```

Use `df_to_gl()` to convert all the individual allele columns into a
single [GL String](https://glstring.org):

``` r
typing_df_gl <- typing_df |>
  dplyr::group_by(id) |> # convert each typing to its own GL String
  df_to_gl()
typing_df_gl
#> # A tibble: 2 × 2
#>   id    glstring                                                                
#>   <chr> <chr>                                                                   
#> 1 001   hla#2024-12-04#HLA-A*01:01+HLA-A*03:01^HLA-B*07:02+HLA-B*08:01^HLA-C*07…
#> 2 002   hla#2024-12-04#HLA-A*02:01+HLA-A*29:02^HLA-B*07:02^HLA-C*05:01
```

Use `gl_to_df()` to go the opposite way: from a dataframe of GL Strings
to one columns per allele:

``` r
typing_df_gl |>
  dplyr::mutate(gl_df = gl_to_df(glstring)) |>
  tidyr::unnest(gl_df)
#> # A tibble: 2 × 11
#>   id    glstring      glstring_index namespace version_or_date A_1   A_2   B_1  
#>   <chr> <chr>                  <int> <chr>     <chr>           <chr> <chr> <chr>
#> 1 001   hla#2024-12-…              1 hla       2024-12-04      HLA-… HLA-… HLA-…
#> 2 002   hla#2024-12-…              2 hla       2024-12-04      HLA-… HLA-… HLA-…
#> # ℹ 3 more variables: B_2 <chr>, C_1 <chr>, C_2 <chr>
```

#### Looking up eplets and alleles

Grab database from the [HLA Eplet
registry](https://www.epregistry.com.br) and use it to lookup which
eplets occur on an HLA allele, or vice versa.

``` r
df_eplets <- load_eplet_registry()
#> Loaded Eplet Registry table ( ),
#> released 2024-08-19, downloaded from https://www.epregistry.com.br
lookup_alleles(df_eplets, "17S")
#> $`17S`
#> [1] "A*01:02" "A*30:01" "A*30:02"
```

``` r
lookup_eplets(df_eplets, "A*01:02")
#> $`A*01:02`
#>  [1] "9S"         "17S"        "44KM"       "62QE"       "62QE+56G"  
#>  [6] "65RA"       "65RNA"      "66N"        "66NH"       "66NM"      
#> [11] "71HS"       "76ANT"      "77N[ABC]"   "77NGT"      "79GT"      
#> [16] "79GT+90D"   "80T"        "80TL"       "90D"        "95I"       
#> [21] "97I"        "99Y"        "109F"       "114R"       "116D"      
#> [26] "138MI"      "138MI+79GT" "144K"       "144KR"      "144KR+151H"
#> [31] "149AH"      "151H"       "152A"       "152HA"      "156R"      
#> [36] "163R"       "163RG"      "166DG"      "193PI"      "275EL"
```

A common use case would be to lookup which eplets occur on a set of
(positive) Luminex beads:

``` r
df_eplets <- load_eplet_registry()
#> Loaded Eplet Registry table ( ),
#> released 2024-08-19, downloaded from https://www.epregistry.com.br
luminex_df <- dplyr::tribble(
  ~sampleID, ~allele, ~positive,
  "001", "A*01:01", TRUE,
  "001", "A*02:01", FALSE,
  "002", "A*01:01", FALSE
)
get_positive_eplets(luminex_df, sampleID, allele, positive, df_eplets)
#> # A tibble: 28 × 2
#>    sampleID eplets_pos
#>    <chr>    <chr>     
#>  1 001      44KM      
#>  2 001      62QE      
#>  3 001      62QE+56G  
#>  4 001      65RNA     
#>  5 001      66N       
#>  6 001      66NH      
#>  7 001      66NM      
#>  8 001      76ANT     
#>  9 001      77N[ABC]  
#> 10 001      77NGT     
#> # ℹ 18 more rows
```

### Antibodies

#### Parsing raw Luminex files

Parse `.csv` files from Luminex Single-Antigen Bead (SAB) assays, along
with the lot-specific `.eds` file that came with the kit, to produce
interpreted results identical to Immucor’s MATCH IT!® Antibody Analysis
Software:

``` r
# path to some mock data (originally for testing purposes)
lum_path <- testthat::test_path("luminex")

# parse the files
read_lum_csv(
  csv_filepath = file.path(lum_path, "LSA1-test.csv"),
  lots_path = lum_path # folder with .eds file
)
#> Joining with `by = join_by(Location, Sample, antigen_id)`
#> Joining with `by = join_by(Location, Sample, antigen_id)`
#> Joining with `by = join_by(Sample, antigen_id)`
#> Joining with `by = join_by(antigen_id)`
#> Joining with `by = join_by(antigen_id)`
#> Joining with `by = join_by(Location, Sample)`
#> Joining with `by = join_by(antigen_id)`
#> # A tibble: 12 × 18
#>    Sample Count antigen_id Cutoff  Median mfi_lra assignment bg_adjusted ad_mfi
#>    <chr>  <int> <chr>       <dbl>   <dbl>   <dbl> <chr>            <dbl>  <dbl>
#>  1 S001      67 103          3.75  8002.   127.   Positive        7822.  6713. 
#>  2 S001     104 PC          NA    18206.    36.3  Negative       16706.  4551. 
#>  3 S001      87 NC          NA      161      1    Negative          NA     NA  
#>  4 S001      73 104          4.01    63      1    Negative        -107     52.1
#>  5 S001      89 105          3.79   127      1    Negative         -53     97.4
#>  6 S001      60 106          3.49   501      1    Negative         256    634. 
#>  7 S002     120 PC          NA    17573    158.   Negative       16073   4393. 
#>  8 S002      74 104          4.01    92.5    1.34 Negative         -77.5   76.4
#>  9 S002      92 NC          NA      113      1    Negative          NA     NA  
#> 10 S002      82 103          3.75    69      1    Negative        -111     57.9
#> 11 S002     113 105          3.79   120      1    Negative         -60     92.0
#> 12 S002      78 106          3.49   112.     1    Negative        -134.   141. 
#> # ℹ 9 more variables: ad_bg_adjusted <dbl>, A <chr>, B <chr>, C <chr>,
#> #   Bw <chr>, A_serology <chr>, B_serology <chr>, C_serology <chr>, RAD <dbl>
```

## Other packages

There’s many other implementations with partly overlapping goals.

- R:
  - [hlatools](https://github.com/gschofl/hlatools) provides access to
    the [IPD-IMGT/HLA database](https://www.ebi.ac.uk/ipd/imgt/hla/)
  - [hlaR](https://cran.r-project.org/web/packages/hlaR) does cleaning
    and imputation of HLA typings
  - [immunotation](http://bioconductor.org/packages/release/bioc/html/immunotation.html)
    formats HLA typings and can work with haplotype/allele frequencies
- Python:
  - [py-ard](https://github.com/nmdp-bioinformatics/py-ard) reduces HLA
    typing resolution. hlapro actually depends on `py-ard` for some of
    this functionality.
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

## Attribution

Package logo made with [hexmake](https://connect.thinkr.fr/hexmake/);
Icon by
[Freepik](https://www.freepik.com/icon/test-tube_9583096#fromView=search&term=hla+gene&page=2&position=35).
