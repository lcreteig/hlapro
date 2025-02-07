# hlapro (development version)

* Update `clean_hla()` to also convert 
  [deleted/changed](https://hla.alleles.org/alleles/deleted.html) alleles
* Minor fixes to `clean_hla()` conversion of v2 alleles

# hlapro 0.4.0

* Update `clean_hla()` to also convert [pre-2010 
  typings](https://hla.alleles.org/nomenclature/naming_2010.html) (v2) 
  to the current nomenclature (v3)
* Fix `load_eplet_registry()` for new release of HLA Eplet registry (2024-08-30)
* Wrap the [`py-ard`](https://github.com/nmdp-bioinformatics/py-ard) python
  package to:
    - Encode/lookup Multiple Allele Codes (`mac_lookup()`, `mac_expand()`)
    - Reduce alleles to serological equivalents (`reduce_to_serology()`)
    - Reduce alleles to 2-field (ARD) level (`reduce_to_field2()`)
    
# hlapro 0.3.0

* New `read_lum_csv()` function to parse raw Luminex Single-Antigen Bead assay
  data
* Fixed tests for new release of HLA Eplet registry (2024-08-19)

# hlapro 0.2.3

* Fixed tests for new release of HLA Eplet registry (2024-03-15)

# hlapro 0.2.2

* Updated for new release of the HLA Eplet registry: incorporating papers by 
  [Bezstarosti et al.](https://doi.org/10.3389/fimmu.2021.800946) 
  and [Osoegawa et al.](https://doi.org/10.1111/tan.14662) (#11)

# hlapro 0.2.1

* Fixed check for `input_loci` argument in `upscale_typings()` (#10).

# hlapro 0.2.0

* New functions to lookup eplets/alleles in the 
  [HLA Eplet registry](https://www.epregistry.com.br). Updated for release 
  "2024-01-22: IPD-IMGT/HLA 3.54" of the registry (#9)
  
  - `load_eplet_registry()` scrapes the database from the website, or loads an
    existing version from disk
  - `lookup_alleles()` and `lookup_eplets()` return the eplets that occur on
    an input vector of HLA alleles, or vice versa.
  - `get_positive_eplets()` returns eplets that occur on positive (but not 
    negative) beads for a dataframe of Luminex results.

* New `upscale_typings()` function to impute high-resolution (two-field) 
  genotypes for a low resolution serological input typing, based on [haplotype
  frequencies released by the NMDP](http://frequency.nmdp.org)

* New `gl_to_df()` and `df_to_gl()` functions to convert between a data frame
  with one column per allele and a [GL String](https://glstring.org) containing
  all alleles.

* New `clean_hla()` function to correct common formatting issues in HLA allele
  typing strings.

* New functions to get serological equivalents of an allele from the 
  [ETRL HLA tables](https://etrl.eurotransplant.org). These tables are also
  shipped with the package as the `etrl_hla` data frame:
  
  - `get_serology()` for the split or else the broad-level serology
  - `get_broad()` for the broad-level serology
  - `get_split()` for the split-level serology
  - `get_public()` for the Bw4 or Bw6 epitope

* New `validate_alleles()` to check whether a(n) (list of) HLA allele(s) is
  well-formed

* New `get_resolution()` to determine resolution of a(n) (list of) HLA allele(s)
  as either low, intermediate, or high.

* New `extract_alleles_str()`, `extract_alleles_df()` gets alleles for each 
  locus from HLA typing string, and separates them into named list elements 
  / new columns in a data frame.
  
  - Use `count_alleles()` to inspect the number of alleles per locus. The
  `extract_alleles()_*` functions also use `count_alleles()` to warn whenever 
  they encounter a typing with more than two alleles per locus.

# hlapro 0.1.0

* New `get_mismatches()` determines mismatched antigens from a pair of 
  donor/recipient typings.
