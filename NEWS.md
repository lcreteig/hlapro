# hlapro (development version)

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
