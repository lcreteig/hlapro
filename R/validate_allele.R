validate_allele <- function(allele) {
  pattern <- stringr::regex(r"(
    ^                  # START of allele string
    (?:HLA-)?          # optional "HLA-" prefix
    (?<locus>          # START of loci
    A|Bw?|Cw?|         # class I loci (incl. "Cw" and "Bw")
    DRA|DR(?:B[1-9])?| # DRA, DR/DRB1-9
    DQA1?|DQ(?:B1)?|   # DQA/DQA1, DQ/DQB1
    DPA1?|DP(?:B1)?    # DPA/DPA1, DP/DPB1
    )                  # END of loci
    [\*-]?             # * (A*01) or - (DP-01), optional
    (?<allele>         # START of allele field
    \d{1,4}            # 1 (A1) to 4 (DP0201) digits
    )                  # END of allele field
                       # ALL that follows is optional
    :?                 # optional semicolon
    (?<protein>        # START of optional protein field
    [A-Z]{2,5}$|       # a 2-5 letter MAC code, then end. OR
    \d{2,3}P?          # 2-3 digits, with optional P group
    )?                 # END of protein field
    (?<coding>         # START of optional 3rd field (coding)
    :                  # semicolon
    \d{2,3}G?          # 2-3 digits, with optional G group
    )?                 # end of 3rd field
    (?<noncoding>      # START of optional 4th field (non-coding)
    :                  # semicolon
    \d{2,3}            # 2-3 digits
    )?                 # END of optional 4th field
    (?<suffix>         # START of optional suffix
    [NLSCAQ]           # an N,L,S,C,A, or Q
    )?                 # end of optional suffix
                       # ALL that follows are optional ambiguities (0 or more)
    (\/                # a forward slash
    (?:\1\*)?          # followed by optional locus (must be same as before)
    \d{2,3}            # 2-3 digits
    (?:                # START of protein-level information
    :                  # a semicolon
    \d{2,3}            # 2-3 digits
    [NLSCAQ]?          # optional suffix
    )*                 # END of protein level (repeat 0 or more times)
    )*                 # END of ambiguities (repeat 0 or more times)
    $                  # END of HLA string
  )", comments = TRUE)

  stringr::str_detect(allele, pattern)
}
