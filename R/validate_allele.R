#' Check whether an HLA typing is well-formed
#'
#' `validate_allele()` takes in a character vector or string of HLA alleles, and
#' returns `TRUE` if the allele is well-formed, and `FALSE` if it isn't.
#'
#' N.B. This function does *not* test whether an allele actually exists (e.g.
#' whether it occurs in the most recent version of the IPD-IMGT/HLA database),
#' but only whether it's string representation conforms to certain standards. An
#' allele can be well-formed but not exist (e.g. `"A*99:01:01"`), or can exist
#' but not be well-formed (e.g. `"HLA-A**02;01"`).
#'
#' The following are explicitly considered valid HLAs:
#'
#' - alleles belonging to loci A, B, C, DRA/DRB, DQA/DQB, DPA/DPB (which is
#'   *not* an exhaustive list of class I or class II HLAs, but simply the ones
#'   that are often typed in the context of transplantation research/matching)
#' - serological/antigen notation, such as `"A2"`, `"DP-0201"`
#' - XX codes, e.g. `"A*02:XX"`
#' - prefixing an allele with `"HLA-"` is allowed
#' - ambiguous alleles such as `"C*01:02/C*01:03/C*01:04/C*01:05/C*01:06"`
#' - Multiple Allele Codes, e.g. `"DRB1*07:GC"` (v3) or `"DPB1*04BDVU"` (v2)
#' - P groups, G groups, and expression-related suffixes (N/L/S/C/A/Q)
#'
#' @inheritParams get_resolution
#'
#' @return A Boolean or logical vector with the same lengths as `allele`, with
#'  `TRUE` or `FALSE` for each element.
#' @export
#'
#' @examples
#' validate_allele("A2")
#' validate_allele("A*99:01:01") # well-formed but non-existing
#' validate_allele("HLA-A**02;01") # existing but not well-formed
#'
#' # also works with character vectors, or in a data frame
#' allele_vec <- c("A2", "A*01:AABJE", "A*24:02:01:02L", "not-an-HLA")
#' validate_allele(allele_vec)
#'
#' df <- tidyr::tibble(alleles = allele_vec)
#' dplyr::mutate(df, alleles_check = validate_allele(alleles))
validate_allele <- function(allele) {
  pattern <- stringr::regex(r"(
    ^                    # START of allele string
    (?:HLA-)?            # optional "HLA-" prefix
    (?<locus>            # START of loci
    A|Bw?|Cw?|           # class I loci (incl. "Cw" and "Bw")
    DRA(?=\*)|           # DRA*
    DR(?:B[1-9](?=\*))?| # DR, or DRB1-9*
    DQA(?=-)|DQA1(?=\*)| # DQA- or DQA1*
    DQ(?:B1(?=\*))?|     # DQ, or DQB1*,
    DPA(?=-)|DPA1(?=\*)| # DPA- or DPA1*
    DP(?=-)|DPB1(?=\*)   # DP- or DPB1*
    )                    # END of loci
    [\*-]?               # * (A*01) or - (DP-01), optional
    (?<allele>           # START of allele field
    \d{1,4}              # 1 (A1) to 4 (DP0201) digits
    )                    # END of allele field
    (?=$|[:A-Z\/])       # either end of string (low-res), or a ":"/capital/"/"
                         # ALL that follows is optional

    :?                   # optional colon
    (?<protein>          # START of optional protein field
    [A-Z]{2,5}$|         # a 2-5 letter MAC code, then end. OR
    \d{2,3}P?            # 2-3 digits, with optional P group
    )?                   # END of protein field
    (?<coding>           # START of optional 3rd field (coding)
    :                    # colon
    \d{2,3}G?            # 2-3 digits, with optional G group
    )?                   # end of 3rd field
    (?<noncoding>        # START of optional 4th field (non-coding)
    :                    # colon
    \d{2,3}              # 2-3 digits
    )?                   # END of optional 4th field
    (?<suffix>           # START of optional suffix
    [NLSCAQ]             # an N,L,S,C,A, or Q
    )?                   # end of optional suffix

                         # ALL that follows are optional ambiguities (0 or more)
    (\/                  # a forward slash
    (?:\1\*)?            # followed by optional locus (must be same as before)
    \d{2,3}              # 2-3 digits (first field)
    (?::\d{2,3}){0,3}    # 0-3 more fields with 2-3 digits each, colon-separated
    [NLSCAQ]?            # optional suffix
    )*                   # END of ambiguities (repeat 0 or more times)
    $                    # END of HLA string
  )", comments = TRUE)

  # if there's a ":", there must also be a *
  !(stringr::str_detect(allele, ":") & !stringr::str_detect(allele, "\\*")) &
    stringr::str_detect(allele, pattern)
}
