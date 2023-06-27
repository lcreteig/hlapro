#' Correct the format of HLA alleles
#'
#' `clean_hla()` takes in a string or character vector of HLA alleles, and
#' performs a number of cleaning steps to correct common issues with formatting.
#'
#' `clean_hla()` performs the following operations, in this order:
#'
#' 1. Removes leading or trailing whitespace
#' 2. Adds a leading zero to fields if necessary (`A*1:1` --> `A*01:01`)
#' 3. Removes redundant "versions" of the allele (e.g. a broad when the split is
#' also specified). See [strip_redundant()]
#' 4. Removes punctuation and symbols that are not part of the notation
#' 5. Adds an ":XX" suffix to molecular alleles with only 1 field (`A*01` -->
#' `A*01:XX`)
#' 6. Propagates loci and allele group fields in ambiguities (`A*01:01/02` -->
#' `A*01:01/A*01:02`)
#'
#' @inheritParams get_resolution
#'
#' @return A string or character vector of the same length as `allele`, with the
#' cleaned-up alleles (or the originals, if no cleaning was necessary).
#' @export
#'
#' @examples
#' clean_hla(" A'*1 ")
#' clean_hla("A25(10)")
#' clean_hla("C*03:01/02")
clean_hla <- function(allele) {
  stringr::str_trim(allele) |>
    add_leading_zero() |>
    strip_redundant() |>
    remove_punctuation() |>
    add_xx_suffix() |>
    prefix_ambiguity()
}

#' Strip redundant alleles from typing string if higher resolution is available
#'
#' @description
#' Sometimes a typing will contain multiple "versions" of an allele at different
#' levels of resolution. Most commonly this is both the split and the broad,
#' e.g. `A24(9)` or `A10 A25`. The latter can cause a typing to contain more
#' than 2 alleles for a given locus, which cannot be handled with
#' [extract_alleles_str()]; the former is not accepted by [validate_allele()].
#' Besides, the broads are redundant in this case, and can always be added back
#' with [get_broad()]. Sometimes even three "versions" are included: the split,
#' the broad, and the two-field allele, e.g. `DR5 DR11 DRB1*11:03`.
#'
#' @param typing A string containing the HLA allele or the full space-separated
#' HLA typing.
#'
#' @return A string without the redundant alleles in the input
#' @keywords internal
#' @export
#'
#' @examples
#' strip_redundant("A24(9)")
#' strip_redundant("A9(24)") # also works when the split is in parentheses
#' strip_redundant("A24(9) A10 A25") # removes both A9 and A10
#' strip_redundant("DR5 DR11 DRB1*11:03") # removes both the split and the broad
#' # also works on character vectors
#' strip_redundant(c("A24(9)", "A25(10)"))
strip_redundant <- function(typing) {
  strip_redundant_1 <- function(string) {
    if (is.na(string)) {
      return(string)
    }
    # remove any parentheses and split up the individual alleles
    typing_clean <- stringr::str_replace_all(
      string,
      c(
        "([A-Za-z]+)(\\d+)\\([A-Za-z]*" = "\\1\\2 \\1", # re-insert locus
        "\\)" = ""
      )
    ) |>
      stringr::str_split_1(" ")
    splits <- typing_clean[is_split(typing_clean)] # all splits in typing
    geno <- typing_clean[!is_serology(typing_clean)] # all molecular typings
    # get lower-level equivalents
    broad_with_split <- get_broad(c(splits, geno))
    split_with_higher <- get_split(geno)
    # remove these from the typing
    typing_clean[!(typing_clean %in% c(broad_with_split, split_with_higher))] |>
      stringr::str_flatten(" ") # make into single string again
  }

  purrr::map_chr(typing, strip_redundant_1)
}

#' Add an XX code to single-field molecular HLA alleles
#'
#' If an HLA-allele is of the form `A*01` or `DRB1*05`, add an XX code suffix
#' to it (`A*01:XX`, `DRB1*05:XX`), as the former notation is not valid
#' according to the IPD IMGT/HLA database conventions (and they mean the same in
#' practice). If the allele is not a single-field molecular typing, it is
#' returned as is.
#'
#' @param allele A string or character vector of HLA alleles
#'
#' @return A vector with the same length as `allele`, with all XX-codes added
#' to molecular typings with just one field
#' @noRd
#'
#' @examples
#' add_xx_suffix("A*01")
add_xx_suffix <- function(allele) {
  ifelse(!is_serology(allele) & get_n_fields(allele) == 1,
    make_xx(allele),
    allele
  )
}

#' Propagate locus and allele group in ambiguous HLA typings
#'
#' Sometimes an ambiguous HLA typings is recorded in a short form, where only
#' the first ambiguity contains the locus and/or allele group.
#' `prefix_ambiguity()` adds the locus and/or allele group to the remaining
#' ambiguities as necessary.
#'
#' @param allele A string or character vector of HLA alleles
#'
#' @return A vector with the same length as `allele`, with all ambiguities
#'  written out
#' @noRd
#'
#' @examples
#' prefix_ambiguity("C*01:02/03/04")
prefix_ambiguity <- function(allele) {
  prefix_ambiguity1 <- function(string) {
    pattern <- stringr::regex(r"(
                            (?<locus>\w+\*) # all letters/numbers and "*"
                            (?<allele>\d{1,4}:) # 1-4 digits and ":"
  )", comments = TRUE)

    prefix <- stringr::str_match(string, pattern) # get locus, allele group
    # if allele doesn't follow this format, return as is
    if (any(is.na(prefix))) {
      return(string)
    }
    ambigs <- stringr::str_split_1(string, r"(\/)") # split ambiguities

    # for each ambiguity, store whether it's missing the locus or allele group
    has_no_locus <- stringr::str_detect(ambigs, "\\*", negate = TRUE)
    has_no_allele <- stringr::str_detect(ambigs, ":", negate = TRUE)

    # if necessary, first prefix with allele group, then with locus
    ambigs_fixed <- ifelse(has_no_allele,
      stringr::str_c(prefix[1, "allele"], ambigs),
      ambigs
    )
    ambigs_fixed <- ifelse(has_no_locus,
      stringr::str_c(prefix[1, "locus"], ambigs_fixed),
      ambigs_fixed
    )

    # join all into one string again
    stringr::str_flatten(ambigs_fixed, collapse = "/")
  }
  purrr::map_chr(allele, prefix_ambiguity1)
}

#' Pad single-digit fields in HLA alleles with leading zeros
#'
#' Sometimes the zero in front of a single-digit field drops of, e.g.
#' A*1 instead of A*01. `add_leading_zero()` left-pads all single digits with a
#' 0 (if the digit is preceded by `*` or `:` or `/`, and followed by another `:`
#' or `/` or nothing).
#'
#' @param allele A string or character vector of HLA alleles
#'
#' @return A vector with the same length as `allele`, with zero-padded fields
#' where necessary
#' @noRd
#'
#' @examples
#' add_leading_zero("A*1:01:11:2")
add_leading_zero <- function(allele) {
  stringr::str_replace_all(allele, "(?<=[:\\*\\-/])(\\d)(?=:|/|$)", "0\\1")
}

#' Delete erroneous punctuation and symbols for HLA allele strings
#'
#' Removes symbols (e.g. `$`) and punctuation (e.g. `,`) that should not occur
#' and have no special meaning in HLA allele strings. Makes an exception for
#' punctuation/symbols that can occur in nomenclature (`:` `*` `-` `/`) and
#' GL String syntax (`?` `^` `|` `+` `~` `.` `#`)
#'
#' @param allele A string or character vector of HLA alleles
#'
#' @return A vector with the same length as `allele`, with punctuation removed
#' where necessary
#' @noRd
#'
#' @examples
#' remove_punctuation("HLA-A'*01:;01,")
remove_punctuation <- function(allele) {
  pattern <- stringr::str_glue("(?!{legal_charset}){punc_symbs}",
    legal_charset = r"([:\*\-/\?\^\|\+~\.#])",
    punc_symbs = r"([[\p{P}][\p{S}]])"
  )
  stringr::str_remove_all(allele, pattern)
}
