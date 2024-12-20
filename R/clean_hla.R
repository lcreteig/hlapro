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
#' 7. Converts v2 to v3 (`A*01010102N` --> `A*01:01:01:02N`). See
#' [convert_v2_to_v3()]
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
    prefix_ambiguity() |>
    convert_v2_to_v3()
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
  ifelse(!is_serology(allele) & !is_v2(allele) & get_n_fields(allele) == 1,
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
#' @param return_v3 If `TRUE` will use [convert_v2_to_v3()] to force all
#'   ambiguities to modern nomenclature
#'
#' @return A vector with the same length as `allele`, with all ambiguities
#'  written out
#' @noRd
#'
#' @examples
#' prefix_ambiguity("C*01:02/03/04")
prefix_ambiguity <- function(allele, return_v3 = TRUE) {
  prefix_ambiguity1 <- function(string) {
    # if not ambiguous or empty, return as is
    if (is.na(string) || stringr::str_detect(string, r"(\/)", negate = TRUE)) {
      return(string)
    }
    pattern <- stringr::regex(r"(
                            (?<locus>\w+\*) # all letters/numbers and "*"
                            # 1-4 digits & ":" (v3), or first 2 of 4 digits (v2)
                            (?<allele>\d{1,4}:|\d{2}(?=\d{2}))
    )", comments = TRUE)

    prefix <- stringr::str_match(string, pattern) # get locus, allele group
    # if allele doesn't follow this format, return as is
    if (any(is.na(prefix))) {
      return(string)
    }
    ambigs <- stringr::str_split_1(string, r"(\/)") # split ambiguities
    v2_flag <- is_v2(ambigs[1])

    # for each ambiguity, store whether it's missing the locus or allele group
    has_no_locus <- stringr::str_detect(ambigs, "\\*", negate = TRUE)
    has_no_group <- !v2_flag & stringr::str_detect(ambigs, ":", negate = TRUE)
    # for v2: allele group is missing if the whole ambiguity is just 2-3 digits,
    # (and an optional suffix)
    pattern_v2 <- "^\\d{2,3}[NLSCAQPG]?$"
    has_no_group_v2 <- v2_flag & stringr::str_detect(ambigs, pattern_v2)

    # if necessary, first prefix with allele group, then with locus
    ambigs_fixed <- ifelse(has_no_group | has_no_group_v2,
      stringr::str_c(prefix[1, "allele"], ambigs),
      ambigs
    )
    ambigs_fixed <- ifelse(has_no_locus,
      stringr::str_c(prefix[1, "locus"], ambigs_fixed),
      ambigs_fixed
    )

    if (return_v3) {
      ambigs_fixed <- convert_v2_to_v3(ambigs_fixed)
    }
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

#' Translates HLA alleles in v2 notation to v3
#'
#' The pre-2010 "v2" notation does not include the field delimiters (`:`) that
#' are now mandatory in v3. This function first tests if an allele is in v2
#' format; if an allele is not in v2 format; it's left alone. But if it is, it
#' looks up its v2 equivalent in the [v2_to_v3] lookup table. If it is not in
#' the table, a v3 version is put together heuristically, by inserting `:` after
#' every two digits.
#'
#' N.B. The heuristic prediction will not work in all cases. For example:
#' - `DPB1*87801N` should be `DPB1*878:01N` (but is output as `DPB1*87:801N`)
#' - `DPB1*152401` should be `DPB1*1524:01` (but is output as `DPB1*15:24:01`)
#'
#' In general it's not possible to make this work purely syntactically without
#' imbuing knowledge on which HLA alleles exist and which do not. For example,
#' should `DRB1*1412601` be `14:126:01` or `14:12:601`? Both are theoretically
#' possible. However, alleles with > 2 digits per field are rare, and were not
#' really around before 2010, so in practice one should rarely encounter them in
#' v2 format.
#'
#' @param allele A string or character vector of HLA alleles
#'
#' @return A vector with the same length as `allele`, with all v2 alleles
#' converted to v3
#' @keywords internal
#' @export
#'
#' @examples
#' convert_v2_to_v3("A*01010101") # known v2 allele
#' convert_v2_to_v3("B*0701") # not a known v2 allele, but heuristic works
#' convert_v2_to_v3("B*9526") # known allele where heuristic would not work
convert_v2_to_v3 <- function(allele) {
  # replace v2 with v3 from lookup table
  v3s <- ifelse(is_v2(allele), unname(lookup_v3[allele]), allele)
  # if not in table
  ifelse(is.na(v3s),
    stringr::str_replace_all(
      allele,
      c(
        "Cw" = "C", # replace Cw with C
        # insert ":" after every 2 digits if followed by 2 digits/letters
        r"((\d{2})(?=\d{2}|[A-Z]{2}))" = "\\1:"
      )
    ),
    v3s
  )
}
