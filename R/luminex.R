determine_positive_beads <- function(specificities, mfi,
                                     ratio_cutoff = 15, mfi_cutoff = 500) {
  loci <- get_loci(specificities) # get locus of each bead

  lra_mfis <- split(mfi, loci) |> # make one vector of MFI values per locus
    purrr::map_dbl(min) # get the lowest MFI for each

  # divide MFI of each bead by the MFI of the Lowest Ranked Antigen in its locus
  ratios <- mfi / lra_mfis[loci]

  # a bead is positive if both the ratio and raw values exceed the cutoff
  ifelse(ratios > ratio_cutoff & mfi > mfi_cutoff, "positive", "negative") |>
    unname()
}

determine_dsa <- function(beads, typing,
                          level = c("split", "broad"), sep = "\\s") {
  # TODO: This works for now, but probably makes more sense if typing is
  # also a vector of alleles. Then DSA is simply: beads == typing & positive.
  level <- arg_match(level)
  convert_typing <- switch(level,
    split = get_serology,
    broad = get_broad
  )

  typing_alleles <- stringr::str_split_1(typing[1], pattern = sep) |>
    convert_typing()

  ifelse(convert_typing(beads) %in% typing_alleles, "yes", "no")
}
