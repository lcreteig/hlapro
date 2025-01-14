# TODO: add test for phased genotype frequencies (ground truth from HaploStats)
test_that("serological typings return expected high-res EURCAU genotypes", {
  # Compared to PROCARE 1.0 data

  # NMDP data cannot be made available; users have to ID and accept license
  skip_on_cran()
  skip_on_ci()

  input_typings <- c(
    # N.B. HaploStats for some reason does not accept DR11, but gives same
    # result when DR11 is swapped for DR5 in A24 A28 B35 B61 DR4 DR11
    paste0(
      "hla#2024#",
      # A24
      "A*24:02/A*24:03/A*24:04/A*24:05/A*24:06/A*24:07/A*24:08/A*24:10/",
      "A*24:13/A*24:14/A*24:15/A*24:17/A*24:20/A*24:22/A*24:25/A*24:26/",
      "A*24:29/A*24:30/A*24:31/A*24:32/A*24:35/A*24:37/A*24:38/A*24:43/",
      "A*24:47/A*24:53/A*24:56/A*24:57/A*24:58/A*24:64/A*24:72/A*24:81/",
      "A*24:95+",
      # A28
      "A*68:01/A*68:02/A*68:03/A*68:05/A*68:06/A*68:07/A*68:08/A*68:10/",
      "A*68:12/A*68:15/A*68:17/A*68:20/A*68:22/A*68:24/A*68:25/A*68:27/",
      "A*68:30/A*68:31/A*68:35/A*68:37/A*68:40/A*69:01^",
      # B35
      "B*35:01/B*35:02/B*35:03/B*35:04/B*35:05/B*35:06/B*35:08/B*35:09/",
      "B*35:10/B*35:11/B*35:12/B*35:13/B*35:14/B*35:15/B*35:16/B*35:17/",
      "B*35:18/B*35:19/B*35:20/B*35:21/B*35:22/B*35:23/B*35:24/B*35:25/",
      "B*35:26/B*35:27/B*35:28/B*35:29/B*35:30/B*35:31/B*35:32/B*35:34/",
      "B*35:37/B*35:41/B*35:43/B*35:47/B*35:48/B*35:49/B*35:55/B*35:77+",
      # B61
      "B*40:02/B*40:03/B*40:04/B*40:06/B*40:09/B*40:11/B*40:16/B*40:20/",
      "B*40:27/B*40:50^",
      # DR4
      "DRB1*04:01/DRB1*04:02/DRB1*04:03/DRB1*04:04/DRB1*04:05/DRB1*04:06/",
      "DRB1*04:07/DRB1*04:08/DRB1*04:09/DRB1*04:10/DRB1*04:11/DRB1*04:12/",
      "DRB1*04:13/DRB1*04:14/DRB1*04:15/DRB1*04:16/DRB1*04:17/DRB1*04:18/",
      "DRB1*04:19/DRB1*04:22/DRB1*04:23/DRB1*04:25/DRB1*04:26/DRB1*04:33/",
      "DRB1*04:35/DRB1*04:38/DRB1*04:40/DRB1*04:41/DRB1*04:50/DRB1*04:51/",
      "DRB1*04:54^",
      # DR11
      "DRB1*11:01/DRB1*11:02/DRB1*11:03/DRB1*11:04/DRB1*11:05/DRB1*11:06/",
      "DRB1*11:07/DRB1*11:08/DRB1*11:09/DRB1*11:10/DRB1*11:11/DRB1*11:12/",
      "DRB1*11:13/DRB1*11:14/DRB1*11:15/DRB1*11:16/DRB1*11:17/DRB1*11:18/",
      "DRB1*11:19/DRB1*11:20/DRB1*11:24/DRB1*11:25/DRB1*11:28/DRB1*11:29/",
      "DRB1*11:32/DRB1*11:34/DRB1*11:36/DRB1*11:37/DRB1*11:39/DRB1*11:42/",
      "DRB1*11:43/DRB1*11:45/DRB1*11:47/DRB1*11:56/DRB1*11:66/DRB1*11:69"
    ),
    paste0(
      "hla#2024#",
      # A2
      "A*02:01/A*02:02/A*02:03/A*02:04/A*02:05/A*02:06/A*02:07/A*02:08/",
      "A*02:10/A*02:11/A*02:119/A*02:12/A*02:122/A*02:123/A*02:13/",
      "A*02:137/A*02:14/A*02:16/A*02:17/A*02:18/A*02:19/A*02:20/A*02:21/",
      "A*02:22/A*02:24/A*02:25/A*02:27/A*02:29/A*02:30/A*02:33/A*02:34/",
      "A*02:35/A*02:36/A*02:38/A*02:39/A*02:44/A*02:45/A*02:49/A*02:58/",
      "A*02:60/A*02:63/A*02:64/A*02:74/A*02:84/A*02:85/A*02:86/A*02:87/",
      "A*02:93/A*02:96+",
      # A3
      "A*03:01/A*03:02/A*03:05/A*03:06/A*03:07/A*03:08/A*03:09/A*03:108/",
      "A*03:17/A*03:22/A*03:50^",
      # B52
      "B*52:01/B*52:02/B*52:03/B*52:04/B*52:21+",
      # B35
      "B*35:01/B*35:02/B*35:03/B*35:04/B*35:05/B*35:06/B*35:08/B*35:09/",
      "B*35:10/B*35:11/B*35:12/B*35:13/B*35:14/B*35:15/B*35:16/B*35:17/",
      "B*35:18/B*35:19/B*35:20/B*35:21/B*35:22/B*35:23/B*35:24/B*35:25/",
      "B*35:26/B*35:27/B*35:28/B*35:29/B*35:30/B*35:31/B*35:32/B*35:34/",
      "B*35:37/B*35:41/B*35:43/B*35:47/B*35:48/B*35:49/B*35:55/B*35:77^",
      # Cw4
      "C*04:01/C*04:03/C*04:04/C*04:05/C*04:06/C*04:07/C*04:08/C*04:10/",
      "C*04:13/C*04:19/C*04:27^",
      # DR11
      "DRB1*11:01/DRB1*11:02/DRB1*11:03/DRB1*11:04/DRB1*11:05/DRB1*11:06/",
      "DRB1*11:07/DRB1*11:08/DRB1*11:09/DRB1*11:10/DRB1*11:11/DRB1*11:12/",
      "DRB1*11:13/DRB1*11:14/DRB1*11:15/DRB1*11:16/DRB1*11:17/DRB1*11:18/",
      "DRB1*11:19/DRB1*11:20/DRB1*11:24/DRB1*11:25/DRB1*11:28/DRB1*11:29/",
      "DRB1*11:32/DRB1*11:34/DRB1*11:36/DRB1*11:37/DRB1*11:39/DRB1*11:42/",
      "DRB1*11:43/DRB1*11:45/DRB1*11:47/DRB1*11:56/DRB1*11:66/DRB1*11:69^",
      # DR52
      "DRB3*01:01/DRB3*01:02/DRB3*01:03/DRB3*02:01/DRB3*02:02/DRB3*02:03/",
      "DRB3*02:06/DRB3*02:10/DRB3*02:11/DRB3*02:17/DRB3*03:01^",
      # DQ3
      "DQB1*03:01/DQB1*03:02/DQB1*03:03/DQB1*03:04/DQB1*03:05/DQB1*03:13"
    )
  )

  res <- upscale_typings(
    filepath = "~/Downloads/A~C~B~DRB3-4-5~DRB1~DQB1.xlsx",
    typing = input_typings,
    loci_input = c("A", "B", "DRB1", "DQB1"),
    n_haplos = 5000
  ) |>
    dplyr::group_by(id_input_typing) |>
    dplyr::slice_max(unphased_prob, n = 1, with_ties = FALSE) |>
    dplyr::mutate(unphased_geno = stringr::str_remove(
      unphased_geno,
      "^hla#[\\d-]+#"
    ))

  # frequency
  expect_equal(dplyr::pull(res, unphased_freq),
    c(2.542e-07, 5.13531e-08),
    tolerance = 1e-03
  )
  # genotype
  expect_equal(
    dplyr::pull(res, unphased_geno),
    c(
      paste0(
        "HLA-A*24:02g+HLA-A*68:01g^HLA-B*35:03g+HLA-B*40:02g^",
        "HLA-C*02:02g+HLA-C*04:01g^HLA-DQB1*03:01g+HLA-DQB1*03:02g^",
        "HLA-DRB1*04:03+HLA-DRB1*11:01g^HLA-DRB3*02:02g+HLA-DRB4*01:01g"
      ),
      paste0(
        "HLA-A*02:01g+HLA-A*03:01g^HLA-B*35:01g+HLA-B*52:01g^",
        "HLA-C*04:01g+HLA-C*12:02^HLA-DQB1*03:01g+HLA-DQB1*03:01g^",
        "HLA-DRB1*11:01g+HLA-DRB1*11:04^HLA-DRB3*02:02g+HLA-DRB3*02:02g"
      )
    )
  )
})

test_that("homozygous low-res genotype returns homozygous high-res genotypes", {
  # Compared to HaploStats output

  # NMDP data cannot be made available; users have to ID and accept license
  skip_on_cran()
  skip_on_ci()

  input_typings <- c(
    paste0(
      "hla#2024#",
      # A2
      "A*02:01/A*02:02/A*02:03/A*02:04/A*02:05/A*02:06/A*02:07/A*02:08/",
      "A*02:10/A*02:11/A*02:119/A*02:12/A*02:122/A*02:123/A*02:13/",
      "A*02:137/A*02:14/A*02:16/A*02:17/A*02:18/A*02:19/A*02:20/A*02:21/",
      "A*02:22/A*02:24/A*02:25/A*02:27/A*02:29/A*02:30/A*02:33/A*02:34/",
      "A*02:35/A*02:36/A*02:38/A*02:39/A*02:44/A*02:45/A*02:49/A*02:58/",
      "A*02:60/A*02:63/A*02:64/A*02:74/A*02:84/A*02:85/A*02:86/A*02:87/",
      "A*02:93/A*02:96^",
      # B7
      "B*07:02/B*07:04/B*07:05/B*07:07/B*07:08/B*07:09/B*07:10/B*07:12/",
      "B*07:13/B*07:14/B*07:15/B*07:16/B*07:20/B*07:23/B*07:26/B*07:36/",
      "B*07:37/B*07:42/B*07:43/B*07:46/B*07:51^",
      # DR15
      "DRB1*15:01/DRB1*15:02/DRB1*15:03/DRB1*15:04/DRB1*15:06/DRB1*15:07/",
      "DRB1*15:10/DRB1*15:11/DRB1*15:14/DRB1*15:18/DRB1*15:20/DRB1*15:22/",
      "DRB1*15:23/DRB1*15:24/DRB1*15:38^",
      # DQ1
      "DQB1*05:01/DQB1*05:02/DQB1*05:03/DQB1*05:04/DQB1*06:01/DQB1*06:02/",
      "DQB1*06:03/DQB1*06:04/DQB1*06:05/DQB1*06:08/DQB1*06:09/DQB1*06:10/",
      "DQB1*06:11"
    ),
    paste0(
      "hla#2024#",
      # A1
      "A*01:01/A*01:02/A*01:03/A*01:06/A*01:09/A*01:12/A*01:17/A*01:25^",
      # B8
      "B*08:01/B*08:02/B*08:03/B*08:04/B*08:07/B*08:09/B*08:13/B*08:18/",
      "B*08:20/B*08:23/B*08:35^",
      # Cw7
      "C*07:01/C*07:02/C*07:04/C*07:05/C*07:07/C*07:10/C*07:12/C*07:13/",
      "C*07:17/C*07:19/C*07:21/C*07:22/C*07:24/C*07:25/C*07:26/C*07:27/",
      "C*07:29/C*07:35/C*07:40/C*07:43/C*07:46/C*07:56/C*07:60/C*07:67/",
      "C*07:72^",
      # DR17
      "DRB1*03:01/DRB1*03:04",
      # DR52
      "DRB3*01:01/DRB3*01:02/DRB3*01:03/DRB3*02:01/DRB3*02:02/DRB3*02:03/",
      "DRB3*02:06/DRB3*02:10/DRB3*02:11/DRB3*02:17/DRB3*03:01^",
      # DQ2
      "DQB1*02:01/DQB1*02:03"
    )
  )

  res <- upscale_typings(
    filepath = "~/Downloads/A~C~B~DRB3-4-5~DRB1~DQB1.xlsx",
    loci_input = c("A", "B", "C", "DRB1", "DQB1"),
    typing = input_typings
  ) |>
    dplyr::group_by(id_input_typing) |>
    dplyr::slice_max(unphased_prob, n = 1, with_ties = FALSE) |>
    dplyr::mutate(unphased_geno = stringr::str_remove(
      unphased_geno,
      "^hla#[\\d-]+#"
    ))

  # frequency
  expect_equal(dplyr::pull(res, unphased_freq),
    c(3.896e-04, 4.259e-03),
    tolerance = 1e-03
  )
  # genotype
  expect_equal(
    dplyr::pull(res, unphased_geno),
    c(
      paste0(
        "HLA-A*02:01g+HLA-A*02:01g^HLA-B*07:02g+HLA-B*07:02g^",
        "HLA-C*07:02g+HLA-C*07:02g^HLA-DQB1*06:02+HLA-DQB1*06:02^",
        "HLA-DRB1*15:01+HLA-DRB1*15:01^HLA-DRB5*01:01+HLA-DRB5*01:01"
      ),
      paste0(
        "HLA-A*01:01g+HLA-A*01:01g^HLA-B*08:01g+HLA-B*08:01g^",
        "HLA-C*07:01g+HLA-C*07:01g^HLA-DQB1*02:01g+HLA-DQB1*02:01g^",
        "HLA-DRB1*03:01+HLA-DRB1*03:01^HLA-DRB3*01:01+HLA-DRB3*01:01"
      )
    )
  )
})
