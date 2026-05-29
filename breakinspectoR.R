#!/usr/bin/env Rscript

{
  cache_dir <- "/home/ljw/sdc1/roukos/"
  genome <- "BSgenome.Hsapiens.UCSC.hg38"
  # non_target_file <- breakinspectoR::read_targets(
  #   file.path(cache_dir, "breaktag_raw_data", "GSM6995157_undig.bed.gz"),
  #   genome = genome,
  #   standard_chromosomes = TRUE,
  #   strandless = TRUE
  # )
  total_guides <- readr::read_csv(file.path(
    cache_dir,
    "sgRNA.csv"
  ))
  for (target_bed in readr::read_csv("Hiplex1.txt", col_names = "target_bed")[[
    "target_bed"
  ]]) {
    browser()
    stem <- stringr::str_split_1(target_bed, stringr::fixed("."))[1]
    PoolNum <- stringr::str_split_1(
      stem,
      stringr::fixed("_")
    )[2]
    guides <- total_guides |>
      dplyr::filter(Dataset == sprintf("%s_Hiplex1", PoolNum)) |>
      dplyr::select(seq) |>
      _[["seq"]]
    target_file <- breakinspectoR::read_targets(
      file.path(cache_dir, "breaktag_raw_data", target_bed),
      genome = genome,
      standard_chromosomes = TRUE,
      strandless = TRUE
    )
    targets <- parallel::mclapply(
      guides,
      function(guide) {
        tryCatch(
          {
            x <- breakinspectoR::breakinspectoR(
              target = target_file,
              nontarget = non_target_file,
              guide = guide,
              PAM = "NNN",
              bsgenome = genome,
              cutsiteFromPAM = 3,
              min_breaks = 8,
              eFDR = TRUE,
              verbose = FALSE
            )
            if (length(x) > 0) {
              x$qval <- ifelse(is.na(x$qval), x$FDR, x$qval)
              breakinspectoR::reduceOT(x, verbose = FALSE)
            }
          },
          error = function(e) x
        )
      },
      mc.cores = 3
    )

    dir_create(file.path(cache_dir, "result", "breakinspectoR"))
    readr::write_rds(
      targets,
      file.path(
        cache_dir,
        "result",
        "breakinspectoR",
        sprintf("%s.targets", stem)
      )
    )
  }
}
