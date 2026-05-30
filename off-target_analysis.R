#!/usr/bin/env Rscript

{
  options(error = function() {
    traceback(3)
    quit(save = "no", status = 1, runLast = FALSE)
  })

  cache_dir <- "/home/ljw/sdc1/roukos/"
  genome <- "BSgenome.Hsapiens.UCSC.hg38"
  save_dir <- file.path(cache_dir, "result", "off-target_analysis")
  load_dir <- file.path(
    cache_dir,
    "breaktag_raw_data"
  )
  fs::dir_create(save_dir)

  non_target_file_str <- file.path(
    load_dir,
    "GSM6995157_undig.bed.gz"
  )
  cat(sprintf("load non target %s\n", non_target_file_str))
  non_target_file <- breakinspectoR::read_targets(
    non_target_file_str,
    genome = genome,
    standard_chromosomes = TRUE,
    strandless = TRUE
  )
  total_guides <- readr::read_csv(file.path(
    cache_dir,
    "sgRNA.csv"
  ))
  for (target_bed in readr::read_csv("Hiplex1.txt", col_names = "target_bed")[[
    "target_bed"
  ]]) {
    stem <- stringr::str_split_1(target_bed, stringr::fixed("."))[1]
    PoolNum <- stringr::str_split_1(
      stem,
      stringr::fixed("_")
    )[2]
    guides <- total_guides |>
      dplyr::filter(Dataset == sprintf("%s_Hiplex1", PoolNum)) |>
      _[["seq"]]

    target_file_str <- file.path(load_dir, target_bed)
    cat(sprintf("load target %s\n", target_file_str))
    target_file <- breakinspectoR::read_targets(
      target_file_str,
      genome = genome,
      standard_chromosomes = TRUE,
      strandless = TRUE
    )

    cat(sprintf("process %s\n", stem))
    targets <- list()
    for (idx in seq_along(guides)) {
      guide <- guides[idx]
      cat(sprintf("process %d: %s\n", idx, guide))
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
        x$qval <- ifelse(is.na(x$qval), x$fdr, x$qval)
        breakinspectoR::reduceOT(x, verbose = FALSE)
      }
      targets[[idx]] <- x
    }

    cat(sprintf("save %s\n", stem))
    readr::write_rds(
      targets,
      file.path(
        save_dir,
        sprintf("%s.targets", stem)
      )
    )
  }
}
