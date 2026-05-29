#!/usr/bin/env Rscript

{
  cache_dir <- "/home/ljw/sdc1/roukos/"
  genome <- "BSgenome.Hsapiens.UCSC.hg38"
  non_target_file <- breakinspectoR::read_targets(
    file.path(cache_dir, "breaktag_raw_data", "GSM6995157_undig.bed.gz"),
    genome = genome,
    standard_chromosomes = TRUE,
    strandless = FALSE
  )
  total_guides <- readr::read_csv(file.path(
    cache_dir,
    "sgRNA.csv"
  ))
  for (target_bed in readr::read_csv("Hiplex1.txt", col_names = "target_bed")) {
    stem <- stringr::str_split_1(target_bed, stringr::fixed("."))[1]
    targets <- readr::read_rds(file.path(
      cache_dir,
      "result",
      "off-target_analysis",
      sprintf("%s.targets", stem)
    ))
    targets <- targets[!sapply(targets, is.null)]
    all_targets <- unlist(as(targets, "GRangesList"))
    target_file <- breakinspectoR::read_targets(
      file.path(cache_dir, "breaktag_raw_data", target_bed),
      genome = genome,
      standard_chromosomes = TRUE,
      strandless = FALSE
    )
    scission_profile <- breakinspectoR::scission_profile_analysis(
      x = all_targets,
      target = target_file,
      nontarget = non_target_file,
      region = 3,
      bsgenome = genome
    )
    fs::dir_create(file.path(
      cache_dir,
      "result",
      "scission_profile_analysis"
    ))
    readr::write_rds(
      scission_profile,
      file.path(
        cache_dir,
        "result",
        "scission_profile_analysis",
        sprintf("%s.profile", stem)
      )
    )
  }
}
