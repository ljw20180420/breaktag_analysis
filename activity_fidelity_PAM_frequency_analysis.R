{
  options(error = function() {
    traceback(3)
    quit(save = "no", status = 1, runLast = FALSE)
  })

  cache_dir <- "/home/ljw/sdc1/roukos/"
  save_dir <- file.path(
    cache_dir,
    "result",
    "activity_fidelity_PAM_frequency_analysis"
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
    ref <- total_guides |>
      dplyr::filter(Dataset == sprintf("%s_Hiplex1", PoolNum)) |>
      _[[1, "sgRNA"]]
    guides <- total_guides |>
      dplyr::filter(Dataset == sprintf("%s_Hiplex1", PoolNum)) |>
      _[["seq"]]
    targets <- readr::read_rds(file.path(
      cache_dir,
      "result",
      "off-target_analysis",
      sprintf("%s.targets", stem)
    ))
    targets <- targets[!sapply(targets, is.null)]
    fs::dir_create(file.path(save_dir, "relative_activity"))
    ggplot2::ggsave(
      file.path(
        save_dir,
        "relative_activity",
        sprintf("%s.pdf", stem)
      ),
      plot = breakinspectoR::plot_relative_activity(
        targets,
        ref = ref,
        what = "all"
      )
    )
    fs::dir_create(file.path(save_dir, "specificity"))
    ggplot2::ggsave(
      file.path(
        save_dir,
        "specificity",
        sprintf("%s.pdf", stem)
      ),
      plot = breakinspectoR::plot_specificity(targets, ref = ref)
    )

    all_targets <- unlist(as(targets, "GRangesList"))
    fs::dir_create(file.path(save_dir, "offtargets_by_pam"))
    ggplot2::ggsave(
      file.path(
        save_dir,
        "offtargets_by_pam",
        sprintf("%s.pdf", stem)
      ),
      plot = breakinspectoR::plot_offtargets_by_pam(
        all_targets,
        fraction = TRUE
      )
    )
    fs::dir_create(file.path(save_dir, "pam_composition"))
    ggplot2::ggsave(
      file.path(
        save_dir,
        "pam_composition",
        sprintf("%s.pdf", stem)
      ),
      plot = breakinspectoR::plot_pam_composition(all_targets)
    )
    fs::dir_create(file.path(save_dir, "pam_logo"))
    ggplot2::ggsave(
      file.path(
        save_dir,
        "pam_logo",
        sprintf("%s.pdf", stem)
      ),
      plot = breakinspectoR::plot_pam_logo(all_targets)
    )

    fs::dir_create(file.path(save_dir, "position_cutsite"))
    ggplot2::ggsave(
      file.path(
        save_dir,
        "position_cutsite",
        sprintf("%s.pdf", stem)
      ),
      plot = breakinspectoR::plot_position_cutsite(
        targets[[1]],
        guide = guides[[1]],
        pam = "NGG"
      )
    )
    fs::dir_create(file.path(save_dir, "sequence_composition"))
    ggplot2::ggsave(
      file.path(
        save_dir,
        "sequence_composition",
        sprintf("%s.pdf", stem)
      ),
      plot = breakinspectoR::plot_sequence_composition(
        targets[[1]],
        guide = guides[[1]],
        pam = "NGG"
      )
    )
    fs::dir_create(file.path(save_dir, "guide_fidelity"))
    ggplot2::ggsave(
      file.path(
        save_dir,
        "guide_fidelity",
        sprintf("%s.pdf", stem)
      ),
      plot = breakinspectoR::plot_guide_fidelity(
        targets[[1]],
        guide = guides[[1]],
        pam = "NGG"
      )
    )
  }
}
