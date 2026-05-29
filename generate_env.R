library(rix)

rix(
  r_ver = "frozen-edge",
  r_pkgs = c(
    "languageserver",
    "devtools",
    "tidyverse",
    "BSgenome.Hsapiens.UCSC.hg38"
  ),
  system_pkgs = NULL,
  git_pkgs = list(
    list(
      package_name = "breakinspectoR",
      repo_url = "https://github.com/roukoslab/breakinspectoR/",
      commit = "c27574244d2eee2bd79a41a817b91636253734af"
    )
  ),
  ide = "none",
  project_path = ".",
  overwrite = TRUE,
  print = TRUE
)
