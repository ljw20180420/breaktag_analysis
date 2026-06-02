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
      commit = "17854efd0be31e28077c042d42580f85ffed9586"
    )
  ),
  ide = "none",
  project_path = ".",
  overwrite = TRUE,
  print = TRUE
)
