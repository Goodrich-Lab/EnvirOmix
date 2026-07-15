#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom tibble tibble
## usethis namespace: end

.require_intermediate_packages <- function() {
  missing_pkgs <- c("epiomics", "xtune")[
    !vapply(c("epiomics", "xtune"), requireNamespace, logical(1), quietly = TRUE)
  ]

  if (length(missing_pkgs) > 0) {
    stop(
      paste0(
        "Intermediate integration requires optional package(s): ",
        paste(missing_pkgs, collapse = ", "), ". ",
        "Install them with one of the following approaches:\n",
        "  remotes::install_github('JingxuanH/xtune')\n",
        "  install.packages('https://cran.r-project.org/src/contrib/Archive/xtune/xtune_2.0.0.tar.gz', repos = NULL, type = 'source')\n",
        "  install.packages('epiomics')"
      ),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

.require_analysis_packages <- function(packages, feature) {
  missing_pkgs <- packages[
    !vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  ]

  if (length(missing_pkgs) > 0) {
    stop(
      paste0(
        feature, " requires optional package(s): ",
        paste(missing_pkgs, collapse = ", "), ". ",
        "Install them with install.packages() before running this analysis."
      ),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

.resolve_hima_index <- function(tbl, default_col = "lf_num") {
  if ("Index" %in% names(tbl)) {
    tbl[[default_col]] <- tbl$Index
  }

  tbl
}

NULL
