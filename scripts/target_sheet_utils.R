read_target_sheet <- function(path) {
  if (!file.exists(path)) {
    stop("Target sheet not found: ", path)
  }

  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]
  lines <- lines[!grepl("^\\s*#", lines)]

  if (length(lines) == 0) {
    stop("Target sheet is empty: ", path)
  }

  dat <- data.table::fread(
    text = paste(lines, collapse = "\n"),
    header = TRUE,
    blank.lines.skip = TRUE,
    data.table = FALSE
  )

  normalized_names <- tolower(gsub("[^A-Za-z0-9]+", "_", colnames(dat)))
  colnames(dat) <- normalized_names

  rename_if_present <- function(candidates, new_name) {
    idx <- match(TRUE, candidates %in% colnames(dat), nomatch = 0)
    if (idx > 0) {
      colnames(dat)[colnames(dat) == candidates[idx]] <<- new_name
    }
  }

  rename_if_present(c("target_name", "target", "name", "id"), "target_name")
  rename_if_present(c("target_type", "type", "kind"), "target_type")

  required_cols <- c("target_name", "target_type")
  missing_cols <- setdiff(required_cols, colnames(dat))
  if (length(missing_cols) > 0) {
    stop("Target sheet is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  parse_flag <- function(values) {
    lowered <- tolower(trimws(as.character(values)))
    lowered %in% c("1", "true", "t", "yes", "y", "x")
  }

  for (flag_col in c("ssgsea", "heatmap", "gene_plot", "volcano_highlight")) {
    if (flag_col %in% colnames(dat)) {
      dat[[flag_col]] <- parse_flag(dat[[flag_col]])
    } else {
      dat[[flag_col]] <- FALSE
    }
  }

  dat$target_name <- trimws(as.character(dat$target_name))
  dat$target_type <- tolower(trimws(as.character(dat$target_type)))
  dat <- dat[nzchar(dat$target_name), , drop = FALSE]

  dat
}

get_targets_by_flag <- function(target_sheet, target_types, flag_col) {
  if (is.null(target_sheet) || nrow(target_sheet) == 0 || !(flag_col %in% colnames(target_sheet))) {
    return(character())
  }

  idx <- target_sheet$target_type %in% tolower(target_types) & target_sheet[[flag_col]]
  unique(target_sheet$target_name[idx])
}