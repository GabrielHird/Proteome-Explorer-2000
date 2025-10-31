`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

to_named_string <- function(x) {
  if (is.null(x) || !length(x)) {
    return("")
  }
  if (is.null(names(x)) || all(names(x) == "")) {
    return(paste(x, collapse = ", "))
  }
  paste(sprintf("%s = %s", names(x), x), collapse = "\n")
}

parse_comma_list <- function(x) {
  vals <- trimws(unlist(strsplit(x %||% "", ",")))
  vals <- vals[nzchar(vals)]
  if (length(vals)) vals else character()
}

parse_named_lines <- function(x) {
  lines <- strsplit(x %||% "", "\n", fixed = TRUE)[[1]]
  lines <- trimws(lines)
  lines <- lines[nzchar(lines)]
  if (!length(lines)) {
    return(character())
  }
  split <- strsplit(lines, "=", fixed = TRUE)
  values <- vapply(split, function(parts) {
    parts <- trimws(parts)
    if (length(parts) < 2) "" else parts[2]
  }, character(1))
  names(values) <- vapply(split, function(parts) {
    parts <- trimws(parts)
    parts[1]
  }, character(1))
  values <- values[nzchar(values) & nzchar(names(values))]
  values
}

empty_to_na <- function(x) {
  x <- trimws(x %||% "")
  if (!nzchar(x)) {
    return(NA)
  }
  x
}

safe_numeric <- function(x) {
  if (is.null(x)) {
    return(NA_real_)
  }
  if (is.numeric(x)) {
    return(as.numeric(x))
  }
  x <- trimws(as.character(x))
  if (!nzchar(x)) {
    return(NA_real_)
  }
  suppressWarnings(as.numeric(x))
}
