log_text <- function(lines) {
  if (!length(lines)) {
    "Pipeline output will appear here after the next run."
  } else {
    paste(lines, collapse = "\n")
  }
}
