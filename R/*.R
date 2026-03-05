#' Index functions defined in R/ by file (static scan)
#'
#' @param dir Directory to scan (default "R")
#' @param recursive Scan subfolders
#' @return data.frame with columns: name, file
index_functions_in_dir <- function(dir = "R", recursive = FALSE) {
  stopifnot(dir.exists(dir))

  files <- list.files(
    dir,
    pattern = "\\.R$",
    full.names = TRUE,
    recursive = recursive
  )

  extract_defs <- function(f) {
    txt <- readLines(f, warn = FALSE)
    # matches: foo <- function(   OR   foo = function(
    m <- regmatches(
      txt,
      regexpr(
        "^\\s*([.A-Za-z][.A-Za-z0-9_]*)\\s*(<-|=)\\s*function\\s*\\(",
        txt,
        perl = TRUE
      )
    )
    if (length(m) == 0) {
      return(character())
    }
    sub("^\\s*([.A-Za-z][.A-Za-z0-9_]*)\\s*(<-|=).*", "\\1", m, perl = TRUE)
  }

  defs <- lapply(files, extract_defs)
  out <- data.frame(
    name = unlist(defs, use.names = FALSE),
    file = rep(files, lengths(defs)),
    stringsAsFactors = FALSE
  )

  out[order(out$name, out$file), , drop = FALSE]
}
