#' @keywords internal
check_character <- function(x, name = deparse(substitute(x))) {
  if (!is.character(x)) {
    cli::cli_abort(c(
      "x" = "{.arg {name}} must be a character vector",
      ">" = "Currently {.arg {name}} is of {.cls {class(x)}}"
    ))
  }
}
