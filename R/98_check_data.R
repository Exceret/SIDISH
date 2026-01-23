#' @keywords internal
check_cm_genes <- function(sc_data, bulk) {
  common_genes <- intersect(
    rownames(bulk),
    rownames(sc_data)
  )
  if (length(common_genes) == 0) {
    cli::cli_abort(c(
      "x" = "No common genes between {.arg bulk} and {.arg sc_data}",
      ">" = "Please check the gene names and try again."
    ))
  }
}
