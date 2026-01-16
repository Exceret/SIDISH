sidish <- function(
  matched_bulk,
  sc_data,
  phenotype,
  label_type = NULL,
  phenotype_class = "survival",
  assay = "RNA",
  ...
) {
  dots <- rlang::list2(...)
  verbose <- dots$verbose %||%
    SigBridgeRUtils::getFuncOption("verbose") %||%
    TRUE
  seed <- dots$seed %||% SigBridgeRUtils::getFuncOption("seed") %||% 123L

  if (!inherits(sc_data, "Seurat")) {
    cli::cli_abort(c(
      "x" = "{.arg sc_data} must be a {.cls SeuratObject}",
      ">" = "but it is currently of {.cls {class(sc_data)})"
    ))
  }

  adata <- anndataR::as_AnnData(
    x = sc_data,
    x_mapping = "counts",
    layers_mapping = TRUE,
    obs_mapping = TRUE,
    var_mapping = TRUE,
    obsm_mapping = TRUE,
    varm_mapping = TRUE,
    obsp_mapping = TRUE,
    varp_mapping = TRUE,
    uns_mapping = TRUE,
    assay_name = 'assay',
    output_class = "ReticulateAnnData"
  )
}
