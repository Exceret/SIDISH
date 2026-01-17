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
      ">" = "Currently {.arg sc_data} is of {.cls {class(sc_data)}}"
    ))
  }

  # * filter genes & make it in order
  common_genes <- intersect(
    rownames(matched_bulk),
    rownames(sc_data)
  )
  if (length(common_genes) == 0) {
    cli::cli_abort(c(
      "x" = "No common genes between {.arg matched_bulk} and {.arg sc_data}",
      ">" = "Please check the gene names and try again."
    ))
  }
  sc_data <- sc_data[common_genes, ]
  matched_bulk <- matched_bulk[common_genes, ]

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
    assay_name = "assay",
    output_class = "ReticulateAnnData"
  )
  py <- reticulate::py
  reticulate::py_require("sidish")

  # * activate binding
  reticulate::py_run_string("print('\n')")
  py$adata <- adata
  py$bulk <- as.data.frame(t(matched_bulk)) # matrix will become np.array

  colnames(phenotype) <- c("Overall_survival_days", "Sample_Status")
  py$survival <- as.data.frame(phenotype)
}
