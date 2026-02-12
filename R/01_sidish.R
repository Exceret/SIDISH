#' @export
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

  set.seed(seed)

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

  # * samples matching
  if (!all(colnames(matched_bulk) == rownames(phenotype))) {
    cli::cli_warn(
      "Samples in {.arg matched_bulk} and {.arg phenotype} do not match, \\
      auto-correcting samples"
    )
    cm_samples <- intersect(
      colnames(matched_bulk),
      rownames(phenotype)
    )
    matched_bulk <- matched_bulk[, cm_samples]
    phenotype <- phenotype[cm_samples, ]
  }

  adata <- anndataR::as_AnnData(
    x = sc_data,
    x_mapping = "data",
    layers_mapping = TRUE,
    obs_mapping = TRUE,
    var_mapping = TRUE,
    obsm_mapping = TRUE,
    varm_mapping = TRUE,
    obsp_mapping = TRUE,
    varp_mapping = TRUE,
    uns_mapping = TRUE,
    assay_name = assay,
    output_class = "ReticulateAnnData"
  )

  py <- reticulate::py
  # avtivate python connection
  reticulate::py_run_string("import os")

  # * params pass to python
  # matrix will become np.array, we need pd.DataFrame
  py$bulk <- reticulate::r_to_py(as.data.frame(t(matched_bulk)))
  py$adata <- adata
  py$survival_df <- reticulate::r_to_py(phenotype)
  py$seed <- reticulate::r_to_py(seed)

  reticulate::py_run_file(
    "inst/python/01_training_SIDISH.py",
    local = TRUE
  )

  res
}
