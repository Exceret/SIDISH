#' @title Run SIDISH for cell-level survival risk estimation
#'
#' @description
#' Wrapper for the SIDISH (Semi-supervised Iterative Deep Learning for Identifying Single-cell High-Risk Populations) algorithm that identifies high-risk cells associated with patient survival outcomes.
#'
#'
#' @param matched_bulk Numeric matrix (genes × samples) of bulk RNA-seq expression.
#' Row names must be gene identifiers matching those in \code{sc_data}. Values should be
#' log-transformed (e.g., log2CPM or log2(TPM+1)).
#' @param sc_data Fully preprocessed Seurat object
#' @param phenotype Data frame containing survival information with samples as rows.
#' For \code{phenotype_class = "survival"}, must contain at least two columns:
#' survival time (numeric) and event indicator (binary 0/1). Row names must match
#' column names of \code{matched_bulk}.
#' @param label_type Character string specifying label of cell
#' @param phenotype_class Character string indicating phenotype type. Currently only
#' \code{"survival"} is supported (expects time/event columns).
#' @param python python executable path. Passed to \code{reticulate::use_python()}.
#' @param sidish_tools Python script for SIDISH training
#' @param assay Character string specifying the assay name for AnnData conversion
#' (default \code{"RNA"}). Used when converting \code{sc_data} to AnnData format.
#' @param ... Additional parameters passed to the Python SIDISH implementation. Key supported
#'   parameters include:
#'   - \code{device}: PyTorch device (\code{"cuda"} or \code{"cpu"}); default \code{"cuda"}
#'   - \code{phase1_epochs}: VAE training epochs; default \code{225L}
#'   - \code{phase1_latent_size}: VAE latent dimension; default \code{32L}
#'   - \code{phase2_epochs}: Deep Cox training epochs; default \code{500L}
#'   - \code{train_percentile}: Percentile threshold for high-risk cells (e.g., 0.95 = top 5%);
#'     default \code{0.95}
#'   - \code{train_iterations}: Number of training iterations; default \code{5L}
#'   - \code{patient_id}: Column name for patient IDs in \code{phenotype}; default \code{"Sample"}
#'   - \code{processed}: Whether input data are pre-normalized; default \code{TRUE}
#'   - \code{verbose}: Show training progress; default inherits from
#'     \code{options(SigBridgeR.verbose)} or \code{TRUE}
#'   - \code{seed}: Random seed for reproducibility; default inherits from
#'     \code{options(SigBridgeR.seed)} or \code{123L}
#'   Unrecognized parameters are passed to Python but ignored by SIDISH. Partial matching is
#'   supported (e.g., \code{phase1_ep} for \code{phase1_epochs}). See the 'Default Parameters'
#'   section for a complete parameter reference.
#'
#' @details
#' Memory considerations:
#'   - GPU acceleration (\code{device = "cuda"}) significantly speeds up training
#'   - For large datasets (>50k cells), reduce \code{phase1_batch_size} or
#'     \code{phase2_batch_size_bulk}
#'   - Set \code{train_num_workers = 0} (default) for reproducibility; increase only if
#'     deterministic results aren't required
#'
#' @note
#' Requires Python environment with:
#'   - \code{SIDISH} package (install via \code{pip install SIDISH})
#'   - PyTorch (\code{torch}), Scanpy (\code{scanpy}), and dependencies
#'   - Proper CUDA setup for GPU usage (\code{device = "cuda"})
#'   The Python script \code{inst/python/01_training_SIDISH.py} must be available in the
#'   package installation directory.
#'
#' @export
sidish <- function(
  matched_bulk,
  sc_data,
  phenotype,
  label_type = NULL,
  phenotype_class = "survival",
  python = NULL,
  sidish_tools = system.file(
    "python/01_training_SIDISH.py",
    package = "rSIDISH"
  ),
  assay = "RNA",
  ...
) {
  dots <- rlang::list2(...)
  verbose <- dots$verbose %||%
    SigBridgeRUtils::getFuncOption("verbose") %||%
    TRUE
  seed <- dots$seed %||% SigBridgeRUtils::getFuncOption("seed") %||% 123L

  params <- modify_sidish_params(usr_list = dots)

  set.seed(seed)

  if (verbose) {
    ts_cli$cli_alert_info(cli::col_green("Start SIDISH screening"))
    ts_cli$cli_alert_info("Checking inputs")
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

  meta_cols <- colnames(sc_data[[]])

  # * samples matching
  if (!all(colnames(matched_bulk) == rownames(phenotype))) {
    cli::cli_warn(
      "Samples in {.arg matched_bulk} and {.arg phenotype} do not match, \\
      try auto-correcting samples"
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

  # * other params
  if (any("" %in% names(params))) {
    cli::cli_abort(c("x" = "Variable name cannot be empty"))
  }
  for (var_name in names(params)) {
    py[[var_name]] <- reticulate::r_to_py(params[[var_name]])
  }

  py$seed <- reticulate::r_to_py(seed)
  py$verbose <- reticulate::r_to_py(verbose)

  reticulate::py_run_file(sidish_tools)

  obs <- reticulate::py_to_r(py$obs) # pd.DataFrame -> R data.frame
  obs <- obs[, -meta_cols]
  colnames(obs) <- paste0("SIDISH_", colnames(obs))

  sc_data <- SeuratObject::AddMetaData(
    object = sc_data,
    metadata =
  )
  sc_data <- SigBridgeRUtils::AddMisc(
    seurat_obj = sc_data,
    SIDISH = list(
      label_type = label_type,
      phenotype_class = phenotype_class
    )
  )

  if (verbose) {
    ts_cli$cli_alert_info(cli::col_green("SIDISH Done"))
  }

  sc_data
}
