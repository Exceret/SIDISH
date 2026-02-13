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
#' @param ... Additional parameters passed to the Python SIDISH implementation. Parameters
#' are grouped by functionality below. Unrecognized parameters are passed to Python but may be ignored.
#'
#' ### Preprocessing parameters (used in `preprocess()`)
#'
#' - **`patient_id`** (default: `"Sample"`): Column name for patient IDs in the `phenotype` data frame.
#' - **`celltype_name`** (default: `"celltype_major"`): Metadata column name containing cell type #' annotations in the Seurat object.
#' - **`processed`** (default: `TRUE`): Whether the scRNA-seq data has been preprocessed (log-normalized). #' If `FALSE`, activates QC filtering steps below.
#' - **`n_genes_by_counts`** (default: `5000`): Minimum number of genes per cell for QC filtering (used #' only when `processed=FALSE`).
#' - **`pct_counts_mt`** (default: `10`): Maximum percentage of mitochondrial counts per cell for QC #' filtering (used only when `processed=FALSE`).
#' - **`batch_correction`** (default: `FALSE`): Enable batch correction during preprocessing.
#' - **`survival_`** (default: `"time"`): Column name for survival time in the `phenotype` data frame.
#' - **`status`** (default: `"status"`): Column name for event status (0 = censored, 1 = event) in the #' `phenotype` data frame.
#'
#' ### Model initialization parameters (used in `sidish()` constructor)
#'
#' - **`device`** (default: `"cuda"`): PyTorch device for computation (`"cuda"` for GPU or `"cpu"`).
#' - **`use_spatial_graph`** (default: `FALSE`): Enable spatial graph integration (relevant for spatial #' transcriptomics data).
#' - **`k_neighbors`** (default: `NULL`): Number of neighbors for spatial graph construction; `NULL` uses #' automatic selection.
#'
#' ### Phase 1 parameters (VAE training via `init_Phase1()`)
#'
#' - **`phase1_epochs`** (default: `225`): Number of epochs for initial VAE training.
#' - **`phase1_i_epochs`** (default: `20`): Number of epochs for VAE retraining in each iteration.
#' - **`phase1_latent_size`** (default: `32`): Dimensionality of the VAE latent space.
#' - **`phase1_layer_dims`** (default: `c(512, 128)`): Hidden layer dimensions for the VAE encoder/decoder.
#' - **`phase1_batch_size`** (default: `256`): Batch size for VAE training.
#' - **`phase1_optimizer`** (default: `"Adam"`): Optimizer used for VAE training.
#' - **`phase1_lr`** (default: `1e-4`): Initial learning rate for VAE training.
#' - **`phase1_lr_3`** (default: `1e-4`): Learning rate used in later VAE training iterations.
#' - **`phase1_dropout`** (default: `0`): Dropout rate applied in VAE layers.
#' - **`phase1_type`** (default: `"Normal"`): Representation type (`"Normal"` or `"Dense"`).
#'
#' ### Phase 2 parameters (Deep Cox training via `init_Phase2()`)
#'
#' - **`phase2_epochs`** (default: `500`): Number of training epochs for the Deep Cox model.
#' - **`phase2_hidden`** (default: `128`): Number of neurons in the hidden layer of the Deep Cox model.
#' - **`phase2_lr`** (default: `1e-4`): Learning rate for Deep Cox training.
#' - **`phase2_dropout`** (default: `0`): Dropout rate applied in the Deep Cox model.
#' - **`phase2_test_size`** (default: `0.2`): Proportion of data allocated to the test set.
#' - **`phase2_batch_size_bulk`** (default: `256`): Batch size for bulk RNA-seq data during Deep Cox #' training.
#'
#' ### Training parameters (used in `train()`)
#'
#' - **`train_iterations`** (default: `5`): Total number of iterative training cycles.
#' - **`train_percentile`** (default: `0.95`): Percentile threshold for defining high-risk cells (e.g., `0.#' 95` selects the top 5% highest-risk cells).
#' - **`train_steepness`** (default: `30`): Scaling factor controlling the steepness of weight updates #' during risk propagation.
#' - **`train_path`** (default: `"./SIDISH_res/"`): Directory path for saving model checkpoints and #' intermediate results.
#' - **`train_num_workers`** (default: `0`): Number of parallel workers for data loading; set to `0` for #' reproducible results.
#' - **`train_distribution_fit`** (default: `"fitted"`): Method for survival distribution fitting #' (`"fitted"` or `"default"`).
#'
#'
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

  obs <- reticulate::py_to_r(py$sidish_obs) # pd.DataFrame -> R data.frame
  colnames(obs)[3] <- paste0("SIDISH_", colnames(obs)[3])
  obs$SIDISH <- ifelse(obs$SIDISH == "h", "Positive", "Other")

  sc_data <- SeuratObject::AddMetaData(
    object = sc_data,
    metadata = obs
  )
  sc_data <- SigBridgeRUtils::AddMisc(
    seurat_obj = sc_data,
    SIDISH = list(
      label_type = label_type,
      phenotype_class = phenotype_class,
      params = params
    )
  )

  if (verbose) {
    ts_cli$cli_alert_info(cli::col_green("SIDISH Done"))
  }

  sc_data
}
