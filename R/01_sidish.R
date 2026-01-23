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

  colnames(phenotype) <- c("duration", "event")
  # matrix will become np.array
  bulk <- reticulate::r_to_py(as.data.frame(t(matched_bulk)))


  # import scanpy as sc
  # import pandas as pd
  # import numpy as np
  # import torch
  # import random
  # import os
  # import matplotlib.pyplot as plt

  sc <- reticulate::import("scanpy", as = "sc")
  pd <- reticulate::import("pandas", as = "pd")
  np <- reticulate::import("numpy", as = "np")
  torch <- reticulate::import("torch")
  random <- reticulate::import("random")
  os <- reticulate::import("os")
  plt <- reticulate::import("matplotlib.pyplot", as = "plt")
  sidish <- reticulate::import("SIDISH")

  torch$manual_seed(seed)
  torch$cuda$manual_seed(seed)
  torch$cuda$manual_seed_all(seed) # if you are using multi-GPU.
  random$seed(seed)
  torch$backends$cudnn$deterministic <- TRUE
  torch$backends$cudnn$benchmark <- FALSE
  ite <- 0

  set_seed(seed = seed, np = np, torch = torch, random = random)

  sdh <- sidish$SIDISH(adata, bulk, "cuda", seed = ite)
  sdh$init_Phase1(225, 20, 32, c(512, 128), 256, "Adam", 1.0e-4, 1e-4, 0)
  sdh$init_Phase2(500, 128, 1e-4, 0, 0.2, 256)
  train_adata <- sdh$train(5, 0.95, 30, "../data/LUNG/", distribution_fit = "fitted")
}
