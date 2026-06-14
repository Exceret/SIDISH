#' @title Modify SIDISH parameters with user-provided values
#' @description
#' This function provides a default set of parameters for the SIDISH method and merges them
#' with any user-provided parameters. The parameters are organized into several categories:
#' preprocessing, VAE training, Deep Cox training, training/risk definition, and execution environment.
#'
#' @param usr_list A named list of parameters to override the defaults. Only parameters that exist in the default list will be modified.
#'
#' @return A list containing the merged parameters with user-provided values overriding the defaults.
#'
#' @section Default Parameters:
#' The following parameters can be passed via \code{...}. Defaults are inherited from
#' the Python SIDISH implementation.
#'
#' Preprocessing parameters:
#'   - \code{patient_id}: Column name for patient/sample identifiers in phenotype data
#'     (default: \code{"Sample"})
#'   - \code{celltype_name}: Column name for cell type annotations in single-cell metadata.
#'     Mapped from R's \code{label_type} parameter (default: \code{"celltype_major"})
#'   - \code{processed}: Whether input expression matrices are already normalized/log-transformed
#'     (default: \code{TRUE})
#'   - \code{use_spatial_graph}: Whether to incorporate spatial neighborhood information
#'     (requires spatial coordinates) (default: \code{FALSE})
#'   - \code{k_neighbors}: Number of neighbors for spatial graph construction.
#'     Ignored when \code{use_spatial_graph = FALSE} (default: \code{NULL})
#'
#' VAE Training (Phase 1) parameters:
#'   - \code{phase1_epochs}: Total number of epochs for initial VAE training (default: \code{225L})
#'   - \code{phase1_i_epochs}: Number of epochs for VAE retraining iterations (default: \code{20L})
#'   - \code{phase1_latent_size}: Dimensionality of the latent space representation (default: \code{32L})
#'   - \code{phase1_layer_dims}: Hidden layer dimensions for encoder/decoder networks
#'     (default: \code{c(512L, 128L)})
#'   - \code{phase1_batch_size}: Mini-batch size for VAE training (default: \code{256L})
#'   - \code{phase1_optimizer}: Optimizer algorithm for VAE training (default: \code{"Adam"})
#'   - \code{phase1_lr}: Initial learning rate for VAE training (default: \code{1.0e-4})
#'   - \code{phase1_lr_3}: Learning rate for later VAE training iterations (default: \code{1.0e-4})
#'   - \code{phase1_dropout}: Dropout rate for VAE layers (default: \code{0L})
#'   - \code{phase1_type}: VAE architecture type: \code{"Normal"} (standard) or \code{"Dense"}
#'     (dense representation) (default: \code{"Normal"})
#'
#' Deep Cox Training (Phase 2) parameters:
#'   - \code{phase2_epochs}: Number of training epochs for Deep Cox proportional hazards model
#'     (default: \code{500L})
#'   - \code{phase2_hidden}: Number of neurons in hidden layer of Deep Cox network
#'     (default: \code{128L})
#'   - \code{phase2_lr}: Learning rate for Deep Cox model training (default: \code{1.0e-4})
#'   - \code{phase2_dropout}: Dropout rate for Deep Cox layers (default: \code{0L})
#'   - \code{phase2_test_size}: Proportion of samples reserved for validation/testing
#'     (default: \code{0.2})
#'   - \code{phase2_batch_size_bulk}: Mini-batch size for bulk data during Deep Cox training
#'     (default: \code{256L})
#'
#' Training & Risk Definition parameters:
#'   - \code{train_iterations}: Number of iterative training cycles to refine cell risk scores
#'     (default: \code{5L})
#'   - \code{train_percentile}: Percentile threshold for defining high-risk cells
#'     (e.g., 0.95 = top 5\% riskiest cells) (default: \code{0.95})
#'   - \code{train_steepness}: Scaling factor controlling steepness of risk score updates
#'     during iterations (default: \code{30L})
#'   - \code{train_path}: Directory path for saving model checkpoints and intermediate results
#'     (default: \code{"./"})
#'   - \code{train_num_workers}: Number of parallel workers for data loading
#'     (set to 0 for reproducibility) (default: \code{0L})
#'   - \code{train_distribution_fit}: Method for fitting risk score distribution:
#'     \code{"fitted"} (empirical) or \code{"default"} (theoretical) (default: \code{"fitted"})
#'
#' Execution Environment parameters:
#'   - \code{device}: Computation device: \code{"cuda"} (GPU) or \code{"cpu"}
#'     (default: \code{"cuda"})
#'   - \code{seed}: Random seed for reproducibility (synchronized across R/Python/PyTorch)
#'     (default: \code{123L})
#'   - \code{verbose}: Whether to display training progress and diagnostics (default: \code{TRUE})
#'
#' @export
modify_sidish_params <- function(usr_list = list()) {
  default <- list(
    # Preprocessing parameters
    patient_id = "Sample",
    celltype_name = "celltype_major",
    processed = TRUE,
    use_spatial_graph = FALSE,
    k_neighbors = NULL,

    # Phase 1: VAE training
    phase1_epochs = 225L,
    phase1_i_epochs = 20L,
    phase1_latent_size = 32L,
    phase1_layer_dims = c(512L, 128L),
    phase1_batch_size = 256L,
    phase1_optimizer = "Adam",
    phase1_lr = 1.0e-4,
    phase1_lr_3 = 1.0e-4,
    phase1_dropout = 0L,
    phase1_type = "Normal",

    # Phase 2: Deep Cox training
    phase2_epochs = 500L,
    phase2_hidden = 128L,
    phase2_lr = 1.0e-4,
    phase2_dropout = 0,
    phase2_test_size = 0.2,
    phase2_batch_size_bulk = 256L,

    # Training & risk definition
    train_iterations = 5L,
    train_percentile = 0.95,
    train_steepness = 30L,
    train_path = "./SIDISH_res/",
    train_num_workers = 0L,
    train_distribution_fit = "fitted",

    # Execution environment
    device = "cuda",
    seed = 123L,
    verbose = TRUE
  )
  utils::modifyList(default, usr_list)
}
