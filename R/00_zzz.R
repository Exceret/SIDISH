# ? Package startup messages
.onAttach <- function(libname, pkgname) {
  pkg_version <- utils::packageVersion(pkgname)

  msg <- cli::cli_fmt(cli::cli_alert_success(
    "{.pkg {pkgname}} v{pkg_version} loaded"
  ))
  packageStartupMessage(msg)
  invisible()
}

.onLoad <- function(libname, pkgname) {
  # Declare Python dependencies (required by inst/python/ scripts)
  reticulate::py_require(c(
    "numpy",
    "pandas",
    "torch",
    "torchvision",
    "scanpy",
    "scikit-learn",
    "tqdm",
    "lifelines",
    "pyro-ppl",
    "seaborn",
    "matplotlib",
    "scipy",
    "statsmodels",
    "imbalanced-learn",
    "shap",
    "torch-geometric",
    "joblib"
  ))
  invisible()
}

#' @keywords internal
ts_cli <- SigBridgeRUtils::CreateTimeStampCliEnv()
