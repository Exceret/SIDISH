#' @title Install PyTorch
#' @description Automatically detect system and hardware architecture,
#' then install PyTorch using pip via reticulate.
#' @param env_name Character, the name of the Python environment to install PyTorch into.
#' @param build Character, "stable" or "preview" build of PyTorch.
#' @param os Character, operating system (auto-detected if not specified).
#' @param compute_platform Character, compute platform (CUDA version or "cpu", auto-detected if NULL).
#' @param ... Additional arguments passed to reticulate functions.
#' @return Invisible NULL
#' @examples
#' \dontrun{
#' pip_install_pytorch()
#' }
#' @references https://pytorch.org/get-started/locally/
#' @export
pip_install_pytorch <- function(
  env_name = "r-reticulate-sidish", # nolint: object_name_linter.
  build = c("stable", "preview"),
  os = c("linux", "macos", "windows"),
  compute_platform = NULL,
  ...
) {
  # Check arguments
  check_character(env_name, "env_name")
  check_character(build, "build")
  check_character(os, "os")
  check_character(compute_platform, "compute_platform")

  dots <- rlang::list2(...)

  # Match arguments
  build <- SigBridgeRUtils::MatchArg(build, c("stable", "preview"))
  os <- SigBridgeRUtils::MatchArg(
    os,
    c("linux", "macos", "windows"),
    tolower(Sys.info()["sysname"])
  )

  # Detect CPU architecture
  cpu_arch <- Sys.info()["machine"]
  if (grepl("x86_64|amd64", cpu_arch)) {
    cpu_arch <- "x86_64"
  } else if (grepl("arm64|aarch64", cpu_arch)) {
    cpu_arch <- "arm64"
  }

  # Detect CUDA availability and version
  compute_platform <- "cpu"

  if (is.null(compute_platform) && os %in% c("linux", "windows")) {
    # Check if nvidia-smi is available and get CUDA version
    nvidia_smi <- rlang::try_fetch(system(
      "which nvidia-smi",
      intern = TRUE,
      ignore.stderr = TRUE
    ), silent = TRUE)

    if (!inherits(nvidia_smi, "try-error") && length(nvidia_smi) > 0) {
      cuda_info <- rlang::try_fetch(system(
        "nvidia-smi --query-gpu=driver_version --format=csv,noheader",
        intern = TRUE, ignore.stderr = TRUE
      ), silent = TRUE)

      if (!inherits(cuda_info, "try-error") && length(cuda_info) > 0) {
        # Assume CUDA 13.0 as default for now
        compute_platform <- "cu130"
      }
    }
  }

  # Generate PyTorch install command based on detected information
  base_url <- "https://download.pytorch.org/whl/"
  index_url <- paste0(base_url, compute_platform)


  # Build install command
  install_cmd <- paste(
    "pip install torch torchvision",
    paste0("--index-url ", index_url)
  )

  # Print information using cli
  cli::cli_h1("PyTorch Installation Information")
  cli::cli_div(theme = list(span.emph = list(color = "blue")))

  cli::cli_alert_info("System Architecture:")
  cli::cli_li("OS: {.emph {detected_os}}")
  cli::cli_li("CPU Architecture: {.emph {cpu_arch}}")
  cli::cli_li("Compute Platform: {.emph {compute_platform}}")

  cli::cli_alert_info("Installation Details:")
  cli::cli_li("Environment Name: {.emph {env_name}}")
  cli::cli_li("Build: {.emph {build}}")
  cli::cli_li("Install Command: {.emph {install_cmd}}")
  cli::cli_end()

  # Interactive confirmation
  confirm <- utils::menu(
    choices = c("Yes", "No"),
    title = "Do you want to proceed with the installation?"
  )

  if (confirm != 1) {
    cli::cli_alert_warning(cli::col_grey("Installation cancelled by user."))
    return(invisible(NULL))
  }

  # Check if environment exists
  existing_envs <- SigBridgeRUtils::ListPyEnv()
  if (!env_name %in% existing_envs$name) {
    cli::cli_abort(c(
      "x" = "Environment {.emph {env_name}} does not exist.",
      ">" = "Available environments: {.emph {existing_envs$name}}"
    ))
  }
  env_type <- existing_envs[existing_envs$name == env_name, "type"]

  # Install PyTorch
  switch(env_type,
    "conda" = rlang::exec(
      reticulate::conda_install,
      env_name,
      c("torch", "torchvision"),
      pip = TRUE,
      pip_options = c(paste0("--index-url=", index_url)),
      !!!SigBridgeRUtils::FilterArgs4Func(dots, reticulate::conda_install)
    ),
    "venv" = rlang::exec(
      reticulate::virtualenv_install,
      envname = env_name,
      packages = c("torch", "torchvision"),
      pip_options = c(paste0("--index-url=", index_url)),
      !!!SigBridgeRUtils::FilterArgs4Func(dots, reticulate::virtualenv_install)
    ),
    cli::cli_abort(c(
      "x" = "Environment {.emph {env_name}} is not \\
             a conda or virtual environment.",
      ">" = "Available environments: {.emph {existing_envs$name}}"
    ))
  )

  return(invisible(NULL))
}
