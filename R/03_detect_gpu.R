#' @title Detect available GPU devices
#' @description
#' Checks for the presence of GPU hardware by detecting NVIDIA, AMD, macOS integrated,
#' or Windows GPUs using system commands.
#'
#' @param verbose Logical indicating whether to print detection messages (default: TRUE)
#'
#' @return Logical value indicating whether any GPU was detected (TRUE) or not (FALSE)
#'
#' @details This function attempts to detect GPU hardware by:
#' \itemize{
#'   \item Checking for NVIDIA GPUs using `nvidia-smi`
#'   \item Checking for AMD GPUs using `rocminfo`
#'   \item Checking for macOS integrated GPUs using `system_profiler`
#'   \item Checking for Windows GPUs using `wmic`
#' }
#' The function provides visual feedback through cli alerts about detection results.
#' @export
detect_gpu <- function(verbose = TRUE) {
  nvidia <- tryCatch(
    {
      out <- system(
        "nvidia-smi --query-gpu=name --format=csv,noheader 2>&1",
        intern = TRUE,
        ignore.stderr = TRUE
      )
      !is.null(out) &&
        length(out) > 0 &&
        !grepl("not found|command not found", out[1], ignore.case = TRUE)
    },
    error = function(e) FALSE
  )

  if (nvidia) {
    if (verbose) {
      cli::cli_alert_success("Detected NVIDIA GPU(s)")
    }
    return(TRUE)
  }

  amd <- tryCatch(
    {
      out <- system(
        "rocminfo 2>&1 | grep -i 'Name' | head -1",
        intern = TRUE,
        ignore.stderr = TRUE
      )
      !is.null(out) &&
        length(out) > 0 &&
        !grepl("not found", out[1], ignore.case = TRUE)
    },
    error = function(e) FALSE
  )

  if (amd) {
    if (verbose) {
      cli::cli_alert_success("Detected AMD GPU")
    }
    return(TRUE)
  }

  # macOS集成GPU
  if (.Platform$OS.type == "unix" && Sys.info()["sysname"] == "Darwin") {
    mac_gpu <- tryCatch(
      length(grep(
        "Chipset Model",
        system("system_profiler SPDisplaysDataType 2>&1", intern = TRUE),
        value = TRUE
      )) >
        0,
      error = function(e) FALSE
    )

    if (mac_gpu) {
      if (verbose) {
        cli::cli_alert_success("Detected GPU on macOS")
      }
      return(TRUE)
    }
  }

  # Windows GPU
  if (.Platform$OS.type == "windows") {
    win_gpu <- tryCatch(
      {
        out <- system(
          'wmic path win32_VideoController get name 2>&1',
          intern = TRUE
        )
        length(out) > 1 && any(nzchar(trimws(out[-1])))
      },
      error = function(e) FALSE
    )

    if (win_gpu) {
      if (verbose) {
        cli::cli_alert_success("Detected GPU on Windows")
      }
      return(TRUE)
    }
  }

  if (verbose) {
    cli::cli_alert_danger("No GPU detected")
  }

  FALSE
}
