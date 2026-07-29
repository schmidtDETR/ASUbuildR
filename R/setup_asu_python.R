#' Setup Python Environment for ASUbuildR
#'
#' This function sets up the Python environment needed for ASUbuildR. It
#' checks for an existing conda installation (installing Miniconda if
#' necessary), creates a conda environment using only the \emph{conda-forge}
#' channel for the Python interpreter itself (avoiding channels that require
#' Terms of Service acceptance), and installs the required Python packages
#' (\code{numpy}, \code{pandas}, \code{networkx}, \code{ortools}) via pip
#' inside that environment. OR-Tools is not distributed on conda-forge, so
#' pip is used for the packages rather than \code{conda install}.
#'
#' @param force Logical. If TRUE, recreates the conda environment even if it exists.
#' @return Invisible NULL. Called for side effects.
#' @export
#' @examples
#' \dontrun{
#' # First-time setup
#' setup_asu_python()
#'
#' # Force reinstall
#' setup_asu_python(force = TRUE)
#' }
setup_asu_python <- function(force = FALSE) {

  message("Setting up Python environment for ASUbuildR...")

  env_name <- "asu-cpsat"

  # Run a conda subcommand and stop with an informative error if it fails,
  # instead of silently continuing as if the step had succeeded.
  run_conda_step <- function(args, step_label) {
    status <- system2(conda_bin, args)
    if (!identical(status, 0L)) {
      stop(
        "Failed to ", step_label, " (conda exit status: ", status, ").\n",
        "Command: ", conda_bin, " ", paste(args, collapse = " "), "\n",
        "See the conda output above for details. Common causes: no internet ",
        "connection, a proxy blocking conda-forge, or insufficient disk space.",
        call. = FALSE
      )
    }
    invisible(status)
  }

  # Determine whether conda is available; if not, install Miniconda
  conda_bin <- tryCatch(reticulate::conda_binary(), error = function(e) "")
  have_conda <- length(conda_bin) == 1L && !is.na(conda_bin) && nzchar(conda_bin)
  just_installed <- FALSE
  if (!have_conda) {
    message("Conda not found. Installing Miniconda...")
    message("This is a one-time installation and may take a few minutes...")
    reticulate::install_miniconda(update = FALSE)
    have_conda <- TRUE
    just_installed <- TRUE
    message("Miniconda installed successfully!")
  }

  if (!have_conda) {
    stop("Conda installation failed; cannot set up Python environment.")
  }

  # Path to conda executable (re-resolved in case Miniconda was just installed)
  conda_bin <- reticulate::conda_binary()

  # Remove existing environment if force = TRUE
  envs <- reticulate::conda_list()$name
  if (force && env_name %in% envs) {
    message("Removing existing conda environment...")
    tryCatch(
      reticulate::conda_remove(envname = env_name, packages = "--all"),
      error = function(e) {
        stop("Failed to remove the existing '", env_name, "' environment: ",
             conditionMessage(e), call. = FALSE)
      }
    )
    envs <- setdiff(envs, env_name)
  }

  # Ensure base conda uses only conda-forge when freshly installed
  if (just_installed) {
    run_conda_step(
      c("update", "--name", "base", "conda", "--yes",
        "--quiet", "--override-channels", "-c", "conda-forge"),
      "update the base conda installation"
    )
  }

  # Create conda environment if needed (using only the conda-forge channel)
  if (!(env_name %in% envs)) {
    message("Creating conda environment '", env_name, "' via conda-forge...")
    run_conda_step(
      c("create", "--yes", "--name", env_name, "python=3.11",
        "--quiet", "--override-channels", "-c", "conda-forge"),
      paste0("create the '", env_name, "' conda environment")
    )
    message("Conda environment created successfully!")
  } else {
    message("Conda environment already exists. Use force=TRUE to recreate.")
  }

  # Activate the env so we can check what's already importable. Checking
  # first (rather than always re-installing) matters in practice: asking
  # conda to re-solve numpy+pandas+networkx+ortools together against the
  # *current* conda-forge metadata can fail (e.g. no mutually satisfiable
  # build combination is resolvable that day) even when every package is
  # already installed and working - re-running setup_asu_python() should
  # not break an environment that's already fine.
  reticulate::use_condaenv(env_name, required = TRUE)

  pkgs <- c("numpy", "pandas", "networkx", "ortools")
  already_available <- vapply(pkgs, reticulate::py_module_available, logical(1))
  missing_pkgs <- pkgs[!already_available]

  if (length(missing_pkgs) == 0) {
    message("All required Python packages are already installed.")
  } else {
    # OR-Tools has no conda-forge distribution, so install via pip (inside
    # the conda environment) rather than `conda install`, which would fail
    # to resolve 'ortools' from conda-forge no matter the Python version.
    message("Installing missing Python packages via pip: ",
            paste(missing_pkgs, collapse = ", "), "...")
    tryCatch(
      reticulate::py_install(
        packages = missing_pkgs,
        envname  = env_name,
        method   = "conda",
        conda    = conda_bin,
        pip      = TRUE
      ),
      error = function(e) {
        stop(
          "Failed to install the required Python packages via pip: ",
          conditionMessage(e), "\n",
          "Common causes: no internet connection or a proxy blocking PyPI.",
          call. = FALSE
        )
      }
    )
  }

  # Verify the packages actually import before declaring success - a conda
  # "solve" can succeed while a package still fails to import (e.g. an ABI
  # mismatch), so don't report success on install exit status alone.
  available <- vapply(pkgs, reticulate::py_module_available, logical(1))
  if (!all(available)) {
    stop(
      "Required Python packages are still not importable: ",
      paste(pkgs[!available], collapse = ", "),
      ". Try setup_asu_python(force = TRUE) to recreate the environment.",
      call. = FALSE
    )
  }

  message("\nPython setup complete! All required packages are importable.")
  message("You can now run launch_ASUbuildR()")

  invisible(NULL)
}

#' Check Python Setup for ASUbuildR
#'
#' Checks if the Python environment is properly configured for ASUbuildR.
#'
#' @return Logical. TRUE if properly configured, FALSE otherwise.
#' @export
check_asu_python <- function() {
  if (!asu_use_python(required = FALSE)) {
    message("Conda environment not found.")
    message("Run setup_asu_python() to set up the Python environment.")
    return(FALSE)
  }

  required <- c("numpy", "pandas", "networkx", "ortools")
  available <- sapply(required, reticulate::py_module_available)

  if (all(available)) {
    message("Python environment is properly configured")
    return(TRUE)
  } else {
    missing <- required[!available]
    message("Missing packages: ", paste(missing, collapse = ", "))
    message("Run setup_asu_python() to install missing packages.")
    return(FALSE)
  }
}
