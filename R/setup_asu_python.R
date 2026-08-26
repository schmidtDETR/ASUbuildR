#' Setup Python Environment for ASUbuildR
#'
#' This function sets up the Python environment needed for ASUbuildR. It
#' checks for an existing conda installation (installing Miniconda if
#' necessary), creates a conda environment using only the \emph{conda-forge}
#' channel for the Python interpreter itself (avoiding channels that require
#' Terms of Service acceptance), and installs the required Python packages
#' (\code{numpy}, \code{pandas}, \code{networkx}, \code{ortools},
#' \code{openpyxl}) via pip inside that environment. OR-Tools is not
#' distributed on conda-forge, so
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
  have_conda <- length(conda_bin) == 1L && !is.na(conda_bin) &&
    nzchar(conda_bin) && file.exists(conda_bin)
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
  conda_bin <- tryCatch(reticulate::conda_binary(), error = function(e) "")
  if (length(conda_bin) != 1L || is.na(conda_bin) || !nzchar(conda_bin) ||
      !file.exists(conda_bin)) {
    stop("Miniconda installation did not provide a usable conda executable.",
         call. = FALSE)
  }

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

  pkgs <- c("numpy", "pandas", "networkx", "ortools", "openpyxl")

  # Try to attach the env so we can skip already-installed packages.
  # If Python is already initialised to a *different* interpreter in this
  # session (reticulate one-time-init rule), fall back to installing all
  # Locate the Python and pip executables by full path to avoid ambiguity
  # when multiple conda installations are present (Anaconda + r-miniconda).
  # This also bypasses reticulate's one-time-init constraint entirely.
  envs_df     <- reticulate::conda_list()
  env_row     <- envs_df[envs_df$name == env_name, ]
  if (nrow(env_row) == 0)
    stop("Conda environment '", env_name, "' not found after creation.", call. = FALSE)
  python_path <- env_row$python[1]
  pip_path    <- if (.Platform$OS.type == "windows")
    file.path(dirname(python_path), "Scripts", "pip.exe")
  else
    file.path(dirname(python_path), "pip")

  pkgs <- c("numpy", "pandas", "networkx", "ortools", "openpyxl")
  pkg_list_py <- paste0("[", paste0('"', pkgs, '"', collapse = ", "), "]")

  # Write Python snippets to temp files to avoid Windows cmd.exe quoting issues.
  check_file  <- tempfile(fileext = ".py")
  verify_file <- tempfile(fileext = ".py")
  on.exit(unlink(c(check_file, verify_file), force = TRUE), add = TRUE)

  writeLines(c(
    "import importlib.util",
    paste0("pkgs = ", pkg_list_py),
    'missing = [p for p in pkgs if importlib.util.find_spec(p) is None]',
    'print("\\n".join(missing))'
  ), check_file)

  raw          <- system2(python_path, args = check_file, stdout = TRUE, stderr = FALSE)
  missing_pkgs <- raw[nzchar(trimws(raw))]

  if (length(missing_pkgs) == 0) {
    message("All required Python packages are already installed.")
  } else {
    message("Installing Python packages via pip: ",
            paste(missing_pkgs, collapse = ", "), "...")
    status <- system2(pip_path, args = c("install", "--quiet", "--isolated", missing_pkgs))
    if (!identical(status, 0L))
      stop("pip install failed (exit status: ", status, "). ",
           "See the output above for details.", call. = FALSE)
  }

  writeLines(c(
    "import importlib.util",
    paste0("pkgs = ", pkg_list_py),
    'bad = [p for p in pkgs if importlib.util.find_spec(p) is None]',
    'print(",".join(bad) if bad else "OK")'
  ), verify_file)

  verify <- trimws(system2(python_path, args = verify_file, stdout = TRUE, stderr = FALSE))
  if (!identical(verify, "OK"))
    stop(
      "Required Python packages are still not importable: ", verify,
      ". Try setup_asu_python(force = TRUE) to recreate the environment.",
      call. = FALSE
    )

  message("\nPython setup complete! You can now run launch_ASUbuildR()")


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

  required <- c("numpy", "pandas", "networkx", "ortools", "openpyxl")
  available <- vapply(required, reticulate::py_module_available, logical(1))

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
