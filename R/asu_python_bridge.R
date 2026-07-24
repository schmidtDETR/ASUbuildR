#' Select the ASU CP-SAT Python environment
#' @keywords internal
asu_use_python <- function(required = FALSE) {
  conda_bin <- tryCatch(reticulate::conda_binary(), error = function(e) "")
  have_conda <- length(conda_bin) == 1L && !is.na(conda_bin) && nzchar(conda_bin)
  if (!have_conda) {
    if (required) stop("Conda was not found. Run ASUbuildR::setup_asu_python() first.")
    return(FALSE)
  }

  envs <- tryCatch(
    reticulate::conda_list(conda = conda_bin)$name,
    error = function(e) character()
  )
  if (!("asu-cpsat" %in% envs)) {
    if (required) stop("Python environment 'asu-cpsat' was not found. Run ASUbuildR::setup_asu_python() first.")
    return(FALSE)
  }

  reticulate::use_condaenv("asu-cpsat", conda = conda_bin, required = required)
  TRUE
}

#' Load the ASU CP-SAT Python module
#' @keywords internal
asu_load_py <- function() {
  path <- system.file("python", "asu_cpsat.py", package = "ASUbuildR")
  if (path == "") stop("Couldn't find inst/python/asu_cpsat.py in the installed package.")

  module_env <- new.env(parent = baseenv())
  reticulate::source_python(path, envir = module_env)
  list(
    build_many_asus_cpsat = get(
      "build_many_asus_cpsat",
      envir = module_env,
      inherits = FALSE
    )
  )
}
