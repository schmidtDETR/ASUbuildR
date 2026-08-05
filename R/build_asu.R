#' Build ASUs using OR-Tools CP-SAT (Python) from an R data frame
#'
#' @param df data.frame with columns: geoid, tract_ASU_unemp, tract_ASU_emp, tract_pop2024
#' @param neighbors list of integer vectors (0-based or 1-based ok). If NULL, you must use CLI with --geometry.
#' @param tau numeric unemployment rate threshold (e.g. 0.0645)
#' @param pop_thresh integer population threshold
#' @param max_asus max number of ASUs to carve
#' @param time_limit seconds per window
#' @param workers CP-SAT threads
#' @param rel_gap optional relative MIP gap (e.g. 0.01)
#' @param verbose logical; print CP-SAT logs
#' @return df with added `asu_id` column (integer; -1 means unassigned)
#' @export
build_asu <- function(
    df,
    neighbors,
    tau = 0.0645,
    pop_thresh = 10000,
    max_asus = 30,
    time_limit = 1200,
    workers = max(1L, parallel::detectCores(logical = TRUE) - 1L),
    rel_gap = NA_real_,
    verbose = interactive()
) {
  asu_use_python(required = TRUE)
  if (!reticulate::py_module_available("ortools"))
    stop("Python env missing 'ortools'. Run ASUbuildR::setup_asu_python() first.")

  # Normalize neighbor indexing to 0-based
  if (is.null(neighbors)) stop("Provide `neighbors` as a list of integer vectors (contiguity).")
  n <- nrow(df)
  nb <- lapply(neighbors, as.integer)
  all_indices <- unlist(nb, use.names = FALSE)
  one_based <- length(all_indices) > 0L && !any(all_indices == 0L) &&
    all(all_indices >= 1L & all_indices <= n)
  if (one_based) nb <- lapply(nb, function(v) v - 1L)
  nb <- lapply(nb, function(v) {
    as.list(v[v >= 0L & v < n])
  })

  # Load python module and call
  mod <- asu_load_py()
  # reticulate converts data.frame -> pandas.DataFrame and list(list(int)) -> Python list of lists
  out <- mod$build_many_asus_cpsat(
    df = df,
    nb = nb,
    tau = tau,
    pop_thresh = as.integer(pop_thresh),
    max_asus = as.integer(max_asus),
    time_limit = as.integer(time_limit),
    workers = as.integer(workers),
    rel_gap = if (is.na(rel_gap)) NULL else as.numeric(rel_gap),
    verbose = isTRUE(verbose)
  )

  df$asu_id <- as.integer(reticulate::py_to_r(out[["asu_id"]]))
  attr(df, "n_asu") <- as.integer(out[["n_asu"]])
  df
}
