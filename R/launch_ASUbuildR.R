# Script to run ASU Flexdashboard RMD file
# This script renders and launches the Shiny-based flexdashboard

# Load required libraries
library(rmarkdown)
library(flexdashboard)
library(shiny)

#' launch_ASUbuildR
#'
#' This function renders and launches a Shiny-based R Markdown Flexdashboard.
#' The dashboard will automatically open in the default web browser.
#'
#' @param rmd_file Character string specifying the path to the R Markdown file.
#'   Default is "ASU Flexdashboard.Rmd".
#'
#' @return This function is called for its side effects (launching the dashboard).
#'   It does not return a value.
#'
#' @examples
#' \dontrun{
#' # Run dashboard with default file name
#' run_asu_dashboard()
#'
#' # Run dashboard with custom file name
#' launch_ASUbuildR("MyCustomDashboard.Rmd")
#' }
#'
#' @seealso \code{\link[rmarkdown]{run}} for more details on running Shiny documents
#'
#' @export
launch_ASUbuildR <- function(rmd_file = "R/ASU Flexdashboard - mapgl.Rmd") {

  # Check if the RMD file exists
  if (!file.exists(rmd_file)) {
    stop(paste("File", rmd_file, "not found in current directory"))
  }

  # Print status message
  cat("Starting ASU Flexdashboard...\n")
  cat("File:", rmd_file, "\n")

  # Run the flexdashboard with shiny
  # This will render and launch the dashboard in the browser
  rmarkdown::run(
    file = rmd_file,
    shiny_args = list(
      host = "127.0.0.1",  # localhost
      port = NULL,         # Let R choose an available port
      launch.browser = TRUE # Automatically open in browser
    )
  )
}
