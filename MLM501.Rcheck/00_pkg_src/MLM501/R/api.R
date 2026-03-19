#' Read the multilevel effects table
#' @return data.frame with columns review_id, analysis_id, study_id, TE, seTE, etc.
#' @export
read_mlm_effects <- function() {
  # Try local path first (for development/scripts)
  p <- "inst/extdata/mlm_effects.csv"
  if (!file.exists(p)) {
    # Fallback to system file for installed package
    p <- system.file("extdata", "mlm_effects.csv", package = "MLM501")
  }
  
  if (!nzchar(p) || !file.exists(p)) {
    stop("mlm_effects.csv not found. Run data-raw/import_cochrane501.R.")
  }
  
  utils::read.csv(p, stringsAsFactors = FALSE)
}
