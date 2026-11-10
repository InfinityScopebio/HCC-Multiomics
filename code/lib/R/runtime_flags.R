resume_enabled <- function(name) {
  value <- Sys.getenv(name, unset = "0")
  if (!value %in% c("0", "1")) {
    stop(name, " must be 0 or 1.")
  }
  identical(value, "1")
}
