R_SEED_MAX <- 2147483647
R_SEED_MODULUS <- 2147483648

validate_r_seed <- function(seed, name = "seed") {
  if (is.character(seed)) {
    if (
      length(seed) != 1L ||
        !grepl("^(0|[1-9][0-9]*)$", seed)
    ) {
      stop(name, " must use canonical non-negative decimal notation.")
    }
  } else if (
    !is.numeric(seed) ||
      is.object(seed) ||
      is.logical(seed)
  ) {
    stop(
      name,
      " must be an unclassed integer/double or canonical character seed."
    )
  }
  numeric_seed <- suppressWarnings(as.numeric(seed))
  if (
    length(numeric_seed) != 1L ||
      is.na(numeric_seed) ||
      !is.finite(numeric_seed) ||
      numeric_seed != floor(numeric_seed) ||
      numeric_seed < 0 ||
      numeric_seed > R_SEED_MAX
  ) {
    stop(
      name,
      " must be an integer in [0, ",
      format(R_SEED_MAX, scientific = FALSE),
      "]."
    )
  }
  as.integer(numeric_seed)
}

derive_r_seed <- function(base_seed, offset = 0, name = "derived seed") {
  validated_base <- validate_r_seed(base_seed, "base seed")
  if (
    !is.numeric(offset) ||
      is.object(offset) ||
      is.logical(offset)
  ) {
    stop(name, " offset must be an unclassed integer or double.")
  }
  numeric_offset <- suppressWarnings(as.numeric(offset))
  if (
    length(numeric_offset) != 1L ||
      is.na(numeric_offset) ||
      !is.finite(numeric_offset) ||
      numeric_offset != floor(numeric_offset) ||
      abs(numeric_offset) > (2^53 - 1)
  ) {
    stop(
      name,
      " offset must be an exactly representable IEEE-754 safe integer."
    )
  }
  offset_modulus <- numeric_offset %% R_SEED_MODULUS
  derived <- (
    as.numeric(validated_base) + offset_modulus
  ) %% R_SEED_MODULUS
  validate_r_seed(derived, name)
}
