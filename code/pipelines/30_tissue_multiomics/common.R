get_project_root <- function() {
  root <- Sys.getenv("LIUFANG_ROOT", unset = "")
  if (nzchar(root)) {
    return(normalizePath(root, mustWork = TRUE))
  }
  normalizePath(getwd(), mustWork = TRUE)
}

project_root <- get_project_root()
source(file.path(project_root, "code", "lib", "R", "deterministic_rng.R"))
result_root <- file.path(project_root, "results", "30_tissue_multiomics")
data_root <- file.path(
  project_root, "results", "00_preprocess", "multiomics", "dat"
)
reference_root <- file.path(project_root, "data", "reference")
audit_root <- file.path(project_root, "audit")

dir.create(result_root, recursive = TRUE, showWarnings = FALSE)
dir.create(audit_root, recursive = TRUE, showWarnings = FALSE)

master_seed <- validate_r_seed(
  Sys.getenv("LIUF_MASTER_SEED", unset = "20260807"),
  "LIUF_MASTER_SEED"
)
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
set.seed(master_seed)

result_path <- function(...) {
  path <- file.path(result_root, ...)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  path
}

data_path <- function(name) {
  file.path(data_root, name)
}

reference_path <- function(...) {
  file.path(reference_root, ...)
}

assert_files <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0L) {
    stop("Missing required local input(s): ", paste(missing, collapse = ", "))
  }
  invisible(paths)
}

paired_tissue_metadata <- function() {
  sample_sheet <- data.table::fread(data_path("sample_sheet.csv"))
  tissue <- sample_sheet[source == "tissue"]
  patients <- intersect(
    tissue[group == "HCC_tumor", patient_num],
    tissue[group == "Adjacent", patient_num]
  )
  list(
    sample_sheet = sample_sheet,
    tissue = tissue,
    patients = sort(patients)
  )
}

build_gene_matrix <- function(mat, symbols, take = "max") {
  keep <- !is.na(symbols) & symbols != "" & rowMeans(!is.na(mat)) >= 0.5
  mat <- mat[keep, , drop = FALSE]
  symbols <- symbols[keep]
  mat <- log2(mat + 1)
  medians <- apply(mat, 2, median, na.rm = TRUE)
  mat <- sweep(mat, 2, medians - median(medians, na.rm = TRUE), "-")
  for (index in seq_len(nrow(mat))) {
    missing <- is.na(mat[index, ])
    if (any(missing) && !all(missing)) {
      mat[index, missing] <- min(mat[index, ], na.rm = TRUE)
    }
  }
  mat[is.na(mat)] <- 0
  if (identical(take, "max")) {
    order_index <- order(-rowMeans(mat))
  } else {
    order_index <- order(-apply(mat, 1, function(values) {
      max(abs(values), na.rm = TRUE)
    }))
  }
  mat <- mat[order_index, , drop = FALSE]
  symbols <- symbols[order_index]
  unique_index <- !duplicated(symbols)
  mat <- mat[unique_index, , drop = FALSE]
  rownames(mat) <- symbols[unique_index]
  mat
}

paired_activity_test <- function(mat, tissue_metadata) {
  output <- data.table::data.table()
  for (item in rownames(mat)) {
    values <- data.table::data.table(
      sample_id = colnames(mat),
      value = as.numeric(mat[item, ])
    )
    values <- merge(
      values,
      tissue_metadata[, .(sample_id, group_4, pair = paired_id)],
      by = "sample_id"
    )
    wide <- data.table::dcast(
      values, pair ~ group_4, value.var = "value", fun.aggregate = mean
    )
    if (!all(c("HCC", "Adjacent") %in% colnames(wide))) {
      next
    }
    wide <- wide[!is.na(HCC) & !is.na(Adjacent)]
    if (nrow(wide) < 5L) {
      next
    }
    test <- t.test(wide$HCC, wide$Adjacent, paired = TRUE)
    output <- rbind(
      output,
      data.table::data.table(
        item = item,
        n = nrow(wide),
        delta = mean(wide$HCC - wide$Adjacent),
        t = as.numeric(test$statistic),
        p = test$p.value
      )
    )
  }
  output[, p_adj := p.adjust(p, method = "BH")]
  output
}
