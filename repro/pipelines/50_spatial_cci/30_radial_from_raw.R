#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(FNN)
  library(jsonlite)
  library(tidyr)
})

source(file.path(dirname(normalizePath(
  sub(
    "^--file=",
    "",
    commandArgs(trailingOnly = FALSE)[
      grepl("^--file=", commandArgs(trailingOnly = FALSE))
    ]
  ),
  mustWork = TRUE
)), "common.R"))

seed <- configure_reproducibility()
root <- project_root()
result_root <- spatial_result_root()
raw_root <- assert_production_path(
  file.path(root, "repo", "data", "public", "raw", "HCC-ST")
)
weight_root <- assert_production_path(file.path(result_root, "rctd", "weights"))
recovery_path <- assert_production_path(
  file.path(result_root, "cm_recovery", "CM_recover.rds")
)
output_root <- file.path(result_root, "radial")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(recovery_path)) {
  stop("Missing production CM recovery: ", recovery_path)
}
recovery <- readRDS(recovery_path)
components <- cm02_rctd_types()

radial_section <- function(section) {
  weight_path <- file.path(weight_root, paste0(section, "_weights.rds"))
  if (!file.exists(weight_path)) {
    stop("Missing production RCTD weights: ", weight_path)
  }
  weights <- readRDS(weight_path)
  missing_components <- setdiff(components, colnames(weights))
  if (length(missing_components) > 0L) {
    stop(
      section,
      ": missing CM02 RCTD components: ",
      paste(missing_components, collapse = ", ")
    )
  }
  section_recovery <- recovery[
    recovery$section == section,
    c("barcode", "CM02", "maxCM"),
    drop = FALSE
  ]
  joined <- merge(
    section_recovery,
    weights[, c("barcode", components), drop = FALSE],
    by = "barcode",
    sort = TRUE
  )
  if (nrow(joined) != nrow(section_recovery)) {
    stop(section, ": barcode join lost recovered spots.")
  }
  section_dir <- assert_production_path(file.path(raw_root, section))
  positions <- read_tissue_positions(section_dir)
  joined <- merge(
    joined,
    positions[, c("barcode", "pxl_row", "pxl_col"), drop = FALSE],
    by = "barcode",
    sort = TRUE
  )
  if (nrow(joined) != nrow(section_recovery)) {
    stop(section, ": barcode join lost tissue coordinates.")
  }
  scale <- fromJSON(
    file.path(section_dir, "spatial", "scalefactors_json.json")
  )$spot_diameter_fullres
  threshold <- quantile(joined$CM02, 0.75, na.rm = TRUE)
  joined$isCore <- joined$CM02 >= threshold & joined$maxCM == "CM02"
  core <- joined[joined$isCore, c("pxl_row", "pxl_col"), drop = FALSE]
  if (nrow(core) < 5L) {
    stop(section, ": fewer than five CM02 core spots.")
  }
  nearest <- get.knnx(
    as.matrix(core),
    as.matrix(joined[, c("pxl_row", "pxl_col"), drop = FALSE]),
    k = 1L
  )
  joined$dist_spots <- nearest$nn.dist[, 1L] / scale
  joined$section <- section
  joined$is_tumor_scope <- section %in% tumor_sections()
  joined
}

set.seed(seed)
radial <- bind_rows(lapply(all_sections(), radial_section))
radial$distance_bin <- cut(
  radial$dist_spots,
  breaks = c(-0.1, 0.1, 1, 2, 4, 8, Inf),
  labels = c("core(0)", "0-1", "1-2", "2-4", "4-8", ">8")
)
saveRDS(radial, file.path(output_root, "radial_all16_raw_rebuilt.rds"))
fwrite(as.data.table(radial), file.path(output_root, "radial_all16_raw_rebuilt.csv.gz"))

primary <- radial %>%
  filter(
    section %in% tumor_sections(),
    maxCM == "CM02",
    !isCore
  )
if (nrow(primary) == 0L) {
  stop("Primary tumor/maxCM==CM02/non-core radial scope is empty.")
}
profile <- primary %>%
  group_by(distance_bin) %>%
  summarise(across(all_of(components), mean, na.rm = TRUE), .groups = "drop")
long_profile <- profile %>%
  pivot_longer(
    all_of(components),
    names_to = "component",
    values_to = "mean_rctd_weight"
  ) %>%
  group_by(component) %>%
  mutate(z_weight = as.numeric(scale(mean_rctd_weight))) %>%
  ungroup()
fwrite(
  as.data.table(long_profile),
  file.path(output_root, "tumor_cm02_radial_profile.tsv"),
  sep = "\t"
)

distance_midpoint <- c(
  "core(0)" = 0,
  "0-1" = 0.5,
  "1-2" = 1.5,
  "2-4" = 3,
  "4-8" = 6,
  ">8" = 12
)
weighted_distance <- long_profile %>%
  mutate(distance_midpoint = distance_midpoint[as.character(distance_bin)]) %>%
  group_by(component) %>%
  summarise(
    weighted_mean_distance = sum(
      mean_rctd_weight * distance_midpoint,
      na.rm = TRUE
    ) / sum(mean_rctd_weight, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(weighted_mean_distance)
fwrite(
  as.data.table(weighted_distance),
  file.path(output_root, "tumor_cm02_component_order.tsv"),
  sep = "\t"
)
fwrite(
  as.data.table(primary),
  file.path(output_root, "tumor_cm02_noncore_spots.csv.gz")
)
cat(
  "Radial rebuild completed: ",
  nrow(radial),
  " all-section spots; ",
  nrow(primary),
  " primary tumor CM02-positive non-core spots.\n",
  sep = ""
)
