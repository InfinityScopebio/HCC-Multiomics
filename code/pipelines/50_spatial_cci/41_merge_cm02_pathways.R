#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
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
pathway_path <- assert_production_path(
  file.path(result_root, "pathways", "pathway_scores_raw_visium.rds")
)
recovery_path <- assert_production_path(
  file.path(result_root, "cm_recovery", "CM_recover.rds")
)
output_root <- file.path(result_root, "pathways")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
dir.create(spatial_audit_root(), recursive = TRUE, showWarnings = FALSE)
if (!file.exists(pathway_path) || !file.exists(recovery_path)) {
  stop("Pathway scores and CM recovery are both required.")
}
scores <- readRDS(pathway_path)
recovery <- attach_historical_recovery_barcodes(readRDS(recovery_path))
components <- cm02_rctd_types()

merge_section <- function(section) {
  weight_path <- active_rctd_weight_path(section)
  weights <- readRDS(weight_path)
  section_recovery <- recovery[
    recovery$section == section,
    c("barcode", "CM01", "CM02", "CM03", "CM04", "CM05", "maxCM"),
    drop = FALSE
  ]
  section_scores <- scores[scores$section == section, , drop = FALSE]
  merged <- merge(
    section_scores,
    section_recovery,
    by = "barcode",
    sort = TRUE
  )
  merged <- merge(
    merged,
    weights[, c("barcode", components), drop = FALSE],
    by = "barcode",
    sort = TRUE
  )
  merged$section <- section
  merged$is_tumor_scope <- section %in% tumor_sections()
  merged$is_cm02_positive <- merged$maxCM == "CM02"
  merged
}

merged <- bind_rows(lapply(all_sections(), merge_section))
saveRDS(merged, file.path(output_root, "pathway_scores_with_cm02.rds"))
fwrite(
  as.data.table(merged),
  file.path(output_root, "pathway_scores_with_cm02.csv.gz")
)
primary <- merged %>%
  filter(section %in% tumor_sections(), maxCM == "CM02")
fwrite(
  as.data.table(primary),
  file.path(output_root, "tumor_cm02_positive_pathway_scores.csv.gz")
)

programs <- c("HYPOXIA", "MAPK_ERK", "mTOR_signal")
program_labels <- c(
  HYPOXIA = "Hypoxia",
  MAPK_ERK = "MAPK/ERK (ERK1/2)",
  mTOR_signal = "PI3K-AKT-mTOR"
)
section_summary <- rbindlist(
  lapply(
    programs,
    function(program) {
      table <- as.data.table(merged)[section %in% tumor_sections()]
      table[
        ,
        value_z := as.numeric(scale(get(program))),
        by = section
      ]
      summary <- table[
        ,
        .(value = median(value_z, na.rm = TRUE)),
        by = .(
          section,
          group = ifelse(
            maxCM == "CM02",
            "CM02-positive",
            "CM02-negative"
          )
        )
      ]
      summary[, program := program_labels[[program]]]
      summary
    }
  )
)
statistics <- section_summary[
  ,
  .(
    p_value = wilcox.test(value ~ group, exact = FALSE)$p.value,
    n_negative = sum(group == "CM02-negative"),
    n_positive = sum(group == "CM02-positive")
  ),
  by = program
]
fwrite(
  section_summary,
  file.path(output_root, "tumor_cm02_pathway_section_summary.tsv"),
  sep = "\t"
)
fwrite(
  statistics,
  file.path(output_root, "tumor_cm02_pathway_statistics.tsv"),
  sep = "\t"
)

cat(
  "Merged pathway/CM02 output: ",
  nrow(merged),
  " spots; primary tumor CM02-positive n=",
  nrow(primary),
  ".\n",
  sep = ""
)
