# Infer simulation layout and condition folders from main_dir.
# Detects either:
# 1) direct/full layout: main_dir contains condition folders
# 2) nested/cell layout: main_dir contains parameter-cell folders that contain condition folders
read_conditions <- function(main_dir) {
  preferred_conditions <- c("EZ_contaminated", "EZRobust_contaminated", "EZ_clean", "EZRobust_clean")
  condition_label_map <- c(
    "EZ_contaminated" = "EZ x Contaminated",
    "EZRobust_contaminated" = "Robust x Contaminated",
    "EZ_clean" = "EZ x Clean",
    "EZRobust_clean" = "Robust x Clean"
  )
  condition_dir_pattern <- "^(EZ|EZRobust)_(clean|contaminated)$"

  top_level_dirs <- list.dirs(main_dir, recursive = FALSE, full.names = FALSE)
  direct_conditions <- top_level_dirs[grepl(condition_dir_pattern, top_level_dirs)]

  if (length(direct_conditions) > 0) {
    conditions <- read_conditions_order(unique(direct_conditions), preferred_conditions)
    return(list(
      layout = "direct",
      conditions = conditions,
      condition_labels = read_conditions_labels(conditions, condition_label_map),
      parameter_cells = character(0),
      message = "Detected direct simulation layout (condition folders found in main_dir)."
    ))
  }

  # Nested layout: parameter-cell folders that contain condition folders
  cell_name_pattern <- "(Drift|Bound)"
  candidate_cells <- top_level_dirs[grepl(cell_name_pattern, top_level_dirs)]
  inferred_cells <- c()
  inferred_conditions <- c()

  for (cell_name in candidate_cells) {
    cell_path <- file.path(main_dir, cell_name)
    nested_dirs <- list.dirs(cell_path, recursive = FALSE, full.names = FALSE)
    valid_conditions <- nested_dirs[grepl(condition_dir_pattern, nested_dirs)]
    if (length(valid_conditions) > 0) {
      inferred_cells <- c(inferred_cells, cell_name)
      inferred_conditions <- c(inferred_conditions, valid_conditions)
    }
  }

  if (length(inferred_cells) > 0) {
    preferred_cells <- c("lowDrift-lowBound", "lowDrift-highBound", "highDrift-lowBound", "highDrift-highBound")
    parameter_cells <- unique(inferred_cells)
    parameter_cells <- c(intersect(preferred_cells, parameter_cells), setdiff(sort(parameter_cells), preferred_cells))
    conditions <- read_conditions_order(unique(inferred_conditions), preferred_conditions)

    return(list(
      layout = "nested",
      conditions = conditions,
      condition_labels = read_conditions_labels(conditions, condition_label_map),
      parameter_cells = unname(parameter_cells),
      message = "Detected nested simulation layout (parameter-cell folders containing condition folders)."
    ))
  }

  stop(
    paste(
      "Could not detect a valid simulation layout in main_dir.",
      "Expected either direct condition folders (e.g., EZ_clean) or nested parameter-cell folders",
      "(e.g., names containing highBound/lowBound) that contain condition folders."
    )
  )
}

# Secondary helper: map condition IDs to display labels.
read_conditions_labels <- function(conditions, condition_label_map) {
  labels <- condition_label_map[conditions]
  missing_labels <- is.na(labels)
  if (any(missing_labels)) {
    labels[missing_labels] <- gsub("_", " x ", conditions[missing_labels], fixed = TRUE)
  }
  unname(labels)
}

# Secondary helper: keep preferred condition order when available.
read_conditions_order <- function(inferred_conditions, preferred_conditions) {
  out <- intersect(preferred_conditions, inferred_conditions)
  if (length(out) == 0) {
    out <- sort(inferred_conditions)
  }
  unname(out)
}
