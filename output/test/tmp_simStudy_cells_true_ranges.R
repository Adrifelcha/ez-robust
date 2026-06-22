# Temporary helper: scan all .RData under demos/simStudy_cells (four factorial cells)
# and report empirical ranges of true parameters.
#
# Population-level entries in true.values / simStudy_Beta$true: bound_mean, drift_mean,
# nondt_mean, *_sdev, betaweight (one value per simulated dataset).
#
# If *_sdev in true.values shows a range across files, it is NOT from "random SD draws"
# when settings supply fixed scalars (length != 2). More common causes:
#   (1) Pooling .RData from runs with different settings (e.g. old vs new simulation_settings).
#   (2) Missing true_sdevs component: in sample_hierarchical_parameters(), NULL
#       true_sdevs$bound_sdev (etc.) triggers fallback bound_mean/5 (or drift_mean/5,
#       nondt_mean/5), which DOES vary across datasets because the means are random.
#   (3) length-2 vector in true_sdevs$* -> runif between endpoints each dataset.
# Participant-level vectors: bound, drift, nondt (length = nPart). Each element is one
# participant's true parameter drawn from the hierarchical model (e.g. bound is
# truncnorm with mean bound_mean, sd bound_sdev, lower limit 0). Min/max across those
# vectors can therefore run from values arbitrarily close to 0 up to an upper tail
# well above the group mean (e.g. ~6 for boundary), which is not a contradiction with
# bound_mean living in a narrow uniform range for a cell study.
#
# Two on-disk layouts are supported (same as bias / RMSE code paths):
#
# 1) Collated files (e.g. under output/RData/cell-simulation/.../EZ_clean/):
#    Object: simStudy_Beta
#    True values: simStudy_Beta$true (matrix or data.frame), one row per simulated
#    dataset, columns = jags parameter names (bound_mean, drift_mean, nondt_mean,
#    betaweight, ...). Built by process_sim_data_memorySafe.R and store_simStudyResults.R
#    from each cell's x$true.values restricted to param_names.
#
# 2) Raw per-seed files (e.g. demos/simStudy_cells/<cell>/samples/seed-*.RData):
#    Object: output (simStudy_runFullSeed.R)
#    True values: nested under output$fixedEffect, output$noEffect, output$betaEffect.
#    Each block is a matrix/data.frame whose columns are conditions (EZ_clean, ...);
#    each cell holds one list from load_JAGS_cellResults -> jags_localResults (load_runJAGS.R),
#    with ground truth in element $true.values (full parameter_set list from
#    get_simulation_parameters: *_mean, *_sdev, betaweight, and participant vectors
#    bound, drift, nondt).
#
# Run from repo root after library(here). Delete this file when you no longer need it.

library(here)

cells_root <- here("demos", "simStudy_cells")
if (!dir.exists(cells_root)) {
  stop("Not found: ", cells_root)
}

cell_dirs <- list.dirs(cells_root, recursive = FALSE, full.names = TRUE)
cell_dirs <- cell_dirs[dir.exists(cell_dirs)]
cell_dirs <- cell_dirs[grepl("Drift", basename(cell_dirs), fixed = TRUE) &
  grepl("Bound", basename(cell_dirs), fixed = TRUE)]

if (length(cell_dirs) == 0) {
  stop("No cell subfolders matching *Drift*/*Bound* under ", cells_root)
}

merge_range_list <- function(acc, chunk) {
  if (is.null(chunk)) {
    return(acc)
  }
  for (nm in names(chunk)) {
    if (is.null(acc[[nm]])) {
      acc[[nm]] <- chunk[[nm]]
    } else {
      acc[[nm]]["min"] <- min(acc[[nm]]["min"], chunk[[nm]]["min"])
      acc[[nm]]["max"] <- max(acc[[nm]]["max"], chunk[[nm]]["max"])
      acc[[nm]]["n"] <- acc[[nm]]["n"] + chunk[[nm]]["n"]
    }
  }
  acc
}

ranges_to_df <- function(rlist) {
  if (length(rlist) == 0) {
    return(data.frame(parameter = character(0), min = numeric(0), max = numeric(0),
                      n_values = numeric(0), stringsAsFactors = FALSE))
  }
  params <- names(rlist)
  data.frame(
    parameter = params,
    min = vapply(rlist, function(x) x["min"], numeric(1)),
    max = vapply(rlist, function(x) x["max"], numeric(1)),
    n_values = vapply(rlist, function(x) x["n"], numeric(1)),
    stringsAsFactors = FALSE
  )
}

# Console printing: keep full numeric storage in ranges_to_df(); display min/max
# without scientific notation (default print() can switch to e-notation for wide scale).
print_ranges_dataframe <- function(df) {
  if (nrow(df) == 0) {
    cat("  (empty)\n")
    return(invisible(NULL))
  }
  fx <- df
  for (nm in intersect(c("min", "max"), names(fx))) {
    fx[[nm]] <- vapply(fx[[nm]], function(z) {
      format(as.numeric(z), scientific = FALSE, digits = 14, trim = TRUE, nsmall = 0)
    }, character(1))
  }
  print(fx, row.names = FALSE, quote = FALSE)
  invisible(NULL)
}

collect_from_simStudy_Beta <- function(sim_beta) {
  if (is.null(sim_beta$true) ||
        (!is.matrix(sim_beta$true) && !is.data.frame(sim_beta$true))) {
    return(NULL)
  }
  true_df <- as.data.frame(sim_beta$true, stringsAsFactors = FALSE)
  out <- list()
  for (nm in names(true_df)) {
    v <- suppressWarnings(as.numeric(true_df[[nm]]))
    v <- v[!is.na(v)]
    if (length(v) == 0) {
      next
    }
    out[[nm]] <- c(min = min(v), max = max(v), n = length(v))
  }
  out
}

# One JAGS cell result: list with $true.values = parameter_set from get_simulation_parameters.
ranges_from_true_values_list <- function(tv) {
  if (is.null(tv)) {
    return(NULL)
  }
  if (!is.list(tv)) {
    tv <- tryCatch(as.list(tv), error = function(e) NULL)
  }
  if (is.null(tv)) {
    return(NULL)
  }
  out <- list()
  for (nm in names(tv)) {
    val <- tv[[nm]]
    if (!is.numeric(val)) {
      next
    }
    v <- as.numeric(val)
    v <- v[!is.na(v)]
    if (length(v) == 0) {
      next
    }
    out[[nm]] <- c(min = min(v), max = max(v), n = length(v))
  }
  out
}

ranges_from_one_result_cell <- function(cell) {
  if (is.null(cell) || !is.list(cell)) {
    return(NULL)
  }
  ranges_from_true_values_list(cell$true.values)
}

# Walk output$fixedEffect / noEffect / betaEffect (matrix or DF of list-cells).
collect_from_seed_output <- function(output) {
  eff_names <- intersect(c("fixedEffect", "noEffect", "betaEffect"), names(output))
  acc <- list()
  for (en in eff_names) {
    eff <- output[[en]]
    if (is.null(eff)) {
      next
    }
    if (is.data.frame(eff)) {
      for (j in seq_len(ncol(eff))) {
        colv <- eff[[j]]
        for (r in seq_along(colv)) {
          acc <- merge_range_list(acc, ranges_from_one_result_cell(colv[[r]]))
        }
      }
    } else if (is.matrix(eff)) {
      for (i in seq_len(nrow(eff))) {
        for (j in seq_len(ncol(eff))) {
          acc <- merge_range_list(acc, ranges_from_one_result_cell(eff[[i, j]]))
        }
      }
    }
  }
  if (length(acc) == 0) {
    NULL
  } else {
    acc
  }
}

collect_true_ranges_from_file <- function(path) {
  e <- new.env(parent = emptyenv())
  load(path, envir = e)
  if (exists("simStudy_Beta", envir = e, inherits = FALSE)) {
    sim_beta <- get("simStudy_Beta", envir = e, inherits = FALSE)
    return(collect_from_simStudy_Beta(sim_beta))
  }
  if (exists("output", envir = e, inherits = FALSE)) {
    output <- get("output", envir = e, inherits = FALSE)
    return(collect_from_seed_output(output))
  }
  NULL
}

cat("Cell folders found:\n")
print(basename(cell_dirs))
cat("\n")

# --- Unique P and T: from output$settings and from embedded p/t in each JAGS cell result ---
extract_pt_from_output <- function(output) {
  sp <- unique(as.numeric(output$settings$participant_levels))
  st <- unique(as.numeric(output$settings$trial_levels))
  ep <- numeric(0)
  et <- numeric(0)
  for (en in intersect(c("fixedEffect", "noEffect", "betaEffect"), names(output))) {
    eff <- output[[en]]
    if (is.null(eff)) {
      next
    }
    walk_cell <- function(cell) {
      if (!is.null(cell) && is.list(cell)) {
        if (!is.null(cell$p)) {
          ep <<- c(ep, as.numeric(cell$p))
        }
        if (!is.null(cell$t)) {
          et <<- c(et, as.numeric(cell$t))
        }
      }
    }
    if (is.data.frame(eff)) {
      for (j in seq_len(ncol(eff))) {
        colv <- eff[[j]]
        for (r in seq_along(colv)) {
          walk_cell(colv[[r]])
        }
      }
    } else if (is.matrix(eff)) {
      for (i in seq_len(nrow(eff))) {
        for (j in seq_len(ncol(eff))) {
          walk_cell(eff[[i, j]])
        }
      }
    }
  }
  list(
    settings_P = sort(unique(sp[!is.na(sp)])),
    settings_T = sort(unique(st[!is.na(st)])),
    empirical_P = sort(unique(ep[!is.na(ep)])),
    empirical_T = sort(unique(et[!is.na(et)]))
  )
}

load_output_from_rdata <- function(path) {
  e <- new.env(parent = emptyenv())
  load(path, envir = e)
  if (!exists("output", envir = e, inherits = FALSE)) {
    return(NULL)
  }
  get("output", envir = e, inherits = FALSE)
}

# First nested JAGS result list that carries true.values (for one-off sanity print).
first_jags_result_cell <- function(output) {
  for (en in intersect(c("fixedEffect", "noEffect", "betaEffect"), names(output))) {
    eff <- output[[en]]
    if (is.null(eff)) {
      next
    }
    if (is.data.frame(eff)) {
      for (j in seq_len(ncol(eff))) {
        colv <- eff[[j]]
        for (r in seq_along(colv)) {
          cell <- colv[[r]]
          if (is.list(cell) && !is.null(cell$true.values)) {
            return(cell)
          }
        }
      }
    } else if (is.matrix(eff)) {
      for (i in seq_len(nrow(eff))) {
        for (j in seq_len(ncol(eff))) {
          cell <- eff[[i, j]]
          if (is.list(cell) && !is.null(cell$true.values)) {
            return(cell)
          }
        }
      }
    }
  }
  NULL
}

cat("========== P and T by cell (first seed file per samples* folder) ==========\n")
all_settings_P <- numeric(0)
all_settings_T <- numeric(0)
all_emp_P <- numeric(0)
all_emp_T <- numeric(0)

for (cell_path in cell_dirs) {
  cell_name <- basename(cell_path)
  subdirs <- list.dirs(cell_path, recursive = FALSE, full.names = TRUE)
  sample_roots <- subdirs[grepl("^samples", basename(subdirs))]
  if (length(sample_roots) == 0) {
    cat(cell_name, ": no directories starting with \"samples\"\n", sep = "")
    next
  }
  cat(cell_name, ":\n", sep = "")
  cat("  sample roots:", paste(basename(sample_roots), collapse = ", "), "\n")
  for (sr in sample_roots) {
    seed_files <- list.files(sr, pattern = "^seed-.*\\.RData$", full.names = TRUE)
    if (length(seed_files) == 0) {
      seed_files <- list.files(sr, pattern = "\\.RData$", full.names = TRUE)
    }
    if (length(seed_files) == 0) {
      cat("    ", basename(sr), ": (no .RData)\n", sep = "")
      next
    }
    out <- tryCatch(load_output_from_rdata(seed_files[1]), error = function(e) NULL)
    if (is.null(out)) {
      cat("    ", basename(sr), ": could not read output from first file\n", sep = "")
      next
    }
    pt <- extract_pt_from_output(out)
    cat("    ", basename(sr), ":\n", sep = "")
    cat("      settings participant_levels:", paste(pt$settings_P, collapse = ", "), "\n")
    cat("      settings trial_levels:      ", paste(pt$settings_T, collapse = ", "), "\n")
    cat("      empirical p in results:     ", paste(pt$empirical_P, collapse = ", "), "\n")
    cat("      empirical t in results:     ", paste(pt$empirical_T, collapse = ", "), "\n")
    all_settings_P <- c(all_settings_P, pt$settings_P)
    all_settings_T <- c(all_settings_T, pt$settings_T)
    all_emp_P <- c(all_emp_P, pt$empirical_P)
    all_emp_T <- c(all_emp_T, pt$empirical_T)
  }
  cat("\n")
}

cat("---------- Unique P and T pooled across all cells (seed-based scan) ----------\n")
cat("Unique P from settings:", paste(sort(unique(all_settings_P)), collapse = ", "), "\n")
cat("Unique T from settings:", paste(sort(unique(all_settings_T)), collapse = ", "), "\n")
cat("Unique P from result cells:", paste(sort(unique(all_emp_P)), collapse = ", "), "\n")
cat("Unique T from result cells:", paste(sort(unique(all_emp_T)), collapse = ", "), "\n")
cat("\n")

cat("========== Reference: settings$true_sdevs (first readable seed, any cell) ==========\n")
cat(
  "If these are fixed scalars, bound_sdev/drift_sdev/nondt_sdev in true.values should be\n",
  "constant per dataset; a wide min/max in the tables below usually means pooled files\n",
  "from different configs, or NULL components (fallback mean/5 in sample_parameters.R).\n\n",
  sep = ""
)
printed_sdev_ref <- FALSE
for (cell_path in cell_dirs) {
  if (printed_sdev_ref) {
    break
  }
  subdirs <- list.dirs(cell_path, recursive = FALSE, full.names = TRUE)
  sample_roots <- subdirs[grepl("^samples", basename(subdirs))]
  for (sr in sample_roots) {
    seed_files <- list.files(sr, pattern = "^seed-.*\\.RData$", full.names = TRUE)
    if (length(seed_files) == 0) {
      seed_files <- list.files(sr, pattern = "\\.RData$", full.names = TRUE)
    }
    if (length(seed_files) == 0) {
      next
    }
    out <- tryCatch(load_output_from_rdata(seed_files[1]), error = function(e) NULL)
    if (is.null(out)) {
      next
    }
    cat("File: ", seed_files[1], "\n", sep = "")
    ts <- out$settings$true_sdevs
    if (is.null(ts)) {
      cat("  settings$true_sdevs is NULL (fallback SDs = mean/5 can vary with drawn means).\n")
    } else {
      cat("  settings$true_sdevs:\n")
      print(ts)
    }
    printed_sdev_ref <- TRUE
    break
  }
}
if (!printed_sdev_ref) {
  cat("  (no seed output found to print settings$true_sdevs)\n")
}
cat("\n")

overall <- list()
per_cell_tables <- list()

for (cell_path in cell_dirs) {
  cell_name <- basename(cell_path)
  rdata_files <- list.files(cell_path, pattern = "\\.RData$", full.names = TRUE, recursive = TRUE)
  cat("--- ", cell_name, " ---\n", sep = "")
  cat("  .RData files:", length(rdata_files), "\n")

  cell_acc <- list()
  for (fp in rdata_files) {
    chunk <- tryCatch(
      collect_true_ranges_from_file(fp),
      error = function(cond) {
        message("  skip (error): ", fp, " -> ", conditionMessage(cond))
        NULL
      }
    )
    cell_acc <- merge_range_list(cell_acc, chunk)
    overall <- merge_range_list(overall, chunk)
  }

  df_cell <- ranges_to_df(cell_acc)
  per_cell_tables[[cell_name]] <- df_cell
  if (nrow(df_cell) > 0) {
    print_ranges_dataframe(df_cell)
  } else {
    cat("  (no simStudy_Beta or seed output$*Effect with true.values)\n")
  }
  cat("\n")
}

cat("========== Pooled across all cells ==========\n")
df_all <- ranges_to_df(overall)
if (nrow(df_all) > 0) {
  print_ranges_dataframe(df_all)
} else {
  cat("No true values collected.\n")
}

cat("\n========== Sanity: stored *_sdev match generative inputs (one seed, one cell) ==========\n")
cat(
  "get_simulation_parameters() writes the same bound_sdev (etc.) into true.values that\n",
  "rtruncnorm(..., sd = bound_sdev) uses; see sample_parameters.R bound <- rtruncnorm(...).\n",
  "The wide min/max for vector \"bound\" in the big table is the GLOBAL extremum over\n",
  "millions of participant draws (order statistics), not the spread within one dataset.\n\n",
  sep = ""
)
printed_sanity <- FALSE
for (cell_path in cell_dirs) {
  if (printed_sanity) {
    break
  }
  subdirs <- list.dirs(cell_path, recursive = FALSE, full.names = TRUE)
  sample_roots <- subdirs[grepl("^samples", basename(subdirs))]
  for (sr in sample_roots) {
    seed_files <- list.files(sr, pattern = "^seed-.*\\.RData$", full.names = TRUE)
    if (length(seed_files) == 0) {
      seed_files <- list.files(sr, pattern = "\\.RData$", full.names = TRUE)
    }
    if (length(seed_files) == 0) {
      next
    }
    out <- tryCatch(load_output_from_rdata(seed_files[1]), error = function(e) NULL)
    if (is.null(out)) {
      next
    }
    cell <- first_jags_result_cell(out)
    if (is.null(cell)) {
      next
    }
    tv <- cell$true.values
    bm <- as.numeric(tv$bound_mean)
    bs <- as.numeric(tv$bound_sdev)
    bvec <- as.numeric(tv$bound)
    cat("Seed file:", seed_files[1], "\n")
    cat("  true.values$bound_mean (scalar):", format(bm, scientific = FALSE, digits = 10), "\n")
    cat("  true.values$bound_sdev (scalar):", format(bs, scientific = FALSE, digits = 10), "\n")
    cat("  true.values$bound length nPart:", length(bvec), "\n")
    cat("  min / max bound within this one dataset:", format(min(bvec), scientific = FALSE, digits = 10),
        " , ", format(max(bvec), scientific = FALSE, digits = 10), "\n")
    cat("  mean / sd of bound within this one dataset:",
        format(mean(bvec), scientific = FALSE, digits = 10), " , ",
        format(stats::sd(bvec), scientific = FALSE, digits = 10), "\n")
    cat("  rough tail check (Normal approx, conservative mu = lower end of bound_mean range):\n")
    n_pool <- if (nrow(df_all) > 0) {
      m <- df_all$n_values[df_all$parameter == "bound"]
      if (length(m) == 1) m[[1]] else NA_real_
    } else {
      NA_real_
    }
    if (!is.na(n_pool) && n_pool > 0 && length(bs) == 1 && bs > 0) {
      mu_lo <- min(as.numeric(df_all$min[df_all$parameter == "bound_mean"]), bm, na.rm = TRUE)
      cat("  pooled count of participant-level bound values in this scan:",
          format(n_pool, scientific = FALSE), "\n")
      cat("  order-stat floor ~ mu_lo - sqrt(2*log(N))*sigma with mu_lo =",
          format(mu_lo, scientific = FALSE, digits = 8), ", sigma =", format(bs, scientific = FALSE), " -> ",
          format(mu_lo - sqrt(2 * log(n_pool)) * bs, scientific = FALSE, digits = 8),
          " (same order as your global min; truncation at 0 is mild when mu is 3.5 to 4).\n",
          sep = "")
    }
    printed_sanity <- TRUE
    break
  }
}
if (!printed_sanity) {
  cat("  (could not find a seed with output + true.values for sanity block)\n")
}

cell_output <- here("output", "RData", "cell-simulation")
if (dir.exists(cell_output)) {
  cat("\n========== ", cell_output, " (collated sim_P* files) ==========\n", sep = "")
  out_cells <- list.dirs(cell_output, recursive = FALSE, full.names = TRUE)
  out_cells <- out_cells[dir.exists(out_cells)]
  out_cells <- out_cells[grepl("Drift", basename(out_cells), fixed = TRUE) &
    grepl("Bound", basename(out_cells), fixed = TRUE)]
  out_overall <- list()
  for (oc in out_cells) {
    cond_dirs <- list.dirs(oc, recursive = FALSE, full.names = TRUE)
    for (cd in cond_dirs) {
      fs <- list.files(cd, pattern = "\\.RData$", full.names = TRUE)
      for (fp in fs) {
        out_overall <- merge_range_list(out_overall, collect_true_ranges_from_file(fp))
      }
    }
  }
  print_ranges_dataframe(ranges_to_df(out_overall))
}
