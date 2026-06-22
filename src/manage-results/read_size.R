# Infer participant/trial size levels from simulation result filenames.
# Expected filename pattern includes: sim_P{p}T{t}_...RData
read_size <- function(parent_dirs) {
  all_files <- c()
  for (dir_path in unique(parent_dirs)) {
    if (!dir.exists(dir_path)) {
      next
    }
    files_i <- list.files(dir_path, pattern = "\\.RData$", full.names = TRUE)
    all_files <- c(all_files, files_i)
  }
  all_files <- unique(all_files)
  rdata_files <- all_files[grepl("_P\\d+T\\d+_", basename(all_files))]
  
  if (length(rdata_files) == 0) {
    stop("No .RData files with pattern '_P{p}T{t}_' were found in the provided directories.")
  }
  
  filenames <- basename(rdata_files)
  p_values <- as.numeric(sub(".*_P(\\d+)T.*", "\\1", filenames))
  t_values <- as.numeric(sub(".*T(\\d+)_.*", "\\1", filenames))
  
  p_levels <- sort(unique(p_values[!is.na(p_values)]))
  t_levels <- sort(unique(t_values[!is.na(t_values)]))
  
  if (length(p_levels) == 0 || length(t_levels) == 0) {
    stop("Could not infer participant or trial levels from result filenames.")
  }
  
  return(list(
    p_levels = p_levels,
    t_levels = t_levels
  ))
}
