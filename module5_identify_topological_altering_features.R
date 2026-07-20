#!/usr/bin/env Rscript

################################################################################
# Module 5: Topological-altering feature identification
#
# Purpose
#   Identify topological-altering features (TAFs; e.g., genes, CpG sites, DMRs,
#   peaks, proteins, metabolites) by testing whether the observed change in
#   participation score between phenotype2 and phenotype1 is unusually large
#   compared with degree-matched background features.
#
# Expected upstream inputs
#   1) Module 4 delta participation score table or Module 4 manifest
#   2) Module 1 output manifest containing phenotype-specific network matrices
#
# Main outputs
#   <prefix>_module5_network_degree_table.tsv
#   <prefix>_module5_degree_matched_null_result.tsv
#   <prefix>_module5_topological_altering_features.tsv
#   <prefix>_module5_summary.tsv
#   <prefix>_module5_volcano_data.tsv
#   plots/<prefix>_module5_volcano_plot.png and/or .pdf
#   <prefix>_module5_output_manifest.tsv
#   <prefix>_module5_run_info.tsv
#
# Design principles
#   - No hard-coded paths
#   - Generic phenotype1 / phenotype2 labels
#   - Modular functions
#   - Clear input validation and error messages
#   - Degree-matched empirical testing separated from Module 4 score calculation
################################################################################

options(stringsAsFactors = FALSE)

##### 1. Utilities #############################################################

message_info <- function(...) message("[INFO] ", paste0(..., collapse = ""))
message_warn <- function(...) warning(paste0(..., collapse = ""), call. = FALSE)
stop_msg <- function(...) stop(paste0(..., collapse = ""), call. = FALSE)

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) == 0 || any(args %in% c("-h", "--help"))) {
    cat(
"Module 5: Topological-altering feature identification\n\n",
"Usage:\n",
"  Rscript module5_identify_topological_altering_features.R \\\n",
"    --module4-manifest module4_output/demo_module5_input_manifest.tsv \\\n",
"    --module1-manifest module1_output/demo_module1_output_manifest.tsv \\\n",
"    --outdir module5_output \\\n",
"    --prefix demo \\\n",
"    --phenotype1 phenotype1 \\\n",
"    --phenotype2 phenotype2 \\\n",
"    --degree-matrix-type topPct \\\n",
"    --degree-threshold 5\n\n",
"Alternative using explicit delta table:\n",
"  Rscript module5_identify_topological_altering_features.R \\\n",
"    --delta-score module4_output/demo_module4_delta_participation_score.tsv \\\n",
"    --module1-manifest module1_output/demo_module1_output_manifest.tsv \\\n",
"    --outdir module5_output \\\n",
"    --prefix demo \\\n",
"    --phenotype1 phenotype1 \\\n",
"    --phenotype2 phenotype2 \\\n",
"    --degree-matrix-type topPct \\\n",
"    --degree-threshold 5\n\n",
"Required arguments:\n",
"  --module1-manifest       Module 1 output manifest TSV\n",
"  --outdir                 Output directory\n",
"  --prefix                 Output prefix\n",
"  --phenotype1             Baseline phenotype label\n",
"  --phenotype2             Comparison phenotype label\n",
"  --degree-threshold       Threshold for network degree calculation\n",
"                           If --degree-matrix-type topPct, use e.g. 5 for top 5% strongest edges.\n",
"                           If --degree-matrix-type distance, edges are 0 < distance <= threshold.\n\n",
"Input arguments, choose one:\n",
"  --delta-score            Module 4 delta participation score TSV\n",
"  --module4-manifest       Module 4 manifest containing delta_participation_score path\n\n",
"Optional arguments:\n",
"  --degree-matrix-type     topPct or distance [default: topPct]\n",
"  --score-col              Delta score column to test [default: delta_participation_score_lifespan_sum]\n",
"  --degree-reference       avg, phenotype1, or phenotype2 [default: avg]\n",
"  --min-matched            Minimum degree-matched background size [default: 20]\n",
"  --max-window             Maximum degree-matching window [default: 50]\n",
"  --alternative            two.sided, greater, or less [default: two.sided]\n",
"  --q-cutoff               BH q-value cutoff for TAFs [default: 0.05]\n",
"  --p-cutoff               Optional raw p-value cutoff; NA disables [default: NA]\n",
"  --min-abs-delta          Minimum absolute delta score for TAFs [default: 0]\n",
"  --exclude-zero-degree    TRUE/FALSE. Exclude degree-zero features from empirical test [default: FALSE]\n",
"  --pseudocount            Empirical p-value pseudocount [default: 1]\n\n",
"Volcano plot arguments:\n",
"  --make-volcano           TRUE/FALSE [default: TRUE]\n",
"  --volcano-format         png, pdf, or both [default: both]\n",
"  --volcano-p-cutoff       Raw empirical p-value cutoff used for plot colors [default: 0.05]\n",
"  --volcano-delta-cutoff   Absolute delta-score cutoff used for plot colors\n",
"                           [default: the value of --min-abs-delta]\n",
"  --volcano-width          Figure width in inches [default: 8]\n",
"  --volcano-height         Figure height in inches [default: 6]\n",
"  --volcano-dpi            PNG resolution [default: 300]\n",
"  --volcano-point-size     Point-size multiplier [default: 0.65]\n\n",
sep = "")
    quit(status = 0)
  }

  out <- list()
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!grepl("^--", key)) stop_msg("Unexpected argument: ", key)
    name <- sub("^--", "", key)
    if (i == length(args) || grepl("^--", args[[i + 1]])) {
      out[[name]] <- TRUE
      i <- i + 1
    } else {
      out[[name]] <- args[[i + 1]]
      i <- i + 2
    }
  }
  out
}

get_arg <- function(args, name, default = NULL, required = FALSE) {
  value <- args[[name]]
  if (is.null(value)) {
    if (required) stop_msg("Missing required argument: --", name)
    return(default)
  }
  value
}

as_logical_arg <- function(x, default = FALSE) {
  if (is.null(x)) return(default)
  if (is.logical(x)) return(x)
  x <- tolower(as.character(x))
  if (x %in% c("true", "t", "1", "yes", "y")) return(TRUE)
  if (x %in% c("false", "f", "0", "no", "n")) return(FALSE)
  stop_msg("Cannot parse logical value: ", x)
}

normalize_path_safe <- function(path) {
  if (is.null(path)) {
    return(NULL)
  }

  path <- as.character(path)
  result <- path

  valid <- !is.na(path) & nzchar(trimws(path))

  if (any(valid)) {
    result[valid] <- normalizePath(
      path.expand(path[valid]),
      mustWork = FALSE
    )
  }

  result
}

read_tsv <- function(path) {
  if (is.null(path) || !file.exists(path)) stop_msg("File not found: ", path)
  tryCatch(
    read.delim(path, sep = "\t", header = TRUE, check.names = FALSE, quote = "", comment.char = ""),
    error = function(e) stop_msg("Failed to read TSV: ", path, "\n", e$message)
  )
}

write_tsv <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(df, file = path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

require_columns <- function(df, cols, label = "table") {
  missing <- setdiff(cols, colnames(df))
  if (length(missing) > 0) {
    stop_msg("Missing required column(s) in ", label, ": ", paste(missing, collapse = ", "))
  }
}

sanitize_label <- function(x) {
  x <- as.character(x)
  x <- gsub("[^A-Za-z0-9_.-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  ifelse(nchar(x) == 0, "phenotype", x)
}

##### 2. Manifest resolution ##################################################

resolve_delta_score_path <- function(delta_score, module4_manifest) {
  if (!is.null(delta_score)) return(normalize_path_safe(delta_score))
  if (is.null(module4_manifest)) {
    stop_msg("Please provide either --delta-score or --module4-manifest")
  }
  manifest <- read_tsv(module4_manifest)
  if (all(c("item", "path") %in% colnames(manifest))) {
    idx <- which(manifest$item %in% c("delta_participation_score", "module4_delta_participation_score"))
    if (length(idx) == 0) idx <- grep("delta.*participation.*score", manifest$item, ignore.case = TRUE)
    if (length(idx) > 0) return(normalize_path_safe(manifest$path[[idx[[1]]]]))
  }
  path_cols <- colnames(manifest)[grepl("path|file", colnames(manifest), ignore.case = TRUE)]
  for (pc in path_cols) {
    idx <- grep("delta.*participation.*score", manifest[[pc]], ignore.case = TRUE)
    if (length(idx) > 0) return(normalize_path_safe(manifest[[pc]][[idx[[1]]]]))
  }
  stop_msg("Could not identify delta participation score path from Module 4 manifest: ", module4_manifest)
}

resolve_network_paths <- function(module1_manifest, phenotype1, phenotype2, degree_matrix_type) {
  manifest <- read_tsv(module1_manifest)
  require_columns(manifest, c("phenotype"), "Module 1 manifest")

  matrix_col <- switch(
    degree_matrix_type,
    topPct = if ("topPct_matrix" %in% colnames(manifest)) "topPct_matrix" else NA_character_,
    distance = if ("distance_matrix" %in% colnames(manifest)) "distance_matrix" else NA_character_,
    stop_msg("Unsupported degree_matrix_type: ", degree_matrix_type)
  )
  if (is.na(matrix_col)) {
    stop_msg("Module 1 manifest does not contain the required matrix column for type '", degree_matrix_type, "'.")
  }

  get_one <- function(ph) {
    idx <- which(as.character(manifest$phenotype) == ph)
    if (length(idx) == 0 && "safe_phenotype_label" %in% colnames(manifest)) {
      idx <- which(as.character(manifest$safe_phenotype_label) == sanitize_label(ph))
    }
    if (length(idx) == 0) {
      stop_msg("Phenotype not found in Module 1 manifest: ", ph,
               ". Available: ", paste(unique(manifest$phenotype), collapse = ", "))
    }
    path <- manifest[[matrix_col]][[idx[[1]]]]
    if (is.na(path) || path == "") stop_msg("Missing ", matrix_col, " path for phenotype: ", ph)
    normalize_path_safe(path)
  }

  list(
    phenotype1_matrix = get_one(phenotype1),
    phenotype2_matrix = get_one(phenotype2),
    matrix_column = matrix_col
  )
}

##### 3. Matrix and degree calculation ########################################

read_named_square_matrix <- function(path) {
  df <- read_tsv(path)
  if (ncol(df) < 2) stop_msg("Matrix file must contain a feature_id column plus numeric columns: ", path)

  feature_col <- colnames(df)[[1]]
  features <- as.character(df[[feature_col]])
  if (anyDuplicated(features)) stop_msg("Duplicated feature IDs in matrix row names: ", path)

  mat_df <- df[, -1, drop = FALSE]
  col_features <- colnames(mat_df)
  if (length(features) != length(col_features)) {
    stop_msg("Matrix is not square in file: ", path,
             ". nrow=", length(features), ", ncol=", length(col_features))
  }
  if (!identical(features, col_features)) {
    # Try to align columns if they contain the same feature set.
    if (setequal(features, col_features)) {
      mat_df <- mat_df[, features, drop = FALSE]
      col_features <- colnames(mat_df)
    } else {
      stop_msg("Matrix row and column feature IDs do not match in file: ", path)
    }
  }

  mat <- as.matrix(mat_df)
  suppressWarnings(storage.mode(mat) <- "numeric")
  if (any(is.na(mat))) stop_msg("Matrix contains NA or non-numeric values: ", path)
  if (any(!is.finite(mat))) stop_msg("Matrix contains non-finite values: ", path)
  rownames(mat) <- features
  colnames(mat) <- features
  mat
}

compute_degree_from_matrix <- function(mat, threshold, matrix_type, phenotype) {
  if (nrow(mat) != ncol(mat)) stop_msg("Degree matrix must be square for phenotype: ", phenotype)
  if (!identical(rownames(mat), colnames(mat))) {
    stop_msg("Degree matrix row/column names are not identical for phenotype: ", phenotype)
  }
  if (is.na(threshold) || !is.finite(threshold)) stop_msg("degree-threshold must be finite.")

  if (matrix_type == "topPct") {
    if (threshold < 0 || threshold > 100) stop_msg("For topPct degree matrix, --degree-threshold must be between 0 and 100.")
  }
  if (matrix_type == "distance" && threshold < 0) {
    stop_msg("For distance degree matrix, --degree-threshold must be >= 0.")
  }

  upper <- upper.tri(mat, diag = FALSE)
  edge_mask <- (mat > 0) & (mat <= threshold) & is.finite(mat)
  edge_mask[!upper] <- FALSE

  degree <- integer(nrow(mat))
  ij <- which(edge_mask, arr.ind = TRUE)
  if (nrow(ij) > 0) {
    degree <- tabulate(c(ij[, 1], ij[, 2]), nbins = nrow(mat))
  }

  data.frame(
    feature_id = rownames(mat),
    phenotype = phenotype,
    degree = as.integer(degree),
    stringsAsFactors = FALSE
  )
}

make_degree_table <- function(path1, path2, phenotype1, phenotype2, threshold, matrix_type) {
  message_info("Reading ", matrix_type, " matrix for ", phenotype1, ": ", path1)
  mat1 <- read_named_square_matrix(path1)
  message_info("Reading ", matrix_type, " matrix for ", phenotype2, ": ", path2)
  mat2 <- read_named_square_matrix(path2)

  common_features <- intersect(rownames(mat1), rownames(mat2))
  if (length(common_features) == 0) stop_msg("No common features between phenotype matrices.")
  if (length(common_features) < nrow(mat1) || length(common_features) < nrow(mat2)) {
    message_warn("Phenotype matrices do not have identical features. Using common features only: ", length(common_features))
  }
  common_features <- sort(common_features)
  mat1 <- mat1[common_features, common_features, drop = FALSE]
  mat2 <- mat2[common_features, common_features, drop = FALSE]

  d1 <- compute_degree_from_matrix(mat1, threshold, matrix_type, phenotype1)
  d2 <- compute_degree_from_matrix(mat2, threshold, matrix_type, phenotype2)
  names(d1)[names(d1) == "degree"] <- paste0("degree_", phenotype1)
  names(d2)[names(d2) == "degree"] <- paste0("degree_", phenotype2)
  d1$phenotype <- NULL
  d2$phenotype <- NULL

  out <- merge(d1, d2, by = "feature_id", all = TRUE, sort = FALSE)
  c1 <- paste0("degree_", phenotype1)
  c2 <- paste0("degree_", phenotype2)
  out[[c1]][is.na(out[[c1]])] <- 0L
  out[[c2]][is.na(out[[c2]])] <- 0L
  out[[c1]] <- as.integer(out[[c1]])
  out[[c2]] <- as.integer(out[[c2]])
  out$degree_avg <- (out[[c1]] + out[[c2]]) / 2
  out$degree_delta <- out[[c2]] - out[[c1]]
  out[order(out$feature_id), , drop = FALSE]
}

##### 4. Empirical testing ####################################################

prepare_delta_table <- function(delta_df, score_col) {
  require_columns(delta_df, c("feature_id"), "Module 4 delta score table")
  if (!score_col %in% colnames(delta_df)) {
    stop_msg("score-col not found in Module 4 delta score table: ", score_col,
             ". Available columns: ", paste(colnames(delta_df), collapse = ", "))
  }
  delta_df$feature_id <- as.character(delta_df$feature_id)
  delta_df[[score_col]] <- suppressWarnings(as.numeric(delta_df[[score_col]]))
  if (any(is.na(delta_df[[score_col]]))) {
    message_warn("Some features have NA score in ", score_col, "; they will be excluded from empirical testing.")
  }
  delta_df
}

select_degree_reference <- function(df, phenotype1, phenotype2, degree_reference) {
  c1 <- paste0("degree_", phenotype1)
  c2 <- paste0("degree_", phenotype2)
  if (degree_reference == "avg") return("degree_avg")
  if (degree_reference == "phenotype1") return(c1)
  if (degree_reference == "phenotype2") return(c2)
  stop_msg("degree-reference must be one of: avg, phenotype1, phenotype2")
}

empirical_pvalue_one <- function(target_score, matched_scores, alternative, pseudocount) {
  matched_scores <- matched_scores[is.finite(matched_scores)]
  if (length(matched_scores) == 0 || !is.finite(target_score)) return(NA_real_)

  if (alternative == "greater") {
    return((pseudocount + sum(matched_scores >= target_score)) / (pseudocount + length(matched_scores)))
  }
  if (alternative == "less") {
    return((pseudocount + sum(matched_scores <= target_score)) / (pseudocount + length(matched_scores)))
  }
  if (alternative == "two.sided") {
    center <- stats::median(matched_scores, na.rm = TRUE)
    return((pseudocount + sum(abs(matched_scores - center) >= abs(target_score - center))) /
             (pseudocount + length(matched_scores)))
  }
  stop_msg("alternative must be one of: two.sided, greater, less")
}

compute_degree_matched_test <- function(df, score_col, degree_col, min_matched, max_window,
                                        alternative, pseudocount, exclude_zero_degree) {
  require_columns(df, c("feature_id", score_col, degree_col), "merged Module 5 input table")

  test_df <- df[!is.na(df[[score_col]]) & is.finite(df[[score_col]]) &
                  !is.na(df[[degree_col]]) & is.finite(df[[degree_col]]), , drop = FALSE]
  if (exclude_zero_degree) {
    test_df <- test_df[test_df[[degree_col]] > 0, , drop = FALSE]
  }
  if (nrow(test_df) == 0) stop_msg("No valid features remain for empirical testing.")

  feature_ids <- test_df$feature_id
  scores <- test_df[[score_col]]
  degrees <- test_df[[degree_col]]
  max_degree_diff <- max(abs(outer(degrees, degrees, "-")), na.rm = TRUE)

  result <- vector("list", length(feature_ids))
  for (i in seq_along(feature_ids)) {
    target_gene <- feature_ids[[i]]
    target_score <- scores[[i]]
    target_degree <- degrees[[i]]

    window <- 0
    matched_idx <- integer(0)
    reason <- "ok"
    repeat {
      matched_idx <- which(abs(degrees - target_degree) <= window & feature_ids != target_gene)
      if (length(matched_idx) >= min_matched) break
      window <- window + 1
      if (window > max_window) {
        reason <- "insufficient_matched_within_max_window"
        break
      }
      if (window > max_degree_diff) {
        reason <- "insufficient_matched_in_dataset"
        break
      }
    }

    matched_scores <- scores[matched_idx]
    p_value <- if (length(matched_scores) >= min_matched) {
      empirical_pvalue_one(target_score, matched_scores, alternative, pseudocount)
    } else {
      NA_real_
    }

    matched_mean <- if (length(matched_scores) > 0) mean(matched_scores, na.rm = TRUE) else NA_real_
    matched_median <- if (length(matched_scores) > 0) median(matched_scores, na.rm = TRUE) else NA_real_
    matched_sd <- if (length(matched_scores) > 1) stats::sd(matched_scores, na.rm = TRUE) else NA_real_
    z_score <- if (!is.na(matched_sd) && matched_sd > 0) (target_score - matched_mean) / matched_sd else NA_real_

    result[[i]] <- data.frame(
      feature_id = target_gene,
      score_tested = target_score,
      degree_reference = target_degree,
      degree_window_used = window,
      n_degree_matched = length(matched_scores),
      matched_delta_mean = matched_mean,
      matched_delta_median = matched_median,
      matched_delta_sd = matched_sd,
      degree_matched_z = z_score,
      empirical_p_value = p_value,
      test_status = reason,
      stringsAsFactors = FALSE
    )
  }

  out <- do.call(rbind, result)
  out$q_value <- stats::p.adjust(out$empirical_p_value, method = "BH")
  out[order(out$q_value, -abs(out$score_tested), out$feature_id), , drop = FALSE]
}

identify_tafs <- function(result_df, q_cutoff, p_cutoff, min_abs_delta) {
  keep <- !is.na(result_df$q_value) & result_df$q_value <= q_cutoff & abs(result_df$score_tested) >= min_abs_delta
  if (!is.na(p_cutoff)) {
    keep <- keep & !is.na(result_df$empirical_p_value) & result_df$empirical_p_value <= p_cutoff
  }
  taf <- result_df[keep, , drop = FALSE]
  taf$direction <- ifelse(taf$score_tested > 0, "increase", ifelse(taf$score_tested < 0, "decrease", "no_change"))
  taf[order(taf$q_value, -abs(taf$score_tested), taf$feature_id), , drop = FALSE]
}

##### 5. Volcano plot #########################################################

# The volcano plot deliberately uses raw empirical p-values on the y-axis,
# matching the conventional -log10(P) definition. Its color classification is
# therefore controlled by volcano_p_cutoff and volcano_delta_cutoff. Formal TAF
# calls remain controlled separately by q_cutoff, p_cutoff, and min_abs_delta.
prepare_volcano_data <- function(result_df, volcano_p_cutoff,
                                 volcano_delta_cutoff) {
  require_columns(
    result_df,
    c("feature_id", "score_tested", "empirical_p_value", "q_value"),
    "Module 5 test result"
  )

  volcano_df <- result_df[, c(
    "feature_id", "score_tested", "empirical_p_value", "q_value"
  ), drop = FALSE]
  volcano_df$score_tested <- suppressWarnings(as.numeric(volcano_df$score_tested))
  volcano_df$empirical_p_value <- suppressWarnings(as.numeric(volcano_df$empirical_p_value))
  volcano_df$q_value <- suppressWarnings(as.numeric(volcano_df$q_value))

  valid <- is.finite(volcano_df$score_tested) &
    is.finite(volcano_df$empirical_p_value) &
    volcano_df$empirical_p_value >= 0 &
    volcano_df$empirical_p_value <= 1
  volcano_df <- volcano_df[valid, , drop = FALSE]

  if (nrow(volcano_df) == 0) {
    volcano_df$plot_p_value <- numeric()
    volcano_df$minus_log10_p <- numeric()
    volcano_df$volcano_category <- character()
    return(volcano_df)
  }

  positive_p <- volcano_df$empirical_p_value[
    volcano_df$empirical_p_value > 0 & is.finite(volcano_df$empirical_p_value)
  ]
  zero_replacement <- if (length(positive_p) > 0) {
    max(min(positive_p) / 10, .Machine$double.xmin)
  } else {
    .Machine$double.xmin
  }

  volcano_df$plot_p_value <- volcano_df$empirical_p_value
  volcano_df$plot_p_value[volcano_df$plot_p_value == 0] <- zero_replacement
  volcano_df$minus_log10_p <- -log10(volcano_df$plot_p_value)

  significant <- volcano_df$empirical_p_value < volcano_p_cutoff &
    abs(volcano_df$score_tested) >= volcano_delta_cutoff
  volcano_df$volcano_category <- "Not significant"
  volcano_df$volcano_category[
    significant & volcano_df$score_tested > 0
  ] <- "Significant increase"
  volcano_df$volcano_category[
    significant & volcano_df$score_tested < 0
  ] <- "Significant decrease"

  category_levels <- c(
    "Not significant", "Significant decrease", "Significant increase"
  )
  volcano_df$volcano_category <- factor(
    volcano_df$volcano_category,
    levels = category_levels
  )
  volcano_df[order(volcano_df$volcano_category, volcano_df$empirical_p_value), , drop = FALSE]
}

draw_volcano_plot <- function(volcano_df, volcano_p_cutoff,
                              volcano_delta_cutoff, point_size) {
  old_par <- par(no.readonly = TRUE)
  on.exit({
    layout(matrix(1))
    par(old_par)
  }, add = TRUE)

  if (nrow(volcano_df) == 0) {
    par(mar = c(5.1, 5.3, 2.0, 2.0), las = 1)
    plot.new()
    title(main = "Volcano plot")
    text(0.5, 0.5, "No valid empirical p-values available")
    return(invisible(NULL))
  }

  category_colors <- c(
    "Not significant" = grDevices::adjustcolor("black", alpha.f = 0.42),
    "Significant decrease" = grDevices::adjustcolor("#2536E8", alpha.f = 0.82),
    "Significant increase" = grDevices::adjustcolor("#FF2A1A", alpha.f = 0.82)
  )

  # Use separate plotting and legend panels. This prevents the legend from
  # covering data points or being clipped by the device boundary, independent
  # of the score range or the number of digits in category counts.
  layout(
    matrix(c(1, 2), nrow = 2, byrow = TRUE),
    heights = c(4.5, 1.5)
  )
  par(mar = c(4.8, 5.3, 1.0, 1.0), las = 1, xpd = FALSE)

  x_values <- volcano_df$score_tested
  y_values <- volcano_df$minus_log10_p
  x_range <- range(x_values, finite = TRUE)
  finite_y <- y_values[is.finite(y_values)]
  y_max <- if (length(finite_y) > 0) max(finite_y) else NA_real_
  if (!all(is.finite(x_range)) || diff(x_range) == 0) {
    center <- if (all(is.finite(x_range))) x_range[[1]] else 0
    x_range <- center + c(-1, 1)
  }
  if (!is.finite(y_max) || y_max <= 0) y_max <- 1

  x_padding <- 0.04 * diff(x_range)
  plot(
    NA_real_, NA_real_,
    xlim = x_range + c(-x_padding, x_padding),
    ylim = c(0, y_max * 1.05),
    xlab = expression(Delta~score),
    ylab = expression(-log[10](italic(P)~value)),
    bty = "l",
    xaxs = "i",
    yaxs = "i"
  )
  grid(col = "grey90", lty = 1)

  abline(
    h = -log10(volcano_p_cutoff),
    col = "grey40",
    lty = 2,
    lwd = 1
  )
  if (volcano_delta_cutoff > 0) {
    abline(
      v = c(-volcano_delta_cutoff, volcano_delta_cutoff),
      col = "grey40",
      lty = 2,
      lwd = 1
    )
  } else {
    abline(v = 0, col = "grey40", lty = 2, lwd = 1)
  }

  draw_order <- c(
    "Not significant", "Significant decrease", "Significant increase"
  )
  for (category in draw_order) {
    df <- volcano_df[as.character(volcano_df$volcano_category) == category, , drop = FALSE]
    if (nrow(df) == 0) next
    points(
      x = df$score_tested,
      y = df$minus_log10_p,
      pch = 16,
      cex = point_size,
      col = category_colors[[category]]
    )
  }

  category_counts <- table(factor(
    volcano_df$volcano_category,
    levels = draw_order
  ))
  legend_labels <- c(
    paste0(
      "Significant increase (n=",
      format(category_counts[["Significant increase"]], big.mark = ","), ")"
    ),
    paste0(
      "Significant decrease (n=",
      format(category_counts[["Significant decrease"]], big.mark = ","), ")"
    ),
    paste0(
      "Not significant (n=",
      format(category_counts[["Not significant"]], big.mark = ","), ")"
    )
  )

  significance_rule <- paste0(
    "Volcano classification: P < ", format(volcano_p_cutoff),
    if (volcano_delta_cutoff > 0) {
      paste0(" and |Delta score| >= ", format(volcano_delta_cutoff))
    } else {
      ""
    }
  )

  # Draw the legend in a dedicated lower panel rather than inside or outside
  # the scatter-plot coordinates.
  par(mar = c(0.2, 5.3, 0.2, 1.0), xpd = FALSE)
  plot.new()
  legend(
    "center",
    legend = legend_labels,
    title = significance_rule,
    col = category_colors[c(
      "Significant increase", "Significant decrease", "Not significant"
    )],
    pch = 16,
    pt.cex = 1.2,
    bty = "n",
    cex = 0.82,
    title.adj = 0
  )
  invisible(NULL)
}

write_volcano_plots <- function(volcano_df, plot_dir, prefix,
                                volcano_p_cutoff, volcano_delta_cutoff,
                                volcano_format, volcano_width,
                                volcano_height, volcano_dpi,
                                volcano_point_size) {
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  formats <- if (volcano_format == "both") c("png", "pdf") else volcano_format
  plot_paths <- character(length(formats))

  for (i in seq_along(formats)) {
    figure_format <- formats[[i]]
    figure_path <- file.path(
      plot_dir,
      paste0(prefix, "_module5_volcano_plot.", figure_format)
    )

    if (figure_format == "png") {
      grDevices::png(
        filename = figure_path,
        width = volcano_width,
        height = volcano_height,
        units = "in",
        res = volcano_dpi
      )
    } else {
      grDevices::pdf(
        file = figure_path,
        width = volcano_width,
        height = volcano_height,
        useDingbats = FALSE
      )
    }

    tryCatch(
      draw_volcano_plot(
        volcano_df = volcano_df,
        volcano_p_cutoff = volcano_p_cutoff,
        volcano_delta_cutoff = volcano_delta_cutoff,
        point_size = volcano_point_size
      ),
      finally = grDevices::dev.off()
    )

    plot_paths[[i]] <- normalize_path_safe(figure_path)
    message_info("Volcano plot: ", figure_path)
  }

  names(plot_paths) <- paste0("volcano_plot_", formats)
  plot_paths
}

make_summary <- function(delta_df, degree_df, merged_df, result_df, taf_df, phenotype1, phenotype2,
                         degree_matrix_type, degree_threshold, score_col, degree_col,
                         min_matched, max_window, alternative, q_cutoff, p_cutoff, min_abs_delta) {
  data.frame(
    item = c(
      "phenotype1",
      "phenotype2",
      "degree_matrix_type",
      "degree_threshold",
      "score_col",
      "degree_reference_col",
      "min_matched",
      "max_window",
      "alternative",
      "q_cutoff",
      "p_cutoff",
      "min_abs_delta",
      "n_features_in_delta_table",
      "n_features_in_degree_table",
      "n_features_merged",
      "n_features_tested",
      "n_features_with_valid_pvalue",
      "n_topological_altering_features",
      "n_increase_TAFs",
      "n_decrease_TAFs"
    ),
    value = c(
      phenotype1,
      phenotype2,
      degree_matrix_type,
      as.character(degree_threshold),
      score_col,
      degree_col,
      as.character(min_matched),
      as.character(max_window),
      alternative,
      as.character(q_cutoff),
      as.character(p_cutoff),
      as.character(min_abs_delta),
      as.character(nrow(delta_df)),
      as.character(nrow(degree_df)),
      as.character(nrow(merged_df)),
      as.character(nrow(result_df)),
      as.character(sum(!is.na(result_df$empirical_p_value))),
      as.character(nrow(taf_df)),
      as.character(sum(taf_df$direction == "increase", na.rm = TRUE)),
      as.character(sum(taf_df$direction == "decrease", na.rm = TRUE))
    ),
    stringsAsFactors = FALSE
  )
}

##### 6. Main #################################################################

main <- function() {
  args <- parse_args()

  module1_manifest <- normalize_path_safe(get_arg(args, "module1-manifest", required = TRUE))
  outdir <- normalize_path_safe(get_arg(args, "outdir", required = TRUE))
  prefix <- get_arg(args, "prefix", required = TRUE)
  phenotype1 <- get_arg(args, "phenotype1", required = TRUE)
  phenotype2 <- get_arg(args, "phenotype2", required = TRUE)

  delta_score_arg <- get_arg(args, "delta-score", default = NULL)
  module4_manifest_arg <- get_arg(args, "module4-manifest", default = NULL)
  degree_matrix_type <- get_arg(args, "degree-matrix-type", default = "topPct")
  degree_threshold <- as.numeric(get_arg(args, "degree-threshold", required = TRUE))
  score_col <- get_arg(args, "score-col", default = "delta_participation_score_lifespan_sum")
  degree_reference <- get_arg(args, "degree-reference", default = "avg")
  min_matched <- as.integer(get_arg(args, "min-matched", default = "20"))
  max_window <- as.numeric(get_arg(args, "max-window", default = "50"))
  alternative <- get_arg(args, "alternative", default = "two.sided")
  q_cutoff <- as.numeric(get_arg(args, "q-cutoff", default = "0.05"))
  p_cutoff_raw <- get_arg(args, "p-cutoff", default = "NA")
  p_cutoff <- suppressWarnings(as.numeric(p_cutoff_raw))
  min_abs_delta <- as.numeric(get_arg(args, "min-abs-delta", default = "0"))
  exclude_zero_degree <- as_logical_arg(get_arg(args, "exclude-zero-degree", default = "FALSE"), default = FALSE)
  pseudocount <- as.numeric(get_arg(args, "pseudocount", default = "1"))
  make_volcano <- as_logical_arg(
    get_arg(args, "make-volcano", default = "TRUE"),
    default = TRUE
  )
  volcano_format <- tolower(get_arg(args, "volcano-format", default = "both"))
  volcano_p_cutoff <- as.numeric(get_arg(args, "volcano-p-cutoff", default = "0.05"))
  volcano_delta_cutoff_arg <- get_arg(args, "volcano-delta-cutoff", default = NULL)
  volcano_delta_cutoff <- if (is.null(volcano_delta_cutoff_arg)) {
    min_abs_delta
  } else {
    as.numeric(volcano_delta_cutoff_arg)
  }
  volcano_width <- as.numeric(get_arg(args, "volcano-width", default = "8"))
  volcano_height <- as.numeric(get_arg(args, "volcano-height", default = "6"))
  volcano_dpi <- as.integer(get_arg(args, "volcano-dpi", default = "300"))
  volcano_point_size <- as.numeric(get_arg(args, "volcano-point-size", default = "0.65"))

  if (!degree_matrix_type %in% c("topPct", "distance")) stop_msg("--degree-matrix-type must be topPct or distance.")
  if (!degree_reference %in% c("avg", "phenotype1", "phenotype2")) stop_msg("--degree-reference must be avg, phenotype1, or phenotype2.")
  if (!alternative %in% c("two.sided", "greater", "less")) stop_msg("--alternative must be two.sided, greater, or less.")
  if (!is.finite(degree_threshold)) stop_msg("--degree-threshold must be numeric and finite.")
  if (!is.finite(min_matched) || min_matched < 1) stop_msg("--min-matched must be a positive integer.")
  if (!is.finite(max_window) || max_window < 0) stop_msg("--max-window must be >= 0.")
  if (!is.finite(q_cutoff) || q_cutoff < 0 || q_cutoff > 1) stop_msg("--q-cutoff must be between 0 and 1.")
  if (is.na(p_cutoff)) p_cutoff <- NA_real_
  if (!is.na(p_cutoff) && (p_cutoff < 0 || p_cutoff > 1)) stop_msg("--p-cutoff must be NA or between 0 and 1.")
  if (!is.finite(min_abs_delta) || min_abs_delta < 0) stop_msg("--min-abs-delta must be >= 0.")
  if (!is.finite(pseudocount) || pseudocount < 0) stop_msg("--pseudocount must be >= 0.")
  if (!volcano_format %in% c("png", "pdf", "both")) {
    stop_msg("--volcano-format must be png, pdf, or both.")
  }
  if (!is.finite(volcano_p_cutoff) || volcano_p_cutoff <= 0 || volcano_p_cutoff > 1) {
    stop_msg("--volcano-p-cutoff must be > 0 and <= 1.")
  }
  if (!is.finite(volcano_delta_cutoff) || volcano_delta_cutoff < 0) {
    stop_msg("--volcano-delta-cutoff must be >= 0.")
  }
  if (!is.finite(volcano_width) || volcano_width <= 0) stop_msg("--volcano-width must be positive.")
  if (!is.finite(volcano_height) || volcano_height <= 0) stop_msg("--volcano-height must be positive.")
  if (is.na(volcano_dpi) || volcano_dpi <= 0) stop_msg("--volcano-dpi must be a positive integer.")
  if (!is.finite(volcano_point_size) || volcano_point_size <= 0) {
    stop_msg("--volcano-point-size must be positive.")
  }
  if (identical(phenotype1, phenotype2)) stop_msg("phenotype1 and phenotype2 must be different.")

  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  delta_path <- resolve_delta_score_path(delta_score_arg, module4_manifest_arg)
  message_info("Reading Module 4 delta score table: ", delta_path)
  delta_df <- prepare_delta_table(read_tsv(delta_path), score_col)

  network_paths <- resolve_network_paths(module1_manifest, phenotype1, phenotype2, degree_matrix_type)
  degree_df <- make_degree_table(
    path1 = network_paths$phenotype1_matrix,
    path2 = network_paths$phenotype2_matrix,
    phenotype1 = phenotype1,
    phenotype2 = phenotype2,
    threshold = degree_threshold,
    matrix_type = degree_matrix_type
  )

  degree_col <- select_degree_reference(degree_df, phenotype1, phenotype2, degree_reference)

  merged_df <- merge(delta_df, degree_df, by = "feature_id", all.x = TRUE, sort = FALSE)
  if (any(is.na(merged_df[[degree_col]]))) {
    n_missing <- sum(is.na(merged_df[[degree_col]]))
    message_warn(n_missing, " features in Module 4 delta table have no degree information and will be excluded from testing.")
  }

  result_core <- compute_degree_matched_test(
    df = merged_df,
    score_col = score_col,
    degree_col = degree_col,
    min_matched = min_matched,
    max_window = max_window,
    alternative = alternative,
    pseudocount = pseudocount,
    exclude_zero_degree = exclude_zero_degree
  )

  # Add selected upstream columns back into result table.
  result_df <- merge(result_core, merged_df, by = "feature_id", all.x = TRUE, sort = FALSE)
  result_df <- result_df[order(result_df$q_value, -abs(result_df$score_tested), result_df$feature_id), , drop = FALSE]

  taf_df <- identify_tafs(
    result_df = result_df,
    q_cutoff = q_cutoff,
    p_cutoff = p_cutoff,
    min_abs_delta = min_abs_delta
  )

  volcano_df <- prepare_volcano_data(
    result_df = result_df,
    volcano_p_cutoff = volcano_p_cutoff,
    volcano_delta_cutoff = volcano_delta_cutoff
  )

  degree_path <- file.path(outdir, paste0(prefix, "_module5_network_degree_table.tsv"))
  merged_path <- file.path(outdir, paste0(prefix, "_module5_merged_score_degree_table.tsv"))
  result_path <- file.path(outdir, paste0(prefix, "_module5_degree_matched_null_result.tsv"))
  taf_path <- file.path(outdir, paste0(prefix, "_module5_topological_altering_features.tsv"))
  summary_path <- file.path(outdir, paste0(prefix, "_module5_summary.tsv"))
  volcano_data_path <- file.path(outdir, paste0(prefix, "_module5_volcano_data.tsv"))
  plot_dir <- file.path(outdir, "plots")
  manifest_path <- file.path(outdir, paste0(prefix, "_module5_output_manifest.tsv"))
  run_info_path <- file.path(outdir, paste0(prefix, "_module5_run_info.tsv"))

  summary_df <- make_summary(
    delta_df = delta_df,
    degree_df = degree_df,
    merged_df = merged_df,
    result_df = result_df,
    taf_df = taf_df,
    phenotype1 = phenotype1,
    phenotype2 = phenotype2,
    degree_matrix_type = degree_matrix_type,
    degree_threshold = degree_threshold,
    score_col = score_col,
    degree_col = degree_col,
    min_matched = min_matched,
    max_window = max_window,
    alternative = alternative,
    q_cutoff = q_cutoff,
    p_cutoff = p_cutoff,
    min_abs_delta = min_abs_delta
  )

  write_tsv(degree_df, degree_path)
  write_tsv(merged_df, merged_path)
  write_tsv(result_df, result_path)
  write_tsv(taf_df, taf_path)
  write_tsv(summary_df, summary_path)
  write_tsv(volcano_df, volcano_data_path)

  volcano_paths <- if (isTRUE(make_volcano)) {
    write_volcano_plots(
      volcano_df = volcano_df,
      plot_dir = plot_dir,
      prefix = prefix,
      volcano_p_cutoff = volcano_p_cutoff,
      volcano_delta_cutoff = volcano_delta_cutoff,
      volcano_format = volcano_format,
      volcano_width = volcano_width,
      volcano_height = volcano_height,
      volcano_dpi = volcano_dpi,
      volcano_point_size = volcano_point_size
    )
  } else {
    character()
  }

  manifest <- data.frame(
    item = c(
      "network_degree_table",
      "merged_score_degree_table",
      "degree_matched_null_result",
      "topological_altering_features",
      "summary",
      "volcano_data",
      names(volcano_paths)
    ),
    path = vapply(
      c(
        degree_path,
        merged_path,
        result_path,
        taf_path,
        summary_path,
        volcano_data_path,
        unname(volcano_paths)
      ),
      normalize_path_safe,
      character(1)
    ),
    stringsAsFactors = FALSE
  )
  write_tsv(manifest, manifest_path)

  run_info <- data.frame(
    parameter = c(
      "module",
      "module1_manifest",
      "module4_manifest",
      "delta_score",
      "outdir",
      "prefix",
      "phenotype1",
      "phenotype2",
      "degree_matrix_type",
      "degree_matrix_column",
      "phenotype1_degree_matrix",
      "phenotype2_degree_matrix",
      "degree_threshold",
      "score_col",
      "degree_reference",
      "degree_reference_col",
      "min_matched",
      "max_window",
      "alternative",
      "q_cutoff",
      "p_cutoff",
      "min_abs_delta",
      "exclude_zero_degree",
      "pseudocount",
      "make_volcano",
      "volcano_format",
      "volcano_p_cutoff",
      "volcano_delta_cutoff",
      "volcano_width",
      "volcano_height",
      "volcano_dpi",
      "volcano_point_size",
      "run_time"
    ),
    value = c(
      "module5_identify_topological_altering_features",
      module1_manifest,
      ifelse(is.null(module4_manifest_arg), "NA", normalize_path_safe(module4_manifest_arg)),
      delta_path,
      outdir,
      prefix,
      phenotype1,
      phenotype2,
      degree_matrix_type,
      network_paths$matrix_column,
      network_paths$phenotype1_matrix,
      network_paths$phenotype2_matrix,
      as.character(degree_threshold),
      score_col,
      degree_reference,
      degree_col,
      as.character(min_matched),
      as.character(max_window),
      alternative,
      as.character(q_cutoff),
      as.character(p_cutoff),
      as.character(min_abs_delta),
      as.character(exclude_zero_degree),
      as.character(pseudocount),
      as.character(make_volcano),
      volcano_format,
      as.character(volcano_p_cutoff),
      as.character(volcano_delta_cutoff),
      as.character(volcano_width),
      as.character(volcano_height),
      as.character(volcano_dpi),
      as.character(volcano_point_size),
      as.character(Sys.time())
    ),
    stringsAsFactors = FALSE
  )
  write_tsv(run_info, run_info_path)

  message_info("Module 5 completed.")
  message_info("Degree table: ", degree_path)
  message_info("Degree-matched result: ", result_path)
  message_info("Topological-altering features: ", taf_path)
  message_info("Volcano data: ", volcano_data_path)
  message_info("Output manifest: ", manifest_path)
}

main()
