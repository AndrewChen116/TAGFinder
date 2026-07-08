#!/usr/bin/env Rscript

################################################################################
# Module 1: Feature-feature network construction
#
# Purpose:
#   Construct phenotype-specific feature-feature networks from a feature table.
#   Rows are samples and columns are features, such as genes, methylation sites,
#   DMRs, peaks, proteins, metabolites, or other quantitative biological features.
#
# Main outputs:
#   1) phenotype-specific correlation matrices
#   2) phenotype-specific distance matrices
#   3) phenotype-specific top-percent filtration matrices with feature names
#   4) phenotype-specific numeric matrices for Ripser++ / Module 2
#   5) selected feature list and quality-control summary tables
#   6) Module 2 input list
#
# Design principles:
#   - no hard-coded project paths
#   - phenotype labels are user-defined, not cancer-specific
#   - a single input feature table and metadata table
#   - consistent feature set across phenotype1 and phenotype2
#   - explicit input validation and informative error messages
#   - modular functions to avoid copy-and-paste
################################################################################

options(stringsAsFactors = FALSE)

##### 1. Lightweight command-line parser ######################################

parse_args <- function(argv) {
  args <- list()
  i <- 1
  while (i <= length(argv)) {
    key <- argv[[i]]
    if (!startsWith(key, "--")) {
      stop("Unexpected argument: ", key, "\nArguments must use --key value format.", call. = FALSE)
    }
    key <- sub("^--", "", key)
    if (i == length(argv) || startsWith(argv[[i + 1]], "--")) {
      args[[key]] <- TRUE
      i <- i + 1
    } else {
      args[[key]] <- argv[[i + 1]]
      i <- i + 2
    }
  }
  args
}

get_arg <- function(args, name, default = NULL, required = FALSE) {
  value <- args[[name]]
  if (is.null(value)) {
    if (required) {
      stop("Missing required argument: --", name, call. = FALSE)
    }
    return(default)
  }
  value
}

as_logical_arg <- function(x) {
  if (is.logical(x)) return(x)
  x <- tolower(as.character(x))
  if (x %in% c("true", "t", "1", "yes", "y")) return(TRUE)
  if (x %in% c("false", "f", "0", "no", "n")) return(FALSE)
  stop("Cannot parse logical value: ", x, call. = FALSE)
}

print_usage <- function() {
  cat(
"\nModule 1: Feature-feature network construction\n\n",
"Required arguments:\n",
"  --feature-table <tsv/csv>     Rows = samples, columns = features.\n",
"  --metadata <tsv/csv>          Sample metadata containing sample ID and phenotype columns.\n",
"  --phenotype1 <label>          First phenotype label.\n",
"  --phenotype2 <label>          Second phenotype label.\n",
"  --outdir <directory>          Output directory.\n\n",
"Common optional arguments:\n",
"  --prefix <string>             Output prefix. Default: analysis\n",
"  --sample-col <string>         Sample ID column in metadata. Default: sample_id\n",
"  --phenotype-col <string>      Phenotype column in metadata. Default: phenotype\n",
"  --feature-id-col <string>     Feature-table ID column if sample IDs are stored in a column.\n",
"                                If omitted, the first column is used when non-numeric; otherwise row names are used.\n",
"  --cor-method <pearson|spearman>              Default: pearson\n",
"  --distance-transform <one_minus_abs_cor|one_minus_cor>  Default: one_minus_abs_cor\n",
"  --feature-filter-mode <common_valid|none>    Default: common_valid\n",
"  --min-samples-per-phenotype <integer>        Default: 3\n",
"  --min-nonmissing-fraction <numeric>          Default: 0.8\n",
"  --min-nonzero-fraction <numeric>             Default: 0\n",
"  --impute-missing <none|median>               Default: none\n",
"  --round-digits <integer>                     Default: 6\n",
"  --write-correlation <true|false>             Default: true\n\n",
"Example:\n",
"  Rscript module1_network_construction_v2.R \\\n",
"    --feature-table expression.tsv \\\n",
"    --metadata metadata.tsv \\\n",
"    --sample-col sample_id \\\n",
"    --phenotype-col phenotype \\\n",
"    --phenotype1 phenotype1 \\\n",
"    --phenotype2 phenotype2 \\\n",
"    --outdir module1_output \\\n",
"    --prefix demo\n\n",
sep = "")
}

##### 2. General utilities #####################################################

message_info <- function(...) message("[INFO] ", paste0(..., collapse = ""))
message_warn <- function(...) warning(paste0(..., collapse = ""), call. = FALSE)
stop_error <- function(...) stop(paste0(..., collapse = ""), call. = FALSE)

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(path)) stop_error("Failed to create output directory: ", path)
}

check_file_exists <- function(path, label) {
  if (is.null(path) || !file.exists(path)) {
    stop_error(label, " not found: ", path)
  }
}

infer_delim <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "csv") return(",")
  "\t"
}

read_table_auto <- function(path) {
  check_file_exists(path, "Input file")
  delim <- infer_delim(path)
  read.table(
    file = path,
    sep = delim,
    header = TRUE,
    quote = "",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

sanitize_label <- function(x) {
  x <- as.character(x)
  x <- gsub("[^A-Za-z0-9_.-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  ifelse(nchar(x) == 0, "phenotype", x)
}

write_tsv <- function(df, path) {
  write.table(df, file = path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

write_matrix_with_names <- function(mat, path, digits = 6) {
  out <- data.frame(feature_id = rownames(mat), mat, check.names = FALSE)
  write.table(out, file = path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

write_numeric_matrix <- function(mat, path, digits = 6) {
  write.table(
    round(mat, digits = digits),
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
}

##### 3. Input loading and validation #########################################

load_feature_table <- function(path, feature_id_col = NULL) {
  df <- read_table_auto(path)
  if (nrow(df) == 0) stop_error("Feature table has zero rows: ", path)
  if (ncol(df) < 2) stop_error("Feature table must contain sample IDs plus at least one feature column.")

  if (!is.null(feature_id_col) && nzchar(feature_id_col)) {
    if (!feature_id_col %in% colnames(df)) {
      stop_error("--feature-id-col not found in feature table: ", feature_id_col)
    }
    sample_ids <- df[[feature_id_col]]
    feature_df <- df[, setdiff(colnames(df), feature_id_col), drop = FALSE]
  } else {
    first_col <- df[[1]]
    first_col_numeric <- suppressWarnings(!any(is.na(as.numeric(first_col))))
    if (!first_col_numeric) {
      sample_ids <- first_col
      feature_df <- df[, -1, drop = FALSE]
    } else {
      if (is.null(rownames(df)) || any(rownames(df) == as.character(seq_len(nrow(df))))) {
        stop_error(
          "Cannot identify sample IDs in feature table. ",
          "Provide a sample ID column using --feature-id-col."
        )
      }
      sample_ids <- rownames(df)
      feature_df <- df
    }
  }

  sample_ids <- as.character(sample_ids)
  if (any(is.na(sample_ids) | sample_ids == "")) stop_error("Feature table contains missing or empty sample IDs.")
  if (anyDuplicated(sample_ids)) {
    dup <- unique(sample_ids[duplicated(sample_ids)])
    stop_error("Duplicated sample IDs in feature table: ", paste(head(dup, 10), collapse = ", "))
  }

  feature_names <- colnames(feature_df)
  if (any(is.na(feature_names) | feature_names == "")) stop_error("Feature table contains missing or empty feature names.")
  if (anyDuplicated(feature_names)) {
    dup <- unique(feature_names[duplicated(feature_names)])
    stop_error("Duplicated feature names in feature table: ", paste(head(dup, 10), collapse = ", "))
  }

  feature_mat <- as.matrix(feature_df)
  suppressWarnings(storage.mode(feature_mat) <- "numeric")
  if (any(is.na(feature_mat) & !is.na(as.matrix(feature_df)))) {
    stop_error("Feature table contains non-numeric values in feature columns.")
  }
  rownames(feature_mat) <- sample_ids
  colnames(feature_mat) <- feature_names

  feature_mat
}

load_metadata <- function(path, sample_col, phenotype_col) {
  metadata <- read_table_auto(path)
  if (!sample_col %in% colnames(metadata)) stop_error("Sample column not found in metadata: ", sample_col)
  if (!phenotype_col %in% colnames(metadata)) stop_error("Phenotype column not found in metadata: ", phenotype_col)

  metadata[[sample_col]] <- as.character(metadata[[sample_col]])
  metadata[[phenotype_col]] <- as.character(metadata[[phenotype_col]])

  if (any(is.na(metadata[[sample_col]]) | metadata[[sample_col]] == "")) {
    stop_error("Metadata contains missing or empty sample IDs.")
  }
  if (anyDuplicated(metadata[[sample_col]])) {
    dup <- unique(metadata[[sample_col]][duplicated(metadata[[sample_col]])])
    stop_error("Duplicated sample IDs in metadata: ", paste(head(dup, 10), collapse = ", "))
  }
  metadata
}

match_samples <- function(feature_mat, metadata, sample_col, phenotype_col, phenotype1, phenotype2) {
  keep_meta <- metadata[metadata[[phenotype_col]] %in% c(phenotype1, phenotype2), , drop = FALSE]
  if (nrow(keep_meta) == 0) {
    stop_error("No metadata samples found for phenotype1/phenotype2.")
  }

  common_samples <- intersect(rownames(feature_mat), keep_meta[[sample_col]])
  if (length(common_samples) == 0) {
    stop_error("No overlapping sample IDs between feature table and metadata.")
  }

  missing_in_feature <- setdiff(keep_meta[[sample_col]], rownames(feature_mat))
  if (length(missing_in_feature) > 0) {
    message_warn("Metadata samples not found in feature table: ", paste(head(missing_in_feature, 10), collapse = ", "))
  }

  keep_meta <- keep_meta[match(common_samples, keep_meta[[sample_col]]), , drop = FALSE]
  feature_mat <- feature_mat[common_samples, , drop = FALSE]

  list(feature_mat = feature_mat, metadata = keep_meta)
}

##### 4. Feature filtering #####################################################

feature_qc_by_group <- function(mat, group_label) {
  nonmissing_fraction <- colMeans(!is.na(mat))
  nonzero_fraction <- colMeans(!is.na(mat) & mat != 0)
  variance <- apply(mat, 2, function(x) stats::var(x, na.rm = TRUE))
  data.frame(
    phenotype = group_label,
    feature_id = colnames(mat),
    n_samples = nrow(mat),
    nonmissing_fraction = nonmissing_fraction,
    nonzero_fraction = nonzero_fraction,
    variance = variance,
    zero_variance = is.na(variance) | variance == 0,
    stringsAsFactors = FALSE
  )
}

select_common_valid_features <- function(feature_mat, metadata, sample_col, phenotype_col,
                                         phenotype1, phenotype2,
                                         min_nonmissing_fraction,
                                         min_nonzero_fraction,
                                         filter_mode) {
  idx1 <- metadata[[phenotype_col]] == phenotype1
  idx2 <- metadata[[phenotype_col]] == phenotype2

  mat1 <- feature_mat[idx1, , drop = FALSE]
  mat2 <- feature_mat[idx2, , drop = FALSE]

  qc1 <- feature_qc_by_group(mat1, phenotype1)
  qc2 <- feature_qc_by_group(mat2, phenotype2)
  qc <- rbind(qc1, qc2)

  if (filter_mode == "none") {
    selected <- colnames(feature_mat)
  } else if (filter_mode == "common_valid") {
    valid1 <- qc1$feature_id[
      qc1$nonmissing_fraction >= min_nonmissing_fraction &
        qc1$nonzero_fraction >= min_nonzero_fraction &
        !qc1$zero_variance
    ]
    valid2 <- qc2$feature_id[
      qc2$nonmissing_fraction >= min_nonmissing_fraction &
        qc2$nonzero_fraction >= min_nonzero_fraction &
        !qc2$zero_variance
    ]
    selected <- intersect(valid1, valid2)
  } else {
    stop_error("Unsupported --feature-filter-mode: ", filter_mode)
  }

  if (length(selected) < 3) {
    stop_error(
      "Only ", length(selected), " features passed filtering. ",
      "At least 3 features are required for persistent homology. ",
      "Consider relaxing --min-nonmissing-fraction or checking zero-variance features."
    )
  }

  list(selected_features = selected, qc = qc)
}

impute_matrix <- function(mat, method) {
  if (method == "none") return(mat)
  if (method != "median") stop_error("Unsupported --impute-missing: ", method)

  for (j in seq_len(ncol(mat))) {
    miss <- is.na(mat[, j])
    if (any(miss)) {
      med <- stats::median(mat[, j], na.rm = TRUE)
      if (is.na(med)) stop_error("Cannot median-impute feature with all missing values: ", colnames(mat)[j])
      mat[miss, j] <- med
    }
  }
  mat
}

##### 5. Network construction #################################################

compute_correlation_matrix <- function(mat, method) {
  if (!method %in% c("pearson", "spearman")) {
    stop_error("Unsupported --cor-method: ", method, ". Use pearson or spearman.")
  }
  cor_mat <- suppressWarnings(stats::cor(mat, method = method, use = "pairwise.complete.obs"))
  if (any(!is.finite(cor_mat))) {
    bad_count <- sum(!is.finite(cor_mat))
    stop_error(
      "Correlation matrix contains non-finite values (n = ", bad_count, "). ",
      "Check missing values, zero variance features, or use --impute-missing median."
    )
  }
  diag(cor_mat) <- 1
  cor_mat
}

correlation_to_distance <- function(cor_mat, transform) {
  if (transform == "one_minus_abs_cor") {
    dist_mat <- 1 - abs(cor_mat)
  } else if (transform == "one_minus_cor") {
    dist_mat <- (1 - cor_mat) / 2
  } else {
    stop_error("Unsupported --distance-transform: ", transform)
  }
  dist_mat[dist_mat < 0 & dist_mat > -1e-12] <- 0
  dist_mat[dist_mat > 1 & dist_mat < 1 + 1e-12] <- 1
  diag(dist_mat) <- 0
  dist_mat
}

distance_to_top_percent <- function(dist_mat, digits = 6) {
  n <- nrow(dist_mat)
  if (n != ncol(dist_mat)) stop_error("Distance matrix must be square.")
  out <- matrix(0, nrow = n, ncol = n, dimnames = dimnames(dist_mat))

  upper_idx <- which(upper.tri(dist_mat), arr.ind = TRUE)
  vals <- dist_mat[upper_idx]
  if (any(!is.finite(vals))) stop_error("Distance matrix contains non-finite off-diagonal values.")

  # Smaller distance = stronger association = appears earlier in filtration.
  # The strongest edge receives a value near 0; the weakest receives 100.
  ord <- order(vals, method = "radix")
  ranks <- integer(length(vals))
  ranks[ord] <- seq_along(vals)

  if (length(vals) == 1) {
    pct <- 0
  } else {
    pct <- 100 * (ranks - 1) / (length(vals) - 1)
  }
  pct <- round(pct, digits = digits)

  out[upper_idx] <- pct
  out[cbind(upper_idx[, 2], upper_idx[, 1])] <- pct
  diag(out) <- 0
  out
}

summarize_network <- function(cor_mat, dist_mat, phenotype, edge_threshold_distance = NULL) {
  upper <- upper.tri(dist_mat)
  distances <- dist_mat[upper]
  correlations <- cor_mat[upper]

  summary <- data.frame(
    phenotype = phenotype,
    n_features = nrow(dist_mat),
    n_edges_complete_graph = length(distances),
    distance_min = min(distances),
    distance_q1 = unname(stats::quantile(distances, 0.25)),
    distance_median = stats::median(distances),
    distance_mean = mean(distances),
    distance_q3 = unname(stats::quantile(distances, 0.75)),
    distance_max = max(distances),
    abs_correlation_mean = mean(abs(correlations)),
    abs_correlation_median = stats::median(abs(correlations)),
    stringsAsFactors = FALSE
  )

  if (!is.null(edge_threshold_distance) && is.finite(edge_threshold_distance)) {
    edge_count <- sum(distances > 0 & distances <= edge_threshold_distance)
    summary$edge_threshold_distance <- edge_threshold_distance
    summary$n_edges_under_threshold <- edge_count
    summary$edge_density_under_threshold <- edge_count / length(distances)
  }

  summary
}

process_one_phenotype <- function(feature_mat, metadata, sample_col, phenotype_col, phenotype_label,
                                  selected_features, cor_method, distance_transform,
                                  impute_missing, round_digits, write_correlation,
                                  outdir, prefix) {
  safe_label <- sanitize_label(phenotype_label)
  idx <- metadata[[phenotype_col]] == phenotype_label
  mat <- feature_mat[idx, selected_features, drop = FALSE]
  mat <- impute_matrix(mat, impute_missing)

  message_info("Computing correlation matrix for phenotype: ", phenotype_label)
  cor_mat <- compute_correlation_matrix(mat, cor_method)
  dist_mat <- correlation_to_distance(cor_mat, distance_transform)
  top_pct <- distance_to_top_percent(dist_mat, digits = round_digits)

  cor_path <- file.path(outdir, paste0(prefix, "_", safe_label, "_correlation_matrix.tsv"))
  dist_path <- file.path(outdir, paste0(prefix, "_", safe_label, "_distance_matrix.tsv"))
  top_path <- file.path(outdir, paste0(prefix, "_", safe_label, "_topPct_matrix.tsv"))
  numeric_path <- file.path(outdir, paste0(prefix, "_", safe_label, "_topPct_numeric.tsv"))

  if (write_correlation) write_matrix_with_names(round(cor_mat, round_digits), cor_path, digits = round_digits)
  write_matrix_with_names(round(dist_mat, round_digits), dist_path, digits = round_digits)
  write_matrix_with_names(top_pct, top_path, digits = round_digits)
  write_numeric_matrix(top_pct, numeric_path, digits = round_digits)

  list(
    phenotype = phenotype_label,
    safe_label = safe_label,
    n_samples = nrow(mat),
    cor_path = if (write_correlation) cor_path else NA_character_,
    dist_path = dist_path,
    top_path = top_path,
    numeric_path = numeric_path,
    summary = summarize_network(cor_mat, dist_mat, phenotype_label)
  )
}

##### 6. Main workflow #########################################################

main <- function() {
  argv <- commandArgs(trailingOnly = TRUE)
  if (length(argv) == 0 || any(argv %in% c("--help", "-h"))) {
    print_usage()
    quit(status = 0)
  }

  args <- parse_args(argv)

  feature_table <- get_arg(args, "feature-table", required = TRUE)
  metadata_file <- get_arg(args, "metadata", required = TRUE)
  phenotype1 <- get_arg(args, "phenotype1", required = TRUE)
  phenotype2 <- get_arg(args, "phenotype2", required = TRUE)
  outdir <- get_arg(args, "outdir", required = TRUE)

  prefix <- get_arg(args, "prefix", "analysis")
  sample_col <- get_arg(args, "sample-col", "sample_id")
  phenotype_col <- get_arg(args, "phenotype-col", "phenotype")
  feature_id_col <- get_arg(args, "feature-id-col", NULL)
  cor_method <- get_arg(args, "cor-method", "pearson")
  distance_transform <- get_arg(args, "distance-transform", "one_minus_abs_cor")
  feature_filter_mode <- get_arg(args, "feature-filter-mode", "common_valid")
  min_samples_per_phenotype <- as.integer(get_arg(args, "min-samples-per-phenotype", "3"))
  min_nonmissing_fraction <- as.numeric(get_arg(args, "min-nonmissing-fraction", "0.8"))
  min_nonzero_fraction <- as.numeric(get_arg(args, "min-nonzero-fraction", "0"))
  impute_missing <- get_arg(args, "impute-missing", "none")
  round_digits <- as.integer(get_arg(args, "round-digits", "6"))
  write_correlation <- as_logical_arg(get_arg(args, "write-correlation", "true"))

  if (phenotype1 == phenotype2) stop_error("--phenotype1 and --phenotype2 must be different.")
  if (!feature_filter_mode %in% c("common_valid", "none")) stop_error("--feature-filter-mode must be common_valid or none.")
  if (!impute_missing %in% c("none", "median")) stop_error("--impute-missing must be none or median.")
  if (min_nonmissing_fraction < 0 || min_nonmissing_fraction > 1) stop_error("--min-nonmissing-fraction must be between 0 and 1.")
  if (min_nonzero_fraction < 0 || min_nonzero_fraction > 1) stop_error("--min-nonzero-fraction must be between 0 and 1.")

  ensure_dir(outdir)

  message_info("Loading feature table: ", feature_table)
  feature_mat <- load_feature_table(feature_table, feature_id_col)

  message_info("Loading metadata: ", metadata_file)
  metadata <- load_metadata(metadata_file, sample_col, phenotype_col)

  matched <- match_samples(feature_mat, metadata, sample_col, phenotype_col, phenotype1, phenotype2)
  feature_mat <- matched$feature_mat
  metadata <- matched$metadata

  n1 <- sum(metadata[[phenotype_col]] == phenotype1)
  n2 <- sum(metadata[[phenotype_col]] == phenotype2)
  if (n1 < min_samples_per_phenotype) {
    stop_error("Phenotype ", phenotype1, " has only ", n1, " samples. Minimum required: ", min_samples_per_phenotype)
  }
  if (n2 < min_samples_per_phenotype) {
    stop_error("Phenotype ", phenotype2, " has only ", n2, " samples. Minimum required: ", min_samples_per_phenotype)
  }

  message_info("Matched samples: ", nrow(feature_mat), " | features before filtering: ", ncol(feature_mat))
  message_info("Phenotype sample counts: ", phenotype1, "=", n1, ", ", phenotype2, "=", n2)

  feature_selection <- select_common_valid_features(
    feature_mat = feature_mat,
    metadata = metadata,
    sample_col = sample_col,
    phenotype_col = phenotype_col,
    phenotype1 = phenotype1,
    phenotype2 = phenotype2,
    min_nonmissing_fraction = min_nonmissing_fraction,
    min_nonzero_fraction = min_nonzero_fraction,
    filter_mode = feature_filter_mode
  )

  selected_features <- feature_selection$selected_features
  qc <- feature_selection$qc

  message_info("Selected common features: ", length(selected_features))

  selected_feature_df <- data.frame(
    feature_id = selected_features,
    feature_index = seq_along(selected_features),
    stringsAsFactors = FALSE
  )
  selected_feature_path <- file.path(outdir, paste0(prefix, "_module1_selected_features.tsv"))
  qc_path <- file.path(outdir, paste0(prefix, "_module1_feature_qc.tsv"))
  matched_sample_path <- file.path(outdir, paste0(prefix, "_module1_matched_samples.tsv"))

  write_tsv(selected_feature_df, selected_feature_path)
  write_tsv(qc, qc_path)
  write_tsv(metadata[, c(sample_col, phenotype_col), drop = FALSE], matched_sample_path)

  res1 <- process_one_phenotype(
    feature_mat = feature_mat,
    metadata = metadata,
    sample_col = sample_col,
    phenotype_col = phenotype_col,
    phenotype_label = phenotype1,
    selected_features = selected_features,
    cor_method = cor_method,
    distance_transform = distance_transform,
    impute_missing = impute_missing,
    round_digits = round_digits,
    write_correlation = write_correlation,
    outdir = outdir,
    prefix = prefix
  )

  res2 <- process_one_phenotype(
    feature_mat = feature_mat,
    metadata = metadata,
    sample_col = sample_col,
    phenotype_col = phenotype_col,
    phenotype_label = phenotype2,
    selected_features = selected_features,
    cor_method = cor_method,
    distance_transform = distance_transform,
    impute_missing = impute_missing,
    round_digits = round_digits,
    write_correlation = write_correlation,
    outdir = outdir,
    prefix = prefix
  )

  network_summary <- rbind(res1$summary, res2$summary)
  network_summary$n_samples <- c(res1$n_samples, res2$n_samples)
  network_summary$correlation_method <- cor_method
  network_summary$distance_transform <- distance_transform
  summary_path <- file.path(outdir, paste0(prefix, "_module1_network_summary.tsv"))
  write_tsv(network_summary, summary_path)

  output_manifest <- data.frame(
    module = "module1_network_construction",
    phenotype = c(res1$phenotype, res2$phenotype),
    safe_phenotype_label = c(res1$safe_label, res2$safe_label),
    n_samples = c(res1$n_samples, res2$n_samples),
    n_features = length(selected_features),
    correlation_matrix = c(res1$cor_path, res2$cor_path),
    distance_matrix = c(res1$dist_path, res2$dist_path),
    topPct_matrix = c(res1$top_path, res2$top_path),
    topPct_numeric_matrix = c(res1$numeric_path, res2$numeric_path),
    stringsAsFactors = FALSE
  )
  manifest_path <- file.path(outdir, paste0(prefix, "_module1_output_manifest.tsv"))
  write_tsv(output_manifest, manifest_path)

  module2_input_list <- file.path(outdir, paste0(prefix, "_module2_ripser_input_list.txt"))
  writeLines(c(res1$numeric_path, res2$numeric_path), con = module2_input_list)

  run_info <- data.frame(
    parameter = c(
      "feature_table", "metadata", "sample_col", "phenotype_col",
      "phenotype1", "phenotype2", "prefix", "cor_method",
      "distance_transform", "feature_filter_mode", "min_samples_per_phenotype",
      "min_nonmissing_fraction", "min_nonzero_fraction", "impute_missing",
      "round_digits", "write_correlation", "selected_features", "module2_input_list"
    ),
    value = c(
      feature_table, metadata_file, sample_col, phenotype_col,
      phenotype1, phenotype2, prefix, cor_method,
      distance_transform, feature_filter_mode, min_samples_per_phenotype,
      min_nonmissing_fraction, min_nonzero_fraction, impute_missing,
      round_digits, write_correlation, length(selected_features), module2_input_list
    ),
    stringsAsFactors = FALSE
  )
  run_info_path <- file.path(outdir, paste0(prefix, "_module1_run_info.tsv"))
  write_tsv(run_info, run_info_path)

  message_info("Module 1 completed successfully.")
  message_info("Selected features: ", selected_feature_path)
  message_info("Output manifest: ", manifest_path)
  message_info("Module 2 input list: ", module2_input_list)
}

main()
