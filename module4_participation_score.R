#!/usr/bin/env Rscript

################################################################################
# Module 4: Participation score calculation
#
# Purpose
#   Convert topological feature composition tables from Module 3 into
#   feature-level participation scores for each phenotype, then calculate
#   phenotype2 - phenotype1 delta scores.
#
# Expected upstream input
#   Module 3 output table with at least the following columns:
#     phenotype
#     topological_feature_id or feature_id
#     birth
#     death
#     lifespan OR finite birth/death values so lifespan can be calculated
#     total_composition OR birth_composition/death_composition
#
# Main outputs
#   <prefix>_module4_participation_score_long.tsv
#   <prefix>_module4_participation_score_wide.tsv
#   <prefix>_module4_delta_participation_score.tsv
#   <prefix>_module5_input_manifest.tsv
#   <prefix>_module4_run_info.tsv
#
# Design principles
#   - No hard-coded project path
#   - Command-line configurable input/output
#   - Modular functions
#   - Automatic input validation and clear error messages
#   - Generic phenotype1 / phenotype2 naming, not cancer/control-specific
################################################################################

options(stringsAsFactors = FALSE)

##### 1. Utility functions ######################################################

message_info <- function(...) {
  message("[INFO] ", paste0(..., collapse = ""))
}

message_warn <- function(...) {
  warning(paste0(..., collapse = ""), call. = FALSE)
}

stop_msg <- function(...) {
  stop(paste0(..., collapse = ""), call. = FALSE)
}

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) == 0 || any(args %in% c("-h", "--help"))) {
    cat(
"Module 4: Participation score calculation\n\n",
"Usage:\n",
"  Rscript module4_participation_score.R \\\n",
"    --composition module3_output/demo_module3_composition_all.tsv \\\n",
"    --outdir module4_output \\\n",
"    --prefix demo \\\n",
"    --phenotype1 phenotype1 \\\n",
"    --phenotype2 phenotype2\n\n",
"Alternative using Module 3 manifest:\n",
"  Rscript module4_participation_score.R \\\n",
"    --module3-manifest module3_output/demo_module4_input_manifest.tsv \\\n",
"    --outdir module4_output \\\n",
"    --prefix demo\n\n",
"Required arguments:\n",
"  --outdir                  Output directory\n",
"  --prefix                  Output prefix\n\n",
"Input arguments, choose one:\n",
"  --composition             Module 3 composition_all TSV\n",
"  --module3-manifest        Module 3 manifest containing composition_all path\n\n",
"Optional arguments:\n",
"  --phenotype1              Baseline phenotype. If omitted, inferred from table order\n",
"  --phenotype2              Comparison phenotype. If omitted, inferred from table order\n",
"  --composition-col         Column used for score attribution: total_composition, birth_composition, or death_composition [default: total_composition]\n",
"  --feature-list            Optional feature list TSV with a feature/gene column; zero scores will be added for absent features\n",
"  --feature-col             Column name in feature-list [default: feature_id]\n",
"  --min-lifespan            Minimum lifespan retained for scoring [default: 0]\n",
"  --include-infinite-death  Keep intervals with Inf death by setting lifespan to --infinite-lifespan-value [default: FALSE]\n",
"  --infinite-lifespan-value Numeric lifespan value assigned to Inf death if included [default: 0]\n",
"  --separator               Separator used in composition columns [default: ;]\n",
"  --write-zero-scores       TRUE/FALSE. Add zero rows for all listed features and phenotypes [default: TRUE]\n",
"  --top-n                   Number of top absolute delta rows for summary [default: 100]\n\n",
sep = "")
    quit(status = 0)
  }

  out <- list()
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!grepl("^--", key)) {
      stop_msg("Unexpected argument: ", key)
    }
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

read_tsv_base <- function(path) {
  if (!file.exists(path)) stop_msg("File not found: ", path)
  tryCatch(
    read.delim(path, sep = "\t", header = TRUE, check.names = FALSE, quote = "", comment.char = ""),
    error = function(e) stop_msg("Failed to read TSV: ", path, "\n", e$message)
  )
}

write_tsv_base <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(df, file = path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

normalize_path <- function(path) {
  if (is.null(path) || is.na(path) || path == "") return(path)
  normalizePath(path, mustWork = FALSE)
}

require_columns <- function(df, cols, label = "table") {
  missing <- setdiff(cols, colnames(df))
  if (length(missing) > 0) {
    stop_msg("Missing required column(s) in ", label, ": ", paste(missing, collapse = ", "))
  }
}

first_existing_col <- function(df, candidates, label = "column") {
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) == 0) {
    stop_msg("Cannot find ", label, ". Tried: ", paste(candidates, collapse = ", "))
  }
  hit[[1]]
}

##### 2. Input resolution ######################################################

resolve_composition_path <- function(composition, module3_manifest) {
  if (!is.null(composition)) return(normalize_path(composition))
  if (is.null(module3_manifest)) {
    stop_msg("Please provide either --composition or --module3-manifest")
  }

  manifest <- read_tsv_base(module3_manifest)
  if (all(c("item", "path") %in% colnames(manifest))) {
    idx <- which(manifest$item %in% c("composition_all", "module3_composition_all"))
    if (length(idx) == 0) {
      idx <- grep("composition.*all", manifest$item, ignore.case = TRUE)
    }
    if (length(idx) > 0) return(normalize_path(manifest$path[[idx[[1]]]]))
  }

  path_cols <- colnames(manifest)[grepl("path|file", colnames(manifest), ignore.case = TRUE)]
  for (pc in path_cols) {
    idx <- grep("composition.*all|module3.*composition", manifest[[pc]], ignore.case = TRUE)
    if (length(idx) > 0) return(normalize_path(manifest[[pc]][[idx[[1]]]]))
  }

  stop_msg("Could not identify composition_all path from module3 manifest: ", module3_manifest)
}

##### 3. Data preparation ######################################################

standardize_composition_table <- function(df, composition_col, include_inf, inf_lifespan_value) {
  require_columns(df, c("phenotype"), "Module 3 composition table")

  topo_id_col <- first_existing_col(
    df,
    c("topological_feature_id", "feature_id", "barcode_id", "ID"),
    label = "topological feature ID column"
  )

  if (!composition_col %in% colnames(df)) {
    stop_msg(
      "Composition column not found: ", composition_col,
      ". Available columns: ", paste(colnames(df), collapse = ", ")
    )
  }

  if (!"lifespan" %in% colnames(df)) {
    require_columns(df, c("birth", "death"), "Module 3 composition table without lifespan")
    df$lifespan <- suppressWarnings(as.numeric(df$death) - as.numeric(df$birth))
  }

  df$birth <- if ("birth" %in% colnames(df)) suppressWarnings(as.numeric(df$birth)) else NA_real_
  df$death <- if ("death" %in% colnames(df)) suppressWarnings(as.numeric(df$death)) else NA_real_
  df$lifespan <- suppressWarnings(as.numeric(df$lifespan))

  if (include_inf) {
    infinite_death <- is.infinite(df$death) | is.infinite(df$lifespan)
    df$lifespan[infinite_death] <- inf_lifespan_value
  }

  df$.topological_feature_id <- as.character(df[[topo_id_col]])
  df$.composition <- as.character(df[[composition_col]])
  df$.composition[is.na(df$.composition)] <- ""

  df
}

split_composition_to_long <- function(df, separator = ";") {
  rows <- vector("list", nrow(df))
  sep_regex <- separator
  if (separator %in% c("|", ".", "+", "*", "?", "^", "$", "(", ")", "[", "]", "{", "}", "\\")) {
    sep_regex <- paste0("\\", separator)
  }

  for (i in seq_len(nrow(df))) {
    comp <- df$.composition[[i]]
    if (is.na(comp) || trimws(comp) == "") next
    features <- unlist(strsplit(comp, sep_regex, fixed = FALSE), use.names = FALSE)
    features <- trimws(features)
    features <- unique(features[features != "" & !is.na(features) & toupper(features) != "NA"])
    if (length(features) == 0) next

    rows[[i]] <- data.frame(
      phenotype = as.character(df$phenotype[[i]]),
      attributed_feature_id = features,
      topological_feature_id = df$.topological_feature_id[[i]],
      birth = df$birth[[i]],
      death = df$death[[i]],
      lifespan = df$lifespan[[i]],
      stringsAsFactors = FALSE
    )
  }

  out <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
  if (is.null(out)) {
    out <- data.frame(
      phenotype = character(),
      attributed_feature_id = character(),
      topological_feature_id = character(),
      birth = numeric(),
      death = numeric(),
      lifespan = numeric(),
      stringsAsFactors = FALSE
    )
  }
  out
}

calculate_scores <- function(comp_long) {
  if (nrow(comp_long) == 0) {
    return(data.frame(
      phenotype = character(),
      feature_id = character(),
      participation_score_lifespan_sum = numeric(),
      participation_score_feature_count = integer(),
      mean_lifespan = numeric(),
      median_lifespan = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  key <- paste(comp_long$phenotype, comp_long$attributed_feature_id, sep = "\r")
  split_idx <- split(seq_len(nrow(comp_long)), key)

  rows <- lapply(split_idx, function(idx) {
    d <- comp_long[idx, , drop = FALSE]
    data.frame(
      phenotype = d$phenotype[[1]],
      feature_id = d$attributed_feature_id[[1]],
      participation_score_lifespan_sum = sum(d$lifespan, na.rm = TRUE),
      participation_score_feature_count = length(unique(d$topological_feature_id)),
      mean_lifespan = mean(d$lifespan, na.rm = TRUE),
      median_lifespan = median(d$lifespan, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out[order(out$phenotype, out$feature_id), , drop = FALSE]
}

read_feature_list <- function(path, feature_col) {
  if (is.null(path)) return(NULL)
  df <- read_tsv_base(path)
  if (!feature_col %in% colnames(df)) {
    # tolerate common alternatives
    alt <- first_existing_col(df, c(feature_col, "feature_id", "feature", "gene", "Gene_ID", "gene_id", "id"), "feature-list column")
    feature_col <- alt
  }
  unique(as.character(df[[feature_col]]))
}

complete_zero_scores <- function(score_df, all_features, phenotypes) {
  if (is.null(all_features)) return(score_df)
  grid <- expand.grid(
    phenotype = phenotypes,
    feature_id = all_features,
    stringsAsFactors = FALSE
  )
  key_grid <- paste(grid$phenotype, grid$feature_id, sep = "\r")
  key_score <- paste(score_df$phenotype, score_df$feature_id, sep = "\r")
  missing <- !(key_grid %in% key_score)

  if (any(missing)) {
    zero_df <- data.frame(
      phenotype = grid$phenotype[missing],
      feature_id = grid$feature_id[missing],
      participation_score_lifespan_sum = 0,
      participation_score_feature_count = 0L,
      mean_lifespan = NA_real_,
      median_lifespan = NA_real_,
      stringsAsFactors = FALSE
    )
    score_df <- rbind(score_df, zero_df)
  }

  score_df[order(score_df$phenotype, score_df$feature_id), , drop = FALSE]
}

make_score_wide <- function(score_df, phenotypes) {
  metrics <- c(
    "participation_score_lifespan_sum",
    "participation_score_feature_count",
    "mean_lifespan",
    "median_lifespan"
  )

  features <- sort(unique(score_df$feature_id))
  wide <- data.frame(feature_id = features, stringsAsFactors = FALSE)

  for (ph in phenotypes) {
    sub <- score_df[score_df$phenotype == ph, c("feature_id", metrics), drop = FALSE]
    colnames(sub)[match(metrics, colnames(sub))] <- paste(metrics, ph, sep = "__")
    wide <- merge(wide, sub, by = "feature_id", all.x = TRUE, sort = FALSE)
  }

  wide[order(wide$feature_id), , drop = FALSE]
}

calculate_delta_scores <- function(score_df, phenotype1, phenotype2) {
  s1 <- score_df[score_df$phenotype == phenotype1, , drop = FALSE]
  s2 <- score_df[score_df$phenotype == phenotype2, , drop = FALSE]

  keep_cols <- c(
    "feature_id",
    "participation_score_lifespan_sum",
    "participation_score_feature_count",
    "mean_lifespan",
    "median_lifespan"
  )
  s1 <- s1[, keep_cols, drop = FALSE]
  s2 <- s2[, keep_cols, drop = FALSE]

  names(s1)[-1] <- paste0(names(s1)[-1], "_", phenotype1)
  names(s2)[-1] <- paste0(names(s2)[-1], "_", phenotype2)

  merged <- merge(s1, s2, by = "feature_id", all = TRUE, sort = FALSE)

  p1_life <- paste0("participation_score_lifespan_sum_", phenotype1)
  p2_life <- paste0("participation_score_lifespan_sum_", phenotype2)
  p1_count <- paste0("participation_score_feature_count_", phenotype1)
  p2_count <- paste0("participation_score_feature_count_", phenotype2)

  for (cn in c(p1_life, p2_life, p1_count, p2_count)) {
    merged[[cn]][is.na(merged[[cn]])] <- 0
  }

  merged$delta_participation_score_lifespan_sum <- merged[[p2_life]] - merged[[p1_life]]
  merged$delta_participation_score_feature_count <- merged[[p2_count]] - merged[[p1_count]]

  merged$direction_lifespan_sum <- ifelse(
    merged$delta_participation_score_lifespan_sum > 0, "increase",
    ifelse(merged$delta_participation_score_lifespan_sum < 0, "decrease", "no_change")
  )
  merged$direction_feature_count <- ifelse(
    merged$delta_participation_score_feature_count > 0, "increase",
    ifelse(merged$delta_participation_score_feature_count < 0, "decrease", "no_change")
  )

  merged$rank_abs_delta_lifespan_sum <- rank(-abs(merged$delta_participation_score_lifespan_sum), ties.method = "first")
  merged$rank_abs_delta_feature_count <- rank(-abs(merged$delta_participation_score_feature_count), ties.method = "first")

  merged[order(merged$rank_abs_delta_lifespan_sum, merged$rank_abs_delta_feature_count), , drop = FALSE]
}

make_summary <- function(comp_df, comp_long, score_df, delta_df, phenotype1, phenotype2, top_n) {
  data.frame(
    item = c(
      "n_topological_features_input",
      "n_feature_composition_rows",
      "n_scored_features_total",
      paste0("n_scored_features_", phenotype1),
      paste0("n_scored_features_", phenotype2),
      "n_delta_features",
      "n_increase_lifespan_sum",
      "n_decrease_lifespan_sum",
      "n_no_change_lifespan_sum",
      "top_n_abs_delta_reported"
    ),
    value = c(
      nrow(comp_df),
      nrow(comp_long),
      length(unique(score_df$feature_id)),
      sum(score_df$phenotype == phenotype1),
      sum(score_df$phenotype == phenotype2),
      nrow(delta_df),
      sum(delta_df$direction_lifespan_sum == "increase", na.rm = TRUE),
      sum(delta_df$direction_lifespan_sum == "decrease", na.rm = TRUE),
      sum(delta_df$direction_lifespan_sum == "no_change", na.rm = TRUE),
      min(top_n, nrow(delta_df))
    ),
    stringsAsFactors = FALSE
  )
}

##### 4. Main ##################################################################

main <- function() {
  args <- parse_args()

  outdir <- normalize_path(get_arg(args, "outdir", required = TRUE))
  prefix <- get_arg(args, "prefix", required = TRUE)
  composition_arg <- get_arg(args, "composition", default = NULL)
  manifest_arg <- get_arg(args, "module3-manifest", default = NULL)

  composition_col <- get_arg(args, "composition-col", default = "total_composition")
  feature_list_path <- get_arg(args, "feature-list", default = NULL)
  feature_col <- get_arg(args, "feature-col", default = "feature_id")
  min_lifespan <- as.numeric(get_arg(args, "min-lifespan", default = "0"))
  include_inf <- as_logical_arg(get_arg(args, "include-infinite-death", default = "FALSE"), default = FALSE)
  inf_lifespan_value <- as.numeric(get_arg(args, "infinite-lifespan-value", default = "0"))
  separator <- get_arg(args, "separator", default = ";")
  write_zero_scores <- as_logical_arg(get_arg(args, "write-zero-scores", default = "TRUE"), default = TRUE)
  top_n <- as.integer(get_arg(args, "top-n", default = "100"))

  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  composition_path <- resolve_composition_path(composition_arg, manifest_arg)
  message_info("Reading Module 3 composition table: ", composition_path)
  comp_raw <- read_tsv_base(composition_path)

  comp <- standardize_composition_table(
    df = comp_raw,
    composition_col = composition_col,
    include_inf = include_inf,
    inf_lifespan_value = inf_lifespan_value
  )

  n_before <- nrow(comp)
  comp <- comp[!is.na(comp$lifespan) & is.finite(comp$lifespan) & comp$lifespan >= min_lifespan, , drop = FALSE]
  message_info("Topological features retained after lifespan filtering: ", nrow(comp), " / ", n_before)

  phenotypes <- unique(as.character(comp$phenotype))
  if (length(phenotypes) < 2) {
    stop_msg("At least two phenotypes are required in Module 3 composition table. Found: ", paste(phenotypes, collapse = ", "))
  }

  phenotype1 <- get_arg(args, "phenotype1", default = phenotypes[[1]])
  phenotype2 <- get_arg(args, "phenotype2", default = phenotypes[[2]])

  if (!phenotype1 %in% phenotypes) stop_msg("phenotype1 not found in composition table: ", phenotype1)
  if (!phenotype2 %in% phenotypes) stop_msg("phenotype2 not found in composition table: ", phenotype2)
  if (identical(phenotype1, phenotype2)) stop_msg("phenotype1 and phenotype2 must be different.")

  target_phenotypes <- c(phenotype1, phenotype2)
  comp <- comp[comp$phenotype %in% target_phenotypes, , drop = FALSE]

  message_info("Using phenotype1 baseline: ", phenotype1)
  message_info("Using phenotype2 comparison: ", phenotype2)
  message_info("Using composition column: ", composition_col)

  comp_long <- split_composition_to_long(comp, separator = separator)
  if (nrow(comp_long) == 0) {
    stop_msg("No feature composition could be expanded from column: ", composition_col)
  }

  score_df <- calculate_scores(comp_long)

  all_features <- read_feature_list(feature_list_path, feature_col)
  if (is.null(all_features)) {
    all_features <- sort(unique(comp_long$attributed_feature_id))
  }
  if (write_zero_scores) {
    score_df <- complete_zero_scores(score_df, all_features = all_features, phenotypes = target_phenotypes)
  }

  score_wide <- make_score_wide(score_df, phenotypes = target_phenotypes)
  delta_df <- calculate_delta_scores(score_df, phenotype1 = phenotype1, phenotype2 = phenotype2)

  # Output paths
  composition_long_path <- file.path(outdir, paste0(prefix, "_module4_feature_composition_long.tsv"))
  score_long_path <- file.path(outdir, paste0(prefix, "_module4_participation_score_long.tsv"))
  score_wide_path <- file.path(outdir, paste0(prefix, "_module4_participation_score_wide.tsv"))
  delta_path <- file.path(outdir, paste0(prefix, "_module4_delta_participation_score.tsv"))
  top_delta_path <- file.path(outdir, paste0(prefix, "_module4_top_abs_delta_participation_score.tsv"))
  summary_path <- file.path(outdir, paste0(prefix, "_module4_score_summary.tsv"))
  module5_manifest_path <- file.path(outdir, paste0(prefix, "_module5_input_manifest.tsv"))
  run_info_path <- file.path(outdir, paste0(prefix, "_module4_run_info.tsv"))

  write_tsv_base(comp_long, composition_long_path)
  write_tsv_base(score_df, score_long_path)
  write_tsv_base(score_wide, score_wide_path)
  write_tsv_base(delta_df, delta_path)
  write_tsv_base(head(delta_df, top_n), top_delta_path)

  summary_df <- make_summary(comp, comp_long, score_df, delta_df, phenotype1, phenotype2, top_n)
  write_tsv_base(summary_df, summary_path)

  module5_manifest <- data.frame(
    item = c(
      "delta_participation_score",
      "participation_score_long",
      "participation_score_wide",
      "feature_composition_long",
      "score_summary"
    ),
    path = normalize_path(c(
      delta_path,
      score_long_path,
      score_wide_path,
      composition_long_path,
      summary_path
    )),
    stringsAsFactors = FALSE
  )
  write_tsv_base(module5_manifest, module5_manifest_path)

  run_info <- data.frame(
    parameter = c(
      "module",
      "composition_path",
      "outdir",
      "prefix",
      "phenotype1",
      "phenotype2",
      "composition_col",
      "feature_list",
      "feature_col",
      "min_lifespan",
      "include_infinite_death",
      "infinite_lifespan_value",
      "separator",
      "write_zero_scores",
      "top_n",
      "run_time"
    ),
    value = c(
      "module4_participation_score",
      composition_path,
      outdir,
      prefix,
      phenotype1,
      phenotype2,
      composition_col,
      ifelse(is.null(feature_list_path), "NA", normalize_path(feature_list_path)),
      feature_col,
      as.character(min_lifespan),
      as.character(include_inf),
      as.character(inf_lifespan_value),
      separator,
      as.character(write_zero_scores),
      as.character(top_n),
      as.character(Sys.time())
    ),
    stringsAsFactors = FALSE
  )
  write_tsv_base(run_info, run_info_path)

  message_info("Module 4 completed.")
  message_info("Score long table: ", score_long_path)
  message_info("Delta score table: ", delta_path)
  message_info("Module 5 input manifest: ", module5_manifest_path)
}

main()
