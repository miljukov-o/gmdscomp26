# ============================================================
# Simple submission CSV helper
# Version 2026-06-11
# ============================================================
#
# This file helps you create the official submission CSV.
# You can use any statistical method you like. The only requirement is that
# your method returns a few values in the format shown below.
#
# The main idea is simple:
#   - Write one R function that analyzes one dataset.
#   - This file runs that function on all datasets in a folder.
#   - This file writes the final Submission_Results_*.csv file.
#
# If you only want to plug in your own method, read sections 1 and 2 first.
# Section 3 is a small toy example. Sections 4-6 can usually be left unchanged.
#

# ============================================================
# 1. QUICK START
# ============================================================
#
# Step 1: Copy my_method_template() and give the copy a new name.
# For example:
#
#   my_group_method <- function(file_path, data, options = list()) {
#     ... fit your model here ...
#     return(list(... values for the submission CSV ...))
#   }
#
# Step 2: Run the batch helper with your method.
#
#   source("R code/submission_csv_automation_module_2026-06-11.R")
#
#   run_submission_batch(
#     input_folder = "path/to/Blinded_data_sets",
#     output_folder = "path/to/output_folder",
#     method_fun = my_group_method
#   )
#
# The output will be a file named like Submission_Results_2026-...csv.
#

# ============================================================
# 2. METHOD TEMPLATE: copy this function and edit the copy
# ============================================================

my_method_template <- function(file_path, data, options = list()) {
  # This function is called once for each dataset.
  #
  # file_path is the path to the current CSV file.
  # data is the current dataset, already loaded as a data.frame.
  # options is an optional list of tuning choices, if you want to use one.
  #
  # Inside this function, you only need to decide what should go into one row
  # of the submission CSV.
  #
  # The returned list below is the important part. Keep the same field names,
  # and replace the placeholder values with results from your own method.

  x_cols <- get_x_columns(data)

  # Start by giving every X variable a default label.
  # Then change the labels selected by your method.
  # Allowed labels are exactly:
  #   "Neither", "Prognostic", "Predictive", "Both", or NA.
  x_classification <- setNames(rep("Neither", length(x_cols)), x_cols)

  # The values below are placeholders.
  # Replace them with the output of your method.
  list(
    endpoint = guess_endpoint_from_y(data),          # "Binary" or "Continuous".
    heterogeneity = "No",                           # TRUE/FALSE or "Yes"/"No".
    subgroup_proportion = 0,                         # Mean of S, where S=1 is the subgroup.
    treatment_effect_subgroup = NA_real_,            # Estimated treatment effect for S=1.
    treatment_effect_complement = NA_real_,          # Estimated treatment effect for S=0.
    simple_rule = "No",                              # Example: "S=1 if X3 > 0.4". Use "No" if none.
    x_classification = x_classification,
    S = rep(0L, nrow(data)),
    metadata = list(method = "template_method")      # Optional notes; safe to remove.
  )
}

# ============================================================
# 3. SMALL TOY EXAMPLE
# ============================================================
#
# This example is deliberately simple. It is only here to show where the pieces go.
#
# The example uses this very simple rule:
#   S = 1 if X1 is above its median.
# Then it compares the treatment effect in S=1 with the treatment effect in S=0.
#
# You can run this example first just to check that the CSV automation works.
# After that, replace it with your real method.

example_x1_split_method <- function(file_path, data, options = list()) {
  x_cols <- get_x_columns(data)
  endpoint <- guess_endpoint_from_y(data)

  if (!all(c("Y", "W") %in% names(data))) {
    stop("Expected columns named Y and W.")
  }

  if (!("X1" %in% names(data))) {
    # If X1 is not available, return a valid row with no subgroup.
    return(my_method_template(file_path, data, options))
  }

  cutoff <- stats::median(data$X1, na.rm = TRUE)
  S <- as.integer(data$X1 > cutoff)

  te_subgroup <- estimate_treatment_effect(data$Y[S == 1], data$W[S == 1])
  te_complement <- estimate_treatment_effect(data$Y[S == 0], data$W[S == 0])

  # This threshold is only for the toy example.
  # In a real method, choose this decision rule carefully.
  threshold <- options$heterogeneity_threshold %||% 0.1
  has_heterogeneity <- is.finite(te_subgroup) && is.finite(te_complement) &&
    abs(te_subgroup - te_complement) > threshold

  x_classification <- setNames(rep("Neither", length(x_cols)), x_cols)
  if (has_heterogeneity) x_classification["X1"] <- "Predictive"

  list(
    endpoint = endpoint,
    heterogeneity = has_heterogeneity,
    subgroup_proportion = mean(S == 1),
    treatment_effect_subgroup = te_subgroup,
    treatment_effect_complement = te_complement,
    simple_rule = if (has_heterogeneity) paste0("S=1 if X1 > ", signif(cutoff, 3)) else "No",
    x_classification = x_classification,
    S = S,
    metadata = list(
      method = "example_x1_split_method",
      rule_cutoff = cutoff,
      n = nrow(data),
      p = length(x_cols)
    )
  )
}

# ============================================================
# 4. BATCH RUNNER: you usually do not need to edit this
# ============================================================

run_submission_batch <- function(input_folder,
                                 output_folder,
                                 method_fun,
                                 options = list(),
                                 max_cov = 50,
                                 file_pattern = "\\.csv$",
                                 output_prefix = "Submission_Results",
                                 metadata_prefix = "Analysis_Metadata",
                                 write_metadata = TRUE,
                                 # If write_s_files = TRUE, your method must also return S:
                                 #   S = integer vector of length nrow(data), containing only 0 and 1.
                                 write_s_files = FALSE,
                                 s_output_subfolder = "Datasets_with_S",
                                 s_column_name = "S",
                                 include_error_rows = TRUE,
                                 verbose = TRUE) {
  if (!dir.exists(input_folder)) stop("input_folder does not exist: ", input_folder)
  if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)
  if (!is.function(method_fun)) stop("method_fun must be a function.")
  s_output_folder <- NULL
  if (isTRUE(write_s_files)) {
    s_output_folder <- file.path(output_folder, s_output_subfolder)
    if (!dir.exists(s_output_folder)) dir.create(s_output_folder, recursive = TRUE)
  }

  file_list <- list.files(input_folder, pattern = file_pattern, full.names = TRUE)

  # If old result files are in the same folder, skip them.
  file_list <- file_list[!grepl(
    "Submission_Results|Analysis_Metadata|_with_S|_results|_variable_classification|_patient_effects|_abc_constraint",
    basename(file_list),
    ignore.case = TRUE
  )]
  file_list <- sort(file_list)

  if (length(file_list) == 0) stop("No raw dataset CSV files found in: ", input_folder)
  if (verbose) cat("Found", length(file_list), "dataset CSV files.\n")

  submission_rows <- list()
  metadata_rows <- list()

  for (file_path in file_list) {
    if (verbose) cat("\nProcessing:", basename(file_path), "\n")

    one_result <- tryCatch({
      data <- read.csv(file_path, check.names = FALSE)

      # This is the key line: run the chosen method on this one dataset.
      method_output <- method_fun(file_path = file_path, data = data, options = options)

      submission_row <- make_submission_row(method_output, file_path, max_cov = max_cov)

      s_file_path <- NULL
      if (isTRUE(write_s_files)) {
        s_file_path <- write_dataset_with_s(
          data = data,
          method_output = method_output,
          file_path = file_path,
          output_folder = s_output_folder,
          s_column_name = s_column_name
        )
      }
      
      if (verbose) {
        cat("  Endpoint:", submission_row[["Endpoint"]],
            "| Heterogeneity:", submission_row[["Subgroup/Treatment Effect Heterogeneity (Yes/No)"]],
            "| Subgroup proportion:", submission_row[["Subgroup Proportion"]], "\n")
      }

      list(
        ok = TRUE,
        submission_row = submission_row,
        metadata = method_output$metadata %||% NULL,
        s_file_path = s_file_path
      )
    }, error = function(e) {
      if (verbose) cat("  ERROR:", conditionMessage(e), "\n")

      if (!include_error_rows) {
        return(list(ok = FALSE, submission_row = NULL, metadata = NULL))
      }

      list(
        ok = FALSE,
        submission_row = make_error_row(file_path, conditionMessage(e), max_cov = max_cov),
        metadata = NULL
      )
    })

    if (!is.null(one_result$submission_row)) {
      submission_rows[[length(submission_rows) + 1]] <- one_result$submission_row
    }

    if (!is.null(one_result$metadata)) {
      metadata_rows[[length(metadata_rows) + 1]] <- make_metadata_row(one_result$metadata, file_path)
    }
  }

  final_submission <- do.call(rbind, submission_rows)
  final_submission <- final_submission[, submission_column_names(max_cov), drop = FALSE]

  run_timestamp <- format(Sys.time(), "%Y-%m-%d_%H%M%S")
  submission_csv_path <- file.path(output_folder, paste0(output_prefix, "_", run_timestamp, ".csv"))
  write.csv(final_submission, submission_csv_path, row.names = FALSE, na = "")

  metadata_csv_path <- NULL
  if (isTRUE(write_metadata) && length(metadata_rows) > 0) {
    metadata_df <- do.call(rbind, metadata_rows)
    metadata_csv_path <- file.path(output_folder, paste0(metadata_prefix, "_", run_timestamp, ".csv"))
    write.csv(metadata_df, metadata_csv_path, row.names = FALSE, na = "")
  }

  if (verbose) {
    cat("\nProcessing complete.\n")
    cat("Submission CSV saved to:\n", submission_csv_path, "\n", sep = "")
    if (!is.null(metadata_csv_path)) cat("Metadata CSV saved to:\n", metadata_csv_path, "\n", sep = "")
    if (!is.null(s_output_folder)) cat("Datasets with S saved to:\n", s_output_folder, "\n", sep = "")
  }

  invisible(list(
    submission = final_submission,
    submission_csv_path = submission_csv_path,
    metadata_csv_path = metadata_csv_path,
    s_output_folder = s_output_folder
  ))
}

# ============================================================
# 5. OPTIONAL EXAMPLE: using an existing method
# ============================================================
#
# You can ignore this section if you are writing your own method.
#
# This section only shows how an existing analysis function can be connected to
# the same CSV helper. The existing function must analyze one dataset and return
# its own result object. The small wrapper below translates that result object
# into the simple list format used in section 2.
#
# Here the existing function is analyze_one_dataset_full_abc().
# It must be loaded before latest_abc_method() is used.
#
# Example:
#
#   source("R code/submission_csv_automation_module_2026-06-11.R")
#   source("R code/abc_algorithm_core_only.R")
#
#   run_submission_batch(
#     input_folder = "path/to/Blinded_data_sets",
#     output_folder = "path/to/output_folder",
#     method_fun = latest_abc_method
#   )
#

latest_abc_method <- function(file_path, data, options = list()) {
  if (!exists("analyze_one_dataset_full_abc", mode = "function")) {
    stop("analyze_one_dataset_full_abc() is not loaded yet.")
  }

  # Keep the call short. The ABC function already contains its own defaults.
  # options can still be used to override anything if needed.
  abc_options <- modifyList(list(
    file_path = file_path,
    out_dir = tempdir(),
    write_outputs = FALSE,
    verbose = FALSE
  ), options)

  result <- do.call(analyze_one_dataset_full_abc, abc_options)

  simple_rule <- if (isTRUE(result$results$simple_rule_available[1])) {
    result$results$simple_rule_text[1]
  } else {
    "No"
  }

  x_classification <- result$variable_classification$Classification
  names(x_classification) <- result$variable_classification$Variable

  list(
    endpoint = result$results$endpoint_type[1],
    heterogeneity = isTRUE(result$results$heterogeneity_detected[1]),
    subgroup_proportion = result$results$subgroup_proportion[1],
    treatment_effect_subgroup = result$results$treatment_effect_subgroup_S1[1],
    treatment_effect_complement = result$results$treatment_effect_complement_S0[1],
    simple_rule = simple_rule,
    x_classification = x_classification,
    # TODO:
     # If the ABC function returns a subgroup assignment vector, include it here. e.g.
     # S = result$subgroup_assignment,
    metadata = result$results
  )
}

# ============================================================
# 6. HELPER FUNCTIONS: no need to edit below this line
# ============================================================

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

submission_column_names <- function(max_cov = 50) {
  c(
    "Dataset",
    "Endpoint",
    "Subgroup/Treatment Effect Heterogeneity (Yes/No)",
    "Subgroup Proportion",
    "Treatment Effect in Subgroup",
    "Treatment Effect in Complement",
    "Simple Rule for Subgroup Assignment (Rule/No)",
    paste0("X", seq_len(max_cov), " (Neither/Prognostic/Predictive/Both/NA)")
  )
}

submission_dataset_name <- function(file_path) {
  tools::file_path_sans_ext(basename(file_path))
}

empty_submission_row <- function(max_cov = 50) {
  out <- as.data.frame(as.list(rep(NA_character_, length(submission_column_names(max_cov)))),
                       stringsAsFactors = FALSE)
  names(out) <- submission_column_names(max_cov)
  out
}

get_x_columns <- function(data) {
  x_cols <- grep("^[Xx][0-9]+$", names(data), value = TRUE)
  x_numbers <- as.integer(sub("^[Xx]", "", x_cols))
  x_cols[order(x_numbers)]
}

guess_endpoint_from_y <- function(data, y_col = "Y") {
  if (!(y_col %in% names(data))) return(NA_character_)
  y <- data[[y_col]]
  y <- y[!is.na(y)]
  if (length(y) == 0) return(NA_character_)

  unique_y <- sort(unique(y))
  if (length(unique_y) <= 2 && all(unique_y %in% c(0, 1))) "Binary" else "Continuous"
}

estimate_treatment_effect <- function(y, w) {
  y <- as.numeric(y)
  w <- as.numeric(w)
  ok <- is.finite(y) & !is.na(w)
  y <- y[ok]
  w <- w[ok]

  if (!any(w == 1) || !any(w == 0)) return(NA_real_)
  mean(y[w == 1]) - mean(y[w == 0])
}

first_value_or_na <- function(x, default = NA) {
  if (is.null(x) || length(x) == 0) return(default)
  x[[1]]
}

round_or_na <- function(x, digits = 4) {
  x <- suppressWarnings(as.numeric(first_value_or_na(x, NA_real_)))
  if (is.na(x) || !is.finite(x)) return(NA_real_)
  round(x, digits)
}

normalize_endpoint <- function(endpoint) {
  endpoint <- tolower(as.character(first_value_or_na(endpoint, NA_character_)))
  if (is.na(endpoint) || !nzchar(endpoint)) return(NA_character_)
  if (endpoint %in% c("binary", "binomial", "logistic")) return("Binary")
  if (endpoint %in% c("continuous", "gaussian", "linear", "normal")) return("Continuous")
  stop("Endpoint must be 'Binary' or 'Continuous'. Got: ", endpoint)
}

normalize_yes_no <- function(x) {
  if (is.logical(x) && length(x) > 0) return(if (isTRUE(x[[1]])) "Yes" else "No")
  x <- tolower(as.character(first_value_or_na(x, NA_character_)))
  if (is.na(x) || !nzchar(x)) return(NA_character_)
  if (x %in% c("yes", "y", "true", "1", "heterogeneous")) return("Yes")
  if (x %in% c("no", "n", "false", "0", "none", "homogeneous")) return("No")
  stop("Heterogeneity must be TRUE/FALSE or Yes/No. Got: ", x)
}

normalize_simple_rule <- function(rule, heterogeneity) {
  rule <- as.character(first_value_or_na(rule, NA_character_))
  if (is.na(rule) || !nzchar(rule)) return("No")
  if (identical(heterogeneity, "No")) return("No")
  rule
}

normalize_x_classification <- function(x_classification, max_cov = 50) {
  allowed <- c("Neither", "Prognostic", "Predictive", "Both", NA_character_)
  x_out <- rep(NA_character_, max_cov)

  if (is.null(x_classification)) return(x_out)

  if (is.data.frame(x_classification)) {
    lower_names <- tolower(names(x_classification))

    if (all(c("variable", "classification") %in% lower_names)) {
      variable_col <- names(x_classification)[match("variable", lower_names)]
      class_col <- names(x_classification)[match("classification", lower_names)]
      values <- as.character(x_classification[[class_col]])
      names(values) <- as.character(x_classification[[variable_col]])
      return(normalize_x_classification(values, max_cov = max_cov))
    }

    x_cols <- grep("^[Xx][0-9]+$", names(x_classification), value = TRUE)
    if (length(x_cols) > 0) {
      values <- as.character(x_classification[1, x_cols, drop = TRUE])
      names(values) <- x_cols
      return(normalize_x_classification(values, max_cov = max_cov))
    }

    stop("x_classification data.frame must have Variable/Classification columns or X1, X2, ... columns.")
  }

  values <- as.character(x_classification)
  x_names <- names(x_classification)
  if (is.null(x_names) || any(!nzchar(x_names))) {
    stop("x_classification must be named, e.g. c(X1 = 'Predictive', X2 = 'Neither').")
  }

  x_numbers <- suppressWarnings(as.integer(sub("^[Xx]", "", x_names)))
  keep <- !is.na(x_numbers) & x_numbers >= 1 & x_numbers <= max_cov
  x_out[x_numbers[keep]] <- values[keep]

  invalid <- !(x_out %in% allowed)
  if (any(invalid, na.rm = TRUE)) {
    bad <- unique(x_out[invalid & !is.na(x_out)])
    stop("Invalid X classification value(s): ", paste(bad, collapse = ", "),
         ". Use Neither, Prognostic, Predictive, Both, or NA.")
  }

  x_out
}

make_submission_row <- function(method_output, file_path, max_cov = 50) {
  if (!is.list(method_output)) stop("Your method must return a named list.")

  endpoint <- normalize_endpoint(method_output$endpoint %||% method_output$endpoint_type)
  heterogeneity <- normalize_yes_no(method_output$heterogeneity %||% method_output$heterogeneity_detected)
  simple_rule <- normalize_simple_rule(method_output$simple_rule %||% method_output$simple_rule_text, heterogeneity)
  x_classes <- normalize_x_classification(
    method_output$x_classification %||% method_output$variable_classification,
    max_cov = max_cov
  )

  row_df <- data.frame(
    `Dataset` = submission_dataset_name(file_path),
    `Endpoint` = endpoint,
    `Subgroup/Treatment Effect Heterogeneity (Yes/No)` = heterogeneity,
    `Subgroup Proportion` = round_or_na(method_output$subgroup_proportion),
    `Treatment Effect in Subgroup` = round_or_na(method_output$treatment_effect_subgroup %||%
                                                   method_output$treatment_effect_subgroup_S1),
    `Treatment Effect in Complement` = round_or_na(method_output$treatment_effect_complement %||%
                                                     method_output$treatment_effect_complement_S0),
    `Simple Rule for Subgroup Assignment (Rule/No)` = simple_rule,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  x_df <- as.data.frame(as.list(x_classes), stringsAsFactors = FALSE)
  names(x_df) <- paste0("X", seq_len(max_cov), " (Neither/Prognostic/Predictive/Both/NA)")

  cbind(row_df, x_df)
}

make_error_row <- function(file_path, error_message, max_cov = 50) {
  row <- empty_submission_row(max_cov = max_cov)
  row[["Dataset"]] <- submission_dataset_name(file_path)
  row[["Subgroup/Treatment Effect Heterogeneity (Yes/No)"]] <- paste("ERROR:", error_message)
  row
}

make_metadata_row <- function(metadata, file_path) {
  metadata <- as.list(metadata)
  metadata <- lapply(metadata, function(value) {
    if (length(value) == 0) return(NA)
    if (length(value) > 1) return(paste(value, collapse = ";"))
    value
  })
  do.call(data.frame, c(
    list(Dataset = submission_dataset_name(file_path)),
    metadata,
    list(stringsAsFactors = FALSE, check.names = FALSE)
  ))
}

get_s_vector <- function(method_output, n) {
  S <- method_output$S %||% method_output$s %||% method_output$subgroup_assignment
  
  if (is.null(S)) {
    stop("write_s_files=TRUE requires method_output$S, method_output$s, or method_output$subgroup_assignment.")
  }
  
  if (is.data.frame(S)) {
    if ("S" %in% names(S)) {
      S <- S[["S"]]
    } else if ("s" %in% names(S)) {
      S <- S[["s"]]
    } else {
      stop("S data.frame must contain a column named S or s.")
    }
  }
  
  if (length(S) != n) {
    stop("S must have length nrow(data). Got length ", length(S), " for n = ", n, ".")
  }
  
  S <- suppressWarnings(as.integer(S))
  
  if (any(is.na(S)) || any(!(S %in% c(0L, 1L)))) {
    stop("S must contain only 0 and 1, with no missing values.")
  }
  
  S
}

write_dataset_with_s <- function(data,
                                 method_output,
                                 file_path,
                                 output_folder,
                                 s_column_name = "S") {
  S <- get_s_vector(method_output, n = nrow(data))
  
  out <- data
  out[[s_column_name]] <- S
  
  if ("ptid" %in% names(out)) {
    out <- out[order(out$ptid), , drop = FALSE]
  }
  
  out_path <- file.path(
    output_folder,
    paste0(submission_dataset_name(file_path), "_with_S.csv")
  )
  
  write.csv(out, out_path, row.names = FALSE, na = "")
  out_path
}
