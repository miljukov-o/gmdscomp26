format_number <- function(x) {
  formatC(x, digits = 2, format = "f")
}

extract_dataset_id_from_path <- function(path) {
  file_name <- fs::path_file(path)
  id_text <- stringr::str_match(file_name, ".*_\\s*(\\d+)\\s*\\.csv$")[, 2]
  as.integer(id_text)
}

get_x_cols <- function(data) {
  names(data) |>
    stringr::str_subset(stringr::regex("^x\\d+$", ignore_case = TRUE)) |>
    stringr::str_sort(numeric = TRUE)
}

detect_endpoint_type <- function(data) {
  y_values <- sort(unique(data$Y))
  
  if (length(y_values) <= 2 && all(y_values %in% c(0, 1))) {
    return("binary")
  }
  
  "continuous"
}

effect_scale_label <- function(endpoint_type) {
  if (identical(endpoint_type, "binary")) {
    return("Difference in response rates")
  }
  
  "Difference in means"
}

outcome_scale_label <- function(endpoint_type) {
  if (identical(endpoint_type, "binary")) {
    return("Response rate")
  }
  
  "Mean outcome"
}

is_binary_covariate <- function(x) {
  x_values <- sort(unique(stats::na.omit(x)))
  length(x_values) <= 2 && all(x_values %in% c(0, 1))
}

column_is_binary <- function(x_var, data) {
  is_binary_covariate(data[[x_var]])
}

default_n_bins <- function(n) {
  if (n <= 150) {
    return(4L)
  }
  
  if (n <= 500) {
    return(6L)
  }
  
  8L
}

compute_overall_effect <- function(data) {
  treatment_mean <- data |>
    dplyr::filter(.data$W == 1) |>
    dplyr::summarise(value = mean(.data$Y), .groups = "drop") |>
    dplyr::pull(.data$value)
  
  control_mean <- data |>
    dplyr::filter(.data$W == 0) |>
    dplyr::summarise(value = mean(.data$Y), .groups = "drop") |>
    dplyr::pull(.data$value)
  
  treatment_mean - control_mean
}

make_binary_grouped_data <- function(data, x_var) {
  data |>
    dplyr::mutate(
      ptid = .data$ptid,
      Y = .data$Y,
      W = .data$W,
      x_value = as.numeric(.data[[x_var]])
    ) |>
    dplyr::mutate(
      group_id = dplyr::if_else(.data$x_value == 0, 1L, 2L),
      x_center = .data$x_value,
      x_label_short = as.character(.data$x_value),
      x_label_long = as.character(.data$x_value)
    )
}

make_continuous_grouped_data <- function(data, x_var, n_bins) {
  data |>
    dplyr::mutate(
      ptid = .data$ptid,
      Y = .data$Y,
      W = .data$W,
      x_value = .data[[x_var]]
    ) |>
    dplyr::mutate(
      group_id = dplyr::ntile(.data$x_value, n_bins)
    ) |>
    dplyr::group_by(.data$group_id) |>
    dplyr::mutate(
      x_center = stats::median(.data$x_value),
      x_label_short = paste0("B", dplyr::first(.data$group_id)),
      x_label_long = glue::glue(
        "{format_number(min(.data$x_value))} to {format_number(max(.data$x_value))}"
      )
    ) |>
    dplyr::ungroup()
}

make_grouped_data <- function(data, x_var, n_bins) {
  if (column_is_binary(x_var, data)) {
    return(make_binary_grouped_data(data, x_var))
  }
  
  make_continuous_grouped_data(data, x_var, n_bins)
}

summarise_arm_level <- function(grouped_data) {
  grouped_data |>
    dplyr::group_by(
      .data$group_id,
      .data$x_center,
      .data$x_label_short,
      .data$x_label_long,
      .data$W
    ) |>
    dplyr::summarise(
      n = dplyr::n(),
      y_mean = mean(.data$Y),
      y_var = stats::var(.data$Y),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      arm = dplyr::if_else(.data$W == 1, "treatment", "control")
    ) |>
    dplyr::select(-.data$W) |>
    tidyr::pivot_wider(
      names_from = .data$arm,
      values_from = c(.data$n, .data$y_mean, .data$y_var),
      names_sep = "_",
      values_fill = list(n = 0, y_mean = NA_real_, y_var = NA_real_)
    ) |>
    dplyr::arrange(.data$group_id)
}

compute_effect_se <- function(arm_summary, endpoint_type) {
  if (identical(endpoint_type, "binary")) {
    return(
      sqrt(
        (arm_summary$y_mean_treatment * (1 - arm_summary$y_mean_treatment) / arm_summary$n_treatment) +
          (arm_summary$y_mean_control * (1 - arm_summary$y_mean_control) / arm_summary$n_control)
      )
    )
  }
  
  treatment_var <- dplyr::if_else(
    arm_summary$n_treatment > 1,
    arm_summary$y_var_treatment,
    NA_real_
  )
  
  control_var <- dplyr::if_else(
    arm_summary$n_control > 1,
    arm_summary$y_var_control,
    NA_real_
  )
  
  sqrt((treatment_var / arm_summary$n_treatment) + (control_var / arm_summary$n_control))
}

summarise_variable_effect <- function(x_var, data, endpoint_type, n_bins) {
  grouped_data <- make_grouped_data(data = data, x_var = x_var, n_bins = n_bins)
  
  arm_summary <- summarise_arm_level(grouped_data = grouped_data)
  
  arm_summary$se <- compute_effect_se(
    arm_summary = arm_summary,
    endpoint_type = endpoint_type
  )
  
  var_type <- if (column_is_binary(x_var, data)) "binary" else "continuous"
  
  arm_summary |>
    dplyr::mutate(
      x_var = x_var,
      var_type = var_type,
      n_total = .data$n_treatment + .data$n_control,
      effect = .data$y_mean_treatment - .data$y_mean_control,
      ci_low = .data$effect - 1.96 * .data$se,
      ci_high = .data$effect + 1.96 * .data$se
    ) |>
    dplyr::select(
      .data$x_var,
      .data$var_type,
      .data$group_id,
      .data$x_center,
      .data$x_label_short,
      .data$x_label_long,
      .data$n_control,
      .data$n_treatment,
      .data$n_total,
      .data$y_mean_control,
      .data$y_mean_treatment,
      .data$effect,
      .data$se,
      .data$ci_low,
      .data$ci_high
    )
}

rank_to_unit <- function(x) {
  out <- rep(NA_real_, length(x))
  ok <- !is.na(x)
  n_ok <- sum(ok)
  
  if (n_ok == 0L) {
    return(out)
  }
  
  if (n_ok == 1L) {
    out[ok] <- 0.5
    return(out)
  }
  
  out[ok] <- (rank(x[ok], ties.method = "average") - 1) / (n_ok - 1)
  out
}

safe_range <- function(x) {
  if (all(is.na(x))) {
    return(NA_real_)
  }
  
  diff(range(x, na.rm = TRUE))
}

compute_trend_strength <- function(group_id, effect) {
  keep <- is.finite(group_id) & is.finite(effect)
  
  if (sum(keep) < 3L || dplyr::n_distinct(group_id[keep]) < 3L) {
    return(0)
  }
  
  abs(suppressWarnings(
    stats::cor(group_id[keep], effect[keep], method = "spearman")
  ))
}

summarise_variable_priority <- function(effect_data_var, overall_effect, min_arm_n = 10L) {
  arm_sizes <- pmin(effect_data_var$n_control, effect_data_var$n_treatment)
  
  tibble::tibble(
    weighted_abs_deviation = stats::weighted.mean(
      x = abs(effect_data_var$effect - overall_effect),
      w = effect_data_var$n_total,
      na.rm = TRUE
    ),
    effect_range = safe_range(effect_data_var$effect),
    trend_strength = compute_trend_strength(
      group_id = effect_data_var$group_id,
      effect = effect_data_var$effect
    ),
    qualitative_signal = as.integer(
      min(effect_data_var$effect, na.rm = TRUE) < 0 &&
        max(effect_data_var$effect, na.rm = TRUE) > 0
    ),
    min_arm_n_observed = min(arm_sizes, na.rm = TRUE),
    unstable_fraction = mean(arm_sizes < min_arm_n, na.rm = TRUE),
    n_groups = dplyr::n_distinct(effect_data_var$group_id)
  )
}

build_variable_priority_table <- function(effect_data, overall_effect, min_arm_n = 10L) {
  priority_data <- effect_data |>
    dplyr::group_by(.data$x_var, .data$var_type) |>
    dplyr::group_modify(
      ~ summarise_variable_priority(
        effect_data_var = .x,
        overall_effect = overall_effect,
        min_arm_n = min_arm_n
      )
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      deviation_rank = rank_to_unit(.data$weighted_abs_deviation),
      range_rank = rank_to_unit(.data$effect_range),
      trend_rank = rank_to_unit(.data$trend_strength),
      priority_score =
        0.45 * .data$deviation_rank +
        0.25 * .data$range_rank +
        0.15 * .data$trend_rank +
        0.25 * .data$qualitative_signal -
        0.25 * .data$unstable_fraction,
      priority_class = dplyr::case_when(
        .data$priority_score >= 0.75 & .data$unstable_fraction <= 0.25 ~ "high",
        .data$priority_score >= 0.45 ~ "medium",
        TRUE ~ "low"
      ),
      priority_class = factor(
        .data$priority_class,
        levels = c("high", "medium", "low")
      ),
      priority_color = dplyr::case_when(
        .data$priority_class == "high" ~ "#b22222",
        .data$priority_class == "medium" ~ "#c77f00",
        TRUE ~ "grey55"
      )
    ) |>
    dplyr::arrange(.data$priority_class, dplyr::desc(.data$priority_score), .data$x_var) |>
    dplyr::mutate(
      priority_rank = dplyr::row_number()
    )
  
  priority_data
}

priority_table_for_display <- function(priority_data) {
  priority_data |>
    dplyr::mutate(
      x_var = .data$x_var,
      var_type = .data$var_type,
      weighted_abs_deviation = .data$weighted_abs_deviation,
      effect_range = .data$effect_range,
      trend_strength = .data$trend_strength,
      qualitative_signal = .data$qualitative_signal,
      min_arm_n_observed = .data$min_arm_n_observed,
      unstable_fraction = .data$unstable_fraction,
      priority_score = .data$priority_score,
      priority_class = .data$priority_class
    )
}

summarise_dataset_priority <- function(priority_data) {
  counts <- priority_data |>
    dplyr::count(.data$priority_class, name = "n_vars") |>
    tidyr::pivot_wider(
      names_from = .data$priority_class,
      values_from = .data$n_vars,
      values_fill = 0
    )
  
  if (!"high" %in% names(counts)) counts$high <- 0L
  if (!"medium" %in% names(counts)) counts$medium <- 0L
  if (!"low" %in% names(counts)) counts$low <- 0L
  
  top_vars <- priority_data |>
    dplyr::arrange(.data$priority_class, dplyr::desc(.data$priority_score), .data$x_var) |>
    dplyr::slice_head(n = 5) |>
    dplyr::pull(.data$x_var) |>
    as.character() |>
    paste(collapse = ", ")
  
  tibble::tibble(
    n_high = as.integer(counts$high[[1]]),
    n_medium = as.integer(counts$medium[[1]]),
    n_low = as.integer(counts$low[[1]]),
    prop_high = as.integer(counts$high[[1]]) / nrow(priority_data),
    prop_medium_or_higher = (as.integer(counts$high[[1]]) + as.integer(counts$medium[[1]])) / nrow(priority_data),
    max_priority_score = max(priority_data$priority_score, na.rm = TRUE),
    mean_priority_score = mean(priority_data$priority_score, na.rm = TRUE),
    top_priority_vars = top_vars
  )
}

derive_dataset_priority_hint <- function(dataset_priority_summary) {
  dplyr::case_when(
    dataset_priority_summary$n_high >= 3 ~ "multiple_high_priority_signals",
    dataset_priority_summary$n_high >= 1 & dataset_priority_summary$n_medium >= 2 ~ "high_plus_additional_support",
    dataset_priority_summary$n_high >= 1 ~ "single_high_priority_signal",
    dataset_priority_summary$n_medium >= 3 ~ "several_medium_priority_signals",
    dataset_priority_summary$n_medium >= 1 ~ "weak_to_moderate_signal_only",
    TRUE ~ "mostly_low_priority_patterns"
  )
}

build_dataset_priority_summary <- function(priority_data) {
  summary_tbl <- summarise_dataset_priority(priority_data)
  
  summary_tbl |>
    dplyr::mutate(
      dataset_priority_hint = derive_dataset_priority_hint(dplyr::pick(dplyr::everything()))
    )
}

summarise_task1_dataset_from_file <- function(csv_file, n_bins = NULL, min_arm_n = 10L) {
  data <- readr::read_csv(csv_file, show_col_types = FALSE)
  dataset_id <- extract_dataset_id_from_path(csv_file)
  
  visual_object <- build_task1_visual_object(
    data = data,
    dataset_id = dataset_id,
    n_bins = n_bins,
    min_arm_n = min_arm_n
  )
  
  build_dataset_priority_summary(visual_object$priority_data) |>
    dplyr::mutate(
      dataset_id = dataset_id,
      endpoint = visual_object$endpoint_type,
      n_obs = nrow(data),
      n_vars = length(visual_object$x_cols),
      n_bins = visual_object$n_bins_used,
      overall_effect = visual_object$overall_effect,
      .before = 1
    )
}

build_task1_dataset_triage_table <- function(csv_files, n_bins = NULL, min_arm_n = 10L) {
  purrr::map(
    csv_files,
    summarise_task1_dataset_from_file,
    n_bins = n_bins,
    min_arm_n = min_arm_n
  ) |>
    purrr::list_rbind() |>
    dplyr::arrange(
      dplyr::desc(.data$n_high),
      dplyr::desc(.data$n_medium),
      dplyr::desc(.data$max_priority_score),
      .data$dataset_id
    )
}

summarise_dataset_manifest_entry <- function(csv_file, n_bins = NULL) {
  data <- readr::read_csv(csv_file, show_col_types = FALSE)
  x_cols <- get_x_cols(data)
  
  binary_flags <- purrr::map(
    x_cols,
    column_is_binary,
    data = data
  )
  
  binary_flags <- unlist(binary_flags, use.names = FALSE)
  
  tibble::tibble(
    dataset_id = extract_dataset_id_from_path(csv_file),
    endpoint = detect_endpoint_type(data),
    n_obs = nrow(data),
    n_vars = length(x_cols),
    binary_covariates = sum(binary_flags),
    continuous_covariates = sum(!binary_flags),
    n_bins = if (is.null(n_bins)) default_n_bins(nrow(data)) else as.integer(n_bins)
  )
}

build_task1_visual_object <- function(data, dataset_id = NA_integer_, n_bins = NULL, min_arm_n = 10L) {
  endpoint_type <- detect_endpoint_type(data)
  x_cols <- get_x_cols(data)
  n_bins_used <- if (is.null(n_bins)) default_n_bins(nrow(data)) else as.integer(n_bins)
  overall_effect <- compute_overall_effect(data)
  
  effect_tables <- purrr::map(
    x_cols,
    summarise_variable_effect,
    data = data,
    endpoint_type = endpoint_type,
    n_bins = n_bins_used
  )
  
  effect_data <- effect_tables |>
    purrr::list_rbind()
  
  priority_data <- build_variable_priority_table(
    effect_data = effect_data,
    overall_effect = overall_effect,
    min_arm_n = min_arm_n
  )
  
  effect_data <- effect_data |>
    dplyr::left_join(
      priority_data |>
        dplyr::select(
          .data$x_var,
          .data$var_type,
          .data$priority_score,
          .data$priority_class,
          .data$priority_color,
          .data$priority_rank,
          .data$weighted_abs_deviation,
          .data$effect_range,
          .data$trend_strength,
          .data$qualitative_signal,
          .data$min_arm_n_observed,
          .data$unstable_fraction
        ),
      by = c("x_var", "var_type")
    ) |>
    dplyr::mutate(
      x_var = factor(.data$x_var, levels = x_cols)
    )
  
  priority_data <- priority_data |>
    dplyr::mutate(
      x_var = factor(.data$x_var, levels = x_cols)
    )
  
  dataset_priority_summary <- build_dataset_priority_summary(priority_data)
  
  list(
    dataset_id = dataset_id,
    endpoint_type = endpoint_type,
    x_cols = x_cols,
    n_bins_used = n_bins_used,
    overall_effect = overall_effect,
    effect_label = effect_scale_label(endpoint_type),
    outcome_label = outcome_scale_label(endpoint_type),
    effect_data = effect_data,
    priority_data = priority_data,
    dataset_priority_summary = dataset_priority_summary
  )
}

plot_single_effect_curve <- function(effect_data_var, overall_effect, effect_label) {
  title_text <- as.character(effect_data_var$x_var[[1]])
  priority_class <- as.character(effect_data_var$priority_class[[1]])
  priority_score <- effect_data_var$priority_score[[1]]
  priority_color <- effect_data_var$priority_color[[1]]
  
  subtitle_text <- paste0(
    effect_data_var$var_type[[1]],
    " covariate | priority: ",
    priority_class#,
   # " | score: ",
  #  format_number(priority_score)
  )
  
  ggplot2::ggplot(
    effect_data_var,
    ggplot2::aes(x = .data$group_id, y = .data$effect)
  ) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.3) +
    ggplot2::geom_hline(
      yintercept = overall_effect,
      linetype = 2,
      linewidth = 0.3
    ) +
    ggplot2::geom_linerange(
      ggplot2::aes(ymin = .data$ci_low, ymax = .data$ci_high),
      linewidth = 0.3
    ) +
    ggplot2::geom_line(linewidth = 0.4) +
    ggplot2::geom_point(
      ggplot2::aes(size = .data$n_total),
      alpha = 0.9
    ) +
    ggplot2::scale_x_continuous(
      breaks = effect_data_var$group_id,
      labels = effect_data_var$x_label_short
    ) +
    ggplot2::scale_size_continuous(range = c(1.5, 4), guide = "none") +
    ggplot2::labs(
      title = title_text,
      subtitle = subtitle_text,
      x = NULL,
      y = effect_label
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", colour = priority_color),
      plot.subtitle = ggplot2::element_text(colour = priority_color),
      plot.background = ggplot2::element_rect(
        fill = NA,
        colour = priority_color,
        linewidth = 0.9
      ),
      axis.text.x = ggplot2::element_text(size = 8)
    )
}

make_plot_for_variable <- function(var_name, effect_data, overall_effect, effect_label) {
  effect_data_var <- effect_data |>
    dplyr::filter(.data$x_var == .env$var_name)
  
  plot_single_effect_curve(
    effect_data_var = effect_data_var,
    overall_effect = overall_effect,
    effect_label = effect_label
  )
}

split_x_vars_into_pages <- function(x_vars, page_size) {
  split(x_vars, ceiling(seq_along(x_vars) / page_size))
}

make_overview_page_plot <- function(x_vars, effect_data, overall_effect, effect_label, ncol) {
  plot_list <- purrr::map(
    x_vars,
    make_plot_for_variable,
    effect_data = effect_data,
    overall_effect = overall_effect,
    effect_label = effect_label
  )
  
  patchwork::wrap_plots(plot_list, ncol = ncol)
}

plot_task1_overview_pages <- function(effect_data, overall_effect, effect_label, page_size = 9L, ncol = 3L) {
  x_vars <- effect_data |>
    dplyr::distinct(.data$x_var, .data$priority_class, .data$priority_score) |>
    dplyr::mutate(
      priority_class = factor(.data$priority_class, levels = c("high", "medium", "low"))
    ) |>
    dplyr::arrange(.data$priority_class, dplyr::desc(.data$priority_score), .data$x_var) |>
    dplyr::pull(.data$x_var) |>
    as.character()
  
  x_pages <- split_x_vars_into_pages(x_vars = x_vars, page_size = page_size)
  
  purrr::map(
    x_pages,
    make_overview_page_plot,
    effect_data = effect_data,
    overall_effect = overall_effect,
    effect_label = effect_label,
    ncol = ncol
  )
}

plot_single_detail_curve <- function(effect_data_var, overall_effect, effect_label, outcome_label) {
  outcome_data <- effect_data_var |>
    dplyr::select(
      .data$group_id,
      .data$x_label_long,
      .data$y_mean_control,
      .data$y_mean_treatment
    ) |>
    tidyr::pivot_longer(
      cols = c(.data$y_mean_control, .data$y_mean_treatment),
      names_to = "arm",
      values_to = "mean_y"
    ) |>
    dplyr::mutate(
      arm = dplyr::recode(
        .data$arm,
        y_mean_control = "Control",
        y_mean_treatment = "Treatment"
      )
    )
  
  upper_plot <- ggplot2::ggplot(
    outcome_data,
    ggplot2::aes(x = .data$group_id, y = .data$mean_y, group = .data$arm)
  ) +
    ggplot2::geom_line(ggplot2::aes(linetype = .data$arm), linewidth = 0.4) +
    ggplot2::geom_point(ggplot2::aes(shape = .data$arm), size = 2) +
    ggplot2::scale_x_continuous(
      breaks = effect_data_var$group_id,
      labels = effect_data_var$x_label_long
    ) +
    ggplot2::labs(
      title = as.character(effect_data_var$x_var[[1]]),
      subtitle = "Arm-specific outcome summary",
      x = NULL,
      y = outcome_label,
      linetype = NULL,
      shape = NULL
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
  
  lower_plot <- ggplot2::ggplot(
    effect_data_var,
    ggplot2::aes(x = .data$group_id, y = .data$effect)
  ) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.3) +
    ggplot2::geom_hline(
      yintercept = overall_effect,
      linetype = 2,
      linewidth = 0.3
    ) +
    ggplot2::geom_linerange(
      ggplot2::aes(ymin = .data$ci_low, ymax = .data$ci_high),
      linewidth = 0.3
    ) +
    ggplot2::geom_line(linewidth = 0.4) +
    ggplot2::geom_point(ggplot2::aes(size = .data$n_total), alpha = 0.9) +
    ggplot2::scale_x_continuous(
      breaks = effect_data_var$group_id,
      labels = effect_data_var$x_label_long
    ) +
    ggplot2::scale_size_continuous(range = c(1.5, 4), guide = "none") +
    ggplot2::labs(
      subtitle = "Observed treatment effect by covariate group",
      x = NULL,
      y = effect_label
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
  
  upper_plot / lower_plot
}

print_plot_object <- function(plot_object) {
  print(plot_object)
}


render_task1_dataset_section <- function(csv_file, n_bins = params$n_bins, min_arm_n = 10L) {
  data <- readr::read_csv(csv_file, show_col_types = FALSE)
  
  dataset_id <- extract_dataset_id_from_path(csv_file)
  
  visual_object <- build_task1_visual_object(
    data = data,
    dataset_id = dataset_id,
    n_bins = n_bins,
    min_arm_n = min_arm_n
  )
  
  info_tbl <- tibble::tibble(
    dataset_id = visual_object$dataset_id,
    endpoint = visual_object$endpoint_type,
    n_obs = nrow(data),
    n_vars = length(visual_object$x_cols),
    n_bins = visual_object$n_bins_used,
    overall_effect = visual_object$overall_effect
  )
  
  priority_summary_tbl <- visual_object$dataset_priority_summary
  
  #priority_tbl <- priority_table_for_display(visual_object$priority_data)
  
  cat("\n\n")
  cat("## Dataset ", visual_object$dataset_id, "\n\n", sep = "")
  cat(
    "All variables are shown. Numeric prioritisation is used only to order and colour the plots for faster visual review.\n\n"
  )
  
  cat(knitr::kable(info_tbl, format = "pipe", digits = 3), sep = "\n")
  cat("\n\n")
  
  cat("### Dataset priority overview\n\n")
  cat(knitr::kable(priority_summary_tbl, format = "pipe", digits = 3), sep = "\n")
  cat("\n\n")
  
  # cat("### Variable priority table\n\n")
  # cat(knitr::kable(priority_tbl, format = "pipe", digits = 3), sep = "\n")
  # cat("\n\n")
  
  overview_pages <- plot_task1_overview_pages(
    effect_data = visual_object$effect_data,
    overall_effect = visual_object$overall_effect,
    effect_label = visual_object$effect_label,
    page_size = 9L,
    ncol = 3L
  )
  
  purrr::walk(overview_pages, print_plot_object)
}
