priority_class_levels <- c("high", "medium", "low")
priority_colors <- c(
  high = "#b22222",
  medium = "#c77f00",
  low = "grey55"
)

get_csv_file_names <- function(data_dir) {
  fs::dir_ls(data_dir, recurse = TRUE, glob = "*.csv") |>
    tibble::tibble(path = _) |>
    dplyr::mutate(dataset_id = purrr::map_int(.data$path, extract_dataset_id_from_path)) |>
    dplyr::arrange(.data$dataset_id) |>
    dplyr::pull(.data$path)
}

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

get_covariate_type <- function(x_var, data) {
  dplyr::if_else(column_is_binary(x_var, data), "binary", "continuous")
}

get_endpoint_metadata <- function(endpoint_type) {
  if (identical(endpoint_type, "binary")) {
    return(list(
      effect_label = "Difference in response rates",
      outcome_label = "Response rate"
    ))
  }
  
  list(
    effect_label = "Difference in means",
    outcome_label = "Mean outcome"
  )
}

detect_endpoint_type <- function(data) {
  y_values <- sort(unique(data$Y))
  
  if (length(y_values) <= 2 && all(y_values %in% c(0, 1))) {
    return("binary")
  }
  
  "continuous"
}


is_binary_covariate <- function(x) {
  x_values <- sort(unique(stats::na.omit(x)))
  length(x_values) <= 2 && all(x_values %in% c(0, 1))
}

column_is_binary <- function(x_var, data) {
  is_binary_covariate(dplyr::pull(data, dplyr::all_of(x_var)))
}

default_n_bins <- function(n) {
  dplyr::case_when(
    n <= 150 ~ 4L,
    n <= 500 ~ 6L,
    TRUE ~ 8L
  )
}

resolve_n_bins <- function(data, n_bins = NULL) {
  if (!is.null(n_bins)) {
    return(as.integer(n_bins))
  }
  
  default_n_bins(nrow(data))
}

read_task1_dataset <- function(csv_file) {
  readr::read_csv(csv_file, show_col_types = FALSE)
}

compute_overall_effect <- function(data) {
  data |>
    dplyr::group_by(.data$W) |>
    dplyr::summarise(y_mean = mean(.data$Y), .groups = "drop") |>
    dplyr::mutate(arm = dplyr::if_else(.data$W == 1, "treatment", "control")) |>
    dplyr::select(-.data$W) |>
    tidyr::pivot_wider(
      names_from = .data$arm,
      values_from = .data$y_mean,
      names_prefix = "y_mean_"
    ) |>
    dplyr::mutate(overall_effect = .data$y_mean_treatment - .data$y_mean_control) |>
    dplyr::pull(.data$overall_effect)
}

make_binary_grouped_data <- function(data, x_var) {
  data |>
    dplyr::mutate(
      ptid = .data$ptid,
      Y = .data$Y,
      W = .data$W,
      x_value = as.numeric(.data[[x_var]]),
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
      x_value = .data[[x_var]],
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
  covariate_type <- get_covariate_type(x_var, data)
  
  if (identical(covariate_type, "binary")) {
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
    dplyr::mutate(arm = dplyr::if_else(.data$W == 1, "treatment", "control")) |>
    dplyr::select(-.data$W) |>
    tidyr::pivot_wider(
      names_from = .data$arm,
      values_from = c(.data$n, .data$y_mean, .data$y_var),
      names_sep = "_",
      values_fill = list(n = 0, y_mean = NA_real_, y_var = NA_real_)
    ) |>
    dplyr::arrange(.data$group_id)
}

compute_binary_effect_se <- function(arm_summary) {
  sqrt(
    (arm_summary$y_mean_treatment * (1 - arm_summary$y_mean_treatment) / arm_summary$n_treatment) +
      (arm_summary$y_mean_control * (1 - arm_summary$y_mean_control) / arm_summary$n_control)
  )
}

compute_continuous_effect_se <- function(arm_summary) {
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

compute_effect_se <- function(arm_summary, endpoint_type) {
  if (identical(endpoint_type, "binary")) {
    return(compute_binary_effect_se(arm_summary))
  }
  
  compute_continuous_effect_se(arm_summary)
}

summarise_variable_effect <- function(x_var, data, endpoint_type, n_bins) {
  grouped_data <- make_grouped_data(
    data = data,
    x_var = x_var,
    n_bins = n_bins
  )
  
  arm_summary <- summarise_arm_level(grouped_data)
  se <- compute_effect_se(arm_summary, endpoint_type)
  var_type <- get_covariate_type(x_var, data)
  
  arm_summary |>
    dplyr::mutate(
      se = se,
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
  ok <- !is.na(x)
  n_ok <- sum(ok)
  
  if (n_ok == 0L) {
    return(rep(NA_real_, length(x)))
  }
  
  if (n_ok == 1L) {
    out <- rep(NA_real_, length(x))
    out[ok] <- 0.5
    return(out)
  }
  
  out <- rep(NA_real_, length(x))
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
  
  abs(suppressWarnings(stats::cor(group_id[keep], effect[keep], method = "spearman")))
}

has_qualitative_signal <- function(effect) {
  if (all(is.na(effect))) {
    return(0L)
  }
  
  as.integer(min(effect, na.rm = TRUE) < 0 && max(effect, na.rm = TRUE) > 0)
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
    qualitative_signal = has_qualitative_signal(effect_data_var$effect),
    min_arm_n_observed = min(arm_sizes, na.rm = TRUE),
    unstable_fraction = mean(arm_sizes < min_arm_n, na.rm = TRUE),
    n_groups = dplyr::n_distinct(effect_data_var$group_id)
  )
}

summarise_priority_from_nested_data <- function(data, overall_effect, min_arm_n) {
  summarise_variable_priority(
    effect_data_var = data,
    overall_effect = overall_effect,
    min_arm_n = min_arm_n
  )
}

apply_priority_scores <- function(priority_data) {
  priority_data |>
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
      priority_class = factor(.data$priority_class, levels = priority_class_levels),
      priority_color = unname(priority_colors[as.character(.data$priority_class)])
    ) |>
    dplyr::arrange(.data$priority_class, dplyr::desc(.data$priority_score), .data$x_var) |>
    dplyr::mutate(priority_rank = dplyr::row_number())
}

build_variable_priority_table <- function(effect_data, overall_effect, min_arm_n = 10L) {
  effect_data |>
    tidyr::nest(data = -c(.data$x_var, .data$var_type)) |>
    dplyr::mutate(
      priority_metrics = purrr::map(
        .x = .data$data,
        .f = summarise_variable_priority,
        overall_effect = overall_effect,
        min_arm_n = min_arm_n
      )
    ) |>
    dplyr::select(-.data$data) |>
    tidyr::unnest(.data$priority_metrics) |>
    apply_priority_scores()
}

priority_table_for_display <- function(priority_data) {
  priority_data |>
    dplyr::select(
      .data$x_var,
      .data$var_type,
      .data$weighted_abs_deviation,
      .data$effect_range,
      .data$trend_strength,
      .data$qualitative_signal,
      .data$min_arm_n_observed,
      .data$unstable_fraction,
      .data$priority_score,
      .data$priority_class
    )
}

complete_priority_counts <- function(priority_data) {
  tibble::tibble(priority_class = factor(priority_class_levels, levels = priority_class_levels)) |>
    dplyr::left_join(
      priority_data |>
        dplyr::count(.data$priority_class, name = "n_vars"),
      by = "priority_class"
    ) |>
    dplyr::mutate(n_vars = tidyr::replace_na(.data$n_vars, 0L))
}

collapse_top_priority_vars <- function(priority_data, n_top = 5L) {
  priority_data |>
    dplyr::arrange(.data$priority_class, dplyr::desc(.data$priority_score), .data$x_var) |>
    dplyr::slice_head(n = n_top) |>
    dplyr::pull(.data$x_var) |>
    as.character() |>
    paste(collapse = ", ")
}

summarise_dataset_priority <- function(priority_data) {
  counts <- complete_priority_counts(priority_data)
  
  tibble::tibble(
    n_high = counts |>
      dplyr::filter(.data$priority_class == "high") |>
      dplyr::pull(.data$n_vars),
    n_medium = counts |>
      dplyr::filter(.data$priority_class == "medium") |>
      dplyr::pull(.data$n_vars),
    n_low = counts |>
      dplyr::filter(.data$priority_class == "low") |>
      dplyr::pull(.data$n_vars),
    prop_high = .data$n_high / nrow(priority_data),
    prop_medium_or_higher = (.data$n_high + .data$n_medium) / nrow(priority_data),
    max_priority_score = max(priority_data$priority_score, na.rm = TRUE),
    mean_priority_score = mean(priority_data$priority_score, na.rm = TRUE),
    top_priority_vars = collapse_top_priority_vars(priority_data)
  )
}

derive_dataset_priority_hint <- function(n_high, n_medium) {
  dplyr::case_when(
    n_high >= 3 ~ "multiple_high_priority_signals",
    n_high >= 1 && n_medium >= 2 ~ "high_plus_additional_support",
    n_high >= 1 ~ "single_high_priority_signal",
    n_medium >= 3 ~ "several_medium_priority_signals",
    n_medium >= 1 ~ "weak_to_moderate_signal_only",
    TRUE ~ "mostly_low_priority_patterns"
  )
}

build_dataset_priority_summary <- function(priority_data) {
  summarise_dataset_priority(priority_data) |>
    dplyr::mutate(
      dataset_priority_hint = derive_dataset_priority_hint(
        n_high = .data$n_high,
        n_medium = .data$n_medium
      )
    )
}



join_priority_metrics <- function(effect_data, priority_data, x_cols) {
  effect_data |>
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
    dplyr::mutate(x_var = factor(.data$x_var, levels = x_cols))
}

build_task1_visual_object <- function(data, dataset_id = NA_integer_, n_bins = NULL, min_arm_n = 10L) {
  endpoint_type <- detect_endpoint_type(data)
  x_cols <- get_x_cols(data)
  n_bins_used <- resolve_n_bins(data, n_bins)
  endpoint_metadata <- get_endpoint_metadata(endpoint_type)
  overall_effect <- compute_overall_effect(data)
  
  effect_data <- purrr::map(
    .x = x_cols,
    .f = summarise_variable_effect,
    data = data,
    endpoint_type = endpoint_type,
    n_bins = n_bins_used
  ) |>
    purrr::list_rbind()
  
  priority_data <- build_variable_priority_table(
    effect_data = effect_data,
    overall_effect = overall_effect,
    min_arm_n = min_arm_n
  ) |>
    dplyr::mutate(x_var = factor(.data$x_var, levels = x_cols))
  
  list(
    dataset_id = dataset_id,
    endpoint_type = endpoint_type,
    x_cols = x_cols,
    n_bins_used = n_bins_used,
    overall_effect = overall_effect,
    effect_label = endpoint_metadata$effect_label,
    outcome_label = endpoint_metadata$outcome_label,
    effect_data = join_priority_metrics(effect_data, priority_data, x_cols),
    priority_data = priority_data,
    dataset_priority_summary = build_dataset_priority_summary(priority_data)
  )
}

extract_plot_metadata <- function(effect_data_var) {
  first_row <- effect_data_var |>
    dplyr::slice_head(n = 1)
  
  list(
    title_text = as.character(first_row$x_var[[1]]),
    subtitle_text = paste0(
      first_row$var_type[[1]],
      " covariate | priority: ",
      as.character(first_row$priority_class[[1]])
    ),
    priority_color = first_row$priority_color[[1]]
  )
}

plot_single_effect_curve <- function(effect_data_var, overall_effect, effect_label) {
  plot_metadata <- extract_plot_metadata(effect_data_var)
  
  ggplot2::ggplot(effect_data_var, ggplot2::aes(x = .data$group_id, y = .data$effect)) +
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
      labels = effect_data_var$x_label_short
    ) +
    ggplot2::scale_size_continuous(range = c(1.5, 4), guide = "none") +
    ggplot2::labs(
      title = plot_metadata$title_text,
      subtitle = plot_metadata$subtitle_text,
      x = NULL,
      y = effect_label
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", colour = plot_metadata$priority_color),
      plot.subtitle = ggplot2::element_text(colour = plot_metadata$priority_color),
      plot.background = ggplot2::element_rect(
        fill = NA,
        colour = plot_metadata$priority_color,
        linewidth = 0.9
      ),
      axis.text.x = ggplot2::element_text(size = 8)
    )
}

make_plot_for_variable <- function(var_name, effect_data, overall_effect, effect_label) {
  effect_data |>
    dplyr::filter(.data$x_var == .env$var_name) |>
    plot_single_effect_curve(
      overall_effect = overall_effect,
      effect_label = effect_label
    )
}

split_x_vars_into_pages <- function(x_vars, page_size) {
  split(x_vars, ceiling(seq_along(x_vars) / page_size))
}

make_overview_page_plot <- function(x_vars, effect_data, overall_effect, effect_label, ncol) {
  purrr::map(
    .x = x_vars,
    .f = make_plot_for_variable,
    effect_data = effect_data,
    overall_effect = overall_effect,
    effect_label = effect_label
  ) |>
    patchwork::wrap_plots(ncol = ncol)
}

get_priority_ordered_x_vars <- function(effect_data) {
  effect_data |>
    dplyr::distinct(.data$x_var, .data$priority_class, .data$priority_score) |>
    dplyr::mutate(
      priority_class = factor(.data$priority_class, levels = priority_class_levels)
    ) |>
    dplyr::arrange(.data$priority_class, dplyr::desc(.data$priority_score), .data$x_var) |>
    dplyr::pull(.data$x_var) |>
    as.character()
}

plot_task1_overview_pages <- function(effect_data, overall_effect, effect_label, page_size = 9L, ncol = 3L) {
  get_priority_ordered_x_vars(effect_data) |>
    split_x_vars_into_pages(page_size = page_size) |>
    purrr::map(
      .f = make_overview_page_plot,
      effect_data = effect_data,
      overall_effect = overall_effect,
      effect_label = effect_label,
      ncol = ncol
    )
}

reshape_detail_outcome_data <- function(effect_data_var) {
  effect_data_var |>
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
}

plot_single_detail_curve <- function(effect_data_var, overall_effect, effect_label, outcome_label) {
  outcome_data <- reshape_detail_outcome_data(effect_data_var)
  title_text <- effect_data_var |>
    dplyr::slice_head(n = 1) |>
    dplyr::pull(.data$x_var) |>
    as.character()
  
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
      title = title_text,
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



render_task1_dataset_section <- function(data) {
  dataset_id <- extract_dataset_id_from_path(csv_file)
  
  visual_object <- build_task1_visual_object(
    data = data,
    dataset_id = dataset_id,
    n_bins = n_bins,
    min_arm_n = min_arm_n
  )
  
  cat("\n\n")
  cat("## Dataset ", visual_object$dataset_id, "\n\n", sep = "")
  cat(
    "All variables are shown. Numeric prioritisation is used only to order and colour the plots for faster visual review.\n\n"
  )
  
  build_dataset_info_table(visual_object, data) |>
    write_kable_output(digits = 3)
  
  cat("### Dataset priority overview\n\n")
  visual_object$dataset_priority_summary |>
    write_kable_output(digits = 3)
  
  plot_task1_overview_pages(
    effect_data = visual_object$effect_data,
    overall_effect = visual_object$overall_effect,
    effect_label = visual_object$effect_label,
    page_size = 9L,
    ncol = 3L
  ) |>
    purrr::walk(print_plot_object)
}

merge_task1_review_sheet <- function(review_template, review_path) {
  if (!fs::file_exists(review_path)) {
    return(review_template)
  }
  
  existing_review <- readr::read_csv2(review_path, show_col_types = FALSE)
  
  shared_cols <- intersect(
    setdiff(names(review_template), "dataset_id"),
    names(existing_review)
  )
  
  extra_existing_cols <- setdiff(
    names(existing_review),
    names(review_template)
  )
  
  merged_review <- review_template |>
    dplyr::left_join(
      existing_review,
      by = "dataset_id",
      suffix = c("_new", "_old")
    )
  
  for (col_name in shared_cols) {
    merged_review <- coalesce_review_column(merged_review, col_name)
  }
  
  merged_review |>
    dplyr::select(
      dplyr::all_of(names(review_template)),
      dplyr::all_of(extra_existing_cols)
    ) |>
    dplyr::arrange(.data$dataset_id)
}

init_task1_review_sheet <- function(dataset_manifest, review_path) {
  review_sheet <- build_review_template(dataset_manifest) |>
    merge_task1_review_sheet(review_path = review_path)
  
  readr::write_csv2(review_sheet, review_path, na = "")
  review_sheet
}

build_review_template <- function(dataset_manifest) {
  dataset_manifest |>
    dplyr::select(
      .data$dataset_id,
      .data$endpoint,
      .data$n_obs,
      .data$n_vars,
      .data$n_bins
    ) |>
    dplyr::mutate(
      visual_teh_assessment = NA_character_,
      reviewer_note = NA_character_
    )
}

review_sheet_base_cols <- function() {
  c(
    "dataset_id",
    "endpoint",
    "n_obs",
    "n_vars",
    "n_bins",
    "visual_teh_assessment",
    "reviewer_note"
  )
}

coalesce_review_column <- function(data, col_name) {
  new_col <- rlang::sym(paste0(col_name, "_new"))
  old_col <- rlang::sym(paste0(col_name, "_old"))
  out_col <- rlang::sym(col_name)
  
  data |>
    dplyr::mutate(
      !!out_col := dplyr::coalesce(!!old_col, !!new_col)
    )
}



build_task1_dataset_record <- function(csv_file, n_bins = NULL, min_arm_n = 10L) {
  data <- read_task1_dataset(csv_file)
  dataset_id <- extract_dataset_id_from_path(csv_file)
  
  visual_object <- build_task1_visual_object(
    data = data,
    dataset_id = dataset_id,
    n_bins = n_bins,
    min_arm_n = min_arm_n
  )
  
  manifest_row <- tibble::tibble(
    dataset_id = dataset_id,
    endpoint = visual_object$endpoint_type,
    n_obs = nrow(data),
    n_vars = length(visual_object$x_cols),
    binary_covariates = sum(
      purrr::map_chr(
        .x = visual_object$x_cols,
        .f = get_covariate_type,
        data = data
      ) == "binary"
    ),
    continuous_covariates = sum(
      purrr::map_chr(
        .x = visual_object$x_cols,
        .f = get_covariate_type,
        data = data
      ) == "continuous"
    ),
    n_bins = visual_object$n_bins_used
  )
  
  triage_row <- build_dataset_priority_summary(visual_object$priority_data) |>
    dplyr::mutate(
      dataset_id = dataset_id,
      endpoint = visual_object$endpoint_type,
      n_obs = nrow(data),
      n_vars = length(visual_object$x_cols),
      n_bins = visual_object$n_bins_used,
      overall_effect = visual_object$overall_effect,
      .before = 1
    )
  
  list(
    csv_file = csv_file,
    dataset_id = dataset_id,
    data = data,
    visual_object = visual_object,
    manifest_row = manifest_row,
    triage_row = triage_row
  )
}

build_task1_dataset_records <- function(csv_files, n_bins = NULL, min_arm_n = 10L) {
  purrr::map(
    .x = csv_files,
    .f = build_task1_dataset_record,
    n_bins = n_bins,
    min_arm_n = min_arm_n
  )
}

build_dataset_manifest <- function(dataset_records) {
  purrr::map(
    .x = dataset_records,
    .f = pluck_manifest_row
  ) |>
    purrr::list_rbind() |>
    dplyr::arrange(.data$dataset_id)
}

pluck_manifest_row <- function(dataset_record) {
  dataset_record$manifest_row
}

build_task1_dataset_triage_table <- function(dataset_records) {
  purrr::map(
    .x = dataset_records,
    .f = pluck_triage_row
  ) |>
    purrr::list_rbind() |>
    dplyr::arrange(
      dplyr::desc(.data$n_high),
      dplyr::desc(.data$n_medium),
      dplyr::desc(.data$max_priority_score),
      .data$dataset_id
    )
}

pluck_triage_row <- function(dataset_record) {
  dataset_record$triage_row
}

build_dataset_info_table <- function(visual_object, data) {
  tibble::tibble(
    dataset_id = visual_object$dataset_id,
    endpoint = visual_object$endpoint_type,
    n_obs = nrow(data),
    n_vars = length(visual_object$x_cols),
    n_bins = visual_object$n_bins_used,
    overall_effect = visual_object$overall_effect
  )
}

write_kable_output <- function(x, digits = 3) {
  cat(knitr::kable(x, format = "pipe", digits = digits), sep = "\n")
  cat("\n\n")
}

render_task1_dataset_record <- function(dataset_record) {
  data <- dataset_record$data
  visual_object <- dataset_record$visual_object
  
  cat("\n\n")
  cat("## Dataset ", visual_object$dataset_id, "\n\n", sep = "")
  cat(
    "All variables are shown. Numeric prioritisation is used only to order and colour the plots for faster visual review.\n\n"
  )
  
  build_dataset_info_table(visual_object, data) |>
    write_kable_output(digits = 3)
  
  cat("### Dataset priority overview\n\n")
  visual_object$dataset_priority_summary |>
    write_kable_output(digits = 3)
  
  plot_task1_overview_pages(
    effect_data = visual_object$effect_data,
    overall_effect = visual_object$overall_effect,
    effect_label = visual_object$effect_label,
    page_size = 9L,
    ncol = 3L
  ) |>
    purrr::walk(.f = print_plot_object)
}

render_task1_dataset_sections <- function(dataset_records) {
  purrr::walk(
    .x = dataset_records,
    .f = render_task1_dataset_record
  )
}
