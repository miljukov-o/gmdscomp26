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
    dplyr::transmute(
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
    dplyr::transmute(
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
    n = nrow(data),
    p = length(x_cols),
    binary_covariates = sum(binary_flags),
    continuous_covariates = sum(!binary_flags),
    n_bins_used = if (is.null(n_bins)) default_n_bins(nrow(data)) else as.integer(n_bins)
  )
}

build_task1_visual_object <- function(data, dataset_id = NA_integer_, n_bins = NULL) {
  endpoint_type <- detect_endpoint_type(data)
  x_cols <- get_x_cols(data)
  n_bins_used <- if (is.null(n_bins)) default_n_bins(nrow(data)) else as.integer(n_bins)
  
  effect_tables <- purrr::map(
    x_cols,
    summarise_variable_effect,
    data = data,
    endpoint_type = endpoint_type,
    n_bins = n_bins_used
  )
  
  effect_data <- effect_tables |>
    purrr::list_rbind() |>
    dplyr::mutate(
      x_var = factor(.data$x_var, levels = x_cols)
    )
  
  list(
    dataset_id = dataset_id,
    endpoint_type = endpoint_type,
    x_cols = x_cols,
    n_bins_used = n_bins_used,
    overall_effect = compute_overall_effect(data),
    effect_label = effect_scale_label(endpoint_type),
    outcome_label = outcome_scale_label(endpoint_type),
    effect_data = effect_data
  )
}

plot_single_effect_curve <- function(effect_data_var, overall_effect, effect_label) {
  title_text <- as.character(effect_data_var$x_var[[1]])
  subtitle_text <- paste(effect_data_var$var_type[[1]], "covariate")
  
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
      plot.title = ggplot2::element_text(face = "bold"),
      axis.text.x = ggplot2::element_text(size = 8)
    )
}

make_plot_for_variable <- function(x_var, effect_data, overall_effect, effect_label) {
  effect_data_var <- effect_data |>
    dplyr::filter(.data$x_var == x_var)
  
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
  x_vars <- levels(effect_data$x_var)
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


render_task1_dataset_section <- function(csv_file, n_bins = params$n_bins) {
  data <- readr::read_csv(csv_file, show_col_types = FALSE)
  
  dataset_id <- extract_dataset_id_from_path(csv_file)
  
  visual_object <- build_task1_visual_object(
    data = data,
    dataset_id = dataset_id,
    n_bins = n_bins
  )
  
  info_tbl <- tibble::tibble(
    dataset_id = visual_object$dataset_id,
    endpoint = visual_object$endpoint_type,
    n = nrow(data),
    p = length(visual_object$x_cols),
    n_bins_used = visual_object$n_bins_used,
    overall_effect = visual_object$overall_effect
  )
  
  cat("\n\n")
  cat("## Dataset ", visual_object$dataset_id, "\n\n", sep = "")
  cat(
    "Univariate visual inspection only. This section does not perform variable selection, subgroup assignment or an automatic dataset-level yes/no decision.\n\n"
  )
  
  cat(knitr::kable(info_tbl, format = "pipe", digits = 3), sep = "\n")
  cat("\n\n")
  
  overview_pages <- plot_task1_overview_pages(
    effect_data = visual_object$effect_data,
    overall_effect = visual_object$overall_effect,
    effect_label = visual_object$effect_label,
    page_size = 9L,
    ncol = 3L
  )
  
  purrr::walk(overview_pages, print_plot_object)
}
