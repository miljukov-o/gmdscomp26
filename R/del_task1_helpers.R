effect_fun <- function(y, w) {
  mean(y[w == 1], na.rm = TRUE) - mean(y[w == 0], na.rm = TRUE)
}

is_binary_x <- function(x) {
  values <- sort(unique(x))
  length(values) == 2 && all(values %in% c(0, 1))
}

make_groups <- function(x, n_bins) {
  stopifnot(length(n_bins) == 1L, is.numeric(n_bins), n_bins >= 2)
  
  if (is_binary_x(x)) {
    return(factor(x, levels = c(0, 1)))
  }
  
  if (dplyr::n_distinct(x) <= 1L) {
    return(factor(rep("all", length(x)), levels = "all"))
  }
  
  breaks <- stats::quantile(
    x,
    probs = seq(0, 1, length.out = n_bins + 1L),
    na.rm = TRUE,
    names = FALSE,
    type = 2
  ) |>
    unique()
  
  cut(
    x,
    breaks = breaks,
    include.lowest = TRUE,
    ordered_result = TRUE
  )
}

summarise_dataset <- function(file, n_bins, top_k_vars) {
  dat <- readr::read_csv(file, show_col_types = FALSE)
  
  dataset_id <- basename(file) |>
    stringr::str_remove("\\.csv$") |>
    stringr::str_squish()
  
  x_vars <- names(dat) |>
    stringr::str_subset("^x\\d+$")
  
  global_te <- effect_fun(dat$Y, dat$W)
  
  var_summary <- purrr::map_dfr(x_vars, ~ summarise_var(., n_bins = n_bins)) |>
    dplyr::arrange(dplyr::desc(te_range))
  
  top_vars <- var_summary |>
    dplyr::slice_head(n = top_k_vars) |>
    dplyr::pull(variable)
  
  list(
    dataset_id = dataset_id,
    data = dat,
    global_te = global_te,
    var_summary = var_summary,
    top_vars = top_vars
  )
}

summarise_var <- function(var, n_bins) {
  groups <- make_groups(dat[[var]], n_bins = n_bins)
  
  plot_data <- tibble::tibble(
    variable = var,
    group = groups,
    Y = dat$Y,
    W = dat$W
  ) |>
    dplyr::group_by(variable, group) |>
    dplyr::summarise(
      n = dplyr::n(),
      local_te = effect_fun(Y, W),
      .groups = "drop"
    )
  
  tibble::tibble(
    variable = var,
    te_range = max(plot_data$local_te, na.rm = TRUE) - min(plot_data$local_te, na.rm = TRUE)
  )
}

generate_decision_tbl <- function(res) {
  tibble::tibble(
    dataset = res$dataset_id,
    global_te = res$global_te,
    top_var_1 = dplyr::first(res$top_vars),
    top_var_2 = dplyr::nth(res$top_vars, 2),
    top_var_3 = dplyr::nth(res$top_vars, 3),
    visual_decision = NA_character_,
    reason = NA_character_
  )
}

make_plot_data <- function(dat, variable, n_bins) {
  stopifnot(variable %in% names(dat))
  stopifnot(all(c("Y", "W") %in% names(dat)))
  
  groups <- make_groups(dat[[variable]], n_bins = n_bins)
  
  tibble::tibble(
    group = groups,
    Y = dat$Y,
    W = dat$W
  ) |>
    dplyr::group_by(group) |>
    dplyr::summarise(
      n = dplyr::n(),
      n_treat = sum(W == 1),
      n_control = sum(W == 0),
      mean_treat = if (sum(W == 1) > 0) mean(Y[W == 1], na.rm = TRUE) else NA_real_,
      mean_control = if (sum(W == 0) > 0) mean(Y[W == 0], na.rm = TRUE) else NA_real_,
      .groups = "drop"
    ) |>
    dplyr::mutate(
      local_te = mean_treat - mean_control
    )
}

plot_one_variable <- function(variable, res, n_bins) {
  plot_data <- make_plot_data(
    dat = res$data,
    variable = variable,
    n_bins = n_bins
  )
  
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = group, y = local_te, group = 1)
  ) +
    ggplot2::geom_hline(
      yintercept = res$global_te,
      linetype = 2
    ) +
    ggplot2::geom_line(na.rm = TRUE) +
    ggplot2::geom_point(size = 2, na.rm = TRUE) +
    ggplot2::geom_text(
      ggplot2::aes(label = paste0("n=", n)),
      vjust = -0.6,
      size = 3,
      na.rm = TRUE
    ) +
    ggplot2::labs(
      title = variable,
      subtitle = paste0("global TE = ", round(res$global_te, 3)),
      x = NULL,
      y = "local TE"
    ) +
    ggplot2::theme_minimal(base_size = 11)
}

plot_top_variables <- function(res, n_bins) {
  plots <- purrr::map(
    res$top_vars,
    plot_one_variable,
    res = res,
    n_bins = n_bins
  )
  
  patchwork::wrap_plots(plots, ncol = 1)
}

