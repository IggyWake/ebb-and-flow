library(tidyverse)
library(arrow)
library(janitor)
library(ggrepel)
library(data.table)
library(ranger)
library(vip)
library(rpart)
library(rpart.plot)
library(purrr)
library(gt)
library(webshot2)
library(rlang)

# ==== FUNCTIONS ====

# 1. Cleaning
exp_cleaning <- function(data) {
  data |> 
    clean_names() |> 
    rename(any_of(c(
      benefit = "b",
      cost = "c",
      synergy = "s",
      trait_mut = "u",
      set_mut = "v",
      set_std = "std",
      dir_exp = "x",
      baseline = "w0",
      set_mean = "k",
      total_sets = "m",
      population = "n",
      coop_freq = "count_individuals_with_strategy_cooperator_n",
      assortment_c = "ifelse_value_any_individuals_with_strategy_cooperator_mean_local_c_freq_of_individuals_with_strategy_cooperator_0",
      assortment_d = "ifelse_value_any_individuals_with_strategy_defector_mean_local_c_freq_of_individuals_with_strategy_defector_0",
      assortment_std = "ifelse_value_count_individuals_with_strategy_cooperator_1_standard_deviation_local_c_freq_of_individuals_with_strategy_cooperator_0",
      degrees_c = "degree_of_individuals_with_strategy_cooperator",
      degrees_d = "degree_of_individuals_with_strategy_defector",
      sets_c = "count_link_neighbors_of_individuals_with_strategy_cooperator",
      sets_d = "count_link_neighbors_of_individuals_with_strategy_defector"
    ))) |>
    mutate(
      # Extract numbers from the NetLogo string format to create list columns
      degrees_c = map(degrees_c, ~ as.numeric(str_extract_all(.x, "-?\\d+\\.?\\d*")[[1]])),
      degrees_d = map(degrees_d, ~ as.numeric(str_extract_all(.x, "-?\\d+\\.?\\d*")[[1]])),
      sets_c = map(sets_c, ~ as.numeric(str_extract_all(.x, "-?\\d+\\.?\\d*")[[1]])),
      sets_d = map(sets_d, ~ as.numeric(str_extract_all(.x, "-?\\d+\\.?\\d*")[[1]])),
      # Calculate metrics in new columns
      mean_degrees_c = map_dbl(degrees_c, mean, na.rm = TRUE),
      mean_degrees_d = map_dbl(degrees_d, mean, na.rm = TRUE),
      std_degrees_c = map_dbl(degrees_c, sd, na.rm = TRUE),
      std_degrees_d = map_dbl(degrees_d, sd, na.rm = TRUE),
      mean_sets_c = map_dbl(sets_c, mean, na.rm = TRUE),
      mean_sets_d = map_dbl(sets_d, mean, na.rm = TRUE),
      std_sets_c = map_dbl(sets_c, sd, na.rm = TRUE),
      std_sets_d = map_dbl(sets_d, sd, na.rm = TRUE),
      run_number = factor(run_number),
      variance_ratio = case_when(
        coop_freq == 0 | coop_freq == 1 ~ 0,
        TRUE ~ assortment_c - assortment_d
      )
    ) |> 
    select(-any_of(c(
      "ifelse_value_count_individuals_with_strategy_defector_1_standard_deviation_local_c_freq_of_individuals_with_strategy_defector_0"
    )))
}

# 2. Phace space plot
phase_space <- function(
    data, x_var, y_var, facet_row, facet_col, fill_var, size_var, 
    plot_title = "Cooperation Success Phase Space",
    plot_subtitle = "Tile color indicates average success. Circle size indicates volatility.",
    x_label = "Synergy (s)",
    y_label = "Benefit (b)",
    fill_label = "Coop Success",
    size_label = "Volatility"
) {
  p <- ggplot(data, aes(x = factor({{ x_var }}), y = factor({{ y_var }}))) +
    geom_tile(aes(fill = {{ fill_var }}), color = "white", linewidth = 0.5) +
    geom_point(
      aes(size = {{ size_var }}, color = {{ fill_var }}), 
      show.legend = c(size = TRUE, color = FALSE)
    ) +
    scale_fill_viridis_c(option = "mako", name = fill_label, limits = c(0, 1)) +
    scale_color_viridis_c(option = "mako", direction = -1, limits = c(0, 1), guide = "none") +
    scale_size_continuous(name = size_label, range = c(0, 6)) +
    labs(
      title = plot_title,
      subtitle = plot_subtitle,
      x = x_label,
      y = y_label
    ) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(size = 10)
    )
  
  # Conditionally add 1D or 2D faceting
  if (!missing(facet_row) && !missing(facet_col)) {
    p <- p + facet_grid(rows = vars({{ facet_row }}), cols = vars({{ facet_col }}), labeller = label_both)
  } else if (!missing(facet_col)) {
    p <- p + facet_wrap(vars({{ facet_col }}), labeller = label_both)
  } else if (!missing(facet_row)) {
    p <- p + facet_wrap(vars({{ facet_row }}), labeller = label_both, dir = "v")
  }
  
  return(p)
}

# 3. Run summary
run_summary <- function(data) {
  data |>
    arrange(run_number, step) |>
    group_by(run_number) |>
    summarise(
      # First time hitting the thresholds (returns NA if it never hit)
      step_max = if_else(any(coop_freq >= 0.95, na.rm = TRUE), step[which.max(coop_freq >= 0.95)], NA_real_),
      step_min = if_else(any(coop_freq <= 0.05, na.rm = TRUE), step[which.max(coop_freq <= 0.05)], NA_real_),
      
      # Resolution speed: The earliest step it hit either extreme
      resolution_speed = if_else(
        is.na(step_max) & is.na(step_min), 
        NA_real_, 
        suppressWarnings(min(c(step_max, step_min), na.rm = TRUE))
      ),
      
      # Success indicator: 1 if max was hit first, 0 if min was hit first
      # if none were hit, takes the mean coop_freq of the last 1000 steps
      coop_success = case_when(
        is.na(resolution_speed) ~ mean(tail(coop_freq, 10), na.rm = TRUE),
        resolution_speed == step_max ~ 1,
        resolution_speed == step_min ~ 0
      ),
      
      # Calculate mean only for the steps that occurred AFTER the resolution_speed
      mean_after_bounce = if_else(
        is.na(resolution_speed) | !any(step > resolution_speed),
        NA_real_,
        mean(coop_freq[step > resolution_speed], na.rm = TRUE)
      ),
      
      bounce_back = case_when(
        is.na(resolution_speed) ~ NA_real_,
        coop_success == 1 ~ 1 - mean_after_bounce,
        coop_success == 0 ~ mean_after_bounce - 0
      ),
      
      # Keeping overall mean and std for runs that never resolve
      mean_coop = mean(coop_freq, na.rm = TRUE),
      std_coop = sd(coop_freq, na.rm = TRUE),
      
      .groups = 'drop'
    ) |> 
    select(-step_max, -step_min) |> 
    left_join(
      data |>
        group_by(run_number) |>
        slice(1) |> # Keep only the first step of each run
        ungroup() |>
        # Select the key and the parameter columns (adjust column names if needed)
        select(run_number, benefit:population), 
      by = "run_number"
    ) |>
    mutate(bounce_back = replace_na(bounce_back, 0),
           K_M_ratio = set_mean/total_sets,
           b_c_ratio = benefit/cost)
}

# 4. VIP plot
vip_plot <- function(rf_model, 
                               plot_title = "Random Forest: Parameter Importance",
                               plot_subtitle = "Ranking the structural and economic drivers of cooperation success") {
  
  vip(rf_model, geom = "col", mapping = aes(fill = Importance)) +
    scale_fill_viridis_c(option = "mako", begin = 0.3, end = 0.8, guide = "none") +
    labs(
      title = plot_title,
      subtitle = plot_subtitle,
      y = "Importance (Permutation Error Increase)",
      x = "Parameters"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.background = element_rect(fill = "#F4F5F7", color = NA), 
      panel.background = element_rect(fill = "#F4F5F7", color = NA),
      panel.grid.major.x = element_line(color = "#FFFFFF", linewidth = 0.8),
      panel.grid.major.y = element_blank(), 
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", size = 18, color = "#2D3748", margin = margin(b = 6)),
      plot.subtitle = element_text(color = "#4A5568", size = 12, margin = margin(b = 20)),
      axis.title = element_text(face = "bold", color = "#4A5568", margin = margin(t = 12, r = 12)),
      axis.text.y = element_text(color = "#2D3748", face = "bold", size = 12),
      axis.text.x = element_text(color = "#718096"),
      plot.margin = margin(t = 20, r = 20, b = 20, l = 20)
    )
}

# 5. Table
nice_table <- function(data, title, subtitle, footnote, x_label, y_label) {
  
  # Extract the actual column names from the data dynamically
  col_names <- colnames(data)
  
  # Build a dynamic list of labels mapping to the current column names
  label_mapping <- list(
    md(x_label),
    md(y_label),
    "Peak Success Rate"
  )
  names(label_mapping) <- c(col_names[1], col_names[2], "max_predicted_outcome")
  
  data |>
    gt() |>
    tab_header(
      title = md(title),
      subtitle = subtitle
    ) |>
    
    # Apply the dynamic labels
    cols_label(
      .list = label_mapping
    ) |>
    
    # Format the first two columns dynamically by their index position (1 and 2)
    fmt_number(
      columns = 1:2,
      decimals = 3
    ) |>
    
    # Format the 3rd column (which is always named max_predicted_outcome by our previous function)
    fmt_percent(
      columns = max_predicted_outcome, 
      decimals = 1
    ) |>
    cols_align(
      align = "center",
      columns = everything()
    ) |>
    tab_options(
      table.border.top.color = "black",
      table.border.top.width = px(2),
      table.border.bottom.color = "black",
      table.border.bottom.width = px(2),
      heading.border.bottom.color = "black",
      column_labels.border.top.color = "black",
      column_labels.border.top.width = px(2),
      column_labels.border.bottom.color = "black",
      column_labels.border.bottom.width = px(1),
      table.width = pct(60),
      data_row.padding = px(6)
    ) |>
    tab_footnote(
      footnote = footnote,
      locations = cells_title(groups = "title")
    )
}

# 6. Extract LOESS peaks
extract_loess_peaks <- function(data, find_peak, nest_by, outcome) {
  # Extract variable names as strings for the formula and base R functions
  peak_name <- rlang::as_name(rlang::enquo(find_peak))
  outcome_name <- rlang::as_name(rlang::enquo(outcome))
  
  # Build the mathematical formula dynamically
  model_formula <- as.formula(paste(outcome_name, "~", peak_name))
  
  data |>
    select({{ nest_by }}, {{ find_peak }}, {{ outcome }}) |>
    drop_na() |>
    nest(.by = {{ nest_by }}) |>
    mutate(
      model = map(data, ~ loess(model_formula, data = .x, span = 0.75)),
      
      # Use dynamic injection (!! :=) to name the column inside tibble()
      grid = map(data, ~ tibble(
        !!peak_name := seq(min(.x[[peak_name]]), max(.x[[peak_name]]), length.out = 10000)
      )),
      
      predictions = map2(model, grid, ~ {
        pred_grid <- .y
        pred_grid$pred_outcome <- predict(.x, newdata = pred_grid)
        return(pred_grid)
      }),
      
      peak = map(predictions, ~ .x |> arrange(desc(pred_outcome)) |> slice(1))
    ) |>
    unnest(peak) |>
    select(
      {{ nest_by }},
      # Dynamically rename the peak column for the final output
      !!paste0("optimal_", peak_name) := {{ find_peak }}, 
      max_predicted_outcome = pred_outcome
    ) |>
    arrange({{ nest_by }})
}

# 7. Plot LOESS peaks
plot_optimal_peak <- function(data, title, subtitle, x_label, y_label) {
  
  # Dynamically grab the column names from the data
  cols <- colnames(data)
  x_var <- cols[1]
  y_var <- cols[2]
  color_var <- cols[3] # Always "max_predicted_outcome"
  
  ggplot(data, aes(x = .data[[x_var]], y = .data[[y_var]])) +
    geom_line(color = "#A0AEC0", linewidth = 1.2, linetype = "dashed") +
    # Map color dynamically to the 3rd column
    geom_point(aes(color = .data[[color_var]]), size = 5) + 
    scale_color_viridis_c(
      option = "mako", 
      begin = 0.3, 
      end = 0.8, 
      name = "Peak Success Rate" # Generalized legend name
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      x = x_label,
      y = y_label
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.background = element_rect(fill = "#F4F5F7", color = NA), 
      panel.background = element_rect(fill = "#F4F5F7", color = NA),
      panel.grid.major = element_line(color = "#FFFFFF", linewidth = 0.8),
      panel.grid.minor = element_line(color = "#FFFFFF", linewidth = 0.4),
      plot.title = element_text(face = "bold", size = 18, color = "#2D3748", margin = margin(b = 6)),
      plot.subtitle = element_text(color = "#4A5568", size = 12, margin = margin(b = 20)),
      axis.title = element_text(face = "bold", color = "#4A5568", margin = margin(t = 12, r = 12)),
      axis.text = element_text(color = "#718096"),
      legend.position = "right",
      legend.title = element_text(face = "bold", color = "#2D3748", size = 11),
      legend.text = element_text(color = "#4A5568"),
      legend.background = element_rect(fill = NA, color = NA),
      plot.margin = margin(t = 20, r = 20, b = 20, l = 20)
    )
}

# 8. 2D LOESS grid
generate_loess_grid <- function(data, x_var, y_var, outcome, span = 0.75, res = 100) {
  x_name <- as_name(enquo(x_var))
  y_name <- as_name(enquo(y_var))
  outcome_name <- as_name(enquo(outcome))
  
  # Dynamically build the formula: outcome ~ x_var * y_var
  model_formula <- as.formula(paste(outcome_name, "~", x_name, "*", y_name))
  loess_2d <- loess(model_formula, data = data, span = span)
  
  # Create the sequences and grid
  x_seq <- seq(min(data[[x_name]], na.rm = TRUE), max(data[[x_name]], na.rm = TRUE), length.out = res)
  y_seq <- seq(min(data[[y_name]], na.rm = TRUE), max(data[[y_name]], na.rm = TRUE), length.out = res)
  
  grid_2d <- expand.grid(x_seq, y_seq)
  colnames(grid_2d) <- c(x_name, y_name)
  
  # Predict and forcefully convert to 1D numeric vector
  grid_2d$pred_outcome <- as.numeric(predict(loess_2d, newdata = grid_2d))
  
  return(grid_2d)
}

# 9. 2D LOESS plot

plot_2d_frontier <- function(grid_data, x_var, y_var, outcome_var, exchange_rate, title, subtitle, x_label, y_label, break_level = 0.9) {
  x_name <- as_name(enquo(x_var))
  y_name <- as_name(enquo(y_var))
  outcome_name <- as_name(enquo(outcome_var))
  
  # Dynamically append the exchange rate to the subtitle
  dynamic_subtitle <- paste(subtitle, "| Required exchange rate =", round(abs(exchange_rate), 3))
  
  ggplot(grid_data, aes(x = .data[[x_name]], y = .data[[y_name]])) +
    geom_raster(aes(fill = .data[[outcome_name]]), interpolate = TRUE) +
    geom_contour(aes(z = .data[[outcome_name]]), breaks = break_level, color = "white", linewidth = 1.2) +
    scale_fill_viridis_c(option = "mako", direction = 1, name = "Predicted Success Rate") +
    labs(
      title = title,
      subtitle = dynamic_subtitle,
      x = x_label,
      y = y_label
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.background = element_rect(fill = "#F4F5F7", color = NA),
      panel.background = element_rect(fill = "#F4F5F7", color = NA),
      panel.grid.major = element_line(color = "#FFFFFF", linewidth = 0.8),
      panel.grid.minor = element_line(color = "#FFFFFF", linewidth = 0.4),
      plot.title = element_text(face = "bold", size = 18, color = "#2D3748", margin = margin(b = 6)),
      plot.subtitle = element_text(color = "#4A5568", size = 12, margin = margin(b = 20)),
      axis.title = element_text(face = "bold", color = "#4A5568", margin = margin(t = 12, r = 12)),
      axis.text = element_text(color = "#718096"),
      legend.position = "right",
      legend.title = element_text(face = "bold", color = "#2D3748", size = 11),
      legend.text = element_text(color = "#4A5568"),
      legend.background = element_rect(fill = NA, color = NA),
      plot.margin = margin(t = 20, r = 20, b = 20, l = 20)
    )
}

# Example usage to plot:
# plot_2d_frontier(
#   grid_data = grid_2d, 
#   x_var = synergy, 
#   y_var = benefit, 
#   outcome_var = pred_outcome, 
#   exchange_rate = exchange_rate,
#   title = "Synergy vs Benefit: Evolutionary Frontier",
#   subtitle = "White line: 90% success threshold",
#   x_label = "Synergy (s)",
#   y_label = "Benefit (b)"
# )

# Example usage to create the grid:
# grid_2d <- generate_loess_grid(svb_x0, synergy, benefit, coop_success)

# ====
# ==== TOOLS ====

# SLICE
# svb_x0 |> 
#   select(-c(degrees_c, degrees_d, sets_c, sets_d)) |> 
#   head(10) |> write.csv(file = "slice_svb_x0.csv")

# FILTERING TOOL
# svb_x0 |> filter(benefit == 0.2, synergy == 0.66, set_mut == 0.5) |> view()

# COUNT UNIQUE RUNS
# s_0 |>
#   group_by(set_std, set_mut, dir_exp) |>
#   summarise(number_of_runs = n_distinct(run_number), .groups = 'drop') |> 
#   print()

# CHECK PARAMETER VALUES (NEEDS SUMMARIES)
s_0_summ |>
  summarise(across(
    c(set_std, trait_mut, set_mut, dir_exp, baseline, set_mean, total_sets, population, b_c_ratio, synergy ), 
    ~ list(unique(.))
  )) |> as.list()

# testing mean of local c freqs equals global c freq
# test |> slice_sample(n = 1) |> 
#   select(c(ifelse_value_any_individuals_with_strategy_cooperator_mean_local_c_freq_of_individuals_with_strategy_cooperator_0,
#            ifelse_value_any_individuals_with_strategy_defector_mean_local_c_freq_of_individuals_with_strategy_defector_0,
#            count_individuals_with_strategy_cooperator_n,
#            n
#   )) |> 
#   mutate(mean_local_c_freq = count_individuals_with_strategy_cooperator_n * 
#            ifelse_value_any_individuals_with_strategy_cooperator_mean_local_c_freq_of_individuals_with_strategy_cooperator_0 + 
#            (1 - count_individuals_with_strategy_cooperator_n) * 
#            ifelse_value_any_individuals_with_strategy_defector_mean_local_c_freq_of_individuals_with_strategy_defector_0) |> 
#   view()

# ====
# ==== LOADING AND CLEANING ====
# svb_x0 <- read_csv("s_vs_b_x0-table.csv", 
#                                            skip = 6) |> exp_cleaning() |> 
#   filter(population == 200, trait_mut == 0)
# 
# svb_x1 <- read_csv("s_vs_b_x1-table.csv", 
#                                            skip = 6) |> exp_cleaning() |> 
#   filter(population == 200, trait_mut == 0)
# 
# svb_ext <- read_csv("s_vs_b_ext-table.csv", 
#                    skip = 6) |> exp_cleaning()
# 
# mvk <- read_csv("m_vs_k-table.csv", 
#                    skip = 6) |> exp_cleaning()
# 
# 
# # test <- read_csv("set-sd-test.csv", 
# #                  skip = 6) |> exp_cleaning()
# 
# s_0 <- read_csv("s_0-table.csv", 
#                    skip = 6) |> exp_cleaning() |> 
#   filter(population == 200, trait_mut == 0)

s_0.5_summ <- read.csv2("s_0.5_summ.csv") |> 
  filter(population == 200, trait_mut == 0) |> mutate(
    K_M_ratio = set_mean/total_sets,
    b_c_ratio = benefit/cost
  )
  

s_1_summ <- read.csv2("s_1_summ.csv") |> 
  filter(population == 200, trait_mut == 0) |> mutate(
    K_M_ratio = set_mean/total_sets,
    b_c_ratio = benefit/cost
  )

s_0.5_summ <- s_0.5_summ |>
  mutate(b_c_ratio = benefit / cost,
         K_M_ratio = set_mean/total_sets,) |> 
  select(-c(cost, benefit))

svb_x0_summ <- read_csv2("svb_x0_summ.csv")
svb_x1_part_summ <- read_csv2("svb_x1_part.csv")
svb_x1_ext_summ <- read_csv2("svb_x1_ext_summ.csv")
svb_x1_summ <- read_csv2("svb_x1_summ.csv")
svb_x0_v0_summ <- read_csv2("svb_x0_v0_summ.csv")
svb_x0_v0.25_summ <- read_csv2("svb_x0_v0.25_summ.csv")
svb_x0_v0.5_summ <- read_csv2("svb_x0_v0.5_summ.csv")
svb_x1_v0.25_summ <- read_csv2("svb_x1_v0.25_summ.csv")
svb_x1_v0.5_summ <- read_csv2("svb_x1_v0.5_summ.csv")
mvk_summ <- read_csv2("mvk_summ.csv")
s_0_summ <- read_csv2("s_0_summ.csv")
s_0.5_summ <- read_csv2("s_0.5_summ.csv")
s_1_summ <- read_csv2("s_1_summ.csv")
 

# ====
# ==== CREATING SUMMARIES ====

# svb_x0_part_summ <- svb_x0 |> run_summary()
# svb_x0_ext_summ <- svb_ext |> filter(dir_exp == 0) |> run_summary()
# svb_x0_summ <- bind_rows(svb_x0_part_summ, svb_x0_ext_summ)
# 
# svb_x1_part <- svb_x1 |> run_summary()
# svb_x1_ext_summ <- svb_ext |> filter(dir_exp == 1) |> run_summary()
# svb_x1_summ <- bind_rows(svb_x1_part, svb_x1_ext_summ)
# 
svb_x0_v0_summ <- svb_x0_summ |> filter(set_mut == 0)
svb_x0_v0.25_summ <- svb_x0_summ |> filter(set_mut == 0.25)
svb_x0_v0.5_summ <- svb_x0_summ |> filter(set_mut == 0.5)
# 
# svb_x1_v0.25_summ <- svb_x1_summ |> filter(set_mut == 0.25)
# svb_x1_v0.5_summ <- svb_x1_summ |> filter(set_mut == 0.5)

# test_summ <- test |> run_summary()
# 
# mvk_summ <- mvk |> run_summary()
# s_0_summ <- s_0 |> run_summary()
# s_0.5_summ <- s_0.5 |> run_summary()
# s_1_summ <- s_1 |> run_summary()

# Write objects to csv2
# write_csv2(svb_x0_summ, "svb_x0_summ.csv")
# write_csv2(svb_x1_part, "svb_x1_part.csv")
# write_csv2(svb_x1_ext_summ, "svb_x1_ext_summ.csv")
# write_csv2(svb_x1_summ, "svb_x1_summ.csv")
# write_csv2(svb_x0_v0_summ, "svb_x0_v0_summ.csv")
# write_csv2(svb_x0_v0.25_summ, "svb_x0_v0.25_summ.csv")
# write_csv2(svb_x0_v0.5_summ, "svb_x0_v0.5_summ.csv")
# write_csv2(svb_x1_v0.25_summ, "svb_x1_v0.25_summ.csv")
# write_csv2(svb_x1_v0.5_summ, "svb_x1_v0.5_summ.csv")
# write_csv2(mvk_summ, "mvk_summ.csv")
# write_csv2(s_0_summ, "s_0_summ.csv")
# write_csv2(s_0.5_summ, "s_0.5_summ.csv")
# write_csv2(s_1_summ, "s_1_summ.csv")

# ====
# ==== RANDOM FOREST VIP ====
# ==== 1. Prepare Data and Train Models

# s = 0
s_0_rf <- s_0_summ |>
  select(coop_success, set_mut, dir_exp, set_mean, total_sets, K_M_ratio, b_c_ratio) |>
  drop_na()

s_0_rfmodel <- ranger(
  formula = coop_success ~ ., 
  data = s_0_rf, 
  importance = "permutation", 
  num.trees = 500
)

# s = 0.5 (Typo corrected here)
s_0.5_rf <- s_0.5_summ |>
  select(coop_success, set_mut, dir_exp, set_mean, total_sets, K_M_ratio, b_c_ratio) |>
  drop_na()

s_0.5_rfmodel <- ranger(
  formula = coop_success ~ ., 
  data = s_0.5_rf, 
  importance = "permutation", 
  num.trees = 500
)

# s = 1
s_1_rf <- s_1_summ |>
  select(coop_success, set_mut, dir_exp, set_mean, total_sets, K_M_ratio, b_c_ratio) |>
  drop_na()

s_1_rfmodel <- ranger(
  formula = coop_success ~ ., 
  data = s_1_rf, 
  importance = "permutation", 
  num.trees = 500
)


# ==== 2. Extract and Combine Importance Scores

# ranger stores the raw scores inside the $variable.importance list
extract_vip <- function(model, synergy_value) {
  data.frame(
    Variable = names(model$variable.importance),
    Importance = as.numeric(model$variable.importance),
    Synergy = synergy_value
  )
}

vip_data <- bind_rows(
  extract_vip(s_0_rfmodel, 0),
  extract_vip(s_0.5_rfmodel, 0.5),
  extract_vip(s_1_rfmodel, 1)
)

# Reorder the Variable factor based on the *average* importance across all 3 models 
# so the bars appear in a consistent, logical order in every facet
vip_data <- vip_data |>
  group_by(Variable) |>
  mutate(Avg_Importance = mean(Importance)) |>
  ungroup() |>
  arrange(Avg_Importance) |>
  mutate(Variable = factor(Variable, levels = unique(Variable)))


# ==== 3. Build the Faceted Plot

ggplot(vip_data, aes(x = Importance, y = Variable, fill = Importance)) +
  geom_col(show.legend = FALSE) +
  # Map continuous importance to the mako scale. direction = 1 ensures high values = bright
  scale_fill_viridis_c(option = "mako", direction = 1) + 
  
  # The labeller ensures the facets say "Synergy (s) = X" instead of just the number
  facet_wrap(~ Synergy, nrow = 1, labeller = labeller(Synergy = function(x) paste("Synergy (s) =", x))) +
  
  labs(
    title = "Random Forest Variable Importance",
    subtitle = "",
    x = "Permutation Importance",
    y = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.background = element_rect(fill = "#F4F5F7", color = NA),
    panel.background = element_rect(fill = "#F4F5F7", color = NA),
    panel.grid.major.x = element_line(color = "#FFFFFF", linewidth = 0.8),
    panel.grid.minor.x = element_line(color = "#FFFFFF", linewidth = 0.4),
    panel.grid.major.y = element_blank(), # Cleaner look for horizontal bar charts
    
    panel.spacing = unit(1.5, "lines"),
    strip.text = element_text(face = "bold", size = 12, color = "#2D3748", margin = margin(b = 10)),
    
    plot.title = element_text(face = "bold", size = 18, color = "#2D3748", margin = margin(b = 6)),
    plot.subtitle = element_text(color = "#4A5568", size = 12, margin = margin(b = 20)),
    axis.title.x = element_text(face = "bold", color = "#4A5568", margin = margin(t = 12)),
    axis.text = element_text(color = "#2D3748", face = "bold"),
    
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20)
  )

ggsave("faceted_rf_vip.png", width = 14, height = 5)
# ====
# ==== DECISION TREES ====
# ==== S=0 ====
# 1. Fit the Decision Tree
# We use method = "anova" because coop_success is a continuous variable (0 to 1).
# The 'cp' (complexity parameter) prevents the tree from over-fitting to noise. 
# If the tree is too massive, increase cp (e.g., 0.02). If it's too small, decrease it (e.g., 0.005).
s_0_dtmodel <- rpart(
  formula = coop_success ~ ., 
  data = s_0_rf, 
  method = "anova",
  control = rpart.control(cp = 0.01, minsplit = 20)
)

# 2. Plot the Decision Tree
rpart.plot(
  s_0_dtmodel, 
  type = 4,               # Draws clear labels directly on the branches
  extra = 101,            # Shows the mean coop_success AND the percentage of data in that node
  fallen.leaves = TRUE,   # Forces all final outcomes (leaves) to line up cleanly at the bottom
  main = "Decision Tree. Synergy = 0",
  box.palette = "-Blues",
  shadow.col = "gray",
  nn = TRUE               # Displays the node numbers for easy reference in text
)

# ==== S=0.5 ====

s_0.5_dtmodel <- rpart(
  formula = coop_success ~ ., 
  data = s_0.5_rf, 
  method = "anova",
  control = rpart.control(cp = 0.01, minsplit = 20)
)

# 2. Plot the Decision Tree
rpart.plot(
  s_0.5_dtmodel, 
  type = 4,               
  extra = 101,           
  fallen.leaves = TRUE,   
  main = "Decision Tree. Synergy = 0.5",
  box.palette = "-Blues", 
  shadow.col = "gray",
  nn = TRUE           
)

# ==== S=1 ====

s_1_dtmodel <- rpart(
  formula = coop_success ~ ., 
  data = s_1_rf, 
  method = "anova",
  control = rpart.control(cp = 0.01, minsplit = 20)
)

# 2. Plot the Decision Tree
rpart.plot(
  s_1_dtmodel, 
  type = 4,               
  extra = 101,           
  fallen.leaves = TRUE,   
  main = "Decision Tree. Synergy = 1",
  box.palette = "-Blues", 
  shadow.col = "gray",
  nn = TRUE           
)

# ====
# ==== 1D LOESS ====
# ==== S = 0 ====
# KM RATIO VS SET MUT
s_0_KM_vs_set_mut <- s_0_summ |> extract_loess_peaks(K_M_ratio, set_mut, coop_success)

s_0_KM_vs_set_mut |> nice_table(
  title = "**Optimal density ratio ~ set mutation**",
  subtitle = "LOESS-derived peaks",
  footnote = "Peaks extracted via localized polynomial regression (span = 0.75).",
  x_label = "Set Mutation (v)",
  y_label = "Optimal *K/M* Ratio"
) |> 
  gtsave("s_0_KM_vs_set_mut.png")

# PLOT

s_0_KM_vs_set_mut |> 
  plot_optimal_peak(
    title = "Structural Optimization",
    subtitle = "Tracking optimal density",
    x_label = "Set Mutation Rate (v)",
    y_label = "Optimal K/M Ratio"
  )

ggsave("s_0_KM_vs_set_mut.png")

# KM RATIO VS BENEFIT
s_0_KM_vs_benefit <- s_0_summ |> extract_loess_peaks(K_M_ratio, b_c_ratio, coop_success)

s_0_KM_vs_benefit |> nice_table(
  title = "**Optimal density ratio ~ benefit**",
  subtitle = "LOESS-derived peaks",
  footnote = "Peaks extracted via localized polynomial regression (span = 0.75).",
  x_label = "Benefit (b)",
  y_label = "Optimal *K/M* Ratio"
) |> 
  gtsave("s_0_KM_vs_benefit.png")

# PLOT

s_0_KM_vs_benefit |> 
  plot_optimal_peak(
    title = "Structural Optimization",
    subtitle = "Tracking optimal density",
    x_label = "Benefit (b)",
    y_label = "Optimal K/M Ratio"
  )

ggsave("s_0_KM_vs_benefit.png")
# ==== S = 0.5 ====
# KM RATIO VS SET MUT
s_0.5_KM_vs_set_mut <- s_0.5_summ |> extract_loess_peaks(K_M_ratio, set_mut, coop_success)

s_0.5_KM_vs_set_mut |> nice_table(
  title = "**Optimal density ratio ~ set mutation**",
  subtitle = "LOESS-derived peaks",
  footnote = "Peaks extracted via localized polynomial regression (span = 0.75).",
  x_label = "Set Mutation (v)",
  y_label = "Optimal *K/M* Ratio"
) |> 
  gtsave("s_0.5_KM_vs_set_mut.png")

# PLOT

s_0.5_KM_vs_set_mut |> 
  plot_optimal_peak(
    title = "Structural Optimization",
    subtitle = "Tracking optimal density",
    x_label = "Set Mutation Rate (v)",
    y_label = "Optimal K/M Ratio"
  )

ggsave("s_0.5_KM_vs_set_mut.png")

# KM RATIO VS BENEFIT
s_0.5_KM_vs_benefit <- s_0.5_summ |> extract_loess_peaks(K_M_ratio, b_c_ratio, coop_success)

s_0.5_KM_vs_benefit |> nice_table(
  title = "**Optimal density ratio ~ benefit**",
  subtitle = "LOESS-derived peaks",
  footnote = "Peaks extracted via localized polynomial regression (span = 0.75).",
  x_label = "Benefit (b)",
  y_label = "Optimal *K/M* Ratio"
) |> 
  gtsave("s_0.5_KM_vs_benefit.png")

# PLOT

s_0.5_KM_vs_benefit |> 
  plot_optimal_peak(
    title = "Structural Optimization",
    subtitle = "Tracking optimal density",
    x_label = "Benefit (b)",
    y_label = "Optimal K/M Ratio"
  )

ggsave("s_0.5_KM_vs_benefit.png")

# ==== S = 1 ====
# KM RATIO VS SET MUT
s_1_KM_vs_set_mut <- s_1_summ |> extract_loess_peaks(K_M_ratio, set_mut, coop_success)

s_1_KM_vs_set_mut |> nice_table(
  title = "**Optimal density ratio ~ set mutation**",
  subtitle = "LOESS-derived peaks",
  footnote = "Peaks extracted via localized polynomial regression (span = 0.75).",
  x_label = "Set Mutation (v)",
  y_label = "Optimal *K/M* Ratio"
) |> 
  gtsave("s_1_KM_vs_set_mut.png")

# PLOT

s_1_KM_vs_set_mut |> 
  plot_optimal_peak(
    title = "Structural Optimization",
    subtitle = "Tracking optimal density",
    x_label = "Set Mutation Rate (v)",
    y_label = "Optimal K/M Ratio"
  )

ggsave("s_1_KM_vs_set_mut.png")

# KM RATIO VS BENEFIT
s_1_KM_vs_benefit <- s_1_summ |> extract_loess_peaks(K_M_ratio, b_c_ratio, coop_success)

s_1_KM_vs_benefit |> nice_table(
  title = "**Optimal density ratio ~ benefit**",
  subtitle = "LOESS-derived peaks",
  footnote = "Peaks extracted via localized polynomial regression (span = 0.75).",
  x_label = "Benefit (b)",
  y_label = "Optimal *K/M* Ratio"
) |> 
  gtsave("s_1_KM_vs_benefit.png")

# PLOT

s_1_KM_vs_benefit |> 
  plot_optimal_peak(
    title = "Structural Optimization",
    subtitle = "Tracking optimal density",
    x_label = "Benefit (b)",
    y_label = "Optimal K/M Ratio"
  )

ggsave("s_1_KM_vs_benefit.png")


# ====
# ==== 2D LOESS - synergy vs b/c ratio ====
svb_v0_loess2d <- generate_loess_grid(
  data = svb_x0_v0_summ,
  x_var = "synergy",
  y_var = "b_c_ratio",
  outcome = "coop_success"
)

svb_v0.25_loess2d <- generate_loess_grid(
  data = svb_x0_v0.25_summ,
  x_var = "synergy",
  y_var = "b_c_ratio",
  outcome = "coop_success"
)

svb_v0.5_loess2d <- generate_loess_grid(
  data = svb_x0_v0.5_summ,
  x_var = "synergy",
  y_var = "b_c_ratio",
  outcome = "coop_success"
)

# 1. Tag each existing grid with its specific mutation rate
svb_v0_loess2d$set_mut <- 0
svb_v0.25_loess2d$set_mut <- 0.25
svb_v0.5_loess2d$set_mut <- 0.5

# 2. Stack them all into one unified dataframe for faceting
faceted_svb_grid <- bind_rows(svb_v0_loess2d, svb_v0.25_loess2d, svb_v0.5_loess2d)

# 3. Calculate and print the exchange rates for each facet safely
target_muts <- c(0, 0.25, 0.5)

cat("--- Exchange Rates at 50% Success Threshold ---\n")
for(mut in target_muts) {
  temp_grid <- faceted_svb_grid |> filter(set_mut == mut)
  
  s_seq <- unique(temp_grid$synergy)
  b_seq <- unique(temp_grid$b_c_ratio)
  pred_matrix <- matrix(temp_grid$pred_outcome, nrow = length(s_seq), ncol = length(b_seq))
  
  tryCatch({
    isocline <- contourLines(x = s_seq, y = b_seq, z = pred_matrix, levels = 0.5)
    df <- data.frame(synergy = isocline[[1]]$x, b_c_ratio = isocline[[1]]$y)
    rate <- abs(coef(lm(b_c_ratio ~ synergy, data = df))["synergy"])
    
    cat(sprintf("v = %-5s | Synergy vs b/c exchange rate: %.3f\n", mut, rate))
  }, error = function(e) {
    cat(sprintf("v = %-5s | Threshold (0.5) not found in this environment.\n", mut))
  })
}
cat("-----------------------------------------------\n")

# 4. Build the faceted plot
ggplot(faceted_svb_grid, aes(x = synergy, y = b_c_ratio)) +
  geom_raster(aes(fill = pred_outcome), interpolate = TRUE) +
  geom_contour(aes(z = pred_outcome), breaks = 0.5, color = "white", linewidth = 1) +
  scale_fill_viridis_c(option = "mako", direction = 1, name = "Predicted Success") +
  
  # Splits the plot into 3 columns based on set_mut
  facet_wrap(~ set_mut, nrow = 1, labeller = labeller(set_mut = function(x) paste("Mutation (v) =", x))) +
  
  labs(
    title = "Synergy vs b_c_ratio",
    subtitle = "White line: 50% success threshold tracking the required exchange rate",
    x = "Synergy (s)",
    y = "b_c_ratio (b)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.background = element_rect(fill = "#F4F5F7", color = NA),
    panel.background = element_rect(fill = "#F4F5F7", color = NA),
    panel.grid.major = element_line(color = "#FFFFFF", linewidth = 0.8),
    panel.grid.minor = element_line(color = "#FFFFFF", linewidth = 0.4),
    
    # Add breathing room between panels
    panel.spacing = unit(1.5, "lines"),
    strip.text = element_text(face = "bold", size = 12, color = "#2D3748", margin = margin(b = 10)),
    
    plot.title = element_text(face = "bold", size = 18, color = "#2D3748", margin = margin(b = 6)),
    plot.subtitle = element_text(color = "#4A5568", size = 12, margin = margin(b = 20)),
    axis.title = element_text(face = "bold", color = "#4A5568", margin = margin(t = 12, r = 12)),
    axis.text = element_text(color = "#718096"),
    legend.position = "right",
    legend.title = element_text(face = "bold", color = "#2D3748", size = 11),
    legend.text = element_text(color = "#4A5568"),
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20)
  )

# Save with a wide aspect ratio so the 3 panels aren't squished
ggsave("faceted_svb_loess2d.png", width = 14, height = 5)
# ====
# ==== 2D LOESS - k/m ratio vs b/c ratio ====
# ==== FACET PLOT S=0 ====
# 1. Define the specific mutation rates you want to plot
target_muts <- c(0, 0.15, 0.3)

# 2. Generate a separate LOESS grid for each mutation rate and combine them
grid_list <- list()
for(mut in target_muts) {
  # Isolate the data for just this mutation rate
  subset_data <- s_0_summ |> filter(set_mut == mut)
  
  # Generate the grid
  grid <- generate_loess_grid(subset_data, K_M_ratio, b_c_ratio, coop_success)
  
  # Stamp the grid with its mutation rate so ggplot knows how to separate them later
  grid$set_mut <- mut
  
  # Store it in our list
  grid_list[[as.character(mut)]] <- grid
}

# Stack all the grids into one massive dataframe
faceted_grid <- bind_rows(grid_list)

# 3. Calculate and print the exchange rates for each facet
cat("--- Exchange Rates at 50% Success Threshold ---\n")
for(mut in target_muts) {
  # Extract the specific grid
  temp_grid <- faceted_grid |> filter(set_mut == mut)
  
  km_seq <- unique(temp_grid$K_M_ratio)
  bc_seq <- unique(temp_grid$b_c_ratio)
  pred_matrix <- matrix(temp_grid$pred_outcome, nrow = length(km_seq), ncol = length(bc_seq))
  
  # Wrap in tryCatch to prevent the loop from crashing if a threshold doesn't exist for one facet
  tryCatch({
    isocline <- contourLines(x = km_seq, y = bc_seq, z = pred_matrix, levels = 0.5)
    df <- data.frame(K_M_ratio = isocline[[1]]$x, b_c_ratio = isocline[[1]]$y)
    rate <- abs(coef(lm(b_c_ratio ~ K_M_ratio, data = df))["K_M_ratio"])
    
    cat(sprintf("v = %-5s | K/M vs b/c exchange rate: %.3f\n", mut, rate))
  }, error = function(e) {
    cat(sprintf("v = %-5s | Threshold (0.5) not found in this environment.\n", mut))
  })
}
cat("-----------------------------------------------\n")
# remove predictions outside of the (0, 1) range
faceted_grid$pred_outcome <- pmin(pmax(faceted_grid$pred_outcome, 0), 1)
# 4. Build the faceted plot
# We build this manually rather than using the function to easily inject facet_wrap
ggplot(faceted_grid, aes(x = K_M_ratio, y = b_c_ratio)) +
  geom_raster(aes(fill = pred_outcome), interpolate = TRUE) +
  geom_contour(aes(z = pred_outcome), breaks = 0.5, color = "white", linewidth = 1) +
  scale_fill_viridis_c(option = "mako", direction = 1, name = "Predicted Success") +
  
  # THIS IS THE MAGIC LINE: It splits the plot into columns based on set_mut
  facet_wrap(~ set_mut, nrow = 1, labeller = labeller(set_mut = function(x) paste("Mutation (v) =", x))) +
  
  labs(
    title = "K/M ratio vs b/c ratio. Synergy = 0",
    subtitle = "White line marks 50% success rate",
    x = "K/M Ratio",
    y = "b/c Ratio"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.background = element_rect(fill = "#F4F5F7", color = NA),
    panel.background = element_rect(fill = "#F4F5F7", color = NA),
    panel.grid.major = element_line(color = "#FFFFFF", linewidth = 0.8),
    panel.grid.minor = element_line(color = "#FFFFFF", linewidth = 0.4),
    
    # Add breathing room between the individual facet panels
    panel.spacing = unit(1.5, "lines"),
    strip.text = element_text(face = "bold", size = 12, color = "#2D3748", margin = margin(b = 10)),
    
    plot.title = element_text(face = "bold", size = 18, color = "#2D3748", margin = margin(b = 6)),
    plot.subtitle = element_text(color = "#4A5568", size = 12, margin = margin(b = 20)),
    axis.title = element_text(face = "bold", color = "#4A5568", margin = margin(t = 12, r = 12)),
    axis.text = element_text(color = "#718096"),
    legend.position = "right",
    legend.title = element_text(face = "bold", color = "#2D3748", size = 11),
    legend.text = element_text(color = "#4A5568"),
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20)
  )

ggsave("faceted_km_vs_bc_loess2d_s0.png", width = 16, height = 5) # Increased width for horizontal facets
# ==== FACET PLOT S=0.5====
# 1. Define the specific mutation rates you want to plot
target_muts <- c(0, 0.2, 0.4)

# 2. Generate a separate LOESS grid for each mutation rate and combine them
grid_list <- list()
for(mut in target_muts) {
  # Isolate the data for just this mutation rate
  subset_data <- s_0.5_summ |> filter(set_mut == mut)
  
  # Generate the grid
  grid <- generate_loess_grid(subset_data, K_M_ratio, b_c_ratio, coop_success)
  
  # Stamp the grid with its mutation rate so ggplot knows how to separate them later
  grid$set_mut <- mut
  
  # Store it in our list
  grid_list[[as.character(mut)]] <- grid
}

# Stack all the grids into one massive dataframe
faceted_grid <- bind_rows(grid_list)

# 3. Calculate and print the exchange rates for each facet
cat("--- Exchange Rates at 50% Success Threshold ---\n")
for(mut in target_muts) {
  # Extract the specific grid
  temp_grid <- faceted_grid |> filter(set_mut == mut)
  
  km_seq <- unique(temp_grid$K_M_ratio)
  bc_seq <- unique(temp_grid$b_c_ratio)
  pred_matrix <- matrix(temp_grid$pred_outcome, nrow = length(km_seq), ncol = length(bc_seq))
  
  # Wrap in tryCatch to prevent the loop from crashing if a threshold doesn't exist for one facet
  tryCatch({
    isocline <- contourLines(x = km_seq, y = bc_seq, z = pred_matrix, levels = 0.5)
    df <- data.frame(K_M_ratio = isocline[[1]]$x, b_c_ratio = isocline[[1]]$y)
    rate <- abs(coef(lm(b_c_ratio ~ K_M_ratio, data = df))["K_M_ratio"])
    
    cat(sprintf("v = %-5s | K/M vs b/c exchange rate: %.3f\n", mut, rate))
  }, error = function(e) {
    cat(sprintf("v = %-5s | Threshold (0.5) not found in this environment.\n", mut))
  })
}
cat("-----------------------------------------------\n")
# remove predictions outside of the (0, 1) range
faceted_grid$pred_outcome <- pmin(pmax(faceted_grid$pred_outcome, 0), 1)

# 4. Build the faceted plot
# We build this manually rather than using the function to easily inject facet_wrap
ggplot(faceted_grid, aes(x = K_M_ratio, y = b_c_ratio)) +
  geom_raster(aes(fill = pred_outcome), interpolate = TRUE) +
  geom_contour(aes(z = pred_outcome), breaks = 0.5, color = "white", linewidth = 1) +
  scale_fill_viridis_c(option = "mako", direction = 1, name = "Predicted Success") +
  
  # THIS IS THE MAGIC LINE: It splits the plot into columns based on set_mut
  facet_wrap(~ set_mut, nrow = 1, labeller = labeller(set_mut = function(x) paste("Mutation (v) =", x))) +
  
  labs(
    title = "K/M ratio vs b/c ratio. Synergy = 0.5",
    subtitle = "White line marks 50% success rate",
    x = "K/M Ratio",
    y = "b/c Ratio"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.background = element_rect(fill = "#F4F5F7", color = NA),
    panel.background = element_rect(fill = "#F4F5F7", color = NA),
    panel.grid.major = element_line(color = "#FFFFFF", linewidth = 0.8),
    panel.grid.minor = element_line(color = "#FFFFFF", linewidth = 0.4),
    
    # Add breathing room between the individual facet panels
    panel.spacing = unit(1.5, "lines"),
    strip.text = element_text(face = "bold", size = 12, color = "#2D3748", margin = margin(b = 10)),
    
    plot.title = element_text(face = "bold", size = 18, color = "#2D3748", margin = margin(b = 6)),
    plot.subtitle = element_text(color = "#4A5568", size = 12, margin = margin(b = 20)),
    axis.title = element_text(face = "bold", color = "#4A5568", margin = margin(t = 12, r = 12)),
    axis.text = element_text(color = "#718096"),
    legend.position = "right",
    legend.title = element_text(face = "bold", color = "#2D3748", size = 11),
    legend.text = element_text(color = "#4A5568"),
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20)
  )

ggsave("faceted_km_vs_bc_loess2d_s0.5.png", width = 16, height = 5) # Increased width for horizontal facets
# ==== FACET PLOT S=1 ====
# 1. Define the specific mutation rates you want to plot
target_muts <- c(0, 0.2, 0.4)

# 2. Generate a separate LOESS grid for each mutation rate and combine them
grid_list <- list()
for(mut in target_muts) {
  # Isolate the data for just this mutation rate
  subset_data <- s_1_summ |> filter(set_mut == mut)
  
  # Generate the grid
  grid <- generate_loess_grid(subset_data, K_M_ratio, b_c_ratio, coop_success)
  
  # Stamp the grid with its mutation rate so ggplot knows how to separate them later
  grid$set_mut <- mut
  
  # Store it in our list
  grid_list[[as.character(mut)]] <- grid
}

# Stack all the grids into one massive dataframe
faceted_grid <- bind_rows(grid_list)

# 3. Calculate and print the exchange rates for each facet
cat("--- Exchange Rates at 50% Success Threshold ---\n")
# remove predictions outside of the (0, 1) range
faceted_grid$pred_outcome <- pmin(pmax(faceted_grid$pred_outcome, 0), 1)
for(mut in target_muts) {
  # Extract the specific grid
  temp_grid <- faceted_grid |> filter(set_mut == mut)
  
  km_seq <- unique(temp_grid$K_M_ratio)
  bc_seq <- unique(temp_grid$b_c_ratio)
  pred_matrix <- matrix(temp_grid$pred_outcome, nrow = length(km_seq), ncol = length(bc_seq))
  
  # Wrap in tryCatch to prevent the loop from crashing if a threshold doesn't exist for one facet
  tryCatch({
    isocline <- contourLines(x = km_seq, y = bc_seq, z = pred_matrix, levels = 0.5)
    df <- data.frame(K_M_ratio = isocline[[1]]$x, b_c_ratio = isocline[[1]]$y)
    rate <- abs(coef(lm(b_c_ratio ~ K_M_ratio, data = df))["K_M_ratio"])
    
    cat(sprintf("v = %-5s | K/M vs b/c exchange rate: %.3f\n", mut, rate))
  }, error = function(e) {
    cat(sprintf("v = %-5s | Threshold (0.5) not found in this environment.\n", mut))
  })
}
cat("-----------------------------------------------\n")

# 4. Build the faceted plot
# We build this manually rather than using the function to easily inject facet_wrap
ggplot(faceted_grid, aes(x = K_M_ratio, y = b_c_ratio)) +
  geom_raster(aes(fill = pred_outcome), interpolate = TRUE) +
  geom_contour(aes(z = pred_outcome), breaks = 0.5, color = "white", linewidth = 1) +
  scale_fill_viridis_c(option = "mako", direction = 1, name = "Predicted Success") +
  
  # THIS IS THE MAGIC LINE: It splits the plot into columns based on set_mut
  facet_wrap(~ set_mut, nrow = 1, labeller = labeller(set_mut = function(x) paste("Mutation (v) =", x))) +
  
  labs(
    title = "K/M ratio vs b/c ratio. Synergy = 1",
    subtitle = "White line marks 50% success rate",
    x = "K/M Ratio",
    y = "b/c Ratio"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.background = element_rect(fill = "#F4F5F7", color = NA),
    panel.background = element_rect(fill = "#F4F5F7", color = NA),
    panel.grid.major = element_line(color = "#FFFFFF", linewidth = 0.8),
    panel.grid.minor = element_line(color = "#FFFFFF", linewidth = 0.4),
    
    # Add breathing room between the individual facet panels
    panel.spacing = unit(1.5, "lines"),
    strip.text = element_text(face = "bold", size = 12, color = "#2D3748", margin = margin(b = 10)),
    
    plot.title = element_text(face = "bold", size = 18, color = "#2D3748", margin = margin(b = 6)),
    plot.subtitle = element_text(color = "#4A5568", size = 12, margin = margin(b = 20)),
    axis.title = element_text(face = "bold", color = "#4A5568", margin = margin(t = 12, r = 12)),
    axis.text = element_text(color = "#718096"),
    legend.position = "right",
    legend.title = element_text(face = "bold", color = "#2D3748", size = 11),
    legend.text = element_text(color = "#4A5568"),
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20)
  )

ggsave("faceted_km_vs_bc_loess2d_s1.png", width = 16, height = 5) # Increased width for horizontal facets
# ====
# ====
# # ==== SUCCESS HEATMAP ====
# 
# # # # # # # # # # # # # # # # # # # # # # # # # #
# # SYNERGY VS BENEFIT WITHOUT DIRECTED EXPLORATION
# # # # # # # # # # # # # # # # # # # # # # # # # #
# 
# svb_x0_heatmap <- svb_x0_summ |>
#   group_by(synergy, benefit, set_mut) |>
#   summarise(
#     avg_success = mean(coop_success, na.rm = TRUE),
#     avg_bounce = mean(bounce_back, na.rm = TRUE),
#     .groups = "drop"
#   )
# 
# phase_space(svb_heatmap, synergy, benefit, set_mut, avg_success, avg_bounce,
#             plot_title = "Cooperation Success Phase Space",
#             plot_subtitle = "Tile color indicates average success. Circle size indicates volatility.",
#             x_label = "Synergy (s)",
#             y_label = "Benefit (b)",
#             fill_label = "Coop Success",
#             size_label = "Volatility")
# ggsave("Graphs/svb_x0_heatmap.png", width = 8, height = 6, dpi = 300)
# 
# # # # # # # # # # # # # # # # # # # # # # # # #
# # SYNERGY VS BENEFIT WITH DIRECTED EXPLORATION
# # # # # # # # # # # # # # # # # # # # # # # # # 
# 
# svb_x1_heatmap <- svb_x1_summ |>
#   # Average multiple runs for each parameter combination
#   group_by(synergy, benefit, set_mut) |>
#   summarise(
#     avg_success = mean(coop_success, na.rm = TRUE),
#     avg_bounce = mean(bounce_back, na.rm = TRUE),
#     .groups = "drop"
#   )
# 
# phase_space(svb_heatmap, synergy, benefit, avg_success, avg_bounce)
# ggsave("Graphs/svb_x1_heatmap.png", width = 8, height = 6, dpi = 300)
# 
# # # # # # # # # # # # # # # # # # #
# # TOTAL SET NUMBER VERSUS SET MEAN
# # # # # # # # # # # # # # # # # # #
# 
# # NO TRAIT MUTATION, NO DIRECTED EXPLORATION
# mvk_heatmap <- mvk_summ |>
#   filter(trait_mut == 0, dir_exp == 0) |> 
#   # Average multiple runs for each parameter combination
#   group_by(total_sets, set_mean, set_mut) |>
#   summarise(
#     avg_success = mean(coop_success, na.rm = TRUE),
#     avg_bounce = mean(bounce_back, na.rm = TRUE),
#     .groups = "drop"
#   )
# 
# phase_space(mvk_heatmap, total_sets, set_mean, set_mut, , avg_success, avg_bounce,
#             plot_title = "Cooperation Success Phase Space",
#             plot_subtitle = "Tile color indicates average success. Circle size indicates volatility.",
#             x_label = "Total sets (M)",
#             y_label = "Set mean (K)",
#             fill_label = "Coop Success",
#             size_label = "Volatility")
# ggsave("Graphs/mvk_u0_x0_heatmap.png", width = 8, height = 6, dpi = 300)
# 
# # NO TRAIT MUTATION, DIRECTED EXPLORATION
# mvk_heatmap <- mvk_summ |>
#   filter(trait_mut == 0, dir_exp == 1) |> 
#   # Average multiple runs for each parameter combination
#   group_by(total_sets, set_mean, set_mut) |>
#   summarise(
#     avg_success = mean(coop_success, na.rm = TRUE),
#     avg_bounce = mean(bounce_back, na.rm = TRUE),
#     .groups = "drop"
#   )
# 
# phase_space(mvk_heatmap, total_sets, set_mean, set_mut, , avg_success, avg_bounce,
#             plot_title = "Cooperation Success Phase Space",
#             plot_subtitle = "Tile color indicates average success. Circle size indicates volatility.",
#             x_label = "Total sets (M)",
#             y_label = "Set mean (K)",
#             fill_label = "Coop Success",
#             size_label = "Volatility")
# ggsave("Graphs/mvk_u0_x1_heatmap.png", width = 8, height = 6, dpi = 300)
# 
# # TRAIT MUTATION, NO DIRECTED EXPLORATION
# mvk_heatmap <- mvk_summ |>
#   filter(trait_mut == 0.15, dir_exp == 0) |> 
#   # Average multiple runs for each parameter combination
#   group_by(total_sets, set_mean, set_mut) |>
#   summarise(
#     avg_success = mean(coop_success, na.rm = TRUE),
#     avg_bounce = mean(bounce_back, na.rm = TRUE),
#     .groups = "drop"
#   )
# 
# phase_space(mvk_heatmap, total_sets, set_mean, set_mut, , avg_success, avg_bounce,
#             plot_title = "Cooperation Success Phase Space",
#             plot_subtitle = "Tile color indicates average success. Circle size indicates volatility.",
#             x_label = "Total sets (M)",
#             y_label = "Set mean (K)",
#             fill_label = "Coop Success",
#             size_label = "Volatility")
# ggsave("Graphs/mvk_u0_x0_heatmap.png", width = 8, height = 6, dpi = 300)
# 
# # TRAIT MUTATION, DIRECTED EXPLORATION
# mvk_heatmap <- mvk_summ |>
#   filter(trait_mut == 0.30, dir_exp == 1) |> 
#   # Average multiple runs for each parameter combination
#   group_by(total_sets, set_mean, set_mut) |>
#   summarise(
#     avg_success = mean(coop_success, na.rm = TRUE),
#     avg_bounce = mean(bounce_back, na.rm = TRUE),
#     .groups = "drop"
#   )
# 
# phase_space(mvk_heatmap, total_sets, set_mean, set_mut, , avg_success, avg_bounce,
#             plot_title = "Cooperation Success Phase Space",
#             plot_subtitle = "Tile color indicates average success. Circle size indicates volatility.",
#             x_label = "Total sets (M)",
#             y_label = "Set mean (K)",
#             fill_label = "Coop Success",
#             size_label = "Volatility")
# ggsave("Graphs/mvk_u0_x1_heatmap.png", width = 8, height = 6, dpi = 300)
# 
# 
# # # # # # # # # 
# # SET STD TEST
# # # # # # # # # 
# 
# # test_heatmap <- test_summ |>
# #   group_by(set_mut, set_std, dir_exp) |>
# #   summarise(
# #     avg_success = mean(coop_success, na.rm = TRUE),
# #     avg_bounce = mean(bounce_back, na.rm = TRUE),
# #     .groups = "drop"
# #   )
# # 
# # phase_space(test_heatmap, set_std, set_mut, dir_exp, , avg_success, avg_bounce,
# #             plot_title = "Cooperation Success Phase Space",
# #             plot_subtitle = "Tile color indicates average success. Circle size indicates volatility.",
# #             x_label = "",
# #             y_label = "",
# #             fill_label = "Coop Success",
# #             size_label = "Volatility")
# # ggsave("Graphs/svb_x0_heatmap.png", width = 8, height = 6, dpi = 300)
# 
# # # # # # # # # # # # # # # # # # # 
# # KM RATIO VS BC RATIO BY SET X = 0
# # # # # # # # # # # # # # # # # # # 
# 
# s0x0_heatmap <- s_0_summ |>
#   filter(dir_exp == 0, set_mut != 0.5) |> 
#   group_by(K_M_ratio, b_c_ratio, set_mut) |>
#   summarise(
#     avg_success = mean(coop_success, na.rm = TRUE),
#     avg_bounce = mean(bounce_back, na.rm = TRUE),
#     .groups = "drop"
#   )
# 
# phase_space(s0x0_heatmap, K_M_ratio, b_c_ratio, set_mut, , avg_success, avg_bounce,
#             plot_title = "Cooperation Success Phase Space",
#             plot_subtitle = "Tile color indicates average success. Circle size indicates volatility.",
#             x_label = "Synergy (s)",
#             y_label = "Benefit (b)",
#             fill_label = "Coop Success",
#             size_label = "Volatility")
# ggsave("Graphs/s0x0_heatmap.png", width = 8, height = 6, dpi = 300)
# 
# # # # # # # # # # # # # # # # # # # 
# # KM RATIO VS BC RATIO BY SET X = 0
# # # # # # # # # # # # # # # # # # # 
# 
# s0x1_heatmap <- s_0_summ |>
#   filter(dir_exp == 1, set_mut != 0.5) |> 
#   group_by(K_M_ratio, b_c_ratio, set_mut) |>
#   summarise(
#     avg_success = mean(coop_success, na.rm = TRUE),
#     avg_bounce = mean(bounce_back, na.rm = TRUE),
#     .groups = "drop"
#   )
# 
# phase_space(s0x1_heatmap, K_M_ratio, b_c_ratio, set_mut, , avg_success, avg_bounce,
#             plot_title = "Cooperation Success Phase Space",
#             plot_subtitle = "Tile color indicates average success. Circle size indicates volatility.",
#             x_label = "Synergy (s)",
#             y_label = "Benefit (b)",
#             fill_label = "Coop Success",
#             size_label = "Volatility")
# ggsave("Graphs/s0x1_heatmap.png", width = 8, height = 6, dpi = 300)
# 
# # ====
# # ==== DOT GRID ====
# 
# # SET MUTATION 0
# run_metrics <- svb_x0 |>
#   group_by(b_c_ratio, synergy, set_mut) |>
#   summarise(
#     mean_coop = mean(coop_freq, na.rm = TRUE),
#     std_coop = sd(coop_freq, na.rm = TRUE),
#     .groups = 'drop'
#   )
# 
# # 2. Aggregate across all runs for a given parameter combination
# # (Alternatively, you can plot run_metrics directly with position_jitter 
# # to see every individual run instead of the aggregate)
# grid_data <- run_metrics |>
#   group_by(set_mean, total_sets, population) |>
#   summarise(
#     avg_mean_coop = mean(mean_coop, na.rm = TRUE),
#     avg_std_coop = mean(std_coop, na.rm = TRUE),
#     .groups = 'drop'
#   )
# 
# # 3. Plot the dot grid
# ggplot(grid_data, aes(x = factor(set_mean), y = factor(total_sets))) +
#   geom_point(aes(size = avg_std_coop, color = avg_mean_coop)) +
#   facet_wrap(~ population, labeller = as_labeller(function(x) paste("Population:", x))) +
#   scale_color_viridis_c(option = "plasma", limits = c(0, 1)) +
#   scale_size_continuous(range = c(3, 12)) + # Adjusts minimum and maximum bubble sizes
#   labs(
#     title = "Cooperation Dynamics: Mean and Volatility",
#     subtitle = "Color = Average Cooperation | Size = Volatility (Std. Dev.)",
#     x = "Mean Sets per Agent (set_mean)",
#     y = "Total Sets Available (total_sets)",
#     color = "Mean Coop\nFreq",
#     size = "Std Dev\n(Volatility)"
#   ) +
#   theme_minimal(base_size = 14) +
#   theme(
#     plot.background = element_rect(fill = "#FEFCF5", color = NA),
#     panel.background = element_rect(fill = "#FEFCF5", color = NA),
#     panel.grid.major = element_line(color = "grey85", linewidth = 0.5, linetype = "dashed"),
#     strip.background = element_rect(fill = "grey90", color = NA),
#     strip.text = element_text(face = "bold", size = 12)
#   )
# 
# 
# 
# # ====
# # ==== SUCCESS HEATMAP ALTERNATIVE ====
# 
# ggplot(heatmap_data, aes(x = factor(synergy), y = factor(b_c_ratio))) +
#   # The background tiles
#   geom_tile(aes(fill = avg_success), color = "white", linewidth = 0.5) +
#   
#   # The bounce-back circles
#   # We map color to avg_success as well, but we will reverse the palette later!
#   geom_point(
#     aes(size = avg_bounce, color = "#EE4B2B"), # put avg_success here for viridis mode
#     show.legend = c(size = TRUE, color = FALSE) # Hide the point color legend
#   ) +
#   
#   scale_fill_gradient(low = "black", high = "white", name = "Coop Success", limits = c(0, 1)) +
#   scale_size_continuous(name = "Bounce Back", range = c(0, 6)) +
#   # scale_fill_viridis_c(name = "Coop Success", limits = c(0, 1)) + # to change to viridis
#   # scale_color_viridis_c(direction = -1, limits = c(0, 1)) +
#   
#   labs(
#     title = "Cooperation Success Phase Space",
#     subtitle = "Tile color indicates average success. Circle size indicates post-resolution volatility.",
#     x = "Synergy (s)",
#     y = "b_c_ratio (b)"
#   ) +
#   theme_minimal() +
#   theme(
#     panel.grid = element_blank(), # Removes grid lines to make tiles pop
#     axis.text = element_text(size = 10)
#   ) +
#   facet_wrap(~ set_mut, labeller = label_both)
# # ==== 
# 
# 
# # ==== CHECK RUN LENGTH ====
# # 1. Keep only the row with the maximum step per runhttp://127.0.0.1:20622/graphics/plot_zoom_png?width=1184&height=861
# df_max_steps <- svb_x0 |>
#   filter(b_c_ratio == 0.2, synergy == 0.66, set_mut == 0.5) |> 
#   group_by(run_number) |>
#   slice_max(step, n = 1, with_ties = FALSE) |>
#   ungroup()
# 
# # 2. Plot as horizontal columns
# ggplot(df_max_steps, aes(x = step, y = reorder(factor(run_number), step))) +
#   geom_col(fill = "steelblue") +
#   labs(
#     title = "Maximum Step Reached by Run",
#     x = "Maximum Step Reached",
#     y = "Run Number"
#   ) +
#   theme_minimal()
# 
# # ====
# # ==== COOP PLOT ====
# 
# svb_x0 |>
#   filter(run_number == 394) |> 
#   ggplot(aes(x = step, color = run_number)) +
#   geom_line(aes(y = coop_freq, linetype = "Cooperator Frequency"), linewidth = 1.2, alpha = 0.85) +
#   geom_line(aes(y = variance_ratio, linetype = "Variance Ratio"), linewidth = 1.2, alpha = 0.85) +
#   scale_color_viridis_d(option = "mako", end = 0.9) +
#   scale_linetype_manual(values = c("Cooperator Frequency" = "solid", "Variance Ratio" = "dashed")) +
#   scale_y_continuous(limits = c(NA, 1), expand = expansion(mult = c(0, 0.05))) +
#   labs(
#     subtitle ="",
#     x = "Simulation Step",
#     y = "Value",
#     color = NULL,
#     linetype = NULL
#   ) +
#   theme_minimal(base_size = 14) +
#   theme(
#     plot.title = element_text(face = "bold", size = 18, margin = margin(b = 5)),
#     plot.subtitle = element_text(color = "grey40", size = 12, margin = margin(b = 15)),
#     legend.position = "top",
#     legend.justification = "left",
#     legend.text = element_text(size = 11),
#     panel.grid.minor = element_blank(), 
#     panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
#     axis.text = element_text(color = "grey30"),
#     axis.title = element_text(face = "bold", margin = margin(t = 10, r = 10)),
#     plot.background = element_rect(fill = "#FEFCF5", color = NA), 
#     panel.background = element_rect(fill = "#FEFCF5", color = NA)
#   )
# # ====