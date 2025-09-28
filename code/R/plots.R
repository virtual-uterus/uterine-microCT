# plots.R
# This script contains plotting functions

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggsignif))
suppressPackageStartupMessages(library(ggpattern))
suppressPackageStartupMessages(library(rstatix))


# Function to plot the statistical analysis results
plot_anova_results <- function(data, model, metric) {
  # Compute maximum y-values for each comparison to avoid overlap
  y_max <- max(data$Value, na.rm = TRUE)
  y_min <- 0

  # Calculate the means and standard deviations for each phase
  data_summary <- data %>%
    group_by(Phase) %>%
    summarise(
      Mean = mean(Value, na.rm = TRUE),
      STD = sd(Value, na.rm = TRUE) / sqrt(n()),
      Value = Mean,
      .groups = "drop"
    )

  # Plot
  p <- ggplot(data, aes(x = Phase, y = Value)) +
    geom_jitter(
      data = data,
      aes(x = Phase, y = Value),
      color = "black",
      size = 3,
      shape = 16,
      show.legend = FALSE,
      width = 0.1
    ) +
    geom_point(
      data = data_summary,
      aes(x = Phase, y = Mean, color = factor(Phase)),
      size = 4,
      shape = 16,
      show.legend = FALSE
    ) +
    geom_errorbar(
      data = data_summary, # Use data_summary for error bars
      aes(x = Phase, ymin = Mean - STD, ymax = Mean + STD, color = Phase),
      linewidth = 1,
      width = 0.2,
      show.legend = FALSE
    ) +
    theme_classic(base_size = 21) +
    scale_x_discrete(labels = function(x) {
      tools::toTitleCase(x)
    }) +
    labs(
      x = "Estrus phase",
      y = get_label(metric)
    ) +
    scale_colour_brewer(palette = "Set1") +
    coord_cartesian(ylim = c(y_min, NA))

  # Pairwise comparison
  comparison_results <- games_howell_test(Value ~ Phase, data = data)

  # Keep only significant results
  significant_comparisons <- comparison_results %>%
    filter(p.adj < 0.05) %>%
    mutate(stars = "*")

  # Prepare list of comparisons for geom_signif
  comparisons <- significant_comparisons %>%
    select(group1, group2) %>%
    pmap(~ c(..1, ..2)) # list of character vectors

  if (nrow(significant_comparisons) > 0) {
    offset <- y_max * 0.05 * (nrow(significant_comparisons)) # Adjust offset based on the number of bars

    # Generate y_positions for each comparison
    bar_positions <- y_max + seq_len(nrow(significant_comparisons)) * offset

    # Plot with adjusted significance bars
    p <- p + geom_signif(
      comparisons = comparisons,
      annotations = significant_comparisons$stars,
      map_signif_level = TRUE,
      textsize = 5,
      tip_length = 0.02, # Controls the length of the brackets
      y_position = bar_positions,
      margin_top = 0.0,
      size = 1.3
    )
  } else {
    # If there are no significant comparisons, reset y_max
    p <- p + coord_cartesian(ylim = c(y_min, y_max))
  }
}


# Function to get the proper y-axis label based on the metric
get_label <- function(metric) {
  # Define y-axis labels based on the metric
  label <- switch(metric,
    "muscle_thickness" = expression("Normalised muscle thickness (mm mg"^-1 * ")"),
    "endometrium_volume" = expression("Normalised endometrium volume (mm"^3 * " mg"^-1 * ")"),
    "muscle_volume" = expression("Normalised myometrium volume (mm"^3 * " mg"^-1 * ")"),
    "length" = expression("Normalised horn length (mm mg"^-1 * ")"),
    "radius" = expression("Normalised horn radius (mm mg"^-1 * ")"),
    "Unknown Metric" # Fallback for metrics not listed
  )
  return(label)
}
