# statistical_analysis.R
# This is the main script

# Load required packages
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(lme4))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggsignif))
suppressPackageStartupMessages(library(emmeans))
suppressPackageStartupMessages(library(car))
suppressPackageStartupMessages(library(FSA))
suppressPackageStartupMessages(library(rstatix))


# Source other scripts
source("utils.R")
source("plots.R")

# Create parser
parser <- ArgumentParser(description = "Perform statistical analysis on data.")
parser$add_argument("dir",
  type = "character",
  help = "path to data directory from base directory"
)
parser$add_argument("metric", type = "character", help = "metric to analyse")
parser$add_argument("--save-dir",
  type = "character", default = "uterine-microCT/figures",
  help = "path to save directory"
)

# Parse arguments
args <- parser$parse_args()

# Retrieve arguments
data_dir <- file.path(base_dir(), args$dir)
metric <- args$metric

# Load the data
file_path <- file.path(data_dir, paste0(metric, ".csv"))
data <- load_metric(file_path)

# One-away ANOVA test
model <- oneway.test(Value ~ Phase, data = data, var.equal = FALSE)
print(model)

# Post-hoc Games-Howell test
games_howellresult <- games_howell_test(Value ~ Phase, data = data)
print(games_howellresult)

# Plot results
plot_anova_results(data, model, metric)

# Save the plot to a file
save_file <- paste0(metric, ".png")
save_path <- file.path(file.path(base_dir(), args$save_dir), save_file)
ggsave(save_path, plot = last_plot(), dpi = 300)
