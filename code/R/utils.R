# utils.R
# This script contains utility functions

# Function that returns base directory
base_dir <- function() {
  return(path.expand("~/Documents/phd"))
}

# Function to load specific metric data and add a phase column
load_metric <- function(path) {
  # Check if the file exists
  if (file.exists(file_path)) {
    # Read the data, first row as column names (experiment names)
    data <- read.csv(file_path)
    data$Phase <- factor(data$Phase,
      levels = c("proestrus", "estrus", "metestrus", "diestrus") # Set order
    )
    return(data)
  }
}
