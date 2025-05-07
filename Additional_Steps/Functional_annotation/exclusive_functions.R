#!/usr/bin/env Rscript

# --- Install missing packages ---
required_packages <- c("ggplot2", "scales", "stringr")
missing_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if (length(missing_packages)) {
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}

# --- Load libraries ---
suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
  library(stringr)
})

# --- Help message ---
if (length(commandArgs(trailingOnly = TRUE)) == 1 && commandArgs(trailingOnly = TRUE)[1] == "--help") {
  cat("Usage: Rscript exclusive_functions.R <input_file> <dataset_name> <sample_type> <top_functions> <colour_graph>\n")
  cat("\nArguments:\n")
  cat("  <input_file>          : The input file containing data for analysis.\n")
  cat("  <dataset_name>        : The name of the organism to be analyzed (e.g., Human, Ruminants).\n")
  cat("  <sample_type>         : The type of sample (e.g., Mags, Genes, miRNA).\n")
  cat("  <top_functions>       : Number of top functions to display.\n")
  cat("  <colour_graph>        : Color to be used for the bar chart (e.g., blue, '#FF5733').\n")
  quit(status = 0)
}

# --- Argument parsing ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop("Usage: Rscript exclusive_functions.R <input_file> <dataset_name> <sample_type> <top_functions> <colour_graph>")
}

input_file <- args[1]
dataset_name <- args[2]
sample_type <- args[3]
top <- as.integer(args[4])
colour_graph <- args[5]

# --- Load data ---
data <- read.table(input_file, header = TRUE, sep = '\t')
num_rows <- nrow(data)

if (num_rows == 0) {
  stop(paste("There are no exclusive pathways for the organism '", dataset_name, "'.", sep = ""))
}

if (top > num_rows) {
  stop(paste("Requested top functions (", top, ") exceed available (", num_rows, "). Adjust accordingly.", sep = ""))
}

# --- Process data ---
expression_means <- rowMeans(as.data.frame(data[, -1]))
data_function <- data.frame(Function = data[, 1], Expression = expression_means)
top_functions <- data_function[order(-data_function$Expression), ][1:top, ]
top_functions$Function <- gsub("_", " ", top_functions$Function)

# --- Extract level ---
level <- sub(".*(level_[0-9]+).*", "\\1", basename(input_file))

# --- Plot ---
plot <- ggplot(top_functions, aes(x = reorder(str_wrap(Function, width = 50), Expression), y = Expression)) +
  geom_bar(stat = "identity", fill = colour_graph) +
  labs(title = paste("Top", top, "Most abundant exclusive functions in", paste0(dataset_name, "'s"), sample_type),
       x = paste("Functions grouped by", level),
       y = "Average expression") +
  theme_minimal() +
  coord_flip() +
  theme(plot.title = element_text(hjust = 0.5, size = 10),
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10)) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
            color = "black", fill = NA, linetype = "solid", size = 1)

# --- Save ---
n_functions <- nrow(top_functions)
dynamic_height <- max(3, n_functions * 0.4)
ggsave(filename = paste0("Top_", top, "_exclusive_functions_", dataset_name, "_", sample_type, "_", level, ".png"),
       plot = plot, width = 8, height = dynamic_height)

# --- Done ---
message(paste0("Functional analysis of exclusive pathways in ", dataset_name, "'s ", sample_type, " is complete."))