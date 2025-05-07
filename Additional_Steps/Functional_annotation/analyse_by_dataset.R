# analyse_by_dataset.R

# Auto-install required packages if missing
packages_needed <- c("ggplot2", "scales")
packages_missing <- packages_needed[!(packages_needed %in% installed.packages()[,"Package"])]
if (length(packages_missing)) {
  install.packages(packages_missing, repos = "https://cloud.r-project.org")
}

# Load packages
suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
})

# Help message
if (length(commandArgs(trailingOnly = TRUE)) == 1 && commandArgs(trailingOnly = TRUE)[1] == "--help") {
  cat("Usage: Rscript analyse_by_dataset.R <input_file> <dataset_name> <sample_type> <abundance_threshold> <prevalence_threshold> <top_functions> <colour_graph>\n")
  cat("\nArguments:\n")
  cat("  <input_file>           : The input file containing data for analysis.\n")
  cat("  <dataset_name>         : The name of the organism or condition (e.g., Human, decrease_RME).\n")
  cat("  <sample_type>          : The type of sample (e.g., Mags, Genes, miRNA).\n")
  cat("  <abundance_threshold>  : Minimum abundance to include a function.\n")
  cat("  <prevalence_threshold> : Minimum fraction of samples expressing the function (0-1).\n")
  cat("  <top_functions>        : Number of top functions to plot.\n")
  cat("  <colour_graph>         : Color for the bars (e.g., 'blue', '#FF5733').\n")
  quit(status = 0)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7) {
  stop("Usage: Rscript analyse_by_dataset.R <input_file> <dataset_name> <sample_type> <abundance_threshold> <prevalence_threshold> <top_functions> <colour_graph>")
}

input_file <- args[1]
dataset_name <- args[2]
sample_type <- args[3]
abundance_threshold <- as.numeric(args[4])
prevalence_threshold <- as.numeric(args[5])
top <- as.integer(args[6])
colour_graph <- args[7]

if (prevalence_threshold > 1 || prevalence_threshold < 0) {
  stop("Prevalence threshold must be between 0 and 1.")
}
if (abundance_threshold < 0) {
  stop("Abundance threshold must be >= 0.")
}

# Load data
data <- read.table(input_file, header = TRUE, sep = ' ')
data <- data[, -ncol(data)]  # remove last column

func_names <- data[, 1]
quantifications <- data[, -1]
expressed <- quantifications > abundance_threshold
percent_expressed <- rowMeans(expressed)
selected_functions <- quantifications[percent_expressed >= prevalence_threshold, ]
selected_functions <- cbind(func_names[percent_expressed >= prevalence_threshold], selected_functions)
colnames(selected_functions)[1] <- colnames(data)[1]

expression_means <- rowMeans(as.data.frame(selected_functions[, -1]))
selected_functions_means <- data.frame(Function = selected_functions[, 1], Expression = expression_means)

num_rows <- nrow(selected_functions_means)
if (top > num_rows) {
  stop(paste("Requested", top, "functions, but only", num_rows, "are available."))
}

top_functions <- selected_functions_means[order(-selected_functions_means$Expression), ][1:top, ]
top_functions$Function <- gsub("_", " ", top_functions$Function)
level <- sub("^(level_[0-9]+).*", "\\1", basename(input_file))

# Plot
plot <- ggplot(top_functions, aes(x = reorder(Function, Expression), y = Expression)) +
  geom_bar(stat = "identity", fill = colour_graph) +
  labs(title = paste("Top", top, "Functions in", dataset_name, sample_type),
       x = paste("Functions grouped by", level),
       y = "Average abundance") +
  theme_minimal() +
  coord_flip() +
  theme(plot.title = element_text(hjust = 0.5, size = 10),
        panel.grid = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10)) + 
        geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
          color = "black", fill = NA, linetype = "solid", linewidth = 1)


# Save
ggsave(filename = paste0("Top_", top, "_functions_", dataset_name, "_", sample_type, "_", level, ".png"),
       plot = plot, width = 8, height = 8)

message(paste("Functional analysis of", dataset_name, sample_type, "completed."))