# analyse_by_dataset.R

# Check if the first argument is -help
if (length(commandArgs(trailingOnly = TRUE)) == 1 && commandArgs(trailingOnly = TRUE)[1] == "--help") {
  cat("Usage: Rscript analyse_by_dataset.R <input_file> <dataset_name> <sample_type> <abundance_threshold> <prevalence_threshold> <top_functions> <colour_graph>\n")
  cat("\nArguments:\n")
  cat("  <input_file>           : The input file containing data for analysis.\n")
  cat("  <dataset_name>         : The name of the organism to be analyzed (e.g., Human, Ruminants) or the conditions of the same organism (e.g., decrease_RME, increase_RME).\n")
  cat("  <sample_type>          : The type of sample (e.g., Mags, Genes, miRNA).\n")
  cat("  <abundance_threshold>  : The threshold for abundance to filter functions. Abundance threshold must be greater than or equal to zero. \n")
  cat("  <prevalence_threshold> : The threshold for prevalence to filter functions. Prevalence threshold must be between 0 and 1.\n")
  cat("  <top_functions>        : The number of top functions to display in the results.\n")
  cat("  <colour_graph>         : The color to be used for the graph (e.g., blue, '#FF5733').\n")
  quit(status = 0)  # Exit the script with a success status
}

args <- commandArgs(trailingOnly = TRUE)

# Check if the correct number of arguments were passed
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

# Check if prevalence_threshold is within valid range (0-1)
if (prevalence_threshold > 1 || prevalence_threshold < 0) {
  stop("Error: Prevalence threshold must be between 0 and 1. Please provide a valid value.")
}

# Check if the abundance threshold is valid
if (abundance_threshold < 0) {
  stop("Error: The abundance threshold must be greater than or equal to zero. Please provide a valid value.")
}

# Load the data
data <- read.table(input_file, header = TRUE, sep = ' ')
data <- data[, -ncol(data)]

func_names <- data[, 1]
quantifications <- data[, -1]
expressed <- quantifications > abundance_threshold
percent_expressed <- rowMeans(expressed)
selected_functions <- quantifications[percent_expressed >= prevalence_threshold, ]
selected_functions <- cbind(func_names[percent_expressed >= prevalence_threshold], selected_functions)
colnames(selected_functions)[1] <- colnames(data)[1]

# Calculate expression means
expression_means <- rowMeans(as.data.frame(selected_functions[, -1]))  
selected_functions_means <- data.frame(Function = selected_functions[, 1], Expression = expression_means)

num_rows <- nrow(selected_functions_means)

if (top > num_rows) {
  stop(paste("Error: The number of top functions requested (", top, ") exceeds the number of available functions (", num_rows, "). Please select a number less than or equal to ", num_rows, ".", sep = ""))
}

# Select the most expressed functions
top_functions <- selected_functions_means[order(-selected_functions_means$Expression), ][1:top, ]

# Replace "_" with space
top_functions$Function <- gsub("_", " ", top_functions$Function)

# Extract level from input filename
level <- sub("^(level_[0-9]+).*", "\\1", basename(input_file))

# Create the plot
library(ggplot2)
library(scales)

ggplot(top_functions, aes(x = reorder(Function, Expression), y = Expression)) +
  geom_bar(stat = "identity", fill = colour_graph) +
  labs(title = paste("Top", top, "Most abundant functions in", paste0(dataset_name, "'s"), sample_type),
       x = paste("Functions grouped by", level), 
       y = "Average abundance") +
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

# Save the plot
ggsave(filename = paste0("Top_" ,top, "_functions_", dataset_name, "_", sample_type, "_", level, ".png"), width = 8, height = 8)

# Final message
message(paste0("Functional analysis of ", dataset_name, "'s ", sample_type, " is complete."))
