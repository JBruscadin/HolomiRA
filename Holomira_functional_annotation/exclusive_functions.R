# exclusive_functions.R

# Check if the first argument is -help
if (length(commandArgs(trailingOnly = TRUE)) == 1 && commandArgs(trailingOnly = TRUE)[1] == "--help") {
  cat("Usage: Rscript exclusive_functions.R <input_file> <dataset_name> <sample_type> <top_functions> <colour_graph>\n")
  cat("\nArguments:\n")
  cat("  <input_file>          : The input file containing data for analysis.\n")
  cat("  <dataset_name>        : The name of the organism to be analyzed (e.g., Human, Ruminants) or the conditions of the same organism (e.g., decrease_RME, increase_RME).\n")
  cat("  <sample_type>         : The type of sample (e.g., Mags, Genes, miRNA).\n")
  cat("  <top_functions>       : The number of top functions to display in the results.\n")
  cat("  <colour_graph>        : The color to be used for the graph (e.g., blue, '#FF5733').\n")
  quit(status = 0)  # Exit the script with a success status
}

args <- commandArgs(trailingOnly = TRUE)

# Check if the correct number of arguments were passed
if (length(args) != 5) {
 stop("Usage: Rscript exclusive_functions.R <input_file> <dataset_name> <sample_type> <top_functions> <colour_graph>")
}

input_file <- args[1]
dataset_name <- args[2]
sample_type <- args[3]
top <- as.integer(args[4])
colour_graph <- args[5]

# Load the data
data <- read.table(input_file, header = TRUE, sep = '\t')

num_rows <- nrow(data)
# Check if there are no rows in the data
if (num_rows == 0) {
  stop(paste("There are no exclusive pathways for the organism '", dataset_name, "'.", sep = ""))
}

if (top > num_rows) {
  stop(paste("The number of top functions requested (", top, ") exceeds the number of available functions (", num_rows, "). Please select a number less than or equal to ", num_rows, ".", sep = ""))
}

# Calculate expression means
expression_means <- rowMeans(as.data.frame(data[, -1]))  
data_function <- data.frame(Function = data[, 1], Expression = expression_means)

# Select the most expressed functions
top_functions <- data_function[order(-data_function$Expression), ][1:top, ]

# Replace "_" with space
top_functions$Function <- gsub("_", " ", top_functions$Function)

# Extract level from input filename
level <- sub(".*(level_[0-9]+).*", "\\1", basename(input_file))

# Create the plot
library(ggplot2)
library(scales)
library(stringr)

ggplot(top_functions, aes(x = reorder(str_wrap(Function, width = 50), Expression), y = Expression)) +
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
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),color = "black", fill = NA, linetype = "solid", size = 1)

n_functions <- nrow(top_functions)
dynamic_height <- max(3, n_functions * 0.4)

# Save the plot
ggsave(filename = paste0("Top_" ,top, "_exclusive_functions_", dataset_name, "_", sample_type, "_", level, ".png"), width = 8, height = dynamic_height)

# Final message
message(paste0("Functional analysis of exclusive pathways ", dataset_name, "'s ", sample_type, " is complete."))
