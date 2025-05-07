#!/usr/bin/env Rscript

# --- Install missing packages ---
required_packages <- c("ggplot2", "reshape2", "dplyr", "VennDiagram", "grid", "futile.logger")
missing_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if (length(missing_packages)) {
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}

# --- Load libraries ---
suppressPackageStartupMessages({
  library(stats)
  library(ggplot2)
  library(reshape2)
  library(dplyr)
  library(VennDiagram)
  library(grid)
  library(futile.logger)
})

if (length(commandArgs(trailingOnly = TRUE)) == 1 && commandArgs(trailingOnly = TRUE)[1] == "--help") {
  cat("Usage: Rscript comparative_statistical_analysis.R <input_data1> <dataset_name1> <input_data2> <dataset_name2> <sample_type> <abundance_threshold> <prevalence_threshold> <padj> <log_threshold> <color_data1> <color_data2> <heatmap_color>\n")
  cat("\nArguments:\n")
  cat("  <input_data1>           : The input data file for the first dataset.\n")
  cat("  <dataset_name1>         : The name of first the organism to be analyzed (e.g., Human, Ruminants) or the first condition of the same organism (e.g., decrease_RME, increase_RME).\n")
  cat("  <input_data2>           : The input data file for the second dataset.\n")
  cat("  <dataset_name2>         : The name of the second organism to be analyzed (e.g., Human, Ruminants) or the second condition of the same organism (e.g., decrease_RME, increase_RME).\n")
  cat("  <sample_type>           : The type of sample (e.g., Mags, Genes, miRNA).\n")
  cat("  <abundance_threshold>   : The threshold for abundance to filter functions. Abundance threshold must be greater than or equal to zero.\n")
  cat("  <prevalence_threshold>  : The threshold for prevalence to filter functions. Prevalence threshold must be between 0 and 1.\n")
  cat("  <padj>                  : The adjusted p-value threshold for significance (should be between 0 and 1).\n")
  cat("  <log_threshold>         : The log2 fold change threshold for filtering results (can be any numeric value).\n")
  cat("  <color_data1>           : The color to be used for the first dataset in visualizations (e.g., blue, '#FF5733').\n")
  cat("  <color_data2>           : The color to be used for the second dataset in visualizations (e.g., red, '#33FF57').\n")
  cat("  <heatmap_color>         : The color palette to be used for heatmap visualization (e.g., purple, green).\n")
  quit(status = 0)  # Exit the script with a success status
}

temp_file <- file("/dev/null", open = "wt")
sink(temp_file)
suppressPackageStartupMessages(library(stats))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(dplyr))
if (!requireNamespace("VennDiagram", quietly = TRUE)) {
  install.packages("VennDiagram", quiet = TRUE)
}
suppressPackageStartupMessages(library(VennDiagram))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(futile.logger))
sink()
close(temp_file)

# Receive command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if the correct number of arguments were passed
if (length(args) != 12) {
  stop("Usage: Rscript comparative_statistical_analysis.R <input_data1> <dataset_name1> <input_data2> <dataset_name2> <sample_type> <abundance_threshold> <prevalence_threshold> <padj> <log_threshold> <color_data1> <color_data2> <heatmap_color>")
}

# Assign arguments to variables
input_data1 <- args[1]
dataset_name1 <- args[2]
input_data2 <- args[3]
dataset_name2 <- args[4]
sample_type <- args[5]
abundance_threshold <- as.numeric(args[6])
prevalence_threshold <- as.numeric(args[7])
padj <- as.numeric(args[8])
log_threshold <- as.numeric(args[9])
color_data1 <- args[10]
color_data2 <- args[11]
heatmap_color <- args[12]

# Extract level and sample type from file names
level1 <- gsub(".*level_(\\d+).*", "\\1", input_data1)
level2 <- gsub(".*level_(\\d+).*", "\\1", input_data2)
sample_type1 <- gsub(".*_(\\w+).*", "\\1", input_data1)
sample_type2 <- gsub(".*_(\\w+).*", "\\1", input_data2)

# Check if level1 is equal to level2 and sample_type1 is equal to sample_type2
if (level1 != level2) {
  stop("Classification levels (level1 and level2) must be the same. Please provide valid inputs.")
}

if (sample_type1 != sample_type2) {
  stop("Sample types (sample_type1 and sample_type2) must be the same. Please provide valid inputs.")
}

# Check if prevalence_threshold is within valid range (0-1)
if (prevalence_threshold > 1 || prevalence_threshold < 0) {
  stop("Prevalence threshold must be between 0 and 1. Please provide a valid value.")
}

# Check if padj is within valid range (0-1)
if ( padj > 1 || padj < 0) {
  stop("Padj value must be between 0 and 1. Please provide a valid value.")
}

# Read the data
data1 <- read.table(input_data1, header = TRUE, sep = ' ')
data1 <- data1[, -ncol(data1)]

data2 <- read.table(input_data2, header = TRUE, sep = ' ')
data2 <- data2[, -ncol(data2)]

num_columns_data1 <- ncol(data1)
num_columns_data2 <- ncol(data2)

# Set threshold for minimum number of columns
min_columns_threshold <- 31

# Check if the number of columns is less than the threshold
if (num_columns_data1 < min_columns_threshold & num_columns_data2 < min_columns_threshold) {
  warning(paste("The sample size for both", input_data1, "and", input_data2, "is too small (less than 30 samples). This may affect the statistical analyses."))
} else if (num_columns_data1 < min_columns_threshold) {
  warning(paste("The sample size for", input_data1, "is too small (less than 30 samples). This may affect the statistical analyses."))
} else if (num_columns_data2 < min_columns_threshold) {
  warning(paste("The sample size for", input_data2, "is too small (less than 30 samples). This may affect the statistical analyses."))
}

# Process the data
func_names_data1 <- data1[, 1]
quantifications_data1 <- data1[, -1]
expressed_data1 <- quantifications_data1 > abundance_threshold
percent_expressed_data1 <- rowMeans(expressed_data1)
selected_functions_data1 <- quantifications_data1[percent_expressed_data1 >= prevalence_threshold, ]
selected_functions_data1 <- cbind(func_names_data1[percent_expressed_data1 >= prevalence_threshold], selected_functions_data1)
colnames(selected_functions_data1)[1] <- colnames(data1)[1]

func_names_data2 <- data2[, 1]
quantifications_data2 <- data2[, -1]
expressed_data2 <- quantifications_data2 > abundance_threshold
percent_expressed_data2 <- rowMeans(expressed_data2)
selected_functions_data2 <- quantifications_data2[percent_expressed_data2 >= prevalence_threshold, ]
selected_functions_data2 <- cbind(func_names_data2[percent_expressed_data2 >= prevalence_threshold], selected_functions_data2)
colnames(selected_functions_data2)[1] <- colnames(data2)[1]

# Compare functions between the two datasets
data1_functions <- selected_functions_data1[,1]
data2_functions <- selected_functions_data2[,1]

#Venn_diagram
venn.plot <- venn.diagram(
  x = list(
    list_function_data1 = data1_functions,
    list_function_data2 = data2_functions
  ),
  category.names = c(dataset_name1, dataset_name2),
  fill = c(color_data1, color_data2),
  alpha = 0.5,
  cat.col = c(color_data1, color_data2),
  cat.cex = 1.5,
  cex = 2,
  main = paste("Venn Diagram -", dataset_name1, "vs", dataset_name2),
  main.cex = 1.5,
  scaled = FALSE,  
  lwd = 2,  
  label.col = "black",  
  fontface = "bold",
  fontfamily = "sans",
  cat.pos = c(0, 0),
  cat.dist = c(0.05, 0.05),
  cat.just = list(c(0.5, 0.5), c(0.5, 0.5)),
  filename = NULL
)

# Save Venn_diagram as .pdf and .png
pdf(file = paste0("Venn_diagram_", sample_type, "_level_", level1, "_", dataset_name1, "_vs_", dataset_name2, ".pdf"))
grid.draw(venn.plot)
dev.off()

png(filename = paste0("Venn_diagram_", sample_type, "_level_", level1, "_", dataset_name1, "_vs_", dataset_name2, ".png"))
grid.draw(venn.plot)
dev.off()

# Check if there are common functions
common_functions <- intersect(data1_functions, data2_functions)
if (length(common_functions) == 0) {
  stop("No common functions were found between the two groups analyzed.")
}

filtered_data1 <- selected_functions_data1[selected_functions_data1[,1] %in% common_functions, ]
filtered_data2 <- selected_functions_data2[selected_functions_data2[,1] %in% common_functions, ]

# Generate exclusive functions for both datasets
exclusive_data1 <- selected_functions_data1[!selected_functions_data1[,1] %in% common_functions, ]
exclusive_data2 <- selected_functions_data2[!selected_functions_data2[,1] %in% common_functions, ]

# Save filtered_data1 and filtered_data2 with common functions
write.table(filtered_data1, file = paste0("Common_functions", sample_type, "_level_", level1, "_", dataset_name1, ".txt"), sep = '\t', row.names = FALSE, quote = FALSE)
write.table(filtered_data2, file = paste0("Common_functions_", sample_type, "_level_", level1, "_", dataset_name2, ".txt"), sep = '\t', row.names = FALSE, quote = FALSE)

# Save exclusive_data1 and exclusive_data2 with exclusive functions
write.table(exclusive_data1, file = paste0("exclusive_functions_", sample_type, "_level_", level1, "_", dataset_name1, ".txt"), sep = '\t', row.names = FALSE, quote = FALSE)
write.table(exclusive_data2, file = paste0("exclusive_functions_", sample_type, "_level_", level1, "_", dataset_name2, ".txt"), sep = '\t', row.names = FALSE, quote = FALSE)

#Comparative analysis
expression_data1 <- filtered_data1[,-1]
expression_data2 <- filtered_data2[,-1]

# Initialize vectors to store results
p_values <- numeric(nrow(filtered_data1))
function_names <- filtered_data1[,1]

# Loop to perform Wilcoxon test for each function
for (i in 1:nrow(filtered_data1)) {
  vector_data1 <- expression_data1[i, ]
  vector_data2 <- expression_data2[i, ]
  vector_data1 <- sapply(vector_data1, as.numeric)
  vector_data2 <- sapply(vector_data2, as.numeric)
  p_values[i] <- wilcox.test(vector_data1, vector_data2, exact = FALSE)$p.value
}

# Apply p-value correction (FDR - Benjamini-Hochberg)
p_values_adjusted <- p.adjust(p_values, method = "fdr")

#print(p_values)
#print(p_values_adjusted)

# Select functions with adjusted p-value <= padj
significant_functions <- function_names[p_values_adjusted <= padj]

# Check if there are significant functions
if (length(significant_functions) == 0) {
  stop("No functions were found with significantly different abundance levels between the two groups. Consider adjusting the padj value.")
}

significant_data1 <- filtered_data1[p_values_adjusted <= padj, ]
significant_data2 <- filtered_data2[p_values_adjusted <= padj, ]
filtered_p_values <- p_values_adjusted[p_values_adjusted <= padj]

# Calculate mean expression per function
mean_data1 <- rowMeans(filtered_data1[,-1])
mean_data2 <- rowMeans(filtered_data2[,-1])

# Calculate the ratio between means of dataset 1 and dataset 2
ratio <- mean_data1 / mean_data2

# Apply logarithmic transformation base 2
log2_ratio <- log2(ratio)


# Select functions with log2_ratio greater than log_threshold in absolute value
selected_functions_log <- log2_ratio[abs(log2_ratio) >= log_threshold]

# Check if there are selected functions after log threshold
if (length(selected_functions_log) == 0) {
  stop("No functions were found with comparatively different abundance levels between the two groups. Consider adjusting the log_threshold value.")
}

# Volcano Plot
volcano_data <- data.frame(
  Function = common_functions,
  log2Ratio = log2_ratio,
  Pvalue = p_values_adjusted,
  Mean_data1 = mean_data1,
  Mean_data2 = mean_data2
)

volcano_data$Function <- gsub("_", " ", volcano_data$Function)

volcano_data$color <- "darkgray" 
volcano_data$color[volcano_data$Pvalue <= padj & volcano_data$log2Ratio >= log_threshold] <- color_data1
volcano_data$color[volcano_data$Pvalue <= padj & volcano_data$log2Ratio <= -log_threshold] <- color_data2

# Create the Volcano Plot
volcano_plot <- ggplot(volcano_data, aes(x = log2Ratio, y = -log10(Pvalue), color = color)) +
  geom_point() +
  scale_color_identity() +
  geom_hline(yintercept = -log10(padj), linetype = "dashed", color = "black") +
  geom_vline(xintercept = log_threshold, linetype = "dashed", color = "black") +
  geom_vline(xintercept = -log_threshold, linetype = "dashed", color = "black") +
  labs(title = paste("Volcano plot - ", sample_type, " (", dataset_name1, " vs ", dataset_name2, ")", sep = ""),
       x = "log2 Fold Change",
       y = "-log10(P-value)") +
  theme_minimal() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 10),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10)) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), color = "black", fill = NA, linetype = "solid", size = 1)

# Save Volcano Plot as .pdf and .png
ggsave(filename = paste0("Volcano_Plot_", sample_type, "_level_", level1, "_", dataset_name1, "_vs_", dataset_name2, ".pdf"), plot = volcano_plot, width = 10, height = 6)
ggsave(filename = paste0("Volcano_Plot_", sample_type, "_level_", level1, "_", dataset_name1, "_vs_", dataset_name2, ".png"), plot = volcano_plot, width = 10, height = 6)

if (all(volcano_data$color == "darkgray")) {
  stop(paste("There are no statistically significant functions for the given adjusted p-value (Padj = ", padj, ") and log2 fold change (log_threshold = ", log_threshold, "). Please consider adjusting these values and rerun the analysis.",sep = ""))
}

volcano_data_filtered <- volcano_data[volcano_data$color != "darkgray", ]
volcano_data_filtered$organism <- ifelse(volcano_data_filtered$log2Ratio > 0, dataset_name1, dataset_name2)
color_values <- setNames(c(color_data1, color_data2), c(dataset_name1, dataset_name2))
volcano_data_filtered$Function <- factor(volcano_data_filtered$Function, levels = volcano_data_filtered$Function[order(volcano_data_filtered$log2Ratio)])

# Create the bar plot

# Create the bar plot (horizontal bars)
functions_signif <- ggplot(volcano_data_filtered, aes(x = Function, y = log2Ratio, fill = organism)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = paste("Differentially abundant functions - ", sample_type, " (", dataset_name1, " vs ", dataset_name2, ")", sep = ""), 
       x = paste("Functions grouped by level", level1),
         "Function", y = "log2 Fold Change") +
  theme(panel.grid = element_blank(),
        legend.position = "right",
        legend.title.align = 0.5,
        plot.title = element_text(hjust = 0.5, size = 10),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        legend.background = element_rect(color = "black", fill = NA, size = 0.5)) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), color = "black", fill = NA, linetype = "solid", size = 1) +
  scale_fill_manual(values = color_values, name = "Most Abundant In") +
  coord_flip()

n_functions <- nrow(volcano_data_filtered)
dynamic_height <- max(6, n_functions * 0.4)

# Save the bar plot
ggsave(filename = paste0("Functions_significant_", sample_type, "_level_", level1, "_", dataset_name1, "_vs_", dataset_name2, ".pdf"), plot = functions_signif, width = 8, height = dynamic_height)
ggsave(filename = paste0("Functions_significant_", sample_type, "_level_", level1, "_", dataset_name1, "_vs_", dataset_name2, ".png"), plot = functions_signif, width = 8, height = dynamic_height)

# Prepare the heatmap data
volcano_data_filtered <- volcano_data_filtered %>% rename(Data1 = Mean_data1, Data2 = Mean_data2)
heatmap_data <- volcano_data_filtered[, c("Function", "Data1", "Data2")]
colnames(heatmap_data) <- c("Function", dataset_name1, dataset_name2)
heatmap_data_long <- melt(heatmap_data, id.vars = "Function", variable.name = "Dataset", value.name = "Mean")

# Create the heatmap
heatmap_plot <- ggplot(heatmap_data_long, aes(x = Function, y = Dataset, fill = Mean)) +
  geom_tile() +
  coord_flip() +
  scale_fill_gradient2(low = "blue", high = heatmap_color, mid = "white", midpoint = 0,
                       limit = c(min(heatmap_data_long$Mean), max(heatmap_data_long$Mean)), name = "Abundance Values") +
  theme_minimal() +
  labs(title = paste("Heatmap of average abundance values - ", sample_type, " (", dataset_name1, " vs ", dataset_name2, ")", sep = ""),
       x = paste("Functions grouped by level", level1)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10),
        panel.grid = element_blank(),
        legend.position = "right",
        legend.title.align = 0.5,
        legend.title = element_text(hjust = 0.5),
        plot.title = element_text(hjust = 0.5, size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.title.x = element_blank(),
        legend.background = element_rect(color = "black", fill = NA, size = 0.5)) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), color = "black", fill = NA, linetype = "solid", size = 1)

# Save the heatmap
ggsave(filename = paste0("Heatmap_abundance_", sample_type, "_level_", level1, "_", dataset_name1, "_vs_", dataset_name2, ".pdf"), plot = heatmap_plot, width = 8, height = dynamic_height)
ggsave(filename = paste0("Heatmap_abundance_", sample_type, "_level_", level1, "_", dataset_name1, "_vs_", dataset_name2, ".png"), plot = heatmap_plot, width = 8, height = dynamic_height)

# Final message
message(paste0("Comparative analysis between ", dataset_name1, " and ", dataset_name2, " for ", sample_type, " resulted in statistically significant findings for the specified padj (", padj, ") and log_threshold (", log_threshold, ") values. Graphs have been generated and saved."))
