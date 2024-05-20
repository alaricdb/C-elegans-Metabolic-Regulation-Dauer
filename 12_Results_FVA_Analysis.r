# ====================================================================================#
# Document name: Results_FVA_Analysis.r
#
# Version: Results_FVA_Analysis.r
# Date: Apr 17, 2024
# Author: Alaric de Biolley, UM, London
# History:
#   Results_FVA_Analysis.r: Creation
# ====================================================================================#

# Load libraries 
library(readr)  
library(ggplot2)  
library(dplyr)
library(tidyr)
library(patchwork)    

# Base path 
base_path <- "/Users/alaric/Desktop/Master_Thesis_Systems_Biology/Code_local/outputs_local/Results/ASJ_neurons/FVA"

# Load various datasets
#all_cells_averaged_t0_FVA <- read_csv(paste0(base_path, "/all_neurons_averaged_t0_FVA.csv"))
all_cells_combined_KO_FVA <- read_csv(paste0(base_path, "/ASJ_neurons_combined_KO_FVA.csv"))
all_cells_imputed_scaled_count_FVA <- read_csv(paste0(base_path, "/ASJ_neurons_imputed_scaled_count_FVA.csv"))
all_cells_KO_daf_7_FVA <- read_csv(paste0(base_path, "/ASJ_neurons_KO_daf_7_FVA.csv"))
all_cells_KO_daf_9_FVA <- read_csv(paste0(base_path, "/ASJ_neurons_KO_daf_9_FVA.csv"))
all_cells_KO_ins_6_FVA <- read_csv(paste0(base_path, "/ASJ_neurons_KO_ins_6_FVA.csv"))
all_cells_KO_tph_1_FVA <- read_csv(paste0(base_path, "/ASJ_neurons_KO_tph_1_FVA.csv"))
all_cells_rna_mean_normalized_FVA <- read_csv(paste0(base_path, "/ASJ_neurons_rna_mean_normalized_FVA.csv"))

# Function to calculate differences
calculate_differences <- function(data, baseline) {
    data %>%
        inner_join(baseline, by = "ReactionID") %>%
        mutate(Diff_MinFlux = MinFlux.x - MinFlux.y,
               Diff_MaxFlux = MaxFlux.x - MaxFlux.y,
               MeanFlux = (Diff_MinFlux + Diff_MaxFlux) / 2) %>%  
        select(ReactionID, Diff_MinFlux, Diff_MaxFlux, MeanFlux)  
}

# Calculate differences with the  MeanFlux column for each condition
#diff_t0 <- calculate_differences(all_cells_averaged_t0_FVA, all_cells_imputed_scaled_count_FVA)
diff_KO_daf_7 <- calculate_differences(all_cells_KO_daf_7_FVA, all_cells_imputed_scaled_count_FVA)
diff_KO_daf_9 <- calculate_differences(all_cells_KO_daf_9_FVA, all_cells_imputed_scaled_count_FVA)
diff_KO_ins_6 <- calculate_differences(all_cells_KO_ins_6_FVA, all_cells_imputed_scaled_count_FVA)
diff_KO_tph_1 <- calculate_differences(all_cells_KO_tph_1_FVA, all_cells_imputed_scaled_count_FVA)
diff_combined_KO <- calculate_differences(all_cells_combined_KO_FVA, all_cells_imputed_scaled_count_FVA)
diff_rna_mean_normalized <- calculate_differences(all_cells_rna_mean_normalized_FVA, all_cells_imputed_scaled_count_FVA)

# Combine all difference dataframes into one with the new MeanFlux column
all_differences <- bind_rows(
    #  mutate(diff_t0, Condition = "T0"),
    mutate(diff_KO_daf_7, Condition = "KO_daf_7"),
    mutate(diff_KO_daf_9, Condition = "KO_daf_9"),
    mutate(diff_KO_ins_6, Condition = "KO_ins_6"),
    mutate(diff_KO_tph_1, Condition = "KO_tph_1"),
    mutate(diff_combined_KO, Condition = "Combined_KO"),
    mutate(diff_rna_mean_normalized, Condition = "RNA_Mean_Normalized")
)

# Melt the data for MeanFlux only
melted_data <- pivot_longer(
    all_differences,
    cols = "MeanFlux",
    names_to = "variable",
    values_to = "value"
)

# Split the data into positive and negative
positive_data <- filter(melted_data, value > 0)
negative_data <- filter(melted_data, value < 0)

# Apply log scale transformation only to positive values
positive_data$value <- log10(positive_data$value)
negative_data$value <- log10(-negative_data$value)  

# Calculate cumulative distributions for each subset
cdf_positive <- positive_data %>%
    arrange(value) %>%
    group_by(Condition) %>%
    mutate(ecdf = rank(value) / n())

cdf_negative <- negative_data %>%
    arrange(value) %>%
    group_by(Condition) %>%
    mutate(ecdf = rank(value) / n())

# Plotting the CDFs for positive MeanFlux differences
plot_positive <- ggplot(cdf_positive, aes(x = value, y = ecdf, color = Condition)) +
    geom_line() +
    scale_x_continuous(name = "Difference in Mean Flux [Log10]", labels = function(x) 10^x) +
    labs(title = "ASJ Neurons Cumulative Distribution Function for Positive MeanFlux Differences",
         y = "Cumulative Distribution") +
    theme_minimal()

# Plotting the CDFs for negative MeanFlux differences
plot_negative <- ggplot(cdf_negative, aes(x = value, y = ecdf, color = Condition)) +
    geom_line() +
    scale_x_continuous(name = "Difference in Mean Flux [Log10]", labels = function(x) -10^(-x)) +
    labs(title = "ASJ Neurons Cumulative Distribution Function for Negative MeanFlux Differences",
         y = "Cumulative Distribution") +
    theme_minimal()

# Combine plots using patchwork
combined_plot <- plot_positive / plot_negative

# Save the combined plot 
ggsave("/Users/alaric/Desktop/Master_Thesis_Systems_Biology/Code_local/outputs_local/Results/ASJ_neurons_combined_flux_differences.png", combined_plot, device = "png", width = 10, height = 8, dpi = 300)

# Additional box plot
# Pivot data for box plotting
melted_data <- pivot_longer(
    all_differences,
    cols = c("MeanFlux"),
    names_to = "variable",
    values_to = "value"
)

# Generate box plots for MeanFlux differences across all conditions
box_plot <- ggplot(melted_data, aes(x = Condition, y = value, fill = Condition)) +
    geom_boxplot(outlier.colour = "red", outlier.shape = 1) +
    facet_wrap(~ variable, scales = "free_y") +
    labs(title = "Box Plot of Flux Differences by Condition",
         x = "Condition",
         y = "Flux Difference",
         fill = "Condition") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Display the plot
print(box_plot)

# Save the box plot
ggsave("/Users/alaric/Desktop/Master_Thesis_Systems_Biology/Code_local/outputs_local/Results/all_cells_boxplot_flux_differences.png", box_plot, device = "png", width = 14, height = 8, dpi = 300)


