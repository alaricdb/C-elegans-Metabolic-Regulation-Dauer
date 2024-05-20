# ====================================================================================#
# Document name: Result_Analysis.R
#
# Version: Result_Analysis.R
# Date: Mar 27, 2024
# Author: Alaric de Biolley, King's College of London, London, UK
# History:
#   Result_Analysis.R: Creation
# ====================================================================================#

# Loading libraries
library(dplyr)
library(tidyr)
library(readr) # For csv files, use read_csv(), which is faster and more reliable than read.csv()
library(ggplot2)
library(pheatmap)
library(tibble)

# Load all the results
combined_KO_reactions_flux_table <- read_csv("/Users/alaric/Desktop/Master_Thesis_Systems_Biology/Code_local/outputs_local/Results/ASJ_neurons/iMAT_FluxTable_Results/ASJ_neurons_combined_KO_altered_expression_flux_table.csv")
imputed_reactions_flux_table <- read_csv("/Users/alaric/Desktop/Master_Thesis_Systems_Biology/Code_local/outputs_local/Results/ASJ_neurons/iMAT_FluxTable_Results/ASJ_neurons_imputed_scaled_count_flux_table.csv")
daf_7_KO_reactions_flux_table <- read_csv("/Users/alaric/Desktop/Master_Thesis_Systems_Biology/Code_local/outputs_local/Results/ASJ_neurons/iMAT_FluxTable_Results/ASJ_neurons_KO_daf_7_altered_expression_flux_table.csv")
daf_9_KO_reactions_flux_table <- read_csv("/Users/alaric/Desktop/Master_Thesis_Systems_Biology/Code_local/outputs_local/Results/ASJ_neurons/iMAT_FluxTable_Results/ASJ_neurons_KO_daf_9_altered_expression_flux_table.csv")
ins_6_KO_reactions_flux_table <- read_csv("/Users/alaric/Desktop/Master_Thesis_Systems_Biology/Code_local/outputs_local/Results/ASJ_neurons/iMAT_FluxTable_Results/ASJ_neurons_KO_ins_6_altered_expression_flux_table.csv")
tph_1_KO_reactions_flux_table <- read_csv("/Users/alaric/Desktop/Master_Thesis_Systems_Biology/Code_local/outputs_local/Results/ASJ_neurons/iMAT_FluxTable_Results/ASJ_neurons_KO_tph_1_altered_expression_flux_table.csv")
#t0_experimental_flux_table <- read_csv("/Users/alaric/Desktop/Master_Thesis_Systems_Biology/Code_local/outputs_local/Results/All_cells/iMAT_FluxTable_Results/t0_experimental_flux_table.csv")

# ====================================================================================#
# 1. Identification of Common and Unique reactions
# ====================================================================================#
# Combine all tables into one and add a source column to identify each
all_reactions <- bind_rows(
    combined_KO_reactions_flux_table %>% mutate(Source = "Combined KO"),
    daf_7_KO_reactions_flux_table %>% mutate(Source = "DAF-7 KO"),
    daf_9_KO_reactions_flux_table %>% mutate(Source = "DAF-9 KO"),
    imputed_reactions_flux_table %>% mutate(Source = "Imputed"),
    ins_6_KO_reactions_flux_table %>% mutate(Source = "INS-6 KO"),
    tph_1_KO_reactions_flux_table %>% mutate(Source = "TPH-1 KO"),
    #  t0_experimental_flux_table %>% mutate(Source = "T0 Experimental")
)

# Identify ALL unique and common reactions
# For common reactions
all_common_reactions <- all_reactions %>%
    group_by(ReactionID) %>%
    summarise(Count = n()) %>%
    filter(Count == n_distinct(all_reactions$Source))

# For unique reactions
all_unique_reactions <- all_reactions %>%
    group_by(ReactionID, Source) %>%
    summarise(Count = n(), .groups = 'drop') %>%
    group_by(ReactionID) %>%
    filter(n() == 1)

# ====================================================================================#
# 2. Filtering each dataset to include only common reactions
# ====================================================================================#
common_reaction_ids <- all_common_reactions$ReactionID

# Filter each table to include only common reactions
combined_KO_filtered <- combined_KO_reactions_flux_table %>% filter(ReactionID %in% common_reaction_ids)
daf_7_KO_filtered <- daf_7_KO_reactions_flux_table %>% filter(ReactionID %in% common_reaction_ids)
daf_9_KO_filtered <- daf_9_KO_reactions_flux_table %>% filter(ReactionID %in% common_reaction_ids)
imputed_filtered <- imputed_reactions_flux_table %>% filter(ReactionID %in% common_reaction_ids)
ins_6_KO_filtered <- ins_6_KO_reactions_flux_table %>% filter(ReactionID %in% common_reaction_ids)
tph_1_KO_filtered <- tph_1_KO_reactions_flux_table %>% filter(ReactionID %in% common_reaction_ids)
#t0_experimental_filtered <- t0_experimental_flux_table %>% filter(ReactionID %in% common_reaction_ids)

# ====================================================================================#
# 3. Preparing data for the heatmap and plot 
# ====================================================================================#
# Combine filtered datasets again, adding a source label to each
filtered_combined <- bind_rows(
    combined_KO_filtered %>% mutate(Source = "Combined KO"),
    daf_7_KO_filtered %>% mutate(Source = "DAF-7 KO"),
    daf_9_KO_filtered %>% mutate(Source = "DAF-9 KO"),
    imputed_filtered %>% mutate(Source = "Imputed"),
    ins_6_KO_filtered %>% mutate(Source = "INS-6 KO"),
    tph_1_KO_filtered %>% mutate(Source = "TPH-1 KO"),
    #  t0_experimental_filtered %>% mutate(Source = "T0 Experimental")
) %>%
    select(-ExpressionCategory) # Exclude the ExpressionCategory column

# Summarize the data to ensure no duplicates and to compute the mean FluxValue for each Source and ReactionID
heatmap_data <- filtered_combined %>%
    pivot_wider(names_from = Source, values_from = FluxValue, values_fn = list(FluxValue = mean), values_fill = list(FluxValue = 0)) %>%
    group_by(ReactionID) %>%
    summarize(across(everything(), mean, na.rm = TRUE)) %>%
    ungroup()

# Identify unique columns
unique_data <- function(data) {
    cols <- colnames(data)
    col_groups <- list()
    
    # Compare each column with every other column
    for (i in 1:(ncol(data) - 1)) {
        for (j in (i + 1):ncol(data)) {
            if (identical(data[, i], data[, j])) {
                # Group columns that are identical
                found <- FALSE
                for (group in seq_along(col_groups)) {
                    if (cols[i] %in% col_groups[[group]] || cols[j] %in% col_groups[[group]]) {
                        col_groups[[group]] <- unique(c(col_groups[[group]], cols[i], cols[j]))
                        found <- TRUE
                        break
                    }
                }
                if (!found) {
                    col_groups[[length(col_groups) + 1]] <- c(cols[i], cols[j])
                }
            }
        }
    }
    
    # Create a new data frame with merged columns
    new_data <- data.frame(ReactionID = data$ReactionID)
    
    # Process groups of columns
    for (group in col_groups) {
        group_name <- paste(group, collapse=" & ")
        new_data[[group_name]] <- rowMeans(data[group], na.rm = TRUE)
    }
    
    # Add columns that were not part of any group
    non_grouped_cols <- setdiff(cols, unlist(col_groups))
    new_data[non_grouped_cols] <- data[non_grouped_cols]
    
    new_data
}

# Apply the unique column identification function and check column names
heatmap_data <- unique_data(heatmap_data)

# Convert to matrix, ensuring the first column (ReactionID) is used for row names
heatmap_matrix <- as.matrix(heatmap_data[-1])
rownames(heatmap_matrix) <- heatmap_data$ReactionID

# Specify the desired order of columns, ensuring to omit duplicates now
ordered_cols <- c("DAF-7 KO & DAF-9 KO & Imputed & INS-6 KO", "Combined KO & TPH-1 KO")
ordered_cols <- ordered_cols[ordered_cols %in% colnames(heatmap_matrix)]

# Reorder the columns in the matrix according to the specified order
heatmap_matrix <- heatmap_matrix[, ordered_cols]

# Filter out rows where all values are zero
heatmap_matrix <- heatmap_matrix[rowSums(heatmap_matrix != 0) > 0, ]

# Filter out rows that have no variance across all conditions
heatmap_matrix_clean <- heatmap_matrix[apply(heatmap_matrix, 1, function(x) var(x) > 0), ]

# Heatmap plot with all 
png(filename = "/Users/alaric/Desktop/Master_Thesis_Systems_Biology/Code_local/outputs_local/Results/All_neurons/figures/all_reactions_all_neurons_Flux_Variations2_Heatmap.png", width = 6400, height = 4800, res = 4800)
pheatmap(heatmap_matrix_clean, 
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         cluster_cols = FALSE,  # Prevent column clustering
         scale = "row",  # Scale rows to z-scores
         color = colorRampPalette(c("blue", "white", "red"))(255),
         show_rownames = TRUE,
         show_colnames = TRUE,
         main = "Flux Heatmap of Common Reactions Across Conditions",
         fontsize_row = 12,
         fontsize_col = 12,
         margins = c(10, 10),  # Adjusted margins
         cex = 0.1)
dev.off()

# Heatmap visually comprehensive 
png(filename = "/Users/alaric/Desktop/Master_Thesis_Systems_Biology/Code_local/outputs_local/Results/All_neurons/figures/visual_All__Flux_Variations_Heatmap_2.png", width = 1600, height = 1200, res = 300)
pheatmap(heatmap_matrix_clean,  # Use the ordered matrix
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         cluster_cols = FALSE,  # Disables column clustering, respects manual order
         scale = "row",
         color = colorRampPalette(c("blue", "white", "red"))(255),
         show_rownames = FALSE,
         show_colnames = TRUE,
         fontsize_col = 6,
         main = "Heatmap of Reactions Flux for All Cells")
dev.off()