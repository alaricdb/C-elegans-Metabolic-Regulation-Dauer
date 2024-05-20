# ====================================================================================#
# Document name: Result_Enrichment_Analysis.R
#
# Version: Result_Enrichment_Analysis.R
# Date: Apr 16, 2024
# Author: Alaric de Biolley, UM, London
# History:
#   Result_Enrichment_Analysis.R: Creation
# ====================================================================================#

# Loading librairies
library(readr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)

# Function to load datasets based on cell type with consistent naming and file structure
load_datasets <- function(cell_type_prefix, cell_type_directory) {
    file_suffixes <- c(
        "averaged_t0_enrichment.csv",
        "combined_KO_enrichment.csv",
        "imputed_scaled_count_enrichment.csv",
        "KO_daf_7_enrichment.csv",
        "KO_daf_9_enrichment.csv",
        "KO_ins_6_enrichment.csv",
        "KO_tph_1_enrichment.csv"
    )
    
    # Construct full paths and read data
    file_names <- paste(cell_type_prefix, file_suffixes, sep = "_")
    paths <- paste("/Users/alaric/Desktop/Master_Thesis_Systems_Biology/Code_local/outputs_local/Results/", cell_type_directory, "/Enrichment_Pathway_Analysis/", file_names, sep = "")
    datasets <- list()
    loaded_names <- list()
    
    for (i in seq_along(paths)) {
        try({
            dataset <- read_csv(paths[i])
            dataset_name <- sub(".csv", "", file_names[i])
            datasets[[dataset_name]] <- dataset
            loaded_names <- c(loaded_names, dataset_name)
        }, silent = TRUE)  # Ignore errors and continue
    }
    
    names(datasets) <- loaded_names  # Naming datasets for successfully loaded files
    return(datasets)
}

# Loading datasets for different cell types using the corrected naming conventions
AllCells <- load_datasets("all_cells", "All_cells")
AllNeurons <- load_datasets("all_neurons", "All_neurons")
SensoryNeurons <- load_datasets("sensory_neurons", "Sensory_neurons")
ASINeurons <- load_datasets("ASI_neurons", "ASI_neurons")
ASJNeurons <- load_datasets("ASJ_neurons", "ASJ_neurons") 

# Function to summarize and visualize enrichment data
analyze_enrichment <- function(datasets) {
    results <- list()
    
    for (name in names(datasets)) {
        dataset <- datasets[[name]]
        top_pathways <- dataset %>%
            arrange(PValue) %>%
            slice_head(n = 5) %>%
            select(Pathway, PValue, AdjustedPValue, OddsRatio) %>%
            mutate(Layer = name)
        
        results[[name]] <- list(
            TopPathways = top_pathways,
            Data = dataset
        )
    }
    
    combined_top_pathways <- bind_rows(lapply(results, function(x) x$TopPathways))
    
    return(list(AllResults = results, TopPathwaysTable = combined_top_pathways))
}

# Applying the analysis function to each group of datasets
all_results <- list(
    AllCells = analyze_enrichment(AllCells),
    AllNeurons = analyze_enrichment(AllNeurons),
    SensoryNeurons = analyze_enrichment(SensoryNeurons),
    ASINeurons = analyze_enrichment(ASINeurons),
    ASJNeurons = analyze_enrichment(ASJNeurons)
)

# Extracting the combined table of top pathways
top_pathways_table <- do.call(rbind, lapply(all_results, function(x) x$TopPathwaysTable))

# Save the table without including the row names
write.csv(top_pathways_table, "/Users/alaric/Desktop/Master_Thesis_Systems_Biology/Code_local/outputs_local/Results/Enrichment_top_pathways_table.csv", row.names = FALSE)

# Adjusting the plot data with modified layer names and adding facet wrap
plot_data <- top_pathways_table %>%
    mutate(
        AdjustedPValueLog = -log10(AdjustedPValue),
        Pathway = factor(Pathway, levels = unique(Pathway)),
        CellTypePrefix = gsub("^([A-Za-z]+_[A-Za-z]+).*", "\\1", Layer),  # Extract the first two words of Layer
        Layer = factor(gsub("enrichment", "", Layer), levels = unique(gsub("enrichment", "", Layer)))
    ) %>%
    filter(
        !Pathway %in% c("Bacterial digestion", "Valine, leucine and isoleucine degradation")
    )

# Creating the dot plot using ggplot2 and adding facet wrap based on CellTypePrefix
dot_plot <- ggplot(plot_data, aes(x = Layer, y = Pathway, size = OddsRatio, color = AdjustedPValueLog)) +
    geom_point(alpha = 0.7) +
    scale_size_continuous(range = c(3, 10), name = "Odds Ratio") +
    scale_color_viridis_c(option = "C", name = "Adjusted P-Value (-log10)") +
    labs(
        title = "Dot Plot of Pathway Enrichment Analysis",
        x = "Layer",
        y = "Pathway"
    ) +
    theme_minimal() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)
    ) +
    facet_wrap(~ CellTypePrefix, scales = "free_x")

# Print the plot
print(dot_plot)

# Save the plot
ggsave("/Users/alaric/Desktop/Master_Thesis_Systems_Biology/Code_local/outputs_local/Results/pathway_enrichment_dot_plot_faceted_by_cell_type.png", plot = dot_plot, width = 20, height = 10, units = "in")
