% Initialize COBRA Toolbox
cd ('/Users/alaric/cobratoolbox/') % to Alaric on local
initCobraToolbox(false)
changeCobraSolver('gurobi', 'all'); % Make sure the right solver is set up

%% Local Load Transcriptomic counts with their associated gene id
rna_mean_normalized = readtable('/Users/alaric/Desktop/Master_Thesis_Systems_Biology/Code_local/outputs_local/Counts_for_iMAT/ASI_neurons/ASI_neurons_KO_tph_1_altered_expression.txt', 'Delimiter', '\t');
model_iCEL1314 = readCbModel('/Users/alaric/Desktop/Master_Thesis_Systems_Biology/Code/originals/iCEL1314.xml');%% Load the iCEL1314 GEM model

% Correctly assigning variable names right after loading the table
rna_mean_normalized.Properties.VariableNames = {'gene', 'value'};

% Mapping gene expression data to model genes
[gene_id, gene_expr] = findUsedGenesLevels(model_iCEL1314, rna_mean_normalized);

%% Visualise the gene_expr distribution
figure;
% Use a large number of bins for finer resolution
nbins = 10000; 
histogram(gene_expr, nbins);
xlabel('Gene Expression Value');
ylabel('Frequency');
title('Distribution of Gene Expression Values');

%% Calculate the mean and standard deviation of gene expression for NaN values, thresholds that dont work for iMAT
gene_expr_filtered = gene_expr(gene_expr >= 0); % Filter for NaN values

% Define thresholds for low and high expression corresponding to the 25th and 75th percentiles
threshold_lb = prctile(gene_expr_filtered, 25); % Lower threshold for lowly expressed genes
threshold_ub = prctile(gene_expr_filtered, 75); % Upper threshold for highly expressed genes

% Map expression data to reactions in the metabolic model
[expressionRxns_iCEL1314, parsedGPR_iCEL1314] = mapExpressionToReactions(model_iCEL1314, rna_mean_normalized);

% Initialize the processedExpression vector with -1 for all reactions
processedExpression = -ones(length(expressionRxns_iCEL1314), 1);

% Find indices of reactions with valid (non-NaN and non-zero) expression data
validDataIndices = ~isnan(expressionRxns_iCEL1314) & expressionRxns_iCEL1314 ~= 0;

% Update only those indices in processedExpression with the actual expression data
processedExpression(validDataIndices) = expressionRxns_iCEL1314(validDataIndices);

%% Visualise data distribution with thresholds
figure;
% Large number of bins for finer resolution
nbins = 10000; 
% Plot normalized histogram to approximate density
histogram(gene_expr_filtered, nbins, 'Normalization', 'probability');
xlabel('Gene Expression Value');
ylabel('Probability Density');
title('Approximate Density Plot of Gene Expression Values');

% Add lines for thresholds
hold on; 
yLimits = ylim;
plot([threshold_lb threshold_lb], yLimits, '--r', 'LineWidth', 2); % Lower bound threshold
plot([threshold_ub threshold_ub], yLimits, '--g', 'LineWidth', 2); % Upper bound threshold
legend('Gene Expression Density', '25th Percentile', '75th Percentile');
hold off;

%% Retrieve a list of all reactions in the iCEL1314 model
all_rxns_iCEL1314 = model_iCEL1314.rxns;

% Determine reaction indices for each expression category based on thresholds
lowly_expressed_indices_iCEL1314 = find(gene_expr <= threshold_lb);
moderately_expressed_indices_iCEL1314 = find(gene_expr > threshold_lb & gene_expr < threshold_ub);
nan_or_negative_indices_iCEL1314 = find(isnan(gene_expr) | gene_expr == -1);

% Adjust the moderately expressed category to include -1 or NaN
moderately_expressed_indices_iCEL1314 = union(moderately_expressed_indices_iCEL1314, nan_or_negative_indices_iCEL1314);
highly_expressed_indices_iCEL1314 = find(gene_expr >= threshold_ub);

% Retrieve reaction IDs for each expression category
lowly_expressed_rxns_iCEL1314 = all_rxns_iCEL1314(lowly_expressed_indices_iCEL1314);
moderately_expressed_rxns_iCEL1314 = all_rxns_iCEL1314(moderately_expressed_indices_iCEL1314);
highly_expressed_rxns_iCEL1314 = all_rxns_iCEL1314(highly_expressed_indices_iCEL1314);


% Run iMAT
% Be careful to use the gurobi solver
iMAT_iCEL1314_rna_mean_normalized = iMAT(model_iCEL1314, processedExpression, threshold_lb, threshold_ub); 

% Load Pathway and Gene Data
pathway_info = readtable("/Users/alaric/Desktop/Master_Thesis_Systems_Biology/Code_local/Dataset_local/iCEL1314.xlsx")

% Preparation of the data for Fisher test
% Extract active reactions
active_rxn_ids = iMAT_iCEL1314_rna_mean_normalized.rxns;

% Find indices of these IDs in the model
active_rxn_indices = find(ismember(model_iCEL1314.rxns, active_rxn_ids));

% Retrieve reaction IDs using indices
active_rxns = model_iCEL1314.rxns(active_rxn_indices);

% Get the full list of pathways and the associated reactions from pathway_info
pathway_ids = pathway_info.ID;
pathway_categories = pathway_info.Pathway;

% Map active reactions to pathways
is_active_in_pathway = ismember(pathway_ids, active_rxns);

% Fisher test, FDR, and logFC calculation
% Unique pathways
unique_pathways = unique(pathway_categories);
p_values = zeros(length(unique_pathways), 1);
odds_ratios = zeros(length(unique_pathways), 1);
logFC_values = zeros(length(unique_pathways), 1);  % Array to hold logFC for each pathway

% Loop through each unique pathway
for i = 1:length(unique_pathways)
    % Extract reactions in the current pathway
    pathway_rxns = pathway_ids(strcmp(pathway_categories, unique_pathways{i}));
    
    % Calculate a, b, c, d for Fisher's exact test
    a = sum(ismember(pathway_rxns, active_rxns)); 
    b = length(active_rxns) - a; 
    c = length(pathway_rxns); 
    d = length(model_iCEL1314.rxns) - c; 

    % Perform Fisher's exact test
    [h, p, stats] = fishertest([a b; c d]);
    p_values(i) = p;
    odds_ratios(i) = stats.OddsRatio; % FoldChange 
    
end

% Adjust p-values for multiple testing using the Benjamini-Hochberg FDR method
adjusted_p_values = mafdr(p_values, 'BHFDR', true);

% Create a table containing pathways, raw p-values, odds ratios, logFC, and their adjusted p-values
results_table = table(unique_pathways, p_values, odds_ratios, adjusted_p_values, ...
                      'VariableNames', {'Pathway', 'PValue', 'OddsRatio', 'AdjustedPValue'});

% Display the results table
disp(results_table);

% Save the result table 
writetable(results_table, '/Users/alaric/Desktop/Master_Thesis_Systems_Biology/Code_local/outputs_local/Results/ASI_neurons/Enrichment_Pathway_Analysis/ASI_neurons_KO_tph_1_enrichment.csv');

%% FVA
[minFlux, maxFlux] = fluxVariability(iMAT_iCEL1314_rna_mean_normalized, 90);

% Create a table with reaction IDs, minimum and maximum fluxes
fva_results = table(iMAT_iCEL1314_rna_mean_normalized.rxns, minFlux, maxFlux, ...
                    'VariableNames', {'ReactionID', 'MinFlux', 'MaxFlux'});

% Save the table to a CSV file
writetable(fva_results, '/Users/alaric/Desktop/Master_Thesis_Systems_Biology/Code_local/outputs_local/Results/ASI_neurons/FVA/ASI_neurons_imputed_scaled_count_enrichment.csv');

%% Optimize the model using optimizeCbModel
printObjective(model_iCEL1314) % Printing the objective before changing: BIO0010, Biomass assembly 

model_iCEL1314 = changeObjective(model_iCEL1314, {'RCC0005'}); % Changing the objective for ATP maintenance

printObjective(model_iCEL1314) % The new objective has correctly been changed to RCC0005

optimizedSolution_iCEL1314 = optimizeCbModel(iMAT_iCEL1314_rna_mean_normalized);% Optimise the model for the objective reaction

% Find reaction indices in the context-specific model for each expression category
iCEL1314_lowly_indices = find(ismember(iMAT_iCEL1314_rna_mean_normalized.rxns, lowly_expressed_rxns_iCEL1314));
iCEL1314_moderately_indices = find(ismember(iMAT_iCEL1314_rna_mean_normalized.rxns, moderately_expressed_rxns_iCEL1314));
iCEL1314_highly_indices = find(ismember(iMAT_iCEL1314_rna_mean_normalized.rxns, highly_expressed_rxns_iCEL1314));

% Get the flux values from the solution for the corresponding reactions,
% based on the objective function, look if it is the growth rate
fluxes_iCEL1314 = optimizedSolution_iCEL1314.x;

% Calculate the number of upregulated and downregulated reactions in each expression category
upregulated_lowly_iCEL1314 = sum(fluxes_iCEL1314(iCEL1314_lowly_indices) > 0);
downregulated_lowly_iCEL1314 = sum(fluxes_iCEL1314(iCEL1314_lowly_indices) < 0);
upregulated_moderately_iCEL1314 = sum(fluxes_iCEL1314(iCEL1314_moderately_indices) > 0);
downregulated_moderately_iCEL1314 = sum(fluxes_iCEL1314(iCEL1314_moderately_indices) < 0);
upregulated_highly_iCEL1314 = sum(fluxes_iCEL1314(iCEL1314_highly_indices) > 0);
downregulated_highly_iCEL1314 = sum(fluxes_iCEL1314(iCEL1314_highly_indices) < 0);

% Print the number of upregulated and downregulated reactions in each expression category
fprintf('Lowly expressed (iCEL1314): upregulated = %d, downregulated = %d\n', upregulated_lowly_iCEL1314, downregulated_lowly_iCEL1314);
fprintf('Moderately expressed (iCEL1314): upregulated = %d, downregulated = %d\n', upregulated_moderately_iCEL1314, downregulated_moderately_iCEL1314);
fprintf('Highly expressed (iCEL1314): upregulated = %d, downregulated = %d\n', upregulated_highly_iCEL1314, downregulated_highly_iCEL1314)

% Find these quantitative fluxes through all reactions. This will give insights into the metabolic shifts for different cell types and during dauer. 
% Combine reaction IDs with their flux values into a table
rxnIDs = iMAT_iCEL1314_rna_mean_normalized.rxns;
fluxValues = fluxes_iCEL1314;
expressionCategories = cell(length(rxnIDs), 1);

% Assign expression categories to each reaction
expressionCategories(iCEL1314_lowly_indices) = {'Low'};
expressionCategories(iCEL1314_moderately_indices) = {'Moderate'};
expressionCategories(iCEL1314_highly_indices) = {'High'};

% Handle reactions that may not have been categorized due to lack of gene expression data
for i = 1:length(expressionCategories)
    if isempty(expressionCategories{i})
        expressionCategories{i} = 'Uncategorized';
    end
end

% Create the table
reactionsFluxTable = table(rxnIDs, expressionCategories, fluxValues, ...
    'VariableNames', {'ReactionID', 'ExpressionCategory', 'FluxValue'});

% Save the table to a CSV file
writetable(reactionsFluxTable, '/Users/alaric/Desktop/Master_Thesis_Systems_Biology/Code_local/outputs_local/Results/ASJ_neurons/iMAT_FluxTable_Results/all_cells_combined_KO_enrichment.csv'); % Adjust the path