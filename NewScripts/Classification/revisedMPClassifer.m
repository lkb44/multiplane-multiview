clear;close all;clc;

if isunix
    addpath(genpath('/Users/leobao/Documents/MultiPlanePipeline/scripts/general/'));
    addpath(genpath('/Users/leobao/Documents/MultiPlanePipeline/scripts/classification/nFold_cross_validation/'));
    addpath(genpath('/Users/leobao/Documents/MultiPlanePipeline/scripts/classification/training_and_testing_sets/'));
    addpath(genpath('/Users/leobao/Documents/MultiPlanePipeline/scripts/classification/classifiers/'));
    addpath(genpath('/Users/leobao/Documents/MultiPlanePipeline/scripts/feature_selection/mrmr_feature_select/'));
    addpath(genpath('/Users/leobao/Documents/MultiPlanePipeline/scripts/'));
    addpath(genpath('/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/Feature_Selection_Classification_Code/Feature_Classifier/'));
    addpath(genpath('/Users/leobao/Documents/MultiPlanePipeline/Scripts-DataPreProcessing/Util-preProcessing'));
    addpath(genpath('/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/2D_Shape_Feature_Extraction_Code/scripts/feature_selection/mrmr_feature_select/estpab.mexmac'));
    % addpath(genpath('/Users/leobao/Documents/MultiPlanePipeline/'));
    % addpath(genpath('/Users/leobao/Documents/MultiPlanePipeline/scripts/feature_selection/mrmr_feature_select/mrmr_mid_d.m'))
end
region = {'MPMR'}; % Can be 'best_fat' or 'best_tumor'
split = {'missingCollage'};
scheme = {'mrmr_qda'};

results_root = '/Volumes/Crucial X6/data copy/Data/MissingCollageResults/';
experiments = dir(fullfile(results_root));
theFilesvol = {experiments.name};

feature_column_path = '/Users/leobao/Documents/MultiPlanePipeline/Feature_Names/Texture_Feature_Names.xlsx';
feature_column_names = readtable(feature_column_path, 'ReadVariableNames',false);
feature_column_names = table2cell(feature_column_names);

axial_proxfat5_feature_column_names = cell(1, numel(feature_column_names));
coronal_proxfat5_feature_column_names = cell(1, numel(feature_column_names));
axial_proxfat10_feature_column_names = cell(1, numel(feature_column_names));
coronal_proxfat10_feature_column_names = cell(1, numel(feature_column_names));
axial_tumor_feature_column_names = cell(1, numel(feature_column_names));
coronal_tumor_feature_column_names = cell(1, numel(feature_column_names));


% Loop through each cell in feature_column_names and add the prefix
for i = 1:numel(feature_column_names)
    axial_tumor_feature_column_names{i} = ['axial_tumor_', feature_column_names{i}];
    coronal_tumor_feature_column_names{i} = ['coronal_tumor_', feature_column_names{i}];
    axial_proxfat5_feature_column_names{i} = ['axial_proxfat5_', feature_column_names{i}];
    coronal_proxfat5_feature_column_names{i} = ['coronal_proxfat5_', feature_column_names{i}];
    axial_proxfat10_feature_column_names{i} = ['axial_proxfat10_', feature_column_names{i}];
    coronal_proxfat10_feature_column_names{i} = ['coronal_proxfat10_', feature_column_names{i}];
end

for k = 1:length(theFilesvol)
    P = theFilesvol;
    P = P(~startsWith(P, '.'));
    axial_proxfat5_wilcoxon_qda_index = find(contains(P, 'Axial_proxfat5_only_wilcoxon_qda'));
    axial_proxfat5_wilcoxon_lda_index = find(contains(P, 'Axial_proxfat5_only_wilcoxon_lda'));
    axial_proxfat5_wilcoxon_rf_index = find(contains(P, 'Axial_proxfat5_only_wilcoxon_rf'));
    % axial_proxfat5_mrmr_qda_index = find(contains(P, 'Axial_proxfat5_only_mrmr_qda'));
    axial_proxfat5_mrmr_lda_index = find(contains(P, 'Axial_proxfat5_only_mrmr_lda'));
    axial_proxfat5_mrmr_rf_index = find(contains(P, 'Axial_proxfat5_only_mrmr_rf'));
    coronal_proxfat5_wilcoxon_qda_index = find(contains(P, 'Coronal_proxfat5_only_wilcoxon_qda'));
    coronal_proxfat5_wilcoxon_lda_index = find(contains(P, 'Coronal_proxfat5_only_wilcoxon_lda'));
    coronal_proxfat5_wilcoxon_rf_index = find(contains(P, 'Coronal_proxfat5_only_wilcoxon_rf'));
    coronal_proxfat5_mrmr_qda_index = find(contains(P, 'Coronal_proxfat5_only_mrmr_qda'));
    coronal_proxfat5_mrmr_lda_index = find(contains(P, 'Coronal_proxfat5_only_mrmr_lda'));
    coronal_proxfat5_mrmr_rf_index = find(contains(P, 'Coronal_proxfat5_only_mrmr_rf'));
    axial_proxfat10_wilcoxon_qda_index = find(contains(P, 'Axial_proxfat10_only_wilcoxon_qda'));
    axial_proxfat10_wilcoxon_lda_index = find(contains(P, 'Axial_proxfat10_only_wilcoxon_lda'));
    axial_proxfat10_wilcoxon_rf_index = find(contains(P, 'Axial_proxfat10_only_wilcoxon_rf'));
    axial_proxfat10_mrmr_qda_index = find(contains(P, 'Axial_proxfat10_only_mrmr_qda'));
    axial_proxfat10_mrmr_lda_index = find(contains(P, 'Axial_proxfat10_only_mrmr_lda'));
    axial_proxfat10_mrmr_rf_index = find(contains(P, 'Axial_proxfat10_only_mrmr_rf'));
    coronal_proxfat10_wilcoxon_qda_index = find(contains(P, 'Coronal_proxfat10_only_wilcoxon_qda'));
    coronal_proxfat10_wilcoxon_lda_index = find(contains(P, 'Coronal_proxfat10_only_wilcoxon_lda'));
    coronal_proxfat10_wilcoxon_rf_index = find(contains(P, 'Coronal_proxfat10_only_wilcoxon_rf'));
    coronal_proxfat10_mrmr_qda_index = find(contains(P, 'Coronal_proxfat10_only_mrmr_qda'));
    coronal_proxfat10_mrmr_lda_index = find(contains(P, 'Coronal_proxfat10_only_mrmr_lda'));
    coronal_proxfat10_mrmr_rf_index = find(contains(P, 'Coronal_proxfat10_only_mrmr_rf'));
    axial_tumor_wilcoxon_qda_index = find(contains(P, 'Axial_tumor_only_wilcoxon_qda'));
    axial_tumor_wilcoxon_lda_index = find(contains(P, 'Axial_tumor_only_wilcoxon_lda'));
    axial_tumor_wilcoxon_rf_index = find(contains(P, 'Axial_tumor_only_wilcoxon_rf'));
    axial_tumor_mrmr_qda_index = find(contains(P, 'Axial_tumor_only_mrmr_qda'));
    axial_tumor_mrmr_lda_index = find(contains(P, 'Axial_tumor_only_mrmr_lda'));
    axial_tumor_mrmr_rf_index = find(contains(P, 'Axial_tumor_only_mrmr_rf'));
    coronal_tumor_wilcoxon_qda_index = find(contains(P, 'Coronal_tumor_only_wilcoxon_qda'));
    coronal_tumor_wilcoxon_lda_index = find(contains(P, 'Coronal_tumor_only_wilcoxon_lda'));
    coronal_tumor_wilcoxon_rf_index = find(contains(P, 'Coronal_tumor_only_wilcoxon_rf'));
    coronal_tumor_mrmr_qda_index = find(contains(P, 'Coronal_tumor_only_mrmr_qda'));
    coronal_tumor_mrmr_lda_index = find(contains(P, 'Coronal_tumor_only_mrmr_lda'));
    coronal_tumor_mrmr_rf_index = find(contains(P, 'Coronal_tumor_only_mrmr_rf'));

    apf5_wilcoxon_qda = load(fullfile(results_root, P{axial_proxfat5_wilcoxon_qda_index}, 'Top_Features_Information.mat'));
    apf5_wilcoxon_lda = load(fullfile(results_root, P{axial_proxfat5_wilcoxon_lda_index}, 'Top_Features_Information.mat'));
    apf5_wilcoxon_rf = load(fullfile(results_root, P{axial_proxfat5_wilcoxon_rf_index}, 'Top_Features_Information.mat'));
    cpf5_wilcoxon_qda = load(fullfile(results_root, P{coronal_proxfat5_wilcoxon_qda_index}, 'Top_Features_Information.mat'));
    cpf5_wilcoxon_lda = load(fullfile(results_root, P{coronal_proxfat5_wilcoxon_lda_index}, 'Top_Features_Information.mat'));
    cpf5_wilcoxon_rf = load(fullfile(results_root, P{coronal_proxfat5_wilcoxon_rf_index}, 'Top_Features_Information.mat'));
    % apf5_mrmr_qda = load(fullfile(results_root, P{axial_proxfat5_mrmr_qda_index}, 'Top_Features_Information.mat'));
    % apf5_mrmr_lda = load(fullfile(results_root, P{axial_proxfat5_mrmr_lda_index}, 'Top_Features_Information.mat'));
    apf5_mrmr_rf = load(fullfile(results_root, P{axial_proxfat5_mrmr_rf_index}, 'Top_Features_Information.mat'));
    cpf5_mrmr_qda = load(fullfile(results_root, P{coronal_proxfat5_mrmr_qda_index}, 'Top_Features_Information.mat'));
    cpf5_mrmr_lda = load(fullfile(results_root, P{coronal_proxfat5_mrmr_lda_index}, 'Top_Features_Information.mat'));
    cpf5_mrmr_rf = load(fullfile(results_root, P{coronal_proxfat5_mrmr_rf_index}, 'Top_Features_Information.mat'));

    apf10_wilcoxon_qda = load(fullfile(results_root, P{axial_proxfat10_wilcoxon_qda_index}, 'Top_Features_Information.mat'));
    apf10_wilcoxon_lda = load(fullfile(results_root, P{axial_proxfat10_wilcoxon_lda_index}, 'Top_Features_Information.mat'));
    apf10_wilcoxon_rf = load(fullfile(results_root, P{axial_proxfat10_wilcoxon_rf_index}, 'Top_Features_Information.mat'));
    cpf10_wilcoxon_qda = load(fullfile(results_root, P{coronal_proxfat10_wilcoxon_qda_index}, 'Top_Features_Information.mat'));
    cpf10_wilcoxon_lda = load(fullfile(results_root, P{coronal_proxfat10_wilcoxon_lda_index}, 'Top_Features_Information.mat'));
    cpf10_wilcoxon_rf = load(fullfile(results_root, P{coronal_proxfat10_wilcoxon_rf_index}, 'Top_Features_Information.mat'));
    apf10_mrmr_qda = load(fullfile(results_root, P{axial_proxfat10_mrmr_qda_index}, 'Top_Features_Information.mat'));
    apf10_mrmr_lda = load(fullfile(results_root, P{axial_proxfat10_mrmr_lda_index}, 'Top_Features_Information.mat'));
    apf10_mrmr_rf = load(fullfile(results_root, P{axial_proxfat10_mrmr_rf_index}, 'Top_Features_Information.mat'));
    cpf10_mrmr_qda = load(fullfile(results_root, P{coronal_proxfat10_mrmr_qda_index}, 'Top_Features_Information.mat'));
    cpf10_mrmr_lda = load(fullfile(results_root, P{coronal_proxfat10_mrmr_lda_index}, 'Top_Features_Information.mat'));
    cpf10_mrmr_rf = load(fullfile(results_root, P{coronal_proxfat10_mrmr_rf_index}, 'Top_Features_Information.mat'));

    at_wilcoxon_qda = load(fullfile(results_root, P{axial_tumor_wilcoxon_qda_index}, 'Top_Features_Information.mat'));
    at_wilcoxon_lda = load(fullfile(results_root, P{axial_tumor_wilcoxon_lda_index}, 'Top_Features_Information.mat'));
    at_wilcoxon_rf = load(fullfile(results_root, P{axial_tumor_wilcoxon_rf_index}, 'Top_Features_Information.mat'));
    ct_wilcoxon_qda = load(fullfile(results_root, P{coronal_tumor_wilcoxon_qda_index}, 'Top_Features_Information.mat'));
    ct_wilcoxon_lda = load(fullfile(results_root, P{coronal_tumor_wilcoxon_lda_index}, 'Top_Features_Information.mat'));
    ct_wilcoxon_rf = load(fullfile(results_root, P{coronal_tumor_wilcoxon_rf_index}, 'Top_Features_Information.mat'));
    at_mrmr_qda = load(fullfile(results_root, P{axial_tumor_mrmr_qda_index}, 'Top_Features_Information.mat'));
    at_mrmr_lda = load(fullfile(results_root, P{axial_tumor_mrmr_lda_index}, 'Top_Features_Information.mat'));
    at_mrmr_rf = load(fullfile(results_root, P{axial_tumor_mrmr_rf_index}, 'Top_Features_Information.mat'));
    ct_mrmr_qda = load(fullfile(results_root, P{coronal_tumor_mrmr_qda_index}, 'Top_Features_Information.mat'));
    ct_mrmr_lda = load(fullfile(results_root, P{coronal_tumor_mrmr_lda_index}, 'Top_Features_Information.mat'));
    ct_mrmr_rf = load(fullfile(results_root, P{coronal_tumor_mrmr_rf_index}, 'Top_Features_Information.mat'));

    apf5w_indices = [apf5_wilcoxon_qda.topResults.indices; apf5_wilcoxon_lda.topResults.indices; apf5_wilcoxon_rf.topResults.indices];
    apf5m_indices = [apf5_mrmr_rf.topResults.indices];
    apf5w_counts = [apf5_wilcoxon_qda.topResults.freq; apf5_wilcoxon_lda.topResults.freq; apf5_wilcoxon_rf.topResults.freq];
    apf5m_counts = [apf5_mrmr_rf.topResults.freq];
    
    apf5w = [apf5w_indices apf5w_counts];
    apf5m = [apf5m_indices apf5m_counts];

    apf10w_indices = [apf10_wilcoxon_qda.topResults.indices; apf10_wilcoxon_lda.topResults.indices; apf10_wilcoxon_rf.topResults.indices];
    apf10m_indices = [apf10_mrmr_qda.topResults.indices; apf10_mrmr_lda.topResults.indices; apf10_mrmr_rf.topResults.indices];
    apf10w_counts = [apf10_wilcoxon_qda.topResults.freq; apf10_wilcoxon_lda.topResults.freq; apf10_wilcoxon_rf.topResults.freq];
    apf10m_counts = [apf10_mrmr_qda.topResults.freq; apf10_mrmr_lda.topResults.freq; apf10_mrmr_rf.topResults.freq];
    
    apf10w = [apf10w_indices apf10w_counts];
    apf10m = [apf10m_indices apf10m_counts];

    atw_indices = [at_wilcoxon_qda.topResults.indices; at_wilcoxon_lda.topResults.indices; at_wilcoxon_rf.topResults.indices];
    atm_indices = [at_mrmr_qda.topResults.indices; at_mrmr_lda.topResults.indices; at_mrmr_rf.topResults.indices];
    atw_counts = [at_wilcoxon_qda.topResults.freq; at_wilcoxon_lda.topResults.freq; at_wilcoxon_rf.topResults.freq];
    atm_counts = [at_mrmr_qda.topResults.freq; at_mrmr_lda.topResults.freq; at_mrmr_rf.topResults.freq];
    
    atw = [atw_indices atw_counts];
    atm = [atm_indices atm_counts];

    cpf5w_indices = [cpf5_wilcoxon_qda.topResults.indices; cpf5_wilcoxon_lda.topResults.indices; cpf5_wilcoxon_rf.topResults.indices];
    cpf5m_indices = [cpf5_mrmr_qda.topResults.indices; cpf5_mrmr_lda.topResults.indices; cpf5_mrmr_rf.topResults.indices];
    cpf5w_counts = [cpf5_wilcoxon_qda.topResults.freq; cpf5_wilcoxon_lda.topResults.freq; cpf5_wilcoxon_rf.topResults.freq];
    cpf5m_counts = [cpf5_mrmr_qda.topResults.freq; cpf5_mrmr_lda.topResults.freq; cpf5_mrmr_rf.topResults.freq];
    
    cpf5w = [cpf5w_indices cpf5w_counts];
    cpf5m = [cpf5m_indices cpf5m_counts];

    cpf10w_indices = [cpf10_wilcoxon_qda.topResults.indices; cpf10_wilcoxon_lda.topResults.indices; cpf10_wilcoxon_rf.topResults.indices];
    cpf10m_indices = [cpf10_mrmr_qda.topResults.indices; cpf10_mrmr_lda.topResults.indices; cpf10_mrmr_rf.topResults.indices];
    cpf10w_counts = [cpf10_wilcoxon_qda.topResults.freq; cpf10_wilcoxon_lda.topResults.freq; cpf10_wilcoxon_rf.topResults.freq];
    cpf10m_counts = [cpf10_mrmr_qda.topResults.freq; cpf10_mrmr_lda.topResults.freq; cpf10_mrmr_rf.topResults.freq];
    
    cpf10w = [cpf10w_indices cpf10w_counts];
    cpf10m = [cpf10m_indices cpf10m_counts];

    ctw_indices = [ct_wilcoxon_qda.topResults.indices; ct_wilcoxon_lda.topResults.indices; ct_wilcoxon_rf.topResults.indices];
    ctm_indices = [ct_mrmr_qda.topResults.indices; ct_mrmr_lda.topResults.indices; ct_mrmr_rf.topResults.indices];
    ctw_counts = [ct_wilcoxon_qda.topResults.freq; ct_wilcoxon_lda.topResults.freq; ct_wilcoxon_rf.topResults.freq];
    ctm_counts = [ct_mrmr_qda.topResults.freq; ct_mrmr_lda.topResults.freq; ct_mrmr_rf.topResults.freq];
    
    ctw = [ctw_indices ctw_counts];
    ctm = [ctm_indices ctm_counts];
end

%%
atw_feature_data = containers.Map('KeyType', 'double', 'ValueType', 'any');
    
for i = 1 : size(atw, 1)
    feature_index = atw(i, 1);
    count = atw(i, 2);
        
    % If the feature is already in the dictionary, update the count and sum
    if isKey(atw_feature_data, feature_index)
        atw_feature_data(feature_index) = [atw_feature_data(feature_index); count];
    else
        atw_feature_data(feature_index) = count;
    end
end

% Initialize a new array to store the results
atw = zeros(length(keys(atw_feature_data)), 2);

% Iterate through the dictionary to calculate the average and store in result_data
atw_keys_array = keys(atw_feature_data);

for i = 1:length(atw_keys_array)
    key = atw_keys_array{i};
    counts = atw_feature_data(key);
    avg_count = sum(counts) / length(counts);

    atw(i, 1) = key;
    atw(i, 2) = avg_count;
    atw = sortrows(atw, -2);
end

atw_feature_names = cell(length(atw_keys_array), 1);

for i = 1 : length(atw)
    indices = atw(:, 1);
    atw_feature_names{i} = axial_tumor_feature_column_names{indices(i)};
end

atm_feature_data = containers.Map('KeyType', 'double', 'ValueType', 'any');
    
for i = 1 : size(atm, 1)
    feature_index = atm(i, 1);
    count = atm(i, 2);
        
    % If the feature is already in the dictionary, update the count and sum
    if isKey(atm_feature_data, feature_index)
        atm_feature_data(feature_index) = [atm_feature_data(feature_index); count];
    else
        atm_feature_data(feature_index) = count;
    end
end

% Initialize a new array to store the results
atm = zeros(length(keys(atm_feature_data)), 2);

% Iterate through the dictionary to calculate the average and store in result_data
atm_keys_array = keys(atm_feature_data);

for i = 1:length(atm_keys_array)
    key = atm_keys_array{i};
    counts = atm_feature_data(key);
    avg_count = sum(counts) / length(counts);

    atm(i, 1) = key;
    atm(i, 2) = avg_count;
    atm = sortrows(atm, -2);
end

atm_feature_names = cell(length(atm_keys_array), 1);

for i = 1 : length(atm)
    indices = atm(:, 1);
    atm_feature_names{i} = axial_tumor_feature_column_names{indices(i)};
end

%%
apf5w_feature_data = containers.Map('KeyType', 'double', 'ValueType', 'any');
    
for i = 1 : size(apf5w, 1)
    feature_index = apf5w(i, 1);
    count = apf5w(i, 2);
        
    % If the feature is already in the dictionary, update the count and sum
    if isKey(apf5w_feature_data, feature_index)
        apf5w_feature_data(feature_index) = [apf5w_feature_data(feature_index); count];
    else
        apf5w_feature_data(feature_index) = count;
    end
end

% Initialize a new array to store the results
apf5w = zeros(length(keys(apf5w_feature_data)), 2);

% Iterate through the dictionary to calculate the average and store in result_data
apf5w_keys_array = keys(apf5w_feature_data);

for i = 1:length(apf5w_keys_array)
    key = apf5w_keys_array{i};
    counts = apf5w_feature_data(key);
    avg_count = sum(counts) / length(counts);

    apf5w(i, 1) = key;
    apf5w(i, 2) = avg_count;
    apf5w = sortrows(apf5w, -2);
end

apf5w_feature_names = cell(length(apf5w_keys_array), 1);

for i = 1 : length(apf5w)
    indices = apf5w(:, 1);
    apf5w_feature_names{i} = axial_proxfat5_feature_column_names{indices(i)};
end

apf5m_feature_data = containers.Map('KeyType', 'double', 'ValueType', 'any');
    
for i = 1 : size(apf5m, 1)
    feature_index = apf5m(i, 1);
    count = apf5m(i, 2);
        
    % If the feature is already in the dictionary, update the count and sum
    if isKey(apf5m_feature_data, feature_index)
        apf5m_feature_data(feature_index) = [apf5m_feature_data(feature_index); count];
    else
        apf5m_feature_data(feature_index) = count;
    end
end

% Initialize a new array to store the results
apf5m = zeros(length(keys(apf5m_feature_data)), 2);

% Iterate through the dictionary to calculate the average and store in result_data
apf5m_keys_array = keys(apf5m_feature_data);

for i = 1:length(apf5m_keys_array)
    key = apf5m_keys_array{i};
    counts = apf5m_feature_data(key);
    avg_count = sum(counts) / length(counts);

    apf5m(i, 1) = key;
    apf5m(i, 2) = avg_count;
    apf5m = sortrows(apf5m, -2);
end

apf5m_feature_names = cell(length(apf5m_keys_array), 1);

for i = 1 : length(apf5m)
    indices = apf5m(:, 1);
    apf5m_feature_names{i} = axial_proxfat5_feature_column_names{indices(i)};
end

%%
apf10w_feature_data = containers.Map('KeyType', 'double', 'ValueType', 'any');
    
for i = 1 : size(apf10w, 1)
    feature_index = apf10w(i, 1);
    count = apf10w(i, 2);
        
    % If the feature is already in the dictionary, update the count and sum
    if isKey(apf10w_feature_data, feature_index)
        apf10w_feature_data(feature_index) = [apf10w_feature_data(feature_index); count];
    else
        apf10w_feature_data(feature_index) = count;
    end
end

% Initialize a new array to store the results
apf10w = zeros(length(keys(apf10w_feature_data)), 2);

% Iterate through the dictionary to calculate the average and store in result_data
apf10w_keys_array = keys(apf10w_feature_data);

for i = 1:length(apf10w_keys_array)
    key = apf10w_keys_array{i};
    counts = apf10w_feature_data(key);
    avg_count = sum(counts) / length(counts);

    apf10w(i, 1) = key;
    apf10w(i, 2) = avg_count;
    apf10w = sortrows(apf10w, -2);
end

apf10w_feature_names = cell(length(apf10w_keys_array), 1);

for i = 1 : length(apf10w)
    indices = apf10w(:, 1);
    apf10w_feature_names{i} = axial_proxfat10_feature_column_names{indices(i)};
end

apf10m_feature_data = containers.Map('KeyType', 'double', 'ValueType', 'any');
    
for i = 1 : size(apf10m, 1)
    feature_index = apf10m(i, 1);
    count = apf10m(i, 2);
        
    % If the feature is already in the dictionary, update the count and sum
    if isKey(apf10m_feature_data, feature_index)
        apf10m_feature_data(feature_index) = [apf10m_feature_data(feature_index); count];
    else
        apf10m_feature_data(feature_index) = count;
    end
end

% Initialize a new array to store the results
apf10m = zeros(length(keys(apf10m_feature_data)), 2);

% Iterate through the dictionary to calculate the average and store in result_data
apf10m_keys_array = keys(apf10m_feature_data);

for i = 1:length(apf10m_keys_array)
    key = apf10m_keys_array{i};
    counts = apf10m_feature_data(key);
    avg_count = sum(counts) / length(counts);

    apf10m(i, 1) = key;
    apf10m(i, 2) = avg_count;
    apf10m = sortrows(apf10m, -2);
end

apf10m_feature_names = cell(length(apf10m_keys_array), 1);

for i = 1 : length(apf10m)
    indices = apf10m(:, 1);
    apf10m_feature_names{i} = axial_proxfat10_feature_column_names{indices(i)};
end

ctw_feature_data = containers.Map('KeyType', 'double', 'ValueType', 'any');
    
for i = 1 : size(ctw, 1)
    feature_index = ctw(i, 1);
    count = ctw(i, 2);
        
    % If the feature is already in the dictionary, update the count and sum
    if isKey(ctw_feature_data, feature_index)
        ctw_feature_data(feature_index) = [ctw_feature_data(feature_index); count];
    else
        ctw_feature_data(feature_index) = count;
    end
end

% Initialize a new array to store the results
ctw = zeros(length(keys(ctw_feature_data)), 2);

% Iterate through the dictionary to calculate the average and store in result_data
ctw_keys_array = keys(ctw_feature_data);

for i = 1:length(ctw_keys_array)
    key = ctw_keys_array{i};
    counts = ctw_feature_data(key);
    avg_count = sum(counts) / length(counts);

    ctw(i, 1) = key;
    ctw(i, 2) = avg_count;
    ctw = sortrows(ctw, -2);
end

ctw_feature_names = cell(length(ctw_keys_array), 1);

for i = 1 : length(ctw)
    indices = ctw(:, 1);
    ctw_feature_names{i} = coronal_tumor_feature_column_names{indices(i)};
end

ctm_feature_data = containers.Map('KeyType', 'double', 'ValueType', 'any');
    
for i = 1 : size(ctm, 1)
    feature_index = ctm(i, 1);
    count = ctm(i, 2);
        
    % If the feature is already in the dictionary, update the count and sum
    if isKey(ctm_feature_data, feature_index)
        ctm_feature_data(feature_index) = [ctm_feature_data(feature_index); count];
    else
        ctm_feature_data(feature_index) = count;
    end
end

% Initialize a new array to store the results
ctm = zeros(length(keys(ctm_feature_data)), 2);

% Iterate through the dictionary to calculate the average and store in result_data
ctm_keys_array = keys(ctm_feature_data);

for i = 1:length(ctm_keys_array)
    key = ctm_keys_array{i};
    counts = ctm_feature_data(key);
    avg_count = sum(counts) / length(counts);

    ctm(i, 1) = key;
    ctm(i, 2) = avg_count;
    ctm = sortrows(ctm, -2);
end

ctm_feature_names = cell(length(ctm_keys_array), 1);

for i = 1 : length(ctm)
    indices = ctm(:, 1);
    ctm_feature_names{i} = coronal_tumor_feature_column_names{indices(i)};
end

cpf5w_feature_data = containers.Map('KeyType', 'double', 'ValueType', 'any');
    
for i = 1 : size(cpf5w, 1)
    feature_index = cpf5w(i, 1);
    count = cpf5w(i, 2);
        
    % If the feature is already in the dictionary, update the count and sum
    if isKey(cpf5w_feature_data, feature_index)
        cpf5w_feature_data(feature_index) = [cpf5w_feature_data(feature_index); count];
    else
        cpf5w_feature_data(feature_index) = count;
    end
end

% Initialize a new array to store the results
cpf5w = zeros(length(keys(cpf5w_feature_data)), 2);

% Iterate through the dictionary to calculate the average and store in result_data
cpf5w_keys_array = keys(cpf5w_feature_data);

for i = 1:length(cpf5w_keys_array)
    key = cpf5w_keys_array{i};
    counts = cpf5w_feature_data(key);
    avg_count = sum(counts) / length(counts);

    cpf5w(i, 1) = key;
    cpf5w(i, 2) = avg_count;
    cpf5w = sortrows(cpf5w, -2);
end

cpf5w_feature_names = cell(length(cpf5w_keys_array), 1);

for i = 1 : length(cpf5w)
    indices = cpf5w(:, 1);
    cpf5w_feature_names{i} = coronal_proxfat5_feature_column_names{indices(i)};
end

cpf5m_feature_data = containers.Map('KeyType', 'double', 'ValueType', 'any');
    
for i = 1 : size(cpf5m, 1)
    feature_index = cpf5m(i, 1);
    count = cpf5m(i, 2);
        
    % If the feature is already in the dictionary, update the count and sum
    if isKey(cpf5m_feature_data, feature_index)
        cpf5m_feature_data(feature_index) = [cpf5m_feature_data(feature_index); count];
    else
        cpf5m_feature_data(feature_index) = count;
    end
end

% Initialize a new array to store the results
cpf5m = zeros(length(keys(cpf5m_feature_data)), 2);

% Iterate through the dictionary to calculate the average and store in result_data
cpf5m_keys_array = keys(cpf5m_feature_data);

for i = 1:length(cpf5m_keys_array)
    key = cpf5m_keys_array{i};
    counts = cpf5m_feature_data(key);
    avg_count = sum(counts) / length(counts);

    cpf5m(i, 1) = key;
    cpf5m(i, 2) = avg_count;
    cpf5m = sortrows(cpf5m, -2);
end

cpf5m_feature_names = cell(length(cpf5m_keys_array), 1);

for i = 1 : length(cpf5m)
    indices = cpf5m(:, 1);
    cpf5m_feature_names{i} = coronal_proxfat5_feature_column_names{indices(i)};
end

cpf10w_feature_data = containers.Map('KeyType', 'double', 'ValueType', 'any');
    
for i = 1 : size(cpf10w, 1)
    feature_index = cpf10w(i, 1);
    count = cpf10w(i, 2);
        
    % If the feature is already in the dictionary, update the count and sum
    if isKey(cpf10w_feature_data, feature_index)
        cpf10w_feature_data(feature_index) = [cpf10w_feature_data(feature_index); count];
    else
        cpf10w_feature_data(feature_index) = count;
    end
end

% Initialize a new array to store the results
cpf10w = zeros(length(keys(cpf10w_feature_data)), 2);

% Iterate through the dictionary to calculate the average and store in result_data
cpf10w_keys_array = keys(cpf10w_feature_data);

for i = 1:length(cpf10w_keys_array)
    key = cpf10w_keys_array{i};
    counts = cpf10w_feature_data(key);
    avg_count = sum(counts) / length(counts);

    cpf10w(i, 1) = key;
    cpf10w(i, 2) = avg_count;
    cpf10w = sortrows(cpf10w, -2);
end

cpf10w_feature_names = cell(length(cpf10w_keys_array), 1);

for i = 1 : length(cpf10w)
    indices = cpf10w(:, 1);
    cpf10w_feature_names{i} = coronal_proxfat10_feature_column_names{indices(i)};
end

cpf10m_feature_data = containers.Map('KeyType', 'double', 'ValueType', 'any');
    
for i = 1 : size(cpf10m, 1)
    feature_index = cpf10m(i, 1);
    count = cpf10m(i, 2);
        
    % If the feature is already in the dictionary, update the count and sum
    if isKey(cpf10m_feature_data, feature_index)
        cpf10m_feature_data(feature_index) = [cpf10m_feature_data(feature_index); count];
    else
        cpf10m_feature_data(feature_index) = count;
    end
end

% Initialize a new array to store the results
cpf10m = zeros(length(keys(cpf10m_feature_data)), 2);

% Iterate through the dictionary to calculate the average and store in result_data
cpf10m_keys_array = keys(cpf10m_feature_data);

for i = 1:length(cpf10m_keys_array)
    key = cpf10m_keys_array{i};
    counts = cpf10m_feature_data(key);
    avg_count = sum(counts) / length(counts);

    cpf10m(i, 1) = key;
    cpf10m(i, 2) = avg_count;
    cpf10m = sortrows(cpf10m, -2);
end

cpf10m_feature_names = cell(length(cpf10m_keys_array), 1);

for i = 1 : length(cpf10m)
    indices = cpf10m(:, 1);
    cpf10m_feature_names{i} = coronal_proxfat10_feature_column_names{indices(i)};
end

matrix_root_path = '/Users/leobao/Documents/MultiPlanePipeline/Data/MissingCollage/';

axial_proxfat5_train_matrix_path = string(strcat(matrix_root_path, region, '/MPMR_Axial_ProxFat5_Training.mat'));
axial_proxfat5_test1_matrix_path = string(strcat(matrix_root_path, region, '/MPMR_Axial_ProxFat5_Testing1.mat'));
axial_proxfat5_test2_matrix_path = string(strcat(matrix_root_path, region, '/MPMR_Axial_ProxFat5_Testing2.mat'));

coronal_proxfat5_train_matrix_path = string(strcat(matrix_root_path, region, '/MPMR_Coronal_ProxFat5_Training.mat'));
coronal_proxfat5_test1_matrix_path = string(strcat(matrix_root_path, region, '/MPMR_Coronal_ProxFat5_Testing1.mat'));
coronal_proxfat5_test2_matrix_path = string(strcat(matrix_root_path, region, '/MPMR_Coronal_ProxFat5_Testing2.mat'));

axial_proxfat10_train_matrix_path = string(strcat(matrix_root_path, region, '/MPMR_Axial_ProxFat10_Training.mat'));
axial_proxfat10_test1_matrix_path = string(strcat(matrix_root_path, region, '/MPMR_Axial_ProxFat10_Testing1.mat'));
axial_proxfat10_test2_matrix_path = string(strcat(matrix_root_path, region, '/MPMR_Axial_ProxFat10_Testing2.mat'));

coronal_proxfat10_train_matrix_path = string(strcat(matrix_root_path, region, '/MPMR_Coronal_ProxFat10_Training.mat'));
coronal_proxfat10_test1_matrix_path = string(strcat(matrix_root_path, region, '/MPMR_Coronal_ProxFat10_Testing1.mat'));
coronal_proxfat10_test2_matrix_path = string(strcat(matrix_root_path, region, '/MPMR_Coronal_ProxFat10_Testing2.mat'));

axial_tumor_train_matrix_path = string(strcat(matrix_root_path, region, '/MPMR_Axial_Tumor_Training.mat'));
axial_tumor_test1_matrix_path = string(strcat(matrix_root_path, region, '/MPMR_Axial_Tumor_Testing1.mat'));
axial_tumor_test2_matrix_path = string(strcat(matrix_root_path, region, '/MPMR_Axial_Tumor_Testing2.mat'));

coronal_tumor_train_matrix_path = string(strcat(matrix_root_path, region, '/MPMR_Coronal_Tumor_Training.mat'));
coronal_tumor_test1_matrix_path = string(strcat(matrix_root_path, region, '/MPMR_Coronal_Tumor_Testing1.mat'));
coronal_tumor_test2_matrix_path = string(strcat(matrix_root_path, region, '/MPMR_Coronal_Tumor_Testing2.mat'));

axial_proxfat5_path_to_train = string(axial_proxfat5_train_matrix_path);
axial_proxfat5_path_to_test1 = string(axial_proxfat5_test1_matrix_path);
axial_proxfat5_path_to_test2 = string(axial_proxfat5_test2_matrix_path);

coronal_proxfat5_path_to_train = string(coronal_proxfat5_train_matrix_path);
coronal_proxfat5_path_to_test1 = string(coronal_proxfat5_test1_matrix_path);
coronal_proxfat5_path_to_test2 = string(coronal_proxfat5_test2_matrix_path);

axial_proxfat10_path_to_train = string(axial_proxfat10_train_matrix_path);
axial_proxfat10_path_to_test1 = string(axial_proxfat10_test1_matrix_path);
axial_proxfat10_path_to_test2 = string(axial_proxfat10_test2_matrix_path);

coronal_proxfat10_path_to_train = string(coronal_proxfat10_train_matrix_path);
coronal_proxfat10_path_to_test1 = string(coronal_proxfat10_test1_matrix_path);
coronal_proxfat10_path_to_test2 = string(coronal_proxfat10_test2_matrix_path);

axial_tumor_path_to_train = string(axial_tumor_train_matrix_path);
axial_tumor_path_to_test1 = string(axial_tumor_test1_matrix_path);
axial_tumor_path_to_test2 = string(axial_tumor_test2_matrix_path);

coronal_tumor_path_to_train = string(coronal_tumor_train_matrix_path);
coronal_tumor_path_to_test1 = string(coronal_tumor_test1_matrix_path);
coronal_tumor_path_to_test2 = string(coronal_tumor_test2_matrix_path);

label_path_root = '/Users/leobao/Documents/MultiPlanePipeline/Data/train_test_labels/MissingCollageLabels/MPMR';

training_label_path = string(strcat(label_path_root, '/train_labels.csv'));
testing1_label_path = string(strcat(label_path_root, '/test1_labels.csv'));
testing2_label_path = string(strcat(label_path_root, '/test2_labels.csv'));

feature_column_path = '/Users/leobao/Documents/MultiPlanePipeline/Feature_Names/Texture_Feature_Names.xlsx';
feature_column_names = readtable(feature_column_path, 'ReadVariableNames', false);
feature_column_names = table2cell(feature_column_names);

if contains(scheme, 'wilcoxon')    
    axial_proxfat5_top_indices = apf5w(:, 1)';
    coronal_proxfat5_top_indices = cpf5w(:, 1)';
    axial_proxfat10_top_indices = apf10w(:, 1)';
    coronal_proxfat10_top_indices = cpf10w(:, 1)';
    axial_tumor_top_indices = atw(:, 1)';
    coronal_tumor_top_indices = ctw(:, 1)';

    axial_proxfat5_feature_names = apf5w_feature_names';
    axial_proxfat10_feature_names = apf10w_feature_names';
    axial_tumor_feature_names = atw_feature_names';
    coronal_proxfat5_feature_names = cpf5w_feature_names';
    coronal_proxfat10_feature_names = cpf10w_feature_names';
    coronal_tumor_feature_names = ctw_feature_names';
elseif contains(scheme, 'mrmr')
    axial_proxfat5_top_indices = apf5m(:, 1)';
    coronal_proxfat5_top_indices = cpf5m(:, 1)';
    axial_proxfat10_top_indices = apf10m(:, 1)';
    coronal_proxfat10_top_indices = cpf10m(:, 1)';
    axial_tumor_top_indices = atm(:, 1)';
    coronal_tumor_top_indices = ctm(:, 1)';

    axial_proxfat5_feature_names = apf5m_feature_names';
    axial_proxfat10_feature_names = apf10m_feature_names';
    axial_tumor_feature_names = atm_feature_names';
    coronal_proxfat5_feature_names = cpf5m_feature_names';
    coronal_proxfat10_feature_names = cpf10m_feature_names';
    coronal_tumor_feature_names = ctm_feature_names';
end


fprintf("Got directory paths to feature matrices! \n");

output_path = "/Users/leobao/Documents/MultiPlanePipeline/Data/MissingCollageResults/";
experiment_date = strrep(string(datetime("now")), " ", "_");
experiment_date = strrep(experiment_date, ":", "_");
output_path = strcat(output_path, experiment_date, "_", region, "_", scheme, "/");
if(~exist(output_path, "dir"))
    mkdir(output_path);
end

fprintf("Created output directory for this experiment! \n");

datasets = struct;

datasets.(strcat("axial_proxfat5_train_MPMR")) = load(axial_proxfat5_path_to_train);
datasets.(strcat("axial_proxfat5_test1_MPMR")) = load(axial_proxfat5_path_to_test1);
datasets.(strcat("axial_proxfat5_test2_MPMR")) = load(axial_proxfat5_path_to_test2);
datasets.(strcat("axial_proxfat10_train_MPMR")) = load(axial_proxfat10_path_to_train);
datasets.(strcat("axial_proxfat10_test1_MPMR")) = load(axial_proxfat10_path_to_test1);
datasets.(strcat("axial_proxfat10_test2_MPMR")) = load(axial_proxfat10_path_to_test2);
datasets.(strcat("axial_tumor_train_MPMR")) = load(axial_tumor_path_to_train);
datasets.(strcat("axial_tumor_test1_MPMR")) = load(axial_tumor_path_to_test1);
datasets.(strcat("axial_tumor_test2_MPMR")) = load(axial_tumor_path_to_test2);
datasets.(strcat("coronal_proxfat5_train_MPMR")) = load(coronal_proxfat5_path_to_train);
datasets.(strcat("coronal_proxfat5_test1_MPMR")) = load(coronal_proxfat5_path_to_test1);
datasets.(strcat("coronal_proxfat5_test2_MPMR")) = load(coronal_proxfat5_path_to_test2);
datasets.(strcat("coronal_proxfat10_train_MPMR")) = load(coronal_proxfat10_path_to_train);
datasets.(strcat("coronal_proxfat10_test1_MPMR")) = load(coronal_proxfat10_path_to_test1);
datasets.(strcat("coronal_proxfat10_test2_MPMR")) = load(coronal_proxfat10_path_to_test2);
datasets.(strcat("coronal_tumor_train_MPMR")) = load(coronal_tumor_path_to_train);
datasets.(strcat("coronal_tumor_test1_MPMR")) = load(coronal_tumor_path_to_test1);
datasets.(strcat("coronal_tumor_test2_MPMR")) = load(coronal_tumor_path_to_test2);

data_labels_training = readmatrix(training_label_path);
data_labels_holdout1 = readmatrix(testing1_label_path);
data_labels_holdout2 = readmatrix(testing2_label_path);

axial_proxfat5_training_features = datasets.axial_proxfat5_train_MPMR.feature_matrix;
axial_proxfat5_testing1_features = datasets.axial_proxfat5_test1_MPMR.feature_matrix;
axial_proxfat5_testing2_features = datasets.axial_proxfat5_test2_MPMR.feature_matrix;
axial_proxfat10_training_features = datasets.axial_proxfat10_train_MPMR.feature_matrix;
axial_proxfat10_testing1_features = datasets.axial_proxfat10_test1_MPMR.feature_matrix;
axial_proxfat10_testing2_features = datasets.axial_proxfat10_test2_MPMR.feature_matrix;
axial_tumor_training_features = datasets.axial_tumor_train_MPMR.feature_matrix;
axial_tumor_testing1_features = datasets.axial_tumor_test1_MPMR.feature_matrix;
axial_tumor_testing2_features = datasets.axial_tumor_test2_MPMR.feature_matrix;
coronal_proxfat5_training_features = datasets.coronal_proxfat5_train_MPMR.feature_matrix;
coronal_proxfat5_testing1_features = datasets.coronal_proxfat5_test1_MPMR.feature_matrix;
coronal_proxfat5_testing2_features = datasets.coronal_proxfat5_test2_MPMR.feature_matrix;
coronal_proxfat10_training_features = datasets.coronal_proxfat10_train_MPMR.feature_matrix;
coronal_proxfat10_testing1_features = datasets.coronal_proxfat10_test1_MPMR.feature_matrix;
coronal_proxfat10_testing2_features = datasets.coronal_proxfat10_test2_MPMR.feature_matrix;
coronal_tumor_training_features = datasets.coronal_tumor_train_MPMR.feature_matrix;
coronal_tumor_testing1_features = datasets.coronal_tumor_test1_MPMR.feature_matrix;
coronal_tumor_testing2_features = datasets.coronal_tumor_test2_MPMR.feature_matrix;

axial_proxfat5_top_training_feats = zeros(size(axial_proxfat5_training_features, 1), 5);
axial_proxfat5_top_testing1_feats = zeros(size(axial_proxfat5_testing1_features, 1), 5);
axial_proxfat5_top_testing2_feats = zeros(size(axial_proxfat5_testing2_features, 1), 5);

axial_proxfat10_top_training_feats = zeros(size(axial_proxfat10_training_features, 1), 5);
axial_proxfat10_top_testing1_feats = zeros(size(axial_proxfat10_testing1_features, 1), 5);
axial_proxfat10_top_testing2_feats = zeros(size(axial_proxfat10_testing2_features, 1), 5);

axial_tumor_top_training_feats = zeros(size(axial_tumor_training_features, 1), 5);
axial_tumor_top_testing1_feats = zeros(size(axial_tumor_testing1_features, 1), 5);
axial_tumor_top_testing2_feats = zeros(size(axial_tumor_testing2_features, 1), 5);

coronal_proxfat5_top_training_feats = zeros(size(coronal_proxfat5_training_features, 1), 5);
coronal_proxfat5_top_testing1_feats = zeros(size(coronal_proxfat5_testing1_features, 1), 5);
coronal_proxfat5_top_testing2_feats = zeros(size(coronal_proxfat5_testing2_features, 1), 5);

coronal_proxfat10_top_training_feats = zeros(size(coronal_proxfat10_training_features, 1), 5);
coronal_proxfat10_top_testing1_feats = zeros(size(coronal_proxfat10_testing1_features, 1), 5);
coronal_proxfat10_top_testing2_feats = zeros(size(coronal_proxfat10_testing2_features, 1), 5);

coronal_tumor_top_training_feats = zeros(size(coronal_tumor_training_features, 1), 5);
coronal_tumor_top_testing1_feats = zeros(size(coronal_tumor_testing1_features, 1), 5);
coronal_tumor_top_testing2_feats = zeros(size(coronal_tumor_testing2_features, 1), 5);



for i = 1 : 5
    axial_proxfat5_index = axial_proxfat5_top_indices(i);
    axial_proxfat10_index = axial_proxfat10_top_indices(i);
    axial_tumor_index = axial_tumor_top_indices(i);
    coronal_proxfat5_index = coronal_proxfat5_top_indices(i);
    coronal_proxfat10_index = coronal_proxfat10_top_indices(i);
    coronal_tumor_index = coronal_tumor_top_indices(i);

    axial_proxfat5_top_training_feats(:, i) = axial_proxfat5_training_features(:, axial_proxfat5_index);
    axial_proxfat5_top_testing1_feats(:, i) = axial_proxfat5_testing1_features(:, axial_proxfat5_index);
    axial_proxfat5_top_testing2_feats(:, i) = axial_proxfat5_testing2_features(:, axial_proxfat5_index);
    axial_proxfat5_feature_names = axial_proxfat5_feature_names(1 : 5);

    axial_proxfat10_top_training_feats(:, i) = axial_proxfat10_training_features(:, axial_proxfat10_index);
    axial_proxfat10_top_testing1_feats(:, i) = axial_proxfat10_testing1_features(:, axial_proxfat10_index);
    axial_proxfat10_top_testing2_feats(:, i) = axial_proxfat10_testing2_features(:, axial_proxfat10_index);
    axial_proxfat10_feature_names = axial_proxfat10_feature_names(1 : 5);

    axial_tumor_top_training_feats(:, i) = axial_tumor_training_features(:, axial_tumor_index);
    axial_tumor_top_testing1_feats(:, i) = axial_tumor_testing1_features(:, axial_tumor_index);
    axial_tumor_top_testing2_feats(:, i) = axial_tumor_testing2_features(:, axial_tumor_index);
    axial_tumor_feature_names = axial_tumor_feature_names(1 : 5);

    coronal_proxfat5_top_training_feats(:, i) = coronal_proxfat5_training_features(:, coronal_proxfat5_index);
    coronal_proxfat5_top_testing1_feats(:, i) = coronal_proxfat5_testing1_features(:, coronal_proxfat5_index);
    coronal_proxfat5_top_testing2_feats(:, i) = coronal_proxfat5_testing2_features(:, coronal_proxfat5_index);
    coronal_proxfat5_feature_names = coronal_proxfat5_feature_names(1 : 5);

    coronal_proxfat10_top_training_feats(:, i) = coronal_proxfat10_training_features(:, coronal_proxfat10_index);
    coronal_proxfat10_top_testing1_feats(:, i) = coronal_proxfat10_testing1_features(:, coronal_proxfat10_index);
    coronal_proxfat10_top_testing2_feats(:, i) = coronal_proxfat10_testing2_features(:, coronal_proxfat10_index);
    coronal_proxfat10_feature_names = coronal_proxfat10_feature_names(1 : 5);

    coronal_tumor_top_training_feats(:, i) = coronal_tumor_training_features(:, coronal_tumor_index);
    coronal_tumor_top_testing1_feats(:, i) = coronal_tumor_testing1_features(:, coronal_tumor_index);
    coronal_tumor_top_testing2_feats(:, i) = coronal_tumor_testing2_features(:, coronal_tumor_index);
    coronal_tumor_feature_names = coronal_tumor_feature_names(1 : 5);
end

feature_column_names = [axial_proxfat5_feature_names, coronal_proxfat5_feature_names, axial_proxfat10_feature_names, coronal_proxfat10_feature_names,axial_tumor_feature_names, coronal_tumor_feature_names];
features_training = [axial_proxfat5_top_training_feats, coronal_proxfat5_top_training_feats, axial_proxfat10_top_training_feats, coronal_proxfat10_top_training_feats, axial_tumor_top_training_feats, coronal_tumor_top_training_feats];
features_holdout1 = [axial_proxfat5_top_testing1_feats, coronal_proxfat5_top_testing1_feats, axial_proxfat10_top_testing1_feats, coronal_proxfat10_top_testing1_feats, axial_tumor_top_testing1_feats, coronal_tumor_top_testing1_feats];
features_holdout2 = [axial_proxfat5_top_testing2_feats, coronal_proxfat5_top_testing2_feats, axial_proxfat10_top_testing2_feats, coronal_proxfat10_top_testing2_feats, axial_tumor_top_testing2_feats, coronal_tumor_top_testing2_feats];

feature_file_name = string(strcat(output_path, 'top_features.xlsx'));
feature_table = cell2table(feature_column_names);
writetable(feature_table, feature_file_name);

holdout_test_size1 = size(features_holdout1);
holdout_test_size2 = size(features_holdout2);
assert(holdout_test_size1(1) == 14, "The size of the dataset is incorrect!");
assert(holdout_test_size2(1) == 26, "The size of the dataset is incorrect!");

features_holdout2(14, :) = [];
data_labels_holdout2(14, :) = [];

fprintf("Loaded in train and holdout testing datasets! \n");
%% Whiten the data
features_training = whitenData(features_training, 'simple');
features_holdout1 = whitenData(features_holdout1, 'simple');
features_holdout2 = whitenData(features_holdout2, 'simple');

fprintf("Whitened the training and holdout testing datasets successfully! \n");

%% Initialize classifier params for classifier
saveParams = true;
if(isequal(saveParams,true))
    if strcmp(scheme, 'wilcoxon_qda')
        params.classifier='QDA';
        params.fsname='wilcoxon';
    elseif strcmp(scheme, 'wilcoxon_lda')
        params.classifier='LDA';
        params.fsname='wilcoxon';
    elseif strcmp(scheme, 'wilcoxon_rf')
        params.classifier='RANDOMFOREST';
        params.fsname='wilcoxon';
    elseif strcmp(scheme, 'mrmr_qda')
        params.classifier='QDA';
        params.fsname='mrmr';
    elseif strcmp(scheme, 'mrmr_lda')
        params.classifier='LDA';
        params.fsname='mrmr';
    elseif strcmp(scheme, 'mrmr_rf')
        params.classifier='RANDOMFOREST';
        params.fsname='mrmr';
    end
    params.shuffle = 1;
    params.n = 3;
    params.nIter = 50;
    params.num_top_feats = 5;
    params.threshmeth = 'euclidean';
    params.osname = 'SMOTE';
    params.fsprunecorr = 'false';
    params.featnames = feature_column_names;
    save(strcat(output_path, 'Classifier_params.mat'),'params');
else
    params = load('Classifier_params.mat');
    params = params.params;
end

%% Create experimental log
log_file_name = strcat(output_path, "Summary_of_Experiment_", experiment_date, ".txt");
fileid = fopen(log_file_name,'w');
fprintf(fileid, strcat("Experiment Date and Time: ", string(datetime("now")), '\n\r'));
fprintf(fileid, "FEATURE FAMILY: 3D Shape \n");
fprintf(fileid, strcat("Trained on ", int2str(length(features_training(:,1))), " patients \n\r"));
fprintf(fileid, strcat("Tested on ", int2str(length(features_holdout1(:,1)))," patients \n\r"));
fprintf(fileid, strcat("Selection Method: ", params.fsname, '\n\r'));
fprintf(fileid, strcat("Classifier: ", params.classifier, '\n\r'));
fprintf(fileid, strcat("Number of Top Features: ", num2str(params.num_top_feats), '\n\r'));
fprintf(fileid, strcat("Number of Folds: ", num2str(params.n), '\n\r'));
fclose(fileid);

%% Run K-fold Cross Validation
fprintf("Running K-Fold Cross Validation with Feature Selection... \n");
feat_stats = nFoldCV_withFS(features_training,data_labels_training,params);
save(strcat(output_path,'CV_Feature_Results.mat'),'feat_stats');
FSResults = getCVTopFeatures(features_training, params, feat_stats, data_labels_training, output_path);

%% Analyze the output of feature selection
[topResults,meanValues] = analyze_CV_results(output_path, feature_column_names, params.num_top_feats);
fprintf("Found the top ranked features! \n")

%% Set the optimal threshold from analyze_CV_results
classifier_thresh = meanValues.Opt_Thres;

%% Train and test the QDA classifier
train_feats = features_training(:,topResults.indices);
test1_feats = features_holdout1(:,topResults.indices);
test2_feats = features_holdout2(:,topResults.indices);
fprintf(strcat("Training on ",int2str(length(train_feats(1,:)))," features for ", int2str(length(features_training(:,1))), " patients \n"));
fprintf(strcat("Testing on ",int2str(length(test1_feats(1,:)))," features for ", int2str(length(features_holdout1(:,1)))," patients \n"));

[stats1, methodstring] = TrainTestModel(params.classifier, train_feats, test1_feats, data_labels_training, data_labels_holdout1,classifier_thresh);
save(strcat(output_path,'Holdout_Testing1_Dataset_Results.mat'),'stats1');

stats_for_export1 = struct;
stats_for_export1.ACC = round(stats1.acc, 3);
stats_for_export1.AUC = round(stats1.AUC, 3);
stats_for_export1.AUCPRC = round(stats1.AUPRC, 3);
stats_for_export1.SENS = round(stats1.sens, 3);
stats_for_export1.SPEC = round(stats1.spec, 3);
stats_for_export1.FSCORE = round(stats1.Fscore, 3);
stats_for_export1.TP = stats1.tp;
stats_for_export1.FP = stats1.fp;
stats_for_export1.FN = stats1.fn;
stats_for_export1.TN = stats1.tn;
stats_for_export1.KAPPA = round(stats1.kappa, 3);
stats_for_export1.MCC = round(stats1.MCC, 3);
    

QDAResults1 = struct2table(stats_for_export1);
writetable(QDAResults1,strcat(output_path,'Holdout_Testing1_Dataset_Results.xlsx'), 'WriteVariableNames',1, 'Sheet', 'Entire_Test_Cohort_Results');

fprintf(strcat("Testing on ",int2str(length(test2_feats(1,:)))," features for ", int2str(length(features_holdout2(:,1)))," patients \n"));

[stats2, methodstring] = TrainTestModel(params.classifier, train_feats, test2_feats, data_labels_training, data_labels_holdout2,classifier_thresh);
save(strcat(output_path,'Holdout_Testing2_Dataset_Results.mat'),'stats2');

stats_for_export2 = struct;
stats_for_export2.ACC = round(stats2.acc, 3);
stats_for_export2.AUC = round(stats2.AUC, 3);
stats_for_export2.AUCPRC = round(stats2.AUPRC, 3);
stats_for_export2.SENS = round(stats2.sens, 3);
stats_for_export2.SPEC = round(stats2.spec, 3);
stats_for_export2.FSCORE = round(stats2.Fscore, 3);
stats_for_export2.TP = stats2.tp;
stats_for_export2.FP = stats2.fp;
stats_for_export2.FN = stats2.fn;
stats_for_export2.TN = stats2.tn;
stats_for_export2.KAPPA = round(stats2.kappa, 3);
stats_for_export2.MCC = round(stats2.MCC, 3);
    

QDAResults2 = struct2table(stats_for_export2);
writetable(QDAResults2,strcat(output_path,'Holdout_Testing2_Dataset_Results.xlsx'), 'WriteVariableNames',1, 'Sheet', 'Entire_Test_Cohort_Results');