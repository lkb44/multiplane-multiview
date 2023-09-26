clc;
clear;
close all;
results_root = '/Users/leobao/Documents/MultiPlanePipeline/Data/MissingCollageResults/';

experiments = dir(fullfile(results_root));
theFilesvol = {experiments.name};

feature_column_path = '/Users/leobao/Documents/MultiPlanePipeline/Feature_Names/Texture_Feature_Names.xlsx';
feature_column_names = readtable(feature_column_path, 'ReadVariableNames',false);
feature_column_names = table2cell(feature_column_names);

axial_fat_feature_column_names = cell(1, numel(feature_column_names));
coronal_fat_feature_column_names = cell(1, numel(feature_column_names));
axial_proxfat5_feature_column_names = cell(1, numel(feature_column_names));
coronal_proxfat5_feature_column_names = cell(1, numel(feature_column_names));
axial_proxfat10_feature_column_names = cell(1, numel(feature_column_names));
coronal_proxfat10_feature_column_names = cell(1, numel(feature_column_names));
axial_proxfat15_feature_column_names = cell(1, numel(feature_column_names));
coronal_proxfat15_feature_column_names = cell(1, numel(feature_column_names));
axial_tumor_feature_column_names = cell(1, numel(feature_column_names));
coronal_tumor_feature_column_names = cell(1, numel(feature_column_names));


% Loop through each cell in feature_column_names and add the prefix
for i = 1:numel(feature_column_names)
    axial_fat_feature_column_names{i} = ['axial_fat_', feature_column_names{i}];
    coronal_fat_feature_column_names{i} = ['coronal_fat_', feature_column_names{i}];
    axial_tumor_feature_column_names{i} = ['axial_tumor_', feature_column_names{i}];
    coronal_tumor_feature_column_names{i} = ['coronal_tumor_', feature_column_names{i}];
    axial_proxfat5_feature_column_names{i} = ['axial_proxfat5_', feature_column_names{i}];
    coronal_proxfat5_feature_column_names{i} = ['coronal_proxfat5_', feature_column_names{i}];
    axial_proxfat10_feature_column_names{i} = ['axial_proxfat10_', feature_column_names{i}];
    coronal_proxfat10_feature_column_names{i} = ['coronal_proxfat10_', feature_column_names{i}];
    axial_proxfat15_feature_column_names{i} = ['axial_proxfat15_', feature_column_names{i}];
    coronal_proxfat15_feature_column_names{i} = ['coronal_proxfat15_', feature_column_names{i}];
end

for k = 1:length(theFilesvol)
    P = theFilesvol;
    P = P(~startsWith(P, '.'));
    axial_fat_wilcoxon_qda_index = find(contains(P, 'Axial_fat_only_wilcoxon_qda'));
    axial_fat_wilcoxon_lda_index = find(contains(P, 'Axial_fat_only_wilcoxon_lda'));
    axial_fat_wilcoxon_rf_index = find(contains(P, 'Axial_fat_only_wilcoxon_rf'));
    axial_fat_wilcoxon_svm_index = find(contains(P, 'Axial_fat_only_wilcoxon_svm'));
    axial_fat_mrmr_qda_index = find(contains(P, 'Axial_fat_only_mrmr_qda'));
    axial_fat_mrmr_lda_index = find(contains(P, 'Axial_fat_only_mrmr_lda'));
    axial_fat_mrmr_rf_index = find(contains(P, 'Axial_fat_only_mrmr_rf'));
    axial_fat_mrmr_svm_index = find(contains(P, 'Axial_fat_only_mrmr_svm'));
    coronal_fat_wilcoxon_qda_index = find(contains(P, 'Coronal_fat_only_wilcoxon_qda'));
    coronal_fat_wilcoxon_lda_index = find(contains(P, 'Coronal_fat_only_wilcoxon_lda'));
    coronal_fat_wilcoxon_rf_index = find(contains(P, 'Coronal_fat_only_wilcoxon_rf'));
    coronal_fat_wilcoxon_svm_index = find(contains(P, 'Coronal_fat_only_wilcoxon_svm'));
    coronal_fat_mrmr_qda_index = find(contains(P, 'Coronal_fat_only_mrmr_qda'));
    coronal_fat_mrmr_lda_index = find(contains(P, 'Coronal_fat_only_mrmr_lda'));
    coronal_fat_mrmr_rf_index = find(contains(P, 'Coronal_fat_only_mrmr_rf'));
    coronal_fat_mrmr_svm_index = find(contains(P, 'Coronal_fat_only_mrmr_svm'));
    axial_proxfat5_wilcoxon_qda_index = find(contains(P, 'Axial_proxfat5_only_wilcoxon_qda'));
    axial_proxfat5_wilcoxon_lda_index = find(contains(P, 'Axial_proxfat5_only_wilcoxon_lda'));
    axial_proxfat5_wilcoxon_rf_index = find(contains(P, 'Axial_proxfat5_only_wilcoxon_rf'));
    axial_proxfat5_wilcoxon_svm_index = find(contains(P, 'Axial_proxfat5_only_wilcoxon_svm'));
    % axial_proxfat5_mrmr_qda_index = find(contains(P, 'Axial_proxfat5_only_mrmr_qda'));
    axial_proxfat5_mrmr_lda_index = find(contains(P, 'Axial_proxfat5_only_mrmr_lda'));
    axial_proxfat5_mrmr_rf_index = find(contains(P, 'Axial_proxfat5_only_mrmr_rf'));
    axial_proxfat5_mrmr_svm_index = find(contains(P, 'Axial_proxfat5_only_mrmr_svm'));
    coronal_proxfat5_wilcoxon_qda_index = find(contains(P, 'Coronal_proxfat5_only_wilcoxon_qda'));
    coronal_proxfat5_wilcoxon_lda_index = find(contains(P, 'Coronal_proxfat5_only_wilcoxon_lda'));
    coronal_proxfat5_wilcoxon_rf_index = find(contains(P, 'Coronal_proxfat5_only_wilcoxon_rf'));
    coronal_proxfat5_wilcoxon_svm_index = find(contains(P, 'Coronal_proxfat5_only_wilcoxon_svm'));
    coronal_proxfat5_mrmr_qda_index = find(contains(P, 'Coronal_proxfat5_only_mrmr_qda'));
    coronal_proxfat5_mrmr_lda_index = find(contains(P, 'Coronal_proxfat5_only_mrmr_lda'));
    coronal_proxfat5_mrmr_rf_index = find(contains(P, 'Coronal_proxfat5_only_mrmr_rf'));
    coronal_proxfat5_mrmr_svm_index = find(contains(P, 'Coronal_proxfat5_only_mrmr_svm'));
    axial_proxfat10_wilcoxon_qda_index = find(contains(P, 'Axial_proxfat10_only_wilcoxon_qda'));
    axial_proxfat10_wilcoxon_lda_index = find(contains(P, 'Axial_proxfat10_only_wilcoxon_lda'));
    axial_proxfat10_wilcoxon_rf_index = find(contains(P, 'Axial_proxfat10_only_wilcoxon_rf'));
    axial_proxfat10_wilcoxon_svm_index = find(contains(P, 'Axial_proxfat10_only_wilcoxon_svm'));
    axial_proxfat10_mrmr_qda_index = find(contains(P, 'Axial_proxfat10_only_mrmr_qda'));
    axial_proxfat10_mrmr_lda_index = find(contains(P, 'Axial_proxfat10_only_mrmr_lda'));
    axial_proxfat10_mrmr_rf_index = find(contains(P, 'Axial_proxfat10_only_mrmr_rf'));
    axial_proxfat10_mrmr_svm_index = find(contains(P, 'Axial_proxfat10_only_mrmr_svm'));
    coronal_proxfat10_wilcoxon_qda_index = find(contains(P, 'Coronal_proxfat10_only_wilcoxon_qda'));
    coronal_proxfat10_wilcoxon_lda_index = find(contains(P, 'Coronal_proxfat10_only_wilcoxon_lda'));
    coronal_proxfat10_wilcoxon_rf_index = find(contains(P, 'Coronal_proxfat10_only_wilcoxon_rf'));
    coronal_proxfat10_wilcoxon_svm_index = find(contains(P, 'Coronal_proxfat10_only_wilcoxon_svm'));
    coronal_proxfat10_mrmr_qda_index = find(contains(P, 'Coronal_proxfat10_only_mrmr_qda'));
    coronal_proxfat10_mrmr_lda_index = find(contains(P, 'Coronal_proxfat10_only_mrmr_lda'));
    coronal_proxfat10_mrmr_rf_index = find(contains(P, 'Coronal_proxfat10_only_mrmr_rf'));
    coronal_proxfat10_mrmr_svm_index = find(contains(P, 'Coronal_proxfat10_only_mrmr_svm'));
    axial_proxfat15_wilcoxon_qda_index = find(contains(P, 'Axial_proxfat15_only_wilcoxon_qda'));
    axial_proxfat15_wilcoxon_lda_index = find(contains(P, 'Axial_proxfat15_only_wilcoxon_lda'));
    axial_proxfat15_wilcoxon_rf_index = find(contains(P, 'Axial_proxfat15_only_wilcoxon_rf'));
    axial_proxfat15_wilcoxon_svm_index = find(contains(P, 'Axial_proxfat15_only_wilcoxon_svm'));
    axial_proxfat15_mrmr_qda_index = find(contains(P, 'Axial_proxfat15_only_mrmr_qda'));
    axial_proxfat15_mrmr_lda_index = find(contains(P, 'Axial_proxfat15_only_mrmr_lda'));
    axial_proxfat15_mrmr_rf_index = find(contains(P, 'Axial_proxfat15_only_mrmr_rf'));
    axial_proxfat15_mrmr_svm_index = find(contains(P, 'Axial_proxfat15_only_mrmr_svm'));
    coronal_proxfat15_wilcoxon_qda_index = find(contains(P, 'Coronal_proxfat15_only_wilcoxon_qda'));
    coronal_proxfat15_wilcoxon_lda_index = find(contains(P, 'Coronal_proxfat15_only_wilcoxon_lda'));
    coronal_proxfat15_wilcoxon_rf_index = find(contains(P, 'Coronal_proxfat15_only_wilcoxon_rf'));
    coronal_proxfat15_wilcoxon_svm_index = find(contains(P, 'Coronal_proxfat15_only_wilcoxon_svm'));
    coronal_proxfat15_mrmr_qda_index = find(contains(P, 'Coronal_proxfat15_only_mrmr_qda'));
    coronal_proxfat15_mrmr_lda_index = find(contains(P, 'Coronal_proxfat15_only_mrmr_lda'));
    coronal_proxfat15_mrmr_rf_index = find(contains(P, 'Coronal_proxfat15_only_mrmr_rf'));
    coronal_proxfat15_mrmr_svm_index = find(contains(P, 'Coronal_proxfat15_only_mrmr_svm'));
    axial_tumor_wilcoxon_qda_index = find(contains(P, 'Axial_tumor_only_wilcoxon_qda'));
    axial_tumor_wilcoxon_lda_index = find(contains(P, 'Axial_tumor_only_wilcoxon_lda'));
    axial_tumor_wilcoxon_rf_index = find(contains(P, 'Axial_tumor_only_wilcoxon_rf'));
    axial_tumor_wilcoxon_svm_index = find(contains(P, 'Axial_tumor_only_wilcoxon_svm'));
    axial_tumor_mrmr_qda_index = find(contains(P, 'Axial_tumor_only_mrmr_qda'));
    axial_tumor_mrmr_lda_index = find(contains(P, 'Axial_tumor_only_mrmr_lda'));
    axial_tumor_mrmr_rf_index = find(contains(P, 'Axial_tumor_only_mrmr_rf'));
    axial_tumor_mrmr_svm_index = find(contains(P, 'Axial_tumor_only_mrmr_svm'));
    coronal_tumor_wilcoxon_qda_index = find(contains(P, 'Coronal_tumor_only_wilcoxon_qda'));
    coronal_tumor_wilcoxon_lda_index = find(contains(P, 'Coronal_tumor_only_wilcoxon_lda'));
    coronal_tumor_wilcoxon_rf_index = find(contains(P, 'Coronal_tumor_only_wilcoxon_rf'));
    coronal_tumor_wilcoxon_svm_index = find(contains(P, 'Coronal_tumor_only_wilcoxon_svm'));
    coronal_tumor_mrmr_qda_index = find(contains(P, 'Coronal_tumor_only_mrmr_qda'));
    coronal_tumor_mrmr_lda_index = find(contains(P, 'Coronal_tumor_only_mrmr_lda'));
    coronal_tumor_mrmr_rf_index = find(contains(P, 'Coronal_tumor_only_mrmr_rf'));
    coronal_tumor_mrmr_svm_index = find(contains(P, 'Coronal_tumor_only_mrmr_svm'));

    af_wilcoxon_qda = load(fullfile(results_root, P{axial_fat_wilcoxon_qda_index}, 'Top_Features_Information.mat'));
    af_wilcoxon_lda = load(fullfile(results_root, P{axial_fat_wilcoxon_lda_index}, 'Top_Features_Information.mat'));
    af_wilcoxon_rf = load(fullfile(results_root, P{axial_fat_wilcoxon_rf_index}, 'Top_Features_Information.mat'));
    af_wilcoxon_svm = load(fullfile(results_root, P{axial_fat_wilcoxon_svm_index}, 'Top_Features_Information.mat'));
    cf_wilcoxon_qda = load(fullfile(results_root, P{coronal_fat_wilcoxon_qda_index}, 'Top_Features_Information.mat'));
    cf_wilcoxon_lda = load(fullfile(results_root, P{coronal_fat_wilcoxon_lda_index}, 'Top_Features_Information.mat'));
    cf_wilcoxon_rf = load(fullfile(results_root, P{coronal_fat_wilcoxon_rf_index}, 'Top_Features_Information.mat'));
    cf_wilcoxon_svm = load(fullfile(results_root, P{coronal_fat_wilcoxon_svm_index}, 'Top_Features_Information.mat'));
    af_mrmr_qda = load(fullfile(results_root, P{axial_fat_mrmr_qda_index}, 'Top_Features_Information.mat'));
    af_mrmr_lda = load(fullfile(results_root, P{axial_fat_mrmr_lda_index}, 'Top_Features_Information.mat'));
    af_mrmr_rf = load(fullfile(results_root, P{axial_fat_mrmr_rf_index}, 'Top_Features_Information.mat'));
    af_mrmr_svm = load(fullfile(results_root, P{axial_fat_mrmr_svm_index}, 'Top_Features_Information.mat'));
    cf_mrmr_qda = load(fullfile(results_root, P{coronal_fat_mrmr_qda_index}, 'Top_Features_Information.mat'));
    cf_mrmr_lda = load(fullfile(results_root, P{coronal_fat_mrmr_lda_index}, 'Top_Features_Information.mat'));
    cf_mrmr_rf = load(fullfile(results_root, P{coronal_fat_mrmr_rf_index}, 'Top_Features_Information.mat'));
    cf_mrmr_svm = load(fullfile(results_root, P{coronal_fat_mrmr_svm_index}, 'Top_Features_Information.mat'));

    apf5_wilcoxon_qda = load(fullfile(results_root, P{axial_proxfat5_wilcoxon_qda_index}, 'Top_Features_Information.mat'));
    apf5_wilcoxon_lda = load(fullfile(results_root, P{axial_proxfat5_wilcoxon_lda_index}, 'Top_Features_Information.mat'));
    apf5_wilcoxon_rf = load(fullfile(results_root, P{axial_proxfat5_wilcoxon_rf_index}, 'Top_Features_Information.mat'));
    apf5_wilcoxon_svm = load(fullfile(results_root, P{axial_proxfat5_wilcoxon_svm_index}, 'Top_Features_Information.mat'));
    cpf5_wilcoxon_qda = load(fullfile(results_root, P{coronal_proxfat5_wilcoxon_qda_index}, 'Top_Features_Information.mat'));
    cpf5_wilcoxon_lda = load(fullfile(results_root, P{coronal_proxfat5_wilcoxon_lda_index}, 'Top_Features_Information.mat'));
    cpf5_wilcoxon_rf = load(fullfile(results_root, P{coronal_proxfat5_wilcoxon_rf_index}, 'Top_Features_Information.mat'));
    cpf5_wilcoxon_svm = load(fullfile(results_root, P{coronal_proxfat5_wilcoxon_svm_index}, 'Top_Features_Information.mat'));
    % apf5_mrmr_qda = load(fullfile(results_root, P{axial_proxfat5_mrmr_qda_index}, 'Top_Features_Information.mat'));
    % apf5_mrmr_lda = load(fullfile(results_root, P{axial_proxfat5_mrmr_lda_index}, 'Top_Features_Information.mat'));
    apf5_mrmr_rf = load(fullfile(results_root, P{axial_proxfat5_mrmr_rf_index}, 'Top_Features_Information.mat'));
    apf5_mrmr_svm = load(fullfile(results_root, P{axial_proxfat5_mrmr_svm_index}, 'Top_Features_Information.mat'));
    cpf5_mrmr_qda = load(fullfile(results_root, P{coronal_proxfat5_mrmr_qda_index}, 'Top_Features_Information.mat'));
    cpf5_mrmr_lda = load(fullfile(results_root, P{coronal_proxfat5_mrmr_lda_index}, 'Top_Features_Information.mat'));
    cpf5_mrmr_rf = load(fullfile(results_root, P{coronal_proxfat5_mrmr_rf_index}, 'Top_Features_Information.mat'));
    cpf5_mrmr_svm = load(fullfile(results_root, P{coronal_proxfat5_mrmr_svm_index}, 'Top_Features_Information.mat'));

    apf10_wilcoxon_qda = load(fullfile(results_root, P{axial_proxfat10_wilcoxon_qda_index}, 'Top_Features_Information.mat'));
    apf10_wilcoxon_lda = load(fullfile(results_root, P{axial_proxfat10_wilcoxon_lda_index}, 'Top_Features_Information.mat'));
    apf10_wilcoxon_rf = load(fullfile(results_root, P{axial_proxfat10_wilcoxon_rf_index}, 'Top_Features_Information.mat'));
    apf10_wilcoxon_svm = load(fullfile(results_root, P{axial_proxfat10_wilcoxon_svm_index}, 'Top_Features_Information.mat'));
    cpf10_wilcoxon_qda = load(fullfile(results_root, P{coronal_proxfat10_wilcoxon_qda_index}, 'Top_Features_Information.mat'));
    cpf10_wilcoxon_lda = load(fullfile(results_root, P{coronal_proxfat10_wilcoxon_lda_index}, 'Top_Features_Information.mat'));
    cpf10_wilcoxon_rf = load(fullfile(results_root, P{coronal_proxfat10_wilcoxon_rf_index}, 'Top_Features_Information.mat'));
    cpf10_wilcoxon_svm = load(fullfile(results_root, P{coronal_proxfat10_wilcoxon_svm_index}, 'Top_Features_Information.mat'));
    apf10_mrmr_qda = load(fullfile(results_root, P{axial_proxfat10_mrmr_qda_index}, 'Top_Features_Information.mat'));
    apf10_mrmr_lda = load(fullfile(results_root, P{axial_proxfat10_mrmr_lda_index}, 'Top_Features_Information.mat'));
    apf10_mrmr_rf = load(fullfile(results_root, P{axial_proxfat10_mrmr_rf_index}, 'Top_Features_Information.mat'));
    apf10_mrmr_svm = load(fullfile(results_root, P{axial_proxfat10_mrmr_svm_index}, 'Top_Features_Information.mat'));
    cpf10_mrmr_qda = load(fullfile(results_root, P{coronal_proxfat10_mrmr_qda_index}, 'Top_Features_Information.mat'));
    cpf10_mrmr_lda = load(fullfile(results_root, P{coronal_proxfat10_mrmr_lda_index}, 'Top_Features_Information.mat'));
    cpf10_mrmr_rf = load(fullfile(results_root, P{coronal_proxfat10_mrmr_rf_index}, 'Top_Features_Information.mat'));
    cpf10_mrmr_svm = load(fullfile(results_root, P{coronal_proxfat10_mrmr_svm_index}, 'Top_Features_Information.mat'));

    apf15_wilcoxon_qda = load(fullfile(results_root, P{axial_proxfat15_wilcoxon_qda_index}, 'Top_Features_Information.mat'));
    apf15_wilcoxon_lda = load(fullfile(results_root, P{axial_proxfat15_wilcoxon_lda_index}, 'Top_Features_Information.mat'));
    apf15_wilcoxon_rf = load(fullfile(results_root, P{axial_proxfat15_wilcoxon_rf_index}, 'Top_Features_Information.mat'));
    apf15_wilcoxon_svm = load(fullfile(results_root, P{axial_proxfat15_wilcoxon_svm_index}, 'Top_Features_Information.mat'));
    cpf15_wilcoxon_qda = load(fullfile(results_root, P{coronal_proxfat15_wilcoxon_qda_index}, 'Top_Features_Information.mat'));
    cpf15_wilcoxon_lda = load(fullfile(results_root, P{coronal_proxfat15_wilcoxon_lda_index}, 'Top_Features_Information.mat'));
    cpf15_wilcoxon_rf = load(fullfile(results_root, P{coronal_proxfat15_wilcoxon_rf_index}, 'Top_Features_Information.mat'));
    cpf15_wilcoxon_svm = load(fullfile(results_root, P{coronal_proxfat15_wilcoxon_svm_index}, 'Top_Features_Information.mat'));
    apf15_mrmr_qda = load(fullfile(results_root, P{axial_proxfat15_mrmr_qda_index}, 'Top_Features_Information.mat'));
    apf15_mrmr_lda = load(fullfile(results_root, P{axial_proxfat15_mrmr_lda_index}, 'Top_Features_Information.mat'));
    apf15_mrmr_rf = load(fullfile(results_root, P{axial_proxfat15_mrmr_rf_index}, 'Top_Features_Information.mat'));
    apf15_mrmr_svm = load(fullfile(results_root, P{axial_proxfat15_mrmr_svm_index}, 'Top_Features_Information.mat'));
    cpf15_mrmr_qda = load(fullfile(results_root, P{coronal_proxfat15_mrmr_qda_index}, 'Top_Features_Information.mat'));
    cpf15_mrmr_lda = load(fullfile(results_root, P{coronal_proxfat15_mrmr_lda_index}, 'Top_Features_Information.mat'));
    cpf15_mrmr_rf = load(fullfile(results_root, P{coronal_proxfat15_mrmr_rf_index}, 'Top_Features_Information.mat'));
    cpf15_mrmr_svm = load(fullfile(results_root, P{coronal_proxfat15_mrmr_svm_index}, 'Top_Features_Information.mat'));
    
    at_wilcoxon_qda = load(fullfile(results_root, P{axial_tumor_wilcoxon_qda_index}, 'Top_Features_Information.mat'));
    at_wilcoxon_lda = load(fullfile(results_root, P{axial_tumor_wilcoxon_lda_index}, 'Top_Features_Information.mat'));
    at_wilcoxon_rf = load(fullfile(results_root, P{axial_tumor_wilcoxon_rf_index}, 'Top_Features_Information.mat'));
    at_wilcoxon_svm = load(fullfile(results_root, P{axial_tumor_wilcoxon_svm_index}, 'Top_Features_Information.mat'));
    ct_wilcoxon_qda = load(fullfile(results_root, P{coronal_tumor_wilcoxon_qda_index}, 'Top_Features_Information.mat'));
    ct_wilcoxon_lda = load(fullfile(results_root, P{coronal_tumor_wilcoxon_lda_index}, 'Top_Features_Information.mat'));
    ct_wilcoxon_rf = load(fullfile(results_root, P{coronal_tumor_wilcoxon_rf_index}, 'Top_Features_Information.mat'));
    ct_wilcoxon_svm = load(fullfile(results_root, P{coronal_tumor_wilcoxon_svm_index}, 'Top_Features_Information.mat'));
    at_mrmr_qda = load(fullfile(results_root, P{axial_tumor_mrmr_qda_index}, 'Top_Features_Information.mat'));
    at_mrmr_lda = load(fullfile(results_root, P{axial_tumor_mrmr_lda_index}, 'Top_Features_Information.mat'));
    at_mrmr_rf = load(fullfile(results_root, P{axial_tumor_mrmr_rf_index}, 'Top_Features_Information.mat'));
    at_mrmr_svm = load(fullfile(results_root, P{axial_tumor_mrmr_svm_index}, 'Top_Features_Information.mat'));
    ct_mrmr_qda = load(fullfile(results_root, P{coronal_tumor_mrmr_qda_index}, 'Top_Features_Information.mat'));
    ct_mrmr_lda = load(fullfile(results_root, P{coronal_tumor_mrmr_lda_index}, 'Top_Features_Information.mat'));
    ct_mrmr_rf = load(fullfile(results_root, P{coronal_tumor_mrmr_rf_index}, 'Top_Features_Information.mat'));
    ct_mrmr_svm = load(fullfile(results_root, P{coronal_tumor_mrmr_svm_index}, 'Top_Features_Information.mat'));
   
    afw_indices = [af_wilcoxon_qda.topResults.indices; af_wilcoxon_lda.topResults.indices; af_wilcoxon_rf.topResults.indices; af_wilcoxon_svm.topResults.indices];
    afm_indices = [af_mrmr_qda.topResults.indices; af_mrmr_lda.topResults.indices; af_mrmr_rf.topResults.indices; af_mrmr_svm.topResults.indices];
    afw_counts = [af_wilcoxon_qda.topResults.freq; af_wilcoxon_lda.topResults.freq; af_wilcoxon_rf.topResults.freq; af_wilcoxon_svm.topResults.freq];
    afm_counts = [af_mrmr_qda.topResults.freq; af_mrmr_lda.topResults.freq; af_mrmr_rf.topResults.freq; af_mrmr_svm.topResults.freq];
    
    afw = [afw_indices afw_counts];
    afm = [afm_indices afm_counts];

    apf5w_indices = [apf5_wilcoxon_qda.topResults.indices; apf5_wilcoxon_lda.topResults.indices; apf5_wilcoxon_rf.topResults.indices; apf5_wilcoxon_svm.topResults.indices];
    apf5m_indices = [apf5_mrmr_rf.topResults.indices; apf5_mrmr_svm.topResults.indices];
    apf5w_counts = [apf5_wilcoxon_qda.topResults.freq; apf5_wilcoxon_lda.topResults.freq; apf5_wilcoxon_rf.topResults.freq; apf5_wilcoxon_svm.topResults.freq];
    apf5m_counts = [apf5_mrmr_rf.topResults.freq; apf5_mrmr_svm.topResults.freq];
    
    apf5w = [apf5w_indices apf5w_counts];
    apf5m = [apf5m_indices apf5m_counts];

    apf10w_indices = [apf10_wilcoxon_qda.topResults.indices; apf10_wilcoxon_lda.topResults.indices; apf10_wilcoxon_rf.topResults.indices; apf10_wilcoxon_svm.topResults.indices];
    apf10m_indices = [apf10_mrmr_qda.topResults.indices; apf10_mrmr_lda.topResults.indices; apf10_mrmr_rf.topResults.indices; apf10_mrmr_svm.topResults.indices];
    apf10w_counts = [apf10_wilcoxon_qda.topResults.freq; apf10_wilcoxon_lda.topResults.freq; apf10_wilcoxon_rf.topResults.freq; apf10_wilcoxon_svm.topResults.freq];
    apf10m_counts = [apf10_mrmr_qda.topResults.freq; apf10_mrmr_lda.topResults.freq; apf10_mrmr_rf.topResults.freq; apf10_mrmr_svm.topResults.freq];
    
    apf10w = [apf10w_indices apf10w_counts];
    apf10m = [apf10m_indices apf10m_counts];

    apf15w_indices = [apf15_wilcoxon_qda.topResults.indices; apf15_wilcoxon_lda.topResults.indices; apf15_wilcoxon_rf.topResults.indices; apf15_wilcoxon_svm.topResults.indices];
    apf15m_indices = [apf15_mrmr_qda.topResults.indices; apf15_mrmr_lda.topResults.indices; apf15_mrmr_rf.topResults.indices; apf15_mrmr_svm.topResults.indices];
    apf15w_counts = [apf15_wilcoxon_qda.topResults.freq; apf15_wilcoxon_lda.topResults.freq; apf15_wilcoxon_rf.topResults.freq; apf15_wilcoxon_svm.topResults.freq];
    apf15m_counts = [apf15_mrmr_qda.topResults.freq; apf15_mrmr_lda.topResults.freq; apf15_mrmr_rf.topResults.freq; apf15_mrmr_svm.topResults.freq];
    
    apf15w = [apf15w_indices apf15w_counts];
    apf15m = [apf15m_indices apf15m_counts];

    atw_indices = [at_wilcoxon_qda.topResults.indices; at_wilcoxon_lda.topResults.indices; at_wilcoxon_rf.topResults.indices; at_wilcoxon_svm.topResults.indices];
    atm_indices = [at_mrmr_qda.topResults.indices; at_mrmr_lda.topResults.indices; at_mrmr_rf.topResults.indices; at_mrmr_svm.topResults.indices];
    atw_counts = [at_wilcoxon_qda.topResults.freq; at_wilcoxon_lda.topResults.freq; at_wilcoxon_rf.topResults.freq; at_wilcoxon_svm.topResults.freq];
    atm_counts = [at_mrmr_qda.topResults.freq; at_mrmr_lda.topResults.freq; at_mrmr_rf.topResults.freq; at_mrmr_svm.topResults.freq];
    
    atw = [atw_indices atw_counts];
    atm = [atm_indices atm_counts];

    cfw_indices = [cf_wilcoxon_qda.topResults.indices; cf_wilcoxon_lda.topResults.indices; cf_wilcoxon_rf.topResults.indices; cf_wilcoxon_svm.topResults.indices];
    cfm_indices = [cf_mrmr_qda.topResults.indices; cf_mrmr_lda.topResults.indices; cf_mrmr_rf.topResults.indices; cf_mrmr_svm.topResults.indices];
    cfw_counts = [cf_wilcoxon_qda.topResults.freq; cf_wilcoxon_lda.topResults.freq; cf_wilcoxon_rf.topResults.freq; cf_wilcoxon_svm.topResults.freq];
    cfm_counts = [cf_mrmr_qda.topResults.freq; cf_mrmr_lda.topResults.freq; cf_mrmr_rf.topResults.freq; cf_mrmr_svm.topResults.freq];
    
    cfw = [cfw_indices cfw_counts];
    cfm = [cfm_indices cfm_counts];

    cpf5w_indices = [cpf5_wilcoxon_qda.topResults.indices; cpf5_wilcoxon_lda.topResults.indices; cpf5_wilcoxon_rf.topResults.indices; cpf5_wilcoxon_svm.topResults.indices];
    cpf5m_indices = [cpf5_mrmr_qda.topResults.indices; cpf5_mrmr_lda.topResults.indices; cpf5_mrmr_rf.topResults.indices; cpf5_mrmr_svm.topResults.indices];
    cpf5w_counts = [cpf5_wilcoxon_qda.topResults.freq; cpf5_wilcoxon_lda.topResults.freq; cpf5_wilcoxon_rf.topResults.freq; cpf5_wilcoxon_svm.topResults.freq];
    cpf5m_counts = [cpf5_mrmr_qda.topResults.freq; cpf5_mrmr_lda.topResults.freq; cpf5_mrmr_rf.topResults.freq; cpf5_mrmr_svm.topResults.freq];
    
    cpf5w = [cpf5w_indices cpf5w_counts];
    cpf5m = [cpf5m_indices cpf5m_counts];

    cpf10w_indices = [cpf10_wilcoxon_qda.topResults.indices; cpf10_wilcoxon_lda.topResults.indices; cpf10_wilcoxon_rf.topResults.indices; cpf10_wilcoxon_svm.topResults.indices];
    cpf10m_indices = [cpf10_mrmr_qda.topResults.indices; cpf10_mrmr_lda.topResults.indices; cpf10_mrmr_rf.topResults.indices; cpf10_mrmr_svm.topResults.indices];
    cpf10w_counts = [cpf10_wilcoxon_qda.topResults.freq; cpf10_wilcoxon_lda.topResults.freq; cpf10_wilcoxon_rf.topResults.freq; cpf10_wilcoxon_svm.topResults.freq];
    cpf10m_counts = [cpf10_mrmr_qda.topResults.freq; cpf10_mrmr_lda.topResults.freq; cpf10_mrmr_rf.topResults.freq; cpf10_mrmr_svm.topResults.freq];
    
    cpf10w = [cpf10w_indices cpf10w_counts];
    cpf10m = [cpf10m_indices cpf10m_counts];

    cpf15w_indices = [cpf15_wilcoxon_qda.topResults.indices; cpf15_wilcoxon_lda.topResults.indices; cpf15_wilcoxon_rf.topResults.indices; cpf15_wilcoxon_svm.topResults.indices];
    cpf15m_indices = [cpf15_mrmr_qda.topResults.indices; cpf15_mrmr_lda.topResults.indices; cpf15_mrmr_rf.topResults.indices; cpf15_mrmr_svm.topResults.indices];
    cpf15w_counts = [cpf15_wilcoxon_qda.topResults.freq; cpf15_wilcoxon_lda.topResults.freq; cpf15_wilcoxon_rf.topResults.freq; cpf15_wilcoxon_svm.topResults.freq];
    cpf15m_counts = [cpf15_mrmr_qda.topResults.freq; cpf15_mrmr_lda.topResults.freq; cpf15_mrmr_rf.topResults.freq; cpf15_mrmr_svm.topResults.freq];
   
    cpf15w = [cpf15w_indices cpf15w_counts];
    cpf15m = [cpf15m_indices cpf15m_counts];

    ctw_indices = [ct_wilcoxon_qda.topResults.indices; ct_wilcoxon_lda.topResults.indices; ct_wilcoxon_rf.topResults.indices; ct_wilcoxon_svm.topResults.indices];
    ctm_indices = [ct_mrmr_qda.topResults.indices; ct_mrmr_lda.topResults.indices; ct_mrmr_rf.topResults.indices; ct_mrmr_svm.topResults.indices];
    ctw_counts = [ct_wilcoxon_qda.topResults.freq; ct_wilcoxon_lda.topResults.freq; ct_wilcoxon_rf.topResults.freq; ct_wilcoxon_svm.topResults.freq];
    ctm_counts = [ct_mrmr_qda.topResults.freq; ct_mrmr_lda.topResults.freq; ct_mrmr_rf.topResults.freq; ct_mrmr_svm.topResults.freq];
    
    ctw = [ctw_indices ctw_counts];
    ctm = [ctm_indices ctm_counts];
end


%%
afw_feature_data = containers.Map('KeyType', 'double', 'ValueType', 'any');
    
for i = 1 : size(afw, 1)
    feature_index = afw(i, 1);
    count = afw(i, 2);
        
    % If the feature is already in the dictionary, update the count and sum
    if isKey(afw_feature_data, feature_index)
        afw_feature_data(feature_index) = [afw_feature_data(feature_index); count];
    else
        afw_feature_data(feature_index) = count;
    end
end

% Initialize a new array to store the results
afw = zeros(length(keys(afw_feature_data)), 2);

% Iterate through the dictionary to calculate the average and store in result_data
afw_keys_array = keys(afw_feature_data);

for i = 1:length(afw_keys_array)
    key = afw_keys_array{i};
    counts = afw_feature_data(key);
    avg_count = sum(counts) / length(counts);

    afw(i, 1) = key;
    afw(i, 2) = avg_count;
    afw = sortrows(afw, -2);
end

afw_feature_names = cell(length(afw_keys_array), 1);

for i = 1 : length(afw)
    indices = afw(:, 1);
    afw_feature_names{i} = axial_fat_feature_column_names{indices(i)};
end

fprintf('Here are the top axial fat, wilcoxon features: \n')
disp(afw_feature_names)
fprintf('Here are their frequencies: \n')
disp(afw(:, 2))

afm_feature_data = containers.Map('KeyType', 'double', 'ValueType', 'any');
    
for i = 1 : size(afm, 1)
    feature_index = afm(i, 1);
    count = afm(i, 2);
        
    % If the feature is already in the dictionary, update the count and sum
    if isKey(afm_feature_data, feature_index)
        afm_feature_data(feature_index) = [afm_feature_data(feature_index); count];
    else
        afm_feature_data(feature_index) = count;
    end
end

% Initialize a new array to store the results
afm = zeros(length(keys(afm_feature_data)), 2);

% Iterate through the dictionary to calculate the average and store in result_data
afm_keys_array = keys(afm_feature_data);

for i = 1:length(afm_keys_array)
    key = afm_keys_array{i};
    counts = afm_feature_data(key);
    avg_count = sum(counts) / length(counts);

    afm(i, 1) = key;
    afm(i, 2) = avg_count;
    afm = sortrows(afm, -2);
end

afm_feature_names = cell(length(afm_keys_array), 1);

for i = 1 : length(afm)
    indices = afm(:, 1);
    afm_feature_names{i} = axial_fat_feature_column_names{indices(i)};
end

fprintf('Here are the top axial fat, mRMR features: \n')
disp(afm_feature_names)
fprintf('Here are their frequencies: \n')
disp(afm(:, 2))

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

fprintf('Here are the top axial tumor, wilcoxon features: \n')
disp(atw_feature_names)
fprintf('Here are their frequencies: \n')
disp(atw(:, 2))

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

fprintf('Here are the top axial tumor, mRMR features: \n')
disp(atm_feature_names)
fprintf('Here are their frequencies: \n')
disp(atm(:, 2))

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

fprintf('Here are the top axial proxfat5, wilcoxon features: \n')
disp(apf5w_feature_names)
fprintf('Here are their frequencies: \n')
disp(apf5w(:, 2))

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

fprintf('Here are the top axial proxfat5, mRMR features: \n')
disp(apf5m_feature_names)
fprintf('Here are their frequencies: \n')
disp(apf5m(:, 2))

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

fprintf('Here are the top axial proxfat10, wilcoxon features: \n')
disp(apf10w_feature_names)
fprintf('Here are their frequencies: \n')
disp(apf10w(:, 2))

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

fprintf('Here are the top axial proxfat10, mRMR features: \n')
disp(apf10m_feature_names)
fprintf('Here are their frequencies: \n')
disp(apf10m(:, 2))

%%
apf15w_feature_data = containers.Map('KeyType', 'double', 'ValueType', 'any');
    
for i = 1 : size(apf15w, 1)
    feature_index = apf15w(i, 1);
    count = apf15w(i, 2);
        
    % If the feature is already in the dictionary, update the count and sum
    if isKey(apf15w_feature_data, feature_index)
        apf15w_feature_data(feature_index) = [apf15w_feature_data(feature_index); count];
    else
        apf15w_feature_data(feature_index) = count;
    end
end

% Initialize a new array to store the results
apf15w = zeros(length(keys(apf15w_feature_data)), 2);

% Iterate through the dictionary to calculate the average and store in result_data
apf15w_keys_array = keys(apf15w_feature_data);

for i = 1:length(apf15w_keys_array)
    key = apf15w_keys_array{i};
    counts = apf15w_feature_data(key);
    avg_count = sum(counts) / length(counts);

    apf15w(i, 1) = key;
    apf15w(i, 2) = avg_count;
    apf15w = sortrows(apf15w, -2);
end

apf15w_feature_names = cell(length(apf15w_keys_array), 1);

for i = 1 : length(apf15w)
    indices = apf15w(:, 1);
    apf15w_feature_names{i} = axial_proxfat15_feature_column_names{indices(i)};
end

fprintf('Here are the top axial proxfat15, wilcoxon features: \n')
disp(apf15w_feature_names)
fprintf('Here are their frequencies: \n')
disp(apf15w(:, 2))

apf15m_feature_data = containers.Map('KeyType', 'double', 'ValueType', 'any');
    
for i = 1 : size(apf15m, 1)
    feature_index = apf15m(i, 1);
    count = apf15m(i, 2);
        
    % If the feature is already in the dictionary, update the count and sum
    if isKey(apf15m_feature_data, feature_index)
        apf15m_feature_data(feature_index) = [apf15m_feature_data(feature_index); count];
    else
        apf15m_feature_data(feature_index) = count;
    end
end

% Initialize a new array to store the results
apf15m = zeros(length(keys(apf15m_feature_data)), 2);

% Iterate through the dictionary to calculate the average and store in result_data
apf15m_keys_array = keys(apf15m_feature_data);

for i = 1:length(apf15m_keys_array)
    key = apf15m_keys_array{i};
    counts = apf15m_feature_data(key);
    avg_count = sum(counts) / length(counts);

    apf15m(i, 1) = key;
    apf15m(i, 2) = avg_count;
    apf15m = sortrows(apf15m, -2);
end

apf15m_feature_names = cell(length(apf15m_keys_array), 1);

for i = 1 : length(apf15m)
    indices = apf15m(:, 1);
    apf15m_feature_names{i} = axial_proxfat15_feature_column_names{indices(i)};
end

fprintf('Here are the top axial proxfat15, mRMR features: \n')
disp(apf15m_feature_names)
fprintf('Here are their frequencies: \n')
disp(apf15m(:, 2))

cfw_feature_data = containers.Map('KeyType', 'double', 'ValueType', 'any');
    
for i = 1 : size(cfw, 1)
    feature_index = cfw(i, 1);
    count = cfw(i, 2);
        
    % If the feature is already in the dictionary, update the count and sum
    if isKey(cfw_feature_data, feature_index)
        cfw_feature_data(feature_index) = [cfw_feature_data(feature_index); count];
    else
        cfw_feature_data(feature_index) = count;
    end
end

% Initialize a new array to store the results
cfw = zeros(length(keys(cfw_feature_data)), 2);

% Iterate through the dictionary to calculate the average and store in result_data
cfw_keys_array = keys(cfw_feature_data);

for i = 1:length(cfw_keys_array)
    key = cfw_keys_array{i};
    counts = cfw_feature_data(key);
    avg_count = sum(counts) / length(counts);

    cfw(i, 1) = key;
    cfw(i, 2) = avg_count;
    cfw = sortrows(cfw, -2);
end

cfw_feature_names = cell(length(cfw_keys_array), 1);

for i = 1 : length(cfw)
    indices = cfw(:, 1);
    cfw_feature_names{i} = coronal_fat_feature_column_names{indices(i)};
end

fprintf('Here are the top coronal fat, wilcoxon features: \n')
disp(cfw_feature_names)
fprintf('Here are their frequencies: \n')
disp(cfw(:, 2))

cfm_feature_data = containers.Map('KeyType', 'double', 'ValueType', 'any');
    
for i = 1 : size(cfm, 1)
    feature_index = cfm(i, 1);
    count = cfm(i, 2);
        
    % If the feature is already in the dictionary, update the count and sum
    if isKey(cfm_feature_data, feature_index)
        cfm_feature_data(feature_index) = [cfm_feature_data(feature_index); count];
    else
        cfm_feature_data(feature_index) = count;
    end
end

% Initialize a new array to store the results
cfm = zeros(length(keys(cfm_feature_data)), 2);

% Iterate through the dictionary to calculate the average and store in result_data
cfm_keys_array = keys(cfm_feature_data);

for i = 1:length(cfm_keys_array)
    key = cfm_keys_array{i};
    counts = cfm_feature_data(key);
    avg_count = sum(counts) / length(counts);

    cfm(i, 1) = key;
    cfm(i, 2) = avg_count;
    cfm = sortrows(cfm, -2);
end

cfm_feature_names = cell(length(cfm_keys_array), 1);

for i = 1 : length(cfm)
    indices = cfm(:, 1);
    cfm_feature_names{i} = coronal_fat_feature_column_names{indices(i)};
end

fprintf('Here are the top coronal fat, mRMR features: \n')
disp(cfm_feature_names)
fprintf('Here are their frequencies: \n')
disp(cfm(:, 2))



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

fprintf('Here are the top coronal tumor, wilcoxon features: \n')
disp(ctw_feature_names)
fprintf('Here are their frequencies: \n')
disp(ctw(:, 2))

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

fprintf('Here are the top coronal tumor, mRMR features: \n')
disp(ctm_feature_names)
fprintf('Here are their frequencies: \n')
disp(ctm(:, 2))



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

fprintf('Here are the top coronal proxfat5, wilcoxon features: \n')
disp(cpf5w_feature_names)
fprintf('Here are their frequencies: \n')
disp(cpf5w(:, 2))

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

fprintf('Here are the top coronal proxfat5, mRMR features: \n')
disp(cpf5m_feature_names)
fprintf('Here are their frequencies: \n')
disp(cpf5m(:, 2))



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

fprintf('Here are the top coronal proxfat10, wilcoxon features: \n')
disp(cpf10w_feature_names)
fprintf('Here are their frequencies: \n')
disp(cpf10w(:, 2))

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

fprintf('Here are the top coronal proxfat10, mRMR features: \n')
disp(cpf10m_feature_names)
fprintf('Here are their frequencies: \n')
disp(cpf10m(:, 2))

cpf15w_feature_data = containers.Map('KeyType', 'double', 'ValueType', 'any');
    
for i = 1 : size(cpf15w, 1)
    feature_index = cpf15w(i, 1);
    count = cpf15w(i, 2);
        
    % If the feature is already in the dictionary, update the count and sum
    if isKey(cpf15w_feature_data, feature_index)
        cpf15w_feature_data(feature_index) = [cpf15w_feature_data(feature_index); count];
    else
        cpf15w_feature_data(feature_index) = count;
    end
end

% Initialize a new array to store the results
cpf15w = zeros(length(keys(cpf15w_feature_data)), 2);

% Iterate through the dictionary to calculate the average and store in result_data
cpf15w_keys_array = keys(cpf15w_feature_data);

for i = 1:length(cpf15w_keys_array)
    key = cpf15w_keys_array{i};
    counts = cpf15w_feature_data(key);
    avg_count = sum(counts) / length(counts);

    cpf15w(i, 1) = key;
    cpf15w(i, 2) = avg_count;
    cpf15w = sortrows(cpf15w, -2);
end

cpf15w_feature_names = cell(length(cpf15w_keys_array), 1);

for i = 1 : length(cpf15w)
    indices = cpf15w(:, 1);
    cpf15w_feature_names{i} = coronal_proxfat15_feature_column_names{indices(i)};
end

fprintf('Here are the top coronal proxfat15, wilcoxon features: \n')
disp(cpf15w_feature_names)
fprintf('Here are their frequencies: \n')
disp(cpf15w(:, 2))

cpf15m_feature_data = containers.Map('KeyType', 'double', 'ValueType', 'any');
    
for i = 1 : size(cpf15m, 1)
    feature_index = cpf15m(i, 1);
    count = cpf15m(i, 2);
        
    % If the feature is already in the dictionary, update the count and sum
    if isKey(cpf15m_feature_data, feature_index)
        cpf15m_feature_data(feature_index) = [cpf15m_feature_data(feature_index); count];
    else
        cpf15m_feature_data(feature_index) = count;
    end
end

% Initialize a new array to store the results
cpf15m = zeros(length(keys(cpf15m_feature_data)), 2);

% Iterate through the dictionary to calculate the average and store in result_data
cpf15m_keys_array = keys(cpf15m_feature_data);

for i = 1:length(cpf15m_keys_array)
    key = cpf15m_keys_array{i};
    counts = cpf15m_feature_data(key);
    avg_count = sum(counts) / length(counts);

    cpf15m(i, 1) = key;
    cpf15m(i, 2) = avg_count;
    cpf15m = sortrows(cpf15m, -2);
end

cpf15m_feature_names = cell(length(cpf15m_keys_array), 1);

for i = 1 : length(cpf15m)
    indices = cpf15m(:, 1);
    cpf15m_feature_names{i} = coronal_proxfat15_feature_column_names{indices(i)};
end

fprintf('Here are the top coronal proxfat15, mRMR features: \n')
disp(cpf15m_feature_names)
fprintf('Here are their frequencies: \n')
disp(cpf15m(:, 2))