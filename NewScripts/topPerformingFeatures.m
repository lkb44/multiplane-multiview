clc;
clear;
close all;
results_root = '/Users/leobao/Documents/MissingCollageResults/';

experiments = dir(fullfile(results_root));
theFilesvol = {experiments.name};

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
    afw_counts = [af_wilcoxon_qda.topResults.counts; af_wilcoxon_lda.topResults.counts; af_wilcoxon_rf.topResults.counts; af_wilcoxon_svm.topResults.counts];
    afm_counts = [af_mrmr_qda.topResults.counts; af_mrmr_lda.topResults.counts; af_mrmr_rf.topResults.counts; af_mrmr_svm.topResults.counts];
    
    afw = [afw_indices afw_counts];
    afm = [afm_indices afm_counts];

    apf5w_indices = [apf5_wilcoxon_qda.topResults.indices; apf5_wilcoxon_lda.topResults.indices; apf5_wilcoxon_rf.topResults.indices; apf5_wilcoxon_svm.topResults.indices];
    apf5m_indices = [apf5_mrmr_rf.topResults.indices; apf5_mrmr_svm.topResults.indices];
    apf5w_counts = [apf5_wilcoxon_qda.topResults.counts; apf5_wilcoxon_lda.topResults.counts; apf5_wilcoxon_rf.topResults.counts; apf5_wilcoxon_svm.topResults.counts];
    apf5m_counts = [apf5_mrmr_rf.topResults.counts; apf5_mrmr_svm.topResults.counts];
    
    apf5w = [apf5w_indices apf5w_counts];
    apf5m = [apf5m_indices apf5m_counts];

    apf10w_indices = [apf10_wilcoxon_qda.topResults.indices; apf10_wilcoxon_lda.topResults.indices; apf10_wilcoxon_rf.topResults.indices; apf10_wilcoxon_svm.topResults.indices];
    apf10m_indices = [apf10_mrmr_qda.topResults.indices; apf10_mrmr_lda.topResults.indices; apf10_mrmr_rf.topResults.indices; apf10_mrmr_svm.topResults.indices];
    apf10w_counts = [apf10_wilcoxon_qda.topResults.counts; apf10_wilcoxon_lda.topResults.counts; apf10_wilcoxon_rf.topResults.counts; apf10_wilcoxon_svm.topResults.counts];
    apf10m_counts = [apf10_mrmr_qda.topResults.counts; apf10_mrmr_lda.topResults.counts; apf10_mrmr_rf.topResults.counts; apf10_mrmr_svm.topResults.counts];
    
    apf10w = [apf10w_indices apf10w_counts];
    apf10m = [apf10m_indices apf10m_counts];

    apf15w_indices = [apf15_wilcoxon_qda.topResults.indices; apf15_wilcoxon_lda.topResults.indices; apf15_wilcoxon_rf.topResults.indices; apf15_wilcoxon_svm.topResults.indices];
    apf15m_indices = [apf15_mrmr_qda.topResults.indices; apf15_mrmr_lda.topResults.indices; apf15_mrmr_rf.topResults.indices; apf15_mrmr_svm.topResults.indices];
    apf15w_counts = [apf15_wilcoxon_qda.topResults.counts; apf15_wilcoxon_lda.topResults.counts; apf15_wilcoxon_rf.topResults.counts; apf15_wilcoxon_svm.topResults.counts];
    apf15m_counts = [apf15_mrmr_qda.topResults.counts; apf15_mrmr_lda.topResults.counts; apf15_mrmr_rf.topResults.counts; apf15_mrmr_svm.topResults.counts];
    
    apf15w = [apf15w_indices apf15w_counts];
    apf15m = [apf15m_indices apf15m_counts];

    atw_indices = [at_wilcoxon_qda.topResults.indices; at_wilcoxon_lda.topResults.indices; at_wilcoxon_rf.topResults.indices; at_wilcoxon_svm.topResults.indices];
    atm_indices = [at_mrmr_qda.topResults.indices; at_mrmr_lda.topResults.indices; at_mrmr_rf.topResults.indices; at_mrmr_svm.topResults.indices];
    atw_counts = [at_wilcoxon_qda.topResults.counts; at_wilcoxon_lda.topResults.counts; at_wilcoxon_rf.topResults.counts; at_wilcoxon_svm.topResults.counts];
    atm_counts = [at_mrmr_qda.topResults.counts; at_mrmr_lda.topResults.counts; at_mrmr_rf.topResults.counts; at_mrmr_svm.topResults.counts];
    
    atw = [atw_indices atw_counts];
    atm = [atm_indices atm_counts];

    cfw_indices = [cf_wilcoxon_qda.topResults.indices; cf_wilcoxon_lda.topResults.indices; cf_wilcoxon_rf.topResults.indices; cf_wilcoxon_svm.topResults.indices];
    cfm_indices = [cf_mrmr_qda.topResults.indices; cf_mrmr_lda.topResults.indices; cf_mrmr_rf.topResults.indices; cf_mrmr_svm.topResults.indices];
    cfw_counts = [cf_wilcoxon_qda.topResults.counts; cf_wilcoxon_lda.topResults.counts; cf_wilcoxon_rf.topResults.counts; cf_wilcoxon_svm.topResults.counts];
    cfm_counts = [cf_mrmr_qda.topResults.counts; cf_mrmr_lda.topResults.counts; cf_mrmr_rf.topResults.counts; cf_mrmr_svm.topResults.counts];
    
    cfw = [cfw_indices cfw_counts];
    cfm = [cfm_indices cfm_counts];

    cpf5w_indices = [cpf5_wilcoxon_qda.topResults.indices; cpf5_wilcoxon_lda.topResults.indices; cpf5_wilcoxon_rf.topResults.indices; cpf5_wilcoxon_svm.topResults.indices];
    cpf5m_indices = [cpf5_mrmr_qda.topResults.indices; cpf5_mrmr_lda.topResults.indices; cpf5_mrmr_rf.topResults.indices; cpf5_mrmr_svm.topResults.indices];
    cpf5w_counts = [cpf5_wilcoxon_qda.topResults.counts; cpf5_wilcoxon_lda.topResults.counts; cpf5_wilcoxon_rf.topResults.counts; cpf5_wilcoxon_svm.topResults.counts];
    cpf5m_counts = [cpf5_mrmr_qda.topResults.counts; cpf5_mrmr_lda.topResults.counts; cpf5_mrmr_rf.topResults.counts; cpf5_mrmr_svm.topResults.counts];
    
    cpf5w = [cpf5w_indices cpf5w_counts];
    cpf5m = [cpf5m_indices cpf5m_counts];

    cpf10w_indices = [cpf10_wilcoxon_qda.topResults.indices; cpf10_wilcoxon_lda.topResults.indices; cpf10_wilcoxon_rf.topResults.indices; cpf10_wilcoxon_svm.topResults.indices];
    cpf10m_indices = [cpf10_mrmr_qda.topResults.indices; cpf10_mrmr_lda.topResults.indices; cpf10_mrmr_rf.topResults.indices; cpf10_mrmr_svm.topResults.indices];
    cpf10w_counts = [cpf10_wilcoxon_qda.topResults.counts; cpf10_wilcoxon_lda.topResults.counts; cpf10_wilcoxon_rf.topResults.counts; cpf10_wilcoxon_svm.topResults.counts];
    cpf10m_counts = [cpf10_mrmr_qda.topResults.counts; cpf10_mrmr_lda.topResults.counts; cpf10_mrmr_rf.topResults.counts; cpf10_mrmr_svm.topResults.counts];
    
    cpf10w = [cpf10w_indices cpf10w_counts];
    cpf10m = [cpf10m_indices cpf10m_counts];

    cpf15w_indices = [cpf15_wilcoxon_qda.topResults.indices; cpf15_wilcoxon_lda.topResults.indices; cpf15_wilcoxon_rf.topResults.indices; cpf15_wilcoxon_svm.topResults.indices];
    cpf15m_indices = [cpf15_mrmr_qda.topResults.indices; cpf15_mrmr_lda.topResults.indices; cpf15_mrmr_rf.topResults.indices; cpf15_mrmr_svm.topResults.indices];
    cpf15w_counts = [cpf15_wilcoxon_qda.topResults.counts; cpf15_wilcoxon_lda.topResults.counts; cpf15_wilcoxon_rf.topResults.counts; cpf15_wilcoxon_svm.topResults.counts];
    cpf15m_counts = [cpf15_mrmr_qda.topResults.counts; cpf15_mrmr_lda.topResults.counts; cpf15_mrmr_rf.topResults.counts; cpf15_mrmr_svm.topResults.counts];
   
    cpf15w = [cpf15w_indices cpf15w_counts];
    cpf15m = [cpf15m_indices cpf15m_counts];

    ctw_indices = [ct_wilcoxon_qda.topResults.indices; ct_wilcoxon_lda.topResults.indices; ct_wilcoxon_rf.topResults.indices; ct_wilcoxon_svm.topResults.indices];
    ctm_indices = [ct_mrmr_qda.topResults.indices; ct_mrmr_lda.topResults.indices; ct_mrmr_rf.topResults.indices; ct_mrmr_svm.topResults.indices];
    ctw_counts = [ct_wilcoxon_qda.topResults.counts; ct_wilcoxon_lda.topResults.counts; ct_wilcoxon_rf.topResults.counts; ct_wilcoxon_svm.topResults.counts];
    ctm_counts = [ct_mrmr_qda.topResults.counts; ct_mrmr_lda.topResults.counts; ct_mrmr_rf.topResults.counts; ct_mrmr_svm.topResults.counts];
    
    ctw = [ctw_indices ctw_counts];
    ctm = [ctm_indices ctm_counts];
end

features = {cfw, cfm, cpf5w, cpf5m, cpf10w, cpf10m, cpf15w, cpf15m, ctw, ctm};

for k = 1 : length(features)
    % Create a dictionary to store the counts and sums for each feature
    feature_data = containers.Map('KeyType', 'double', 'ValueType', 'any');
    
    % Iterate through the data to update the counts and sums
    for i = 1 : size(features{k}, 1)
        feature_index = features{k}(i, 1);
        count = features{k}(i, 2);
        
        % If the feature is already in the dictionary, update the count and sum
        if isKey(feature_data, feature_index)
            feature_data(feature_index) = [feature_data(feature_index); count];
        else
            feature_data(feature_index) = count;
        end
    end
    
    % Initialize a new array to store the results
    features{k} = zeros(length(keys(feature_data)), 2);
    
    % Iterate through the dictionary to calculate the average and store in result_data
    keys_array = keys(feature_data);
    for i = 1:length(keys_array)
        key = keys_array{i};
        counts = feature_data(key);
        avg_count = sum(counts) / length(counts);
        
        features{k}(i, 1) = key;
        features{k}(i, 2) = avg_count;
    end
    
    % Display the result
    disp(features{k});

end