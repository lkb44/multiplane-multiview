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
end
best_fat_roi = {'Proximal_Fat10'};
rois = {'Proximal_Fat10', 'Tumor'};
plane = {'Multi-Plane'}; % Can be 'Multi-Plane', 'Axial', or 'Coronal
scheme = {'wilcoxon_rf'};

if strcmp(scheme(1), 'wilcoxon_qda')
    selector = {'wilcoxon'};
    classifier = {'qda'};
elseif strcmp(scheme(1), 'wilcoxon_lda')
    selector = {'wilcoxon'};
    classifier = {'lda'};
elseif strcmp(scheme(1), 'wilcoxon_rf')
    selector = {'wilcoxon'};
    classifier = {'rf'};
elseif strcmp(scheme(1), 'mrmr_qda')
    selector = {'mrmr'};
    classifier = {'qda'};
elseif strcmp(scheme(1), 'mrmr_lda')
    selector = {'mrmr'};
    classifier = {'lda'};
elseif strcmp(scheme(1), 'mrmr_rf')
    selector = {'mrmr'};
    classifier = {'rf'};
end

if length(rois) == 2
    region = {'multiregion'};
elseif length(rois) == 1
    if strcmp(rois{1}, 'Proximal_Fat5')
        region = {'proxfat5'};
    elseif strcmp(rois{1}, 'Proximal_Fat10')
        region = {'proxfat10'};
    elseif strcmp(rois{1}, 'Fat')
        region = {'fat'};
    elseif strcmp(rois{1}, 'Tumor')
        region = {'tumor'};
    end
end

if strcmp(plane{1}, 'Multi-Plane')
    view = {'multiplane'};
elseif strcmp(plane{1}, 'Axial')
    view = {'axial'};
elseif strcmp(plane{1}, 'Coronal')
    view = {'coronal'};
end

results_root = '/Users/leobao/Documents/MultiPlanePipeline/AACR2023/Results/';
experiments = dir(fullfile(results_root));
theFilesvol = {experiments.name};

for k = 1:length(theFilesvol)
    P = theFilesvol;
    P = P(~startsWith(P, '.'));
    axial_proxfat5_wilcoxon_qda_index = find(contains(P, 'Axial_proxfat5_only_wilcoxon_qda'));
    axial_proxfat5_wilcoxon_lda_index = find(contains(P, 'Axial_proxfat5_only_wilcoxon_lda'));
    axial_proxfat5_wilcoxon_rf_index = find(contains(P, 'Axial_proxfat5_only_wilcoxon_rf'));
    % axial_proxfat5_mrmr_lda_index = find(contains(P, 'Axial_proxfat5_only_mrmr_lda'));
    % axial_proxfat5_mrmr_rf_index = find(contains(P, 'Axial_proxfat5_only_mrmr_rf'));
    coronal_proxfat5_wilcoxon_qda_index = find(contains(P, 'Coronal_proxfat5_only_wilcoxon_qda'));
    coronal_proxfat5_wilcoxon_lda_index = find(contains(P, 'Coronal_proxfat5_only_wilcoxon_lda'));
    coronal_proxfat5_wilcoxon_rf_index = find(contains(P, 'Coronal_proxfat5_only_wilcoxon_rf'));
    % coronal_proxfat5_mrmr_qda_index = find(contains(P, 'Coronal_proxfat5_only_mrmr_qda'));
    % coronal_proxfat5_mrmr_lda_index = find(contains(P, 'Coronal_proxfat5_only_mrmr_lda'));
    % coronal_proxfat5_mrmr_rf_index = find(contains(P, 'Coronal_proxfat5_only_mrmr_rf'));
    axial_proxfat10_wilcoxon_qda_index = find(contains(P, 'Axial_proxfat10_only_wilcoxon_qda'));
    axial_proxfat10_wilcoxon_lda_index = find(contains(P, 'Axial_proxfat10_only_wilcoxon_lda'));
    axial_proxfat10_wilcoxon_rf_index = find(contains(P, 'Axial_proxfat10_only_wilcoxon_rf'));
    % axial_proxfat10_mrmr_qda_index = find(contains(P, 'Axial_proxfat10_only_mrmr_qda'));
    % axial_proxfat10_mrmr_lda_index = find(contains(P, 'Axial_proxfat10_only_mrmr_lda'));
    % axial_proxfat10_mrmr_rf_index = find(contains(P, 'Axial_proxfat10_only_mrmr_rf'));
    coronal_proxfat10_wilcoxon_qda_index = find(contains(P, 'Coronal_proxfat10_only_wilcoxon_qda'));
    coronal_proxfat10_wilcoxon_lda_index = find(contains(P, 'Coronal_proxfat10_only_wilcoxon_lda'));
    coronal_proxfat10_wilcoxon_rf_index = find(contains(P, 'Coronal_proxfat10_only_wilcoxon_rf'));
    % coronal_proxfat10_mrmr_qda_index = find(contains(P, 'Coronal_proxfat10_only_mrmr_qda'));
    % coronal_proxfat10_mrmr_lda_index = find(contains(P, 'Coronal_proxfat10_only_mrmr_lda'));
    % coronal_proxfat10_mrmr_rf_index = find(contains(P, 'Coronal_proxfat10_only_mrmr_rf'));
    axial_tumor_wilcoxon_qda_index = find(contains(P, 'Axial_tumor_only_wilcoxon_qda'));
    axial_tumor_wilcoxon_lda_index = find(contains(P, 'Axial_tumor_only_wilcoxon_lda'));
    axial_tumor_wilcoxon_rf_index = find(contains(P, 'Axial_tumor_only_wilcoxon_rf'));
    % axial_tumor_mrmr_qda_index = find(contains(P, 'Axial_tumor_only_mrmr_qda'));
    % axial_tumor_mrmr_lda_index = find(contains(P, 'Axial_tumor_only_mrmr_lda'));
    % axial_tumor_mrmr_rf_index = find(contains(P, 'Axial_tumor_only_mrmr_rf'));
    coronal_tumor_wilcoxon_qda_index = find(contains(P, 'Coronal_tumor_only_wilcoxon_qda'));
    coronal_tumor_wilcoxon_lda_index = find(contains(P, 'Coronal_tumor_only_wilcoxon_lda'));
    coronal_tumor_wilcoxon_rf_index = find(contains(P, 'Coronal_tumor_only_wilcoxon_rf'));
    % coronal_tumor_mrmr_qda_index = find(contains(P, 'Coronal_tumor_only_mrmr_qda'));
    % coronal_tumor_mrmr_lda_index = find(contains(P, 'Coronal_tumor_only_mrmr_lda'));
    % coronal_tumor_mrmr_rf_index = find(contains(P, 'Coronal_tumor_only_mrmr_rf'));
end

%%
matrix_root_path = '/Users/leobao/Documents/MultiPlanePipeline/AACR2023/';

if strcmp(selector{1}, 'wilcoxon')
    if strcmp(classifier{1}, 'qda')
        if strcmp(plane{1}, 'Multi-Plane')
            if length(rois) == 2
                if strcmp(best_fat_roi{1}, 'Proximal_Fat5')
                    axial_fat_feats_train = load(fullfile(results_root, P{axial_proxfat5_wilcoxon_qda_index}, 'top5_training_feats.mat')).train_feats;
                    axial_fat_feats_test1 = load(fullfile(results_root, P{axial_proxfat5_wilcoxon_qda_index}, 'top5_test1_feats.mat')).test1_feats;
                    axial_fat_feats_test2 = load(fullfile(results_root, P{axial_proxfat5_wilcoxon_qda_index}, 'top5_test2_feats.mat')).test2_feats;
                    axial_tumor_feats_train = load(fullfile(results_root, P{axial_tumor_wilcoxon_qda_index}, 'top5_training_feats.mat')).train_feats;
                    axial_tumor_feats_test1 = load(fullfile(results_root, P{axial_tumor_wilcoxon_qda_index}, 'top5_test1_feats.mat')).test1_feats;
                    axial_tumor_feats_test2 = load(fullfile(results_root, P{axial_tumor_wilcoxon_qda_index}, 'top5_test2_feats.mat')).test2_feats;
                    coronal_fat_feats_train = load(fullfile(results_root, P{coronal_proxfat5_wilcoxon_qda_index}, 'top5_training_feats.mat')).train_feats;
                    coronal_fat_feats_test1 = load(fullfile(results_root, P{coronal_proxfat5_wilcoxon_qda_index}, 'top5_test1_feats.mat')).test1_feats;
                    coronal_fat_feats_test2 = load(fullfile(results_root, P{coronal_proxfat5_wilcoxon_qda_index}, 'top5_test2_feats.mat')).test2_feats;
                    coronal_tumor_feats_train = load(fullfile(results_root, P{coronal_tumor_wilcoxon_qda_index}, 'top5_training_feats.mat')).train_feats;
                    coronal_tumor_feats_test1 = load(fullfile(results_root, P{coronal_tumor_wilcoxon_qda_index}, 'top5_test1_feats.mat')).test1_feats;
                    coronal_tumor_feats_test2 = load(fullfile(results_root, P{coronal_tumor_wilcoxon_qda_index}, 'top5_test2_feats.mat')).test2_feats;
                    axial_fat_feat_names = load(fullfile(results_root, P{axial_proxfat5_wilcoxon_qda_index}, 'Top_Features_Information.mat')).topResults.names';
                    coronal_fat_feat_names = load(fullfile(results_root, P{coronal_proxfat5_wilcoxon_qda_index}, 'Top_Features_Information.mat')).topResults.names';
                    axial_tumor_feat_names = load(fullfile(results_root, P{axial_tumor_wilcoxon_qda_index}, 'Top_Features_Information.mat')).topResults.names';
                    coronal_tumor_feat_names = load(fullfile(results_root, P{coronal_tumor_wilcoxon_qda_index}, 'Top_Features_Information.mat')).topResults.names';
                    feature_column_names = [axial_fat_feat_names, coronal_fat_feat_names, axial_tumor_feat_names, coronal_tumor_feat_names];
                    features_training = [axial_fat_feats_train, coronal_fat_feats_train, axial_tumor_feats_train, coronal_tumor_feats_train];
                    features_holdout1 = [axial_fat_feats_test1, coronal_fat_feats_test1, axial_tumor_feats_test1, coronal_tumor_feats_test1];
                    features_holdout2 = [axial_fat_feats_test2, coronal_fat_feats_test2, axial_tumor_feats_test2, coronal_tumor_feats_test2];
                elseif strcmp(best_fat_roi{1}, 'Proximal_Fat10')
                    axial_fat_feats_train = load(fullfile(results_root, P{axial_proxfat10_wilcoxon_qda_index}, 'top5_training_feats.mat')).train_feats;
                    axial_fat_feats_test1 = load(fullfile(results_root, P{axial_proxfat10_wilcoxon_qda_index}, 'top5_test1_feats.mat')).test1_feats;
                    axial_fat_feats_test2 = load(fullfile(results_root, P{axial_proxfat10_wilcoxon_qda_index}, 'top5_test2_feats.mat')).test2_feats;
                    axial_tumor_feats_train = load(fullfile(results_root, P{axial_tumor_wilcoxon_qda_index}, 'top5_training_feats.mat')).train_feats;
                    axial_tumor_feats_test1 = load(fullfile(results_root, P{axial_tumor_wilcoxon_qda_index}, 'top5_test1_feats.mat')).test1_feats;
                    axial_tumor_feats_test2 = load(fullfile(results_root, P{axial_tumor_wilcoxon_qda_index}, 'top5_test2_feats.mat')).test2_feats;
                    coronal_fat_feats_train = load(fullfile(results_root, P{coronal_proxfat10_wilcoxon_qda_index}, 'top5_training_feats.mat')).train_feats;
                    coronal_fat_feats_test1 = load(fullfile(results_root, P{coronal_proxfat10_wilcoxon_qda_index}, 'top5_test1_feats.mat')).test1_feats;
                    coronal_fat_feats_test2 = load(fullfile(results_root, P{coronal_proxfat10_wilcoxon_qda_index}, 'top5_test2_feats.mat')).test2_feats;
                    coronal_tumor_feats_train = load(fullfile(results_root, P{coronal_tumor_wilcoxon_qda_index}, 'top5_training_feats.mat')).train_feats;
                    coronal_tumor_feats_test1 = load(fullfile(results_root, P{coronal_tumor_wilcoxon_qda_index}, 'top5_test1_feats.mat')).test1_feats;
                    coronal_tumor_feats_test2 = load(fullfile(results_root, P{coronal_tumor_wilcoxon_qda_index}, 'top5_test2_feats.mat')).test2_feats;
                    axial_fat_feat_names = load(fullfile(results_root, P{axial_proxfat5_wilcoxon_qda_index}, 'Top_Features_Information.mat')).topResults.names';
                    coronal_fat_feat_names = load(fullfile(results_root, P{coronal_proxfat5_wilcoxon_qda_index}, 'Top_Features_Information.mat')).topResults.names';
                    axial_tumor_feat_names = load(fullfile(results_root, P{axial_tumor_wilcoxon_qda_index}, 'Top_Features_Information.mat')).topResults.names';
                    coronal_tumor_feat_names = load(fullfile(results_root, P{coronal_tumor_wilcoxon_qda_index}, 'Top_Features_Information.mat')).topResults.names';
                    feature_column_names = [axial_fat_feat_names, coronal_fat_feat_names, axial_tumor_feat_names, coronal_tumor_feat_names];
                    features_training = [axial_fat_feats_train, coronal_fat_feats_train, axial_tumor_feats_train, coronal_tumor_feats_train];
                    features_holdout1 = [axial_fat_feats_test1, coronal_fat_feats_test1, axial_tumor_feats_test1, coronal_tumor_feats_test1];
                    features_holdout2 = [axial_fat_feats_test2, coronal_fat_feats_test2, axial_tumor_feats_test2, coronal_tumor_feats_test2];
                elseif strcmp(best_fat_roi{1}, 'Fat')
                    axial_fat_feats_train = load(fullfile(results_root, P{axial_fat_wilcoxon_qda_index}, 'top5_training_feats.mat')).train_feats;
                    axial_fat_feats_test1 = load(fullfile(results_root, P{axial_fat_wilcoxon_qda_index}, 'top5_test1_feats.mat')).test1_feats;
                    axial_fat_feats_test2 = load(fullfile(results_root, P{axial_fat_wilcoxon_qda_index}, 'top5_test2_feats.mat')).test2_feats;
                    axial_tumor_feats_train = load(fullfile(results_root, P{axial_tumor_wilcoxon_qda_index}, 'top5_training_feats.mat')).train_feats;
                    axial_tumor_feats_test1 = load(fullfile(results_root, P{axial_tumor_wilcoxon_qda_index}, 'top5_test1_feats.mat')).test1_feats;
                    axial_tumor_feats_test2 = load(fullfile(results_root, P{axial_tumor_wilcoxon_qda_index}, 'top5_test2_feats.mat')).test2_feats;
                    coronal_fat_feats_train = load(fullfile(results_root, P{coronal_fat_wilcoxon_qda_index}, 'top5_training_feats.mat')).train_feats;
                    coronal_fat_feats_test1 = load(fullfile(results_root, P{coronal_fat_wilcoxon_qda_index}, 'top5_test1_feats.mat')).test1_feats;
                    coronal_fat_feats_test2 = load(fullfile(results_root, P{coronal_fat_wilcoxon_qda_index}, 'top5_test2_feats.mat')).test2_feats;
                    coronal_tumor_feats_train = load(fullfile(results_root, P{coronal_tumor_wilcoxon_qda_index}, 'top5_training_feats.mat')).train_feats;
                    coronal_tumor_feats_test1 = load(fullfile(results_root, P{coronal_tumor_wilcoxon_qda_index}, 'top5_test1_feats.mat')).test1_feats;
                    coronal_tumor_feats_test2 = load(fullfile(results_root, P{coronal_tumor_wilcoxon_qda_index}, 'top5_test2_feats.mat')).test2_feats;
                    axial_fat_feat_names = load(fullfile(results_root, P{axial_fat_wilcoxon_qda_index}, 'Top_Features_Information.mat')).topResults.names';
                    coronal_fat_feat_names = load(fullfile(results_root, P{coronal_fat_wilcoxon_qda_index}, 'Top_Features_Information.mat')).topResults.names';
                    axial_tumor_feat_names = load(fullfile(results_root, P{axial_tumor_wilcoxon_qda_index}, 'Top_Features_Information.mat')).topResults.names';
                    coronal_tumor_feat_names = load(fullfile(results_root, P{coronal_tumor_wilcoxon_qda_index}, 'Top_Features_Information.mat')).topResults.names';
                    feature_column_names = [axial_fat_feat_names, coronal_fat_feat_names, axial_tumor_feat_names, coronal_tumor_feat_names];
                    features_training = [axial_fat_feats_train, coronal_fat_feats_train, axial_tumor_feats_train, coronal_tumor_feats_train];
                    features_holdout1 = [axial_fat_feats_test1, coronal_fat_feats_test1, axial_tumor_feats_test1, coronal_tumor_feats_test1];
                    features_holdout2 = [axial_fat_feats_test2, coronal_fat_feats_test2, axial_tumor_feats_test2, coronal_tumor_feats_test2];
                end
            elseif length(rois) == 1
                if strcmp(rois{1}, 'Proximal_Fat5')
                    axial_feats_train = load(fullfile(results_root, P{axial_proxfat5_wilcoxon_qda_index}, 'top5_training_feats.mat')).train_feats;
                    axial_feats_test1 = load(fullfile(results_root, P{axial_proxfat5_wilcoxon_qda_index}, 'top5_test1_feats.mat')).test1_feats;
                    axial_feats_test2 = load(fullfile(results_root, P{axial_proxfat5_wilcoxon_qda_index}, 'top5_test2_feats.mat')).test2_feats;
                    coronal_feats_train = load(fullfile(results_root, P{coronal_proxfat5_wilcoxon_qda_index}, 'top5_training_feats.mat')).train_feats;
                    coronal_feats_test1 = load(fullfile(results_root, P{coronal_proxfat5_wilcoxon_qda_index}, 'top5_test1_feats.mat')).test1_feats;
                    coronal_feats_test2 = load(fullfile(results_root, P{coronal_proxfat5_wilcoxon_qda_index}, 'top5_test2_feats.mat')).test2_feats;
                    axial_feat_names = load(fullfile(results_root, P{axial_proxfat5_wilcoxon_qda_index}, 'Top_Features_Information.mat')).topResults.names';
                    coronal_feat_names = load(fullfile(results_root, P{coronal_proxfat5_wilcoxon_qda_index}, 'Top_Features_Information.mat')).topResults.names';
                    feature_column_names = [axial_feat_names, coronal_feat_names];
                    features_training = [axial_feats_train, coronal_feats_train];
                    features_holdout1 = [axial_feats_test1, coronal_feats_test1];
                    features_holdout2 = [axial_feats_test2, coronal_feats_test2];
                elseif strcmp(rois{1}, 'Proximal_Fat10')
                    axial_feats_train = load(fullfile(results_root, P{axial_proxfat10_wilcoxon_qda_index}, 'top5_training_feats.mat')).train_feats;
                    axial_feats_test1 = load(fullfile(results_root, P{axial_proxfat10_wilcoxon_qda_index}, 'top5_test1_feats.mat')).test1_feats;
                    axial_feats_test2 = load(fullfile(results_root, P{axial_proxfat10_wilcoxon_qda_index}, 'top5_test2_feats.mat')).test2_feats;
                    coronal_feats_train = load(fullfile(results_root, P{coronal_proxfat10_wilcoxon_qda_index}, 'top5_training_feats.mat')).train_feats;
                    coronal_feats_test1 = load(fullfile(results_root, P{coronal_proxfat10_wilcoxon_qda_index}, 'top5_test1_feats.mat')).test1_feats;
                    coronal_feats_test2 = load(fullfile(results_root, P{coronal_proxfat10_wilcoxon_qda_index}, 'top5_test2_feats.mat')).test2_feats;
                    axial_feat_names = load(fullfile(results_root, P{axial_proxfat10_wilcoxon_qda_index}, 'Top_Features_Information.mat')).topResults.names';
                    coronal_feat_names = load(fullfile(results_root, P{coronal_proxfat10_wilcoxon_qda_index}, 'Top_Features_Information.mat')).topResults.names';
                    feature_column_names = [axial_feat_names, coronal_feat_names];
                    features_training = [axial_feats_train, coronal_feats_train];
                    features_holdout1 = [axial_feats_test1, coronal_feats_test1];
                    features_holdout2 = [axial_feats_test2, coronal_feats_test2];
                elseif strcmp(rois{1}, 'Fat')
                    axial_feats_train = load(fullfile(results_root, P{axial_fat_wilcoxon_qda_index}, 'top5_training_feats.mat')).train_feats;
                    axial_feats_test1 = load(fullfile(results_root, P{axial_fat_wilcoxon_qda_index}, 'top5_test1_feats.mat')).test1_feats;
                    axial_feats_test2 = load(fullfile(results_root, P{axial_fat_wilcoxon_qda_index}, 'top5_test2_feats.mat')).test2_feats;
                    coronal_feats_train = load(fullfile(results_root, P{coronal_fat_wilcoxon_qda_index}, 'top5_training_feats.mat')).train_feats;
                    coronal_feats_test1 = load(fullfile(results_root, P{coronal_fat_wilcoxon_qda_index}, 'top5_test1_feats.mat')).test1_feats;
                    coronal_feats_test2 = load(fullfile(results_root, P{coronal_fat_wilcoxon_qda_index}, 'top5_test2_feats.mat')).test2_feats;
                    axial_feat_names = load(fullfile(results_root, P{axial_fat_wilcoxon_qda_index}, 'Top_Features_Information.mat')).topResults.names';
                    coronal_feat_names = load(fullfile(results_root, P{coronal_fat_wilcoxon_qda_index}, 'Top_Features_Information.mat')).topResults.names';
                    feature_column_names = [axial_feat_names, coronal_feat_names];
                    features_training = [axial_feats_train, coronal_feats_train];
                    features_holdout1 = [axial_feats_test1, coronal_feats_test1];
                    features_holdout2 = [axial_feats_test2, coronal_feats_test2];
                elseif strcmp(rois{1}, 'Tumor')
                    axial_feats_train = load(fullfile(results_root, P{axial_tumor_wilcoxon_qda_index}, 'top5_training_feats.mat')).train_feats;
                    axial_feats_test1 = load(fullfile(results_root, P{axial_tumor_wilcoxon_qda_index}, 'top5_test1_feats.mat')).test1_feats;
                    axial_feats_test2 = load(fullfile(results_root, P{axial_tumor_wilcoxon_qda_index}, 'top5_test2_feats.mat')).test2_feats;
                    coronal_feats_train = load(fullfile(results_root, P{coronal_tumor_wilcoxon_qda_index}, 'top5_training_feats.mat')).train_feats;
                    coronal_feats_test1 = load(fullfile(results_root, P{coronal_tumor_wilcoxon_qda_index}, 'top5_test1_feats.mat')).test1_feats;
                    coronal_feats_test2 = load(fullfile(results_root, P{coronal_tumor_wilcoxon_qda_index}, 'top5_test2_feats.mat')).test2_feats;
                    axial_feat_names = load(fullfile(results_root, P{axial_tumor_wilcoxon_qda_index}, 'Top_Features_Information.mat')).topResults.names';
                    coronal_feat_names = load(fullfile(results_root, P{coronal_tumor_wilcoxon_qda_index}, 'Top_Features_Information.mat')).topResults.names';
                    feature_column_names = [axial_feat_names, coronal_feat_names];
                    features_training = [axial_feats_train, coronal_feats_train];
                    features_holdout1 = [axial_feats_test1, coronal_feats_test1];
                    features_holdout2 = [axial_feats_test2, coronal_feats_test2];
                end
            end
        elseif strcmp(plane{1}, 'Axial')
            if strcmp(best_fat_roi{1}, 'Proximal_Fat5')
                axial_fat_feats_train = load(fullfile(results_root, P{axial_proxfat5_wilcoxon_qda_index}, 'top5_training_feats.mat')).train_feats;
                axial_fat_feats_test1 = load(fullfile(results_root, P{axial_proxfat5_wilcoxon_qda_index}, 'top5_test1_feats.mat')).test1_feats;
                axial_fat_feats_test2 = load(fullfile(results_root, P{axial_proxfat5_wilcoxon_qda_index}, 'top5_test2_feats.mat')).test2_feats;
                axial_tumor_feats_train = load(fullfile(results_root, P{axial_tumor_wilcoxon_qda_index}, 'top5_training_feats.mat')).train_feats;
                axial_tumor_feats_test1 = load(fullfile(results_root, P{axial_tumor_wilcoxon_qda_index}, 'top5_test1_feats.mat')).test1_feats;
                axial_tumor_feats_test2 = load(fullfile(results_root, P{axial_tumor_wilcoxon_qda_index}, 'top5_test2_feats.mat')).test2_feats;
                axial_fat_feat_names = load(fullfile(results_root, P{axial_proxfat5_wilcoxon_qda_index}, 'Top_Features_Information.mat')).topResults.names';
                axial_tumor_feat_names = load(fullfile(results_root, P{axial_tumor_wilcoxon_qda_index}, 'Top_Features_Information.mat')).topResults.names';
                feature_column_names = [axial_fat_feat_names, axial_tumor_feat_names];
                features_training = [axial_fat_feats_train, axial_tumor_feats_train];
                features_holdout1 = [axial_fat_feats_test1, axial_tumor_feats_test1];
                features_holdout2 = [axial_fat_feats_test2, axial_tumor_feats_test2];
            elseif strcmp(best_fat_roi{1}, 'Proximal_Fat10')
                axial_fat_feats_train = load(fullfile(results_root, P{axial_proxfat10_wilcoxon_qda_index}, 'top5_training_feats.mat')).train_feats;
                axial_fat_feats_test1 = load(fullfile(results_root, P{axial_proxfat10_wilcoxon_qda_index}, 'top5_test1_feats.mat')).test1_feats;
                axial_fat_feats_test2 = load(fullfile(results_root, P{axial_proxfat10_wilcoxon_qda_index}, 'top5_test2_feats.mat')).test2_feats;
                axial_tumor_feats_train = load(fullfile(results_root, P{axial_tumor_wilcoxon_qda_index}, 'top5_training_feats.mat')).train_feats;
                axial_tumor_feats_test1 = load(fullfile(results_root, P{axial_tumor_wilcoxon_qda_index}, 'top5_test1_feats.mat')).test1_feats;
                axial_tumor_feats_test2 = load(fullfile(results_root, P{axial_tumor_wilcoxon_qda_index}, 'top5_test2_feats.mat')).test2_feats;
                axial_fat_feat_names = load(fullfile(results_root, P{axial_proxfat5_wilcoxon_qda_index}, 'Top_Features_Information.mat')).topResults.names';
                axial_tumor_feat_names = load(fullfile(results_root, P{axial_tumor_wilcoxon_qda_index}, 'Top_Features_Information.mat')).topResults.names';
                feature_column_names = [axial_fat_feat_names, axial_tumor_feat_names];
                features_training = [axial_fat_feats_train, axial_tumor_feats_train];
                features_holdout1 = [axial_fat_feats_test1, axial_tumor_feats_test1];
                features_holdout2 = [axial_fat_feats_test2, axial_tumor_feats_test2];
            elseif strcmp(best_fat_roi{1}, 'Fat')
                axial_fat_feats_train = load(fullfile(results_root, P{axial_fat_wilcoxon_qda_index}, 'top5_training_feats.mat')).train_feats;
                axial_fat_feats_test1 = load(fullfile(results_root, P{axial_fat_wilcoxon_qda_index}, 'top5_test1_feats.mat')).test1_feats;
                axial_fat_feats_test2 = load(fullfile(results_root, P{axial_fat_wilcoxon_qda_index}, 'top5_test2_feats.mat')).test2_feats;
                axial_tumor_feats_train = load(fullfile(results_root, P{axial_tumor_wilcoxon_qda_index}, 'top5_training_feats.mat')).train_feats;
                axial_tumor_feats_test1 = load(fullfile(results_root, P{axial_tumor_wilcoxon_qda_index}, 'top5_test1_feats.mat')).test1_feats;
                axial_tumor_feats_test2 = load(fullfile(results_root, P{axial_tumor_wilcoxon_qda_index}, 'top5_test2_feats.mat')).test2_feats;
                axial_fat_feat_names = load(fullfile(results_root, P{axial_fat_wilcoxon_qda_index}, 'Top_Features_Information.mat')).topResults.names';
                axial_tumor_feat_names = load(fullfile(results_root, P{axial_tumor_wilcoxon_qda_index}, 'Top_Features_Information.mat')).topResults.names';
                feature_column_names = [axial_fat_feat_names, axial_tumor_feat_names];
                features_training = [axial_fat_feats_train, axial_tumor_feats_train];
                features_holdout1 = [axial_fat_feats_test1, axial_tumor_feats_test1];
                features_holdout2 = [axial_fat_feats_test2, axial_tumor_feats_test2];
            end
        elseif strcmp(plane{1}, 'Coronal')
            if strcmp(best_fat_roi{1}, 'Proximal_Fat5')
                coronal_fat_feats_train = load(fullfile(results_root, P{coronal_proxfat5_wilcoxon_qda_index}, 'top5_training_feats.mat')).train_feats;
                coronal_fat_feats_test1 = load(fullfile(results_root, P{coronal_proxfat5_wilcoxon_qda_index}, 'top5_test1_feats.mat')).test1_feats;
                coronal_fat_feats_test2 = load(fullfile(results_root, P{coronal_proxfat5_wilcoxon_qda_index}, 'top5_test2_feats.mat')).test2_feats;
                coronal_tumor_feats_train = load(fullfile(results_root, P{coronal_tumor_wilcoxon_qda_index}, 'top5_training_feats.mat')).train_feats;
                coronal_tumor_feats_test1 = load(fullfile(results_root, P{coronal_tumor_wilcoxon_qda_index}, 'top5_test1_feats.mat')).test1_feats;
                coronal_tumor_feats_test2 = load(fullfile(results_root, P{coronal_tumor_wilcoxon_qda_index}, 'top5_test2_feats.mat')).test2_feats;
                coronal_fat_feat_names = load(fullfile(results_root, P{coronal_proxfat5_wilcoxon_qda_index}, 'Top_Features_Information.mat')).topResults.names';
                coronal_tumor_feat_names = load(fullfile(results_root, P{coronal_tumor_wilcoxon_qda_index}, 'Top_Features_Information.mat')).topResults.names';
                feature_column_names = [coronal_fat_feat_names, coronal_tumor_feat_names];
                features_training = [coronal_fat_feats_train, coronal_tumor_feats_train];
                features_holdout1 = [coronal_fat_feats_test1, coronal_tumor_feats_test1];
                features_holdout2 = [coronal_fat_feats_test2, coronal_tumor_feats_test2];
            elseif strcmp(best_fat_roi{1}, 'Proximal_Fat10')
                coronal_fat_feats_train = load(fullfile(results_root, P{coronal_proxfat10_wilcoxon_qda_index}, 'top5_training_feats.mat')).train_feats;
                coronal_fat_feats_test1 = load(fullfile(results_root, P{coronal_proxfat10_wilcoxon_qda_index}, 'top5_test1_feats.mat')).test1_feats;
                coronal_fat_feats_test2 = load(fullfile(results_root, P{coronal_proxfat10_wilcoxon_qda_index}, 'top5_test2_feats.mat')).test2_feats;
                coronal_tumor_feats_train = load(fullfile(results_root, P{coronal_tumor_wilcoxon_qda_index}, 'top5_training_feats.mat')).train_feats;
                coronal_tumor_feats_test1 = load(fullfile(results_root, P{coronal_tumor_wilcoxon_qda_index}, 'top5_test1_feats.mat')).test1_feats;
                coronal_tumor_feats_test2 = load(fullfile(results_root, P{coronal_tumor_wilcoxon_qda_index}, 'top5_test2_feats.mat')).test2_feats;
                coronal_fat_feat_names = load(fullfile(results_root, P{coronal_proxfat5_wilcoxon_qda_index}, 'Top_Features_Information.mat')).topResults.names';
                coronal_tumor_feat_names = load(fullfile(results_root, P{coronal_tumor_wilcoxon_qda_index}, 'Top_Features_Information.mat')).topResults.names';
                feature_column_names = [coronal_fat_feat_names, coronal_tumor_feat_names];
                features_training = [coronal_fat_feats_train, coronal_tumor_feats_train];
                features_holdout1 = [coronal_fat_feats_test1, coronal_tumor_feats_test1];
                features_holdout2 = [coronal_fat_feats_test2, coronal_tumor_feats_test2];
            elseif strcmp(best_fat_roi{1}, 'Fat')
                coronal_fat_feats_train = load(fullfile(results_root, P{coronal_fat_wilcoxon_qda_index}, 'top5_training_feats.mat')).train_feats;
                coronal_fat_feats_test1 = load(fullfile(results_root, P{coronal_fat_wilcoxon_qda_index}, 'top5_test1_feats.mat')).test1_feats;
                coronal_fat_feats_test2 = load(fullfile(results_root, P{coronal_fat_wilcoxon_qda_index}, 'top5_test2_feats.mat')).test2_feats;
                coronal_tumor_feats_train = load(fullfile(results_root, P{coronal_tumor_wilcoxon_qda_index}, 'top5_training_feats.mat')).train_feats;
                coronal_tumor_feats_test1 = load(fullfile(results_root, P{coronal_tumor_wilcoxon_qda_index}, 'top5_test1_feats.mat')).test1_feats;
                coronal_tumor_feats_test2 = load(fullfile(results_root, P{coronal_tumor_wilcoxon_qda_index}, 'top5_test2_feats.mat')).test2_feats;
                coronal_fat_feat_names = load(fullfile(results_root, P{coronal_fat_wilcoxon_qda_index}, 'Top_Features_Information.mat')).topResults.names';
                coronal_tumor_feat_names = load(fullfile(results_root, P{coronal_tumor_wilcoxon_qda_index}, 'Top_Features_Information.mat')).topResults.names';
                feature_column_names = [coronal_fat_feat_names, coronal_tumor_feat_names];
                features_training = [coronal_fat_feats_train, coronal_tumor_feats_train];
                features_holdout1 = [coronal_fat_feats_test1, coronal_tumor_feats_test1];
                features_holdout2 = [coronal_fat_feats_test2, coronal_tumor_feats_test2];
            end
        end
    elseif strcmp(classifier(1), 'lda')
        if strcmp(plane{1}, 'Multi-Plane')
            if length(rois) == 2
                if strcmp(best_fat_roi{1}, 'Proximal_Fat5')
                    axial_fat_feats_train = load(fullfile(results_root, P{axial_proxfat5_wilcoxon_lda_index}, 'top5_training_feats.mat')).train_feats;
                    axial_fat_feats_test1 = load(fullfile(results_root, P{axial_proxfat5_wilcoxon_lda_index}, 'top5_test1_feats.mat')).test1_feats;
                    axial_fat_feats_test2 = load(fullfile(results_root, P{axial_proxfat5_wilcoxon_lda_index}, 'top5_test2_feats.mat')).test2_feats;
                    axial_tumor_feats_train = load(fullfile(results_root, P{axial_tumor_wilcoxon_lda_index}, 'top5_training_feats.mat')).train_feats;
                    axial_tumor_feats_test1 = load(fullfile(results_root, P{axial_tumor_wilcoxon_lda_index}, 'top5_test1_feats.mat')).test1_feats;
                    axial_tumor_feats_test2 = load(fullfile(results_root, P{axial_tumor_wilcoxon_lda_index}, 'top5_test2_feats.mat')).test2_feats;
                    coronal_fat_feats_train = load(fullfile(results_root, P{coronal_proxfat5_wilcoxon_lda_index}, 'top5_training_feats.mat')).train_feats;
                    coronal_fat_feats_test1 = load(fullfile(results_root, P{coronal_proxfat5_wilcoxon_lda_index}, 'top5_test1_feats.mat')).test1_feats;
                    coronal_fat_feats_test2 = load(fullfile(results_root, P{coronal_proxfat5_wilcoxon_lda_index}, 'top5_test2_feats.mat')).test2_feats;
                    coronal_tumor_feats_train = load(fullfile(results_root, P{coronal_tumor_wilcoxon_lda_index}, 'top5_training_feats.mat')).train_feats;
                    coronal_tumor_feats_test1 = load(fullfile(results_root, P{coronal_tumor_wilcoxon_lda_index}, 'top5_test1_feats.mat')).test1_feats;
                    coronal_tumor_feats_test2 = load(fullfile(results_root, P{coronal_tumor_wilcoxon_lda_index}, 'top5_test2_feats.mat')).test2_feats;
                    axial_fat_feat_names = load(fullfile(results_root, P{axial_proxfat5_wilcoxon_lda_index}, 'Top_Features_Information.mat')).topResults.names';
                    coronal_fat_feat_names = load(fullfile(results_root, P{coronal_proxfat5_wilcoxon_lda_index}, 'Top_Features_Information.mat')).topResults.names';
                    axial_tumor_feat_names = load(fullfile(results_root, P{axial_tumor_wilcoxon_lda_index}, 'Top_Features_Information.mat')).topResults.names';
                    coronal_tumor_feat_names = load(fullfile(results_root, P{coronal_tumor_wilcoxon_lda_index}, 'Top_Features_Information.mat')).topResults.names';
                    feature_column_names = [axial_fat_feat_names, coronal_fat_feat_names, axial_tumor_feat_names, coronal_tumor_feat_names];
                    features_training = [axial_fat_feats_train, coronal_fat_feats_train, axial_tumor_feats_train, coronal_tumor_feats_train];
                    features_holdout1 = [axial_fat_feats_test1, coronal_fat_feats_test1, axial_tumor_feats_test1, coronal_tumor_feats_test1];
                    features_holdout2 = [axial_fat_feats_test2, coronal_fat_feats_test2, axial_tumor_feats_test2, coronal_tumor_feats_test2];
                elseif strcmp(best_fat_roi{1}, 'Proximal_Fat10')
                    axial_fat_feats_train = load(fullfile(results_root, P{axial_proxfat10_wilcoxon_lda_index}, 'top5_training_feats.mat')).train_feats;
                    axial_fat_feats_test1 = load(fullfile(results_root, P{axial_proxfat10_wilcoxon_lda_index}, 'top5_test1_feats.mat')).test1_feats;
                    axial_fat_feats_test2 = load(fullfile(results_root, P{axial_proxfat10_wilcoxon_lda_index}, 'top5_test2_feats.mat')).test2_feats;
                    axial_tumor_feats_train = load(fullfile(results_root, P{axial_tumor_wilcoxon_lda_index}, 'top5_training_feats.mat')).train_feats;
                    axial_tumor_feats_test1 = load(fullfile(results_root, P{axial_tumor_wilcoxon_lda_index}, 'top5_test1_feats.mat')).test1_feats;
                    axial_tumor_feats_test2 = load(fullfile(results_root, P{axial_tumor_wilcoxon_lda_index}, 'top5_test2_feats.mat')).test2_feats;
                    coronal_fat_feats_train = load(fullfile(results_root, P{coronal_proxfat10_wilcoxon_lda_index}, 'top5_training_feats.mat')).train_feats;
                    coronal_fat_feats_test1 = load(fullfile(results_root, P{coronal_proxfat10_wilcoxon_lda_index}, 'top5_test1_feats.mat')).test1_feats;
                    coronal_fat_feats_test2 = load(fullfile(results_root, P{coronal_proxfat10_wilcoxon_lda_index}, 'top5_test2_feats.mat')).test2_feats;
                    coronal_tumor_feats_train = load(fullfile(results_root, P{coronal_tumor_wilcoxon_lda_index}, 'top5_training_feats.mat')).train_feats;
                    coronal_tumor_feats_test1 = load(fullfile(results_root, P{coronal_tumor_wilcoxon_lda_index}, 'top5_test1_feats.mat')).test1_feats;
                    coronal_tumor_feats_test2 = load(fullfile(results_root, P{coronal_tumor_wilcoxon_lda_index}, 'top5_test2_feats.mat')).test2_feats;
                    axial_fat_feat_names = load(fullfile(results_root, P{axial_proxfat5_wilcoxon_lda_index}, 'Top_Features_Information.mat')).topResults.names';
                    coronal_fat_feat_names = load(fullfile(results_root, P{coronal_proxfat5_wilcoxon_lda_index}, 'Top_Features_Information.mat')).topResults.names';
                    axial_tumor_feat_names = load(fullfile(results_root, P{axial_tumor_wilcoxon_lda_index}, 'Top_Features_Information.mat')).topResults.names';
                    coronal_tumor_feat_names = load(fullfile(results_root, P{coronal_tumor_wilcoxon_lda_index}, 'Top_Features_Information.mat')).topResults.names';
                    feature_column_names = [axial_fat_feat_names, coronal_fat_feat_names, axial_tumor_feat_names, coronal_tumor_feat_names];
                    features_training = [axial_fat_feats_train, coronal_fat_feats_train, axial_tumor_feats_train, coronal_tumor_feats_train];
                    features_holdout1 = [axial_fat_feats_test1, coronal_fat_feats_test1, axial_tumor_feats_test1, coronal_tumor_feats_test1];
                    features_holdout2 = [axial_fat_feats_test2, coronal_fat_feats_test2, axial_tumor_feats_test2, coronal_tumor_feats_test2];
                elseif strcmp(best_fat_roi{1}, 'Fat')
                    axial_fat_feats_train = load(fullfile(results_root, P{axial_fat_wilcoxon_lda_index}, 'top5_training_feats.mat')).train_feats;
                    axial_fat_feats_test1 = load(fullfile(results_root, P{axial_fat_wilcoxon_lda_index}, 'top5_test1_feats.mat')).test1_feats;
                    axial_fat_feats_test2 = load(fullfile(results_root, P{axial_fat_wilcoxon_lda_index}, 'top5_test2_feats.mat')).test2_feats;
                    axial_tumor_feats_train = load(fullfile(results_root, P{axial_tumor_wilcoxon_lda_index}, 'top5_training_feats.mat')).train_feats;
                    axial_tumor_feats_test1 = load(fullfile(results_root, P{axial_tumor_wilcoxon_lda_index}, 'top5_test1_feats.mat')).test1_feats;
                    axial_tumor_feats_test2 = load(fullfile(results_root, P{axial_tumor_wilcoxon_lda_index}, 'top5_test2_feats.mat')).test2_feats;
                    coronal_fat_feats_train = load(fullfile(results_root, P{coronal_fat_wilcoxon_lda_index}, 'top5_training_feats.mat')).train_feats;
                    coronal_fat_feats_test1 = load(fullfile(results_root, P{coronal_fat_wilcoxon_lda_index}, 'top5_test1_feats.mat')).test1_feats;
                    coronal_fat_feats_test2 = load(fullfile(results_root, P{coronal_fat_wilcoxon_lda_index}, 'top5_test2_feats.mat')).test2_feats;
                    coronal_tumor_feats_train = load(fullfile(results_root, P{coronal_tumor_wilcoxon_lda_index}, 'top5_training_feats.mat')).train_feats;
                    coronal_tumor_feats_test1 = load(fullfile(results_root, P{coronal_tumor_wilcoxon_lda_index}, 'top5_test1_feats.mat')).test1_feats;
                    coronal_tumor_feats_test2 = load(fullfile(results_root, P{coronal_tumor_wilcoxon_lda_index}, 'top5_test2_feats.mat')).test2_feats;
                    axial_fat_feat_names = load(fullfile(results_root, P{axial_fat_wilcoxon_lda_index}, 'Top_Features_Information.mat')).topResults.names';
                    coronal_fat_feat_names = load(fullfile(results_root, P{coronal_fat_wilcoxon_lda_index}, 'Top_Features_Information.mat')).topResults.names';
                    axial_tumor_feat_names = load(fullfile(results_root, P{axial_tumor_wilcoxon_lda_index}, 'Top_Features_Information.mat')).topResults.names';
                    coronal_tumor_feat_names = load(fullfile(results_root, P{coronal_tumor_wilcoxon_lda_index}, 'Top_Features_Information.mat')).topResults.names';
                    feature_column_names = [axial_fat_feat_names, coronal_fat_feat_names, axial_tumor_feat_names, coronal_tumor_feat_names];
                    features_training = [axial_fat_feats_train, coronal_fat_feats_train, axial_tumor_feats_train, coronal_tumor_feats_train];
                    features_holdout1 = [axial_fat_feats_test1, coronal_fat_feats_test1, axial_tumor_feats_test1, coronal_tumor_feats_test1];
                    features_holdout2 = [axial_fat_feats_test2, coronal_fat_feats_test2, axial_tumor_feats_test2, coronal_tumor_feats_test2];
                end
            elseif length(rois) == 1
                if strcmp(rois{1}, 'Proximal_Fat5')
                    axial_feats_train = load(fullfile(results_root, P{axial_proxfat5_wilcoxon_lda_index}, 'top5_training_feats.mat')).train_feats;
                    axial_feats_test1 = load(fullfile(results_root, P{axial_proxfat5_wilcoxon_lda_index}, 'top5_test1_feats.mat')).test1_feats;
                    axial_feats_test2 = load(fullfile(results_root, P{axial_proxfat5_wilcoxon_lda_index}, 'top5_test2_feats.mat')).test2_feats;
                    coronal_feats_train = load(fullfile(results_root, P{coronal_proxfat5_wilcoxon_lda_index}, 'top5_training_feats.mat')).train_feats;
                    coronal_feats_test1 = load(fullfile(results_root, P{coronal_proxfat5_wilcoxon_lda_index}, 'top5_test1_feats.mat')).test1_feats;
                    coronal_feats_test2 = load(fullfile(results_root, P{coronal_proxfat5_wilcoxon_lda_index}, 'top5_test2_feats.mat')).test2_feats;
                    axial_feat_names = load(fullfile(results_root, P{axial_proxfat5_wilcoxon_lda_index}, 'Top_Features_Information.mat')).topResults.names';
                    coronal_feat_names = load(fullfile(results_root, P{coronal_proxfat5_wilcoxon_lda_index}, 'Top_Features_Information.mat')).topResults.names';
                    feature_column_names = [axial_feat_names, coronal_feat_names];
                    features_training = [axial_feats_train, coronal_feats_train];
                    features_holdout1 = [axial_feats_test1, coronal_feats_test1];
                    features_holdout2 = [axial_feats_test2, coronal_feats_test2];
                elseif strcmp(rois{1}, 'Proximal_Fat10')
                    axial_feats_train = load(fullfile(results_root, P{axial_proxfat10_wilcoxon_lda_index}, 'top5_training_feats.mat')).train_feats;
                    axial_feats_test1 = load(fullfile(results_root, P{axial_proxfat10_wilcoxon_lda_index}, 'top5_test1_feats.mat')).test1_feats;
                    axial_feats_test2 = load(fullfile(results_root, P{axial_proxfat10_wilcoxon_lda_index}, 'top5_test2_feats.mat')).test2_feats;
                    coronal_feats_train = load(fullfile(results_root, P{coronal_proxfat10_wilcoxon_lda_index}, 'top5_training_feats.mat')).train_feats;
                    coronal_feats_test1 = load(fullfile(results_root, P{coronal_proxfat10_wilcoxon_lda_index}, 'top5_test1_feats.mat')).test1_feats;
                    coronal_feats_test2 = load(fullfile(results_root, P{coronal_proxfat10_wilcoxon_lda_index}, 'top5_test2_feats.mat')).test2_feats;
                    axial_feat_names = load(fullfile(results_root, P{axial_proxfat10_wilcoxon_lda_index}, 'Top_Features_Information.mat')).topResults.names';
                    coronal_feat_names = load(fullfile(results_root, P{coronal_proxfat10_wilcoxon_lda_index}, 'Top_Features_Information.mat')).topResults.names';
                    feature_column_names = [axial_feat_names, coronal_feat_names];
                    features_training = [axial_feats_train, coronal_feats_train];
                    features_holdout1 = [axial_feats_test1, coronal_feats_test1];
                    features_holdout2 = [axial_feats_test2, coronal_feats_test2];
                elseif strcmp(rois{1}, 'Fat')
                    axial_feats_train = load(fullfile(results_root, P{axial_fat_wilcoxon_lda_index}, 'top5_training_feats.mat')).train_feats;
                    axial_feats_test1 = load(fullfile(results_root, P{axial_fat_wilcoxon_lda_index}, 'top5_test1_feats.mat')).test1_feats;
                    axial_feats_test2 = load(fullfile(results_root, P{axial_fat_wilcoxon_lda_index}, 'top5_test2_feats.mat')).test2_feats;
                    coronal_feats_train = load(fullfile(results_root, P{coronal_fat_wilcoxon_lda_index}, 'top5_training_feats.mat')).train_feats;
                    coronal_feats_test1 = load(fullfile(results_root, P{coronal_fat_wilcoxon_lda_index}, 'top5_test1_feats.mat')).test1_feats;
                    coronal_feats_test2 = load(fullfile(results_root, P{coronal_fat_wilcoxon_lda_index}, 'top5_test2_feats.mat')).test2_feats;
                    axial_feat_names = load(fullfile(results_root, P{axial_fat_wilcoxon_lda_index}, 'Top_Features_Information.mat')).topResults.names';
                    coronal_feat_names = load(fullfile(results_root, P{coronal_fat_wilcoxon_lda_index}, 'Top_Features_Information.mat')).topResults.names';
                    feature_column_names = [axial_feat_names, coronal_feat_names];
                    features_training = [axial_feats_train, coronal_feats_train];
                    features_holdout1 = [axial_feats_test1, coronal_feats_test1];
                    features_holdout2 = [axial_feats_test2, coronal_feats_test2];
                elseif strcmp(rois{1}, 'Tumor')
                    axial_feats_train = load(fullfile(results_root, P{axial_tumor_wilcoxon_lda_index}, 'top5_training_feats.mat')).train_feats;
                    axial_feats_test1 = load(fullfile(results_root, P{axial_tumor_wilcoxon_lda_index}, 'top5_test1_feats.mat')).test1_feats;
                    axial_feats_test2 = load(fullfile(results_root, P{axial_tumor_wilcoxon_lda_index}, 'top5_test2_feats.mat')).test2_feats;
                    coronal_feats_train = load(fullfile(results_root, P{coronal_tumor_wilcoxon_lda_index}, 'top5_training_feats.mat')).train_feats;
                    coronal_feats_test1 = load(fullfile(results_root, P{coronal_tumor_wilcoxon_lda_index}, 'top5_test1_feats.mat')).test1_feats;
                    coronal_feats_test2 = load(fullfile(results_root, P{coronal_tumor_wilcoxon_lda_index}, 'top5_test2_feats.mat')).test2_feats;
                    axial_feat_names = load(fullfile(results_root, P{axial_tumor_wilcoxon_lda_index}, 'Top_Features_Information.mat')).topResults.names';
                    coronal_feat_names = load(fullfile(results_root, P{coronal_tumor_wilcoxon_lda_index}, 'Top_Features_Information.mat')).topResults.names';
                    feature_column_names = [axial_feat_names, coronal_feat_names];
                    features_training = [axial_feats_train, coronal_feats_train];
                    features_holdout1 = [axial_feats_test1, coronal_feats_test1];
                    features_holdout2 = [axial_feats_test2, coronal_feats_test2];
                end
            end
        elseif strcmp(plane{1}, 'Axial')
            if strcmp(best_fat_roi{1}, 'Proximal_Fat5')
                axial_fat_feats_train = load(fullfile(results_root, P{axial_proxfat5_wilcoxon_lda_index}, 'top5_training_feats.mat')).train_feats;
                axial_fat_feats_test1 = load(fullfile(results_root, P{axial_proxfat5_wilcoxon_lda_index}, 'top5_test1_feats.mat')).test1_feats;
                axial_fat_feats_test2 = load(fullfile(results_root, P{axial_proxfat5_wilcoxon_lda_index}, 'top5_test2_feats.mat')).test2_feats;
                axial_tumor_feats_train = load(fullfile(results_root, P{axial_tumor_wilcoxon_lda_index}, 'top5_training_feats.mat')).train_feats;
                axial_tumor_feats_test1 = load(fullfile(results_root, P{axial_tumor_wilcoxon_lda_index}, 'top5_test1_feats.mat')).test1_feats;
                axial_tumor_feats_test2 = load(fullfile(results_root, P{axial_tumor_wilcoxon_lda_index}, 'top5_test2_feats.mat')).test2_feats;
                axial_fat_feat_names = load(fullfile(results_root, P{axial_proxfat5_wilcoxon_lda_index}, 'Top_Features_Information.mat')).topResults.names';
                axial_tumor_feat_names = load(fullfile(results_root, P{axial_tumor_wilcoxon_lda_index}, 'Top_Features_Information.mat')).topResults.names';
                feature_column_names = [axial_fat_feat_names, axial_tumor_feat_names];
                features_training = [axial_fat_feats_train, axial_tumor_feats_train];
                features_holdout1 = [axial_fat_feats_test1, axial_tumor_feats_test1];
                features_holdout2 = [axial_fat_feats_test2, axial_tumor_feats_test2];
            elseif strcmp(best_fat_roi{1}, 'Proximal_Fat10')
                axial_fat_feats_train = load(fullfile(results_root, P{axial_proxfat10_wilcoxon_lda_index}, 'top5_training_feats.mat')).train_feats;
                axial_fat_feats_test1 = load(fullfile(results_root, P{axial_proxfat10_wilcoxon_lda_index}, 'top5_test1_feats.mat')).test1_feats;
                axial_fat_feats_test2 = load(fullfile(results_root, P{axial_proxfat10_wilcoxon_lda_index}, 'top5_test2_feats.mat')).test2_feats;
                axial_tumor_feats_train = load(fullfile(results_root, P{axial_tumor_wilcoxon_lda_index}, 'top5_training_feats.mat')).train_feats;
                axial_tumor_feats_test1 = load(fullfile(results_root, P{axial_tumor_wilcoxon_lda_index}, 'top5_test1_feats.mat')).test1_feats;
                axial_tumor_feats_test2 = load(fullfile(results_root, P{axial_tumor_wilcoxon_lda_index}, 'top5_test2_feats.mat')).test2_feats;
                axial_fat_feat_names = load(fullfile(results_root, P{axial_proxfat5_wilcoxon_lda_index}, 'Top_Features_Information.mat')).topResults.names';
                axial_tumor_feat_names = load(fullfile(results_root, P{axial_tumor_wilcoxon_lda_index}, 'Top_Features_Information.mat')).topResults.names';
                feature_column_names = [axial_fat_feat_names, axial_tumor_feat_names];
                features_training = [axial_fat_feats_train, axial_tumor_feats_train];
                features_holdout1 = [axial_fat_feats_test1, axial_tumor_feats_test1];
                features_holdout2 = [axial_fat_feats_test2, axial_tumor_feats_test2];
            elseif strcmp(best_fat_roi{1}, 'Fat')
                axial_fat_feats_train = load(fullfile(results_root, P{axial_fat_wilcoxon_lda_index}, 'top5_training_feats.mat')).train_feats;
                axial_fat_feats_test1 = load(fullfile(results_root, P{axial_fat_wilcoxon_lda_index}, 'top5_test1_feats.mat')).test1_feats;
                axial_fat_feats_test2 = load(fullfile(results_root, P{axial_fat_wilcoxon_lda_index}, 'top5_test2_feats.mat')).test2_feats;
                axial_tumor_feats_train = load(fullfile(results_root, P{axial_tumor_wilcoxon_lda_index}, 'top5_training_feats.mat')).train_feats;
                axial_tumor_feats_test1 = load(fullfile(results_root, P{axial_tumor_wilcoxon_lda_index}, 'top5_test1_feats.mat')).test1_feats;
                axial_tumor_feats_test2 = load(fullfile(results_root, P{axial_tumor_wilcoxon_lda_index}, 'top5_test2_feats.mat')).test2_feats;
                axial_fat_feat_names = load(fullfile(results_root, P{axial_fat_wilcoxon_lda_index}, 'Top_Features_Information.mat')).topResults.names';
                axial_tumor_feat_names = load(fullfile(results_root, P{axial_tumor_wilcoxon_lda_index}, 'Top_Features_Information.mat')).topResults.names';
                feature_column_names = [axial_fat_feat_names, axial_tumor_feat_names];
                features_training = [axial_fat_feats_train, axial_tumor_feats_train];
                features_holdout1 = [axial_fat_feats_test1, axial_tumor_feats_test1];
                features_holdout2 = [axial_fat_feats_test2, axial_tumor_feats_test2];
            end
        elseif strcmp(plane{1}, 'Coronal')
            if strcmp(best_fat_roi{1}, 'Proximal_Fat5')
                coronal_fat_feats_train = load(fullfile(results_root, P{coronal_proxfat5_wilcoxon_lda_index}, 'top5_training_feats.mat')).train_feats;
                coronal_fat_feats_test1 = load(fullfile(results_root, P{coronal_proxfat5_wilcoxon_lda_index}, 'top5_test1_feats.mat')).test1_feats;
                coronal_fat_feats_test2 = load(fullfile(results_root, P{coronal_proxfat5_wilcoxon_lda_index}, 'top5_test2_feats.mat')).test2_feats;
                coronal_tumor_feats_train = load(fullfile(results_root, P{coronal_tumor_wilcoxon_lda_index}, 'top5_training_feats.mat')).train_feats;
                coronal_tumor_feats_test1 = load(fullfile(results_root, P{coronal_tumor_wilcoxon_lda_index}, 'top5_test1_feats.mat')).test1_feats;
                coronal_tumor_feats_test2 = load(fullfile(results_root, P{coronal_tumor_wilcoxon_lda_index}, 'top5_test2_feats.mat')).test2_feats;
                coronal_fat_feat_names = load(fullfile(results_root, P{coronal_proxfat5_wilcoxon_lda_index}, 'Top_Features_Information.mat')).topResults.names';
                coronal_tumor_feat_names = load(fullfile(results_root, P{coronal_tumor_wilcoxon_lda_index}, 'Top_Features_Information.mat')).topResults.names';
                feature_column_names = [coronal_fat_feat_names, coronal_tumor_feat_names];
                features_training = [coronal_fat_feats_train, coronal_tumor_feats_train];
                features_holdout1 = [coronal_fat_feats_test1, coronal_tumor_feats_test1];
                features_holdout2 = [coronal_fat_feats_test2, coronal_tumor_feats_test2];
            elseif strcmp(best_fat_roi{1}, 'Proximal_Fat10')
                coronal_fat_feats_train = load(fullfile(results_root, P{coronal_proxfat10_wilcoxon_lda_index}, 'top5_training_feats.mat')).train_feats;
                coronal_fat_feats_test1 = load(fullfile(results_root, P{coronal_proxfat10_wilcoxon_lda_index}, 'top5_test1_feats.mat')).test1_feats;
                coronal_fat_feats_test2 = load(fullfile(results_root, P{coronal_proxfat10_wilcoxon_lda_index}, 'top5_test2_feats.mat')).test2_feats;
                coronal_tumor_feats_train = load(fullfile(results_root, P{coronal_tumor_wilcoxon_lda_index}, 'top5_training_feats.mat')).train_feats;
                coronal_tumor_feats_test1 = load(fullfile(results_root, P{coronal_tumor_wilcoxon_lda_index}, 'top5_test1_feats.mat')).test1_feats;
                coronal_tumor_feats_test2 = load(fullfile(results_root, P{coronal_tumor_wilcoxon_lda_index}, 'top5_test2_feats.mat')).test2_feats;
                coronal_fat_feat_names = load(fullfile(results_root, P{coronal_proxfat5_wilcoxon_lda_index}, 'Top_Features_Information.mat')).topResults.names';
                coronal_tumor_feat_names = load(fullfile(results_root, P{coronal_tumor_wilcoxon_lda_index}, 'Top_Features_Information.mat')).topResults.names';
                feature_column_names = [coronal_fat_feat_names, coronal_tumor_feat_names];
                features_training = [coronal_fat_feats_train, coronal_tumor_feats_train];
                features_holdout1 = [coronal_fat_feats_test1, coronal_tumor_feats_test1];
                features_holdout2 = [coronal_fat_feats_test2, coronal_tumor_feats_test2];
            elseif strcmp(best_fat_roi{1}, 'Fat')
                coronal_fat_feats_train = load(fullfile(results_root, P{coronal_fat_wilcoxon_lda_index}, 'top5_training_feats.mat')).train_feats;
                coronal_fat_feats_test1 = load(fullfile(results_root, P{coronal_fat_wilcoxon_lda_index}, 'top5_test1_feats.mat')).test1_feats;
                coronal_fat_feats_test2 = load(fullfile(results_root, P{coronal_fat_wilcoxon_lda_index}, 'top5_test2_feats.mat')).test2_feats;
                coronal_tumor_feats_train = load(fullfile(results_root, P{coronal_tumor_wilcoxon_lda_index}, 'top5_training_feats.mat')).train_feats;
                coronal_tumor_feats_test1 = load(fullfile(results_root, P{coronal_tumor_wilcoxon_lda_index}, 'top5_test1_feats.mat')).test1_feats;
                coronal_tumor_feats_test2 = load(fullfile(results_root, P{coronal_tumor_wilcoxon_lda_index}, 'top5_test2_feats.mat')).test2_feats;
                coronal_fat_feat_names = load(fullfile(results_root, P{coronal_fat_wilcoxon_lda_index}, 'Top_Features_Information.mat')).topResults.names';
                coronal_tumor_feat_names = load(fullfile(results_root, P{coronal_tumor_wilcoxon_lda_index}, 'Top_Features_Information.mat')).topResults.names';
                feature_column_names = [coronal_fat_feat_names, coronal_tumor_feat_names];
                features_training = [coronal_fat_feats_train, coronal_tumor_feats_train];
                features_holdout1 = [coronal_fat_feats_test1, coronal_tumor_feats_test1];
                features_holdout2 = [coronal_fat_feats_test2, coronal_tumor_feats_test2];
            end
        end
    elseif strcmp(classifier(1), 'rf')
        if strcmp(plane{1}, 'Multi-Plane')
            if length(rois) == 2
                if strcmp(best_fat_roi{1}, 'Proximal_Fat5')
                    axial_fat_feats_train = load(fullfile(results_root, P{axial_proxfat5_wilcoxon_rf_index}, 'top5_training_feats.mat')).train_feats;
                    axial_fat_feats_test1 = load(fullfile(results_root, P{axial_proxfat5_wilcoxon_rf_index}, 'top5_test1_feats.mat')).test1_feats;
                    axial_fat_feats_test2 = load(fullfile(results_root, P{axial_proxfat5_wilcoxon_rf_index}, 'top5_test2_feats.mat')).test2_feats;
                    axial_tumor_feats_train = load(fullfile(results_root, P{axial_tumor_wilcoxon_rf_index}, 'top5_training_feats.mat')).train_feats;
                    axial_tumor_feats_test1 = load(fullfile(results_root, P{axial_tumor_wilcoxon_rf_index}, 'top5_test1_feats.mat')).test1_feats;
                    axial_tumor_feats_test2 = load(fullfile(results_root, P{axial_tumor_wilcoxon_rf_index}, 'top5_test2_feats.mat')).test2_feats;
                    coronal_fat_feats_train = load(fullfile(results_root, P{coronal_proxfat5_wilcoxon_rf_index}, 'top5_training_feats.mat')).train_feats;
                    coronal_fat_feats_test1 = load(fullfile(results_root, P{coronal_proxfat5_wilcoxon_rf_index}, 'top5_test1_feats.mat')).test1_feats;
                    coronal_fat_feats_test2 = load(fullfile(results_root, P{coronal_proxfat5_wilcoxon_rf_index}, 'top5_test2_feats.mat')).test2_feats;
                    coronal_tumor_feats_train = load(fullfile(results_root, P{coronal_tumor_wilcoxon_rf_index}, 'top5_training_feats.mat')).train_feats;
                    coronal_tumor_feats_test1 = load(fullfile(results_root, P{coronal_tumor_wilcoxon_rf_index}, 'top5_test1_feats.mat')).test1_feats;
                    coronal_tumor_feats_test2 = load(fullfile(results_root, P{coronal_tumor_wilcoxon_rf_index}, 'top5_test2_feats.mat')).test2_feats;
                    axial_fat_feat_names = load(fullfile(results_root, P{axial_proxfat5_wilcoxon_rf_index}, 'Top_Features_Information.mat')).topResults.names';
                    coronal_fat_feat_names = load(fullfile(results_root, P{coronal_proxfat5_wilcoxon_rf_index}, 'Top_Features_Information.mat')).topResults.names';
                    axial_tumor_feat_names = load(fullfile(results_root, P{axial_tumor_wilcoxon_rf_index}, 'Top_Features_Information.mat')).topResults.names';
                    coronal_tumor_feat_names = load(fullfile(results_root, P{coronal_tumor_wilcoxon_rf_index}, 'Top_Features_Information.mat')).topResults.names';
                    feature_column_names = [axial_fat_feat_names, coronal_fat_feat_names, axial_tumor_feat_names, coronal_tumor_feat_names];
                    features_training = [axial_fat_feats_train, coronal_fat_feats_train, axial_tumor_feats_train, coronal_tumor_feats_train];
                    features_holdout1 = [axial_fat_feats_test1, coronal_fat_feats_test1, axial_tumor_feats_test1, coronal_tumor_feats_test1];
                    features_holdout2 = [axial_fat_feats_test2, coronal_fat_feats_test2, axial_tumor_feats_test2, coronal_tumor_feats_test2];
                elseif strcmp(best_fat_roi{1}, 'Proximal_Fat10')
                    axial_fat_feats_train = load(fullfile(results_root, P{axial_proxfat10_wilcoxon_rf_index}, 'top5_training_feats.mat')).train_feats;
                    axial_fat_feats_test1 = load(fullfile(results_root, P{axial_proxfat10_wilcoxon_rf_index}, 'top5_test1_feats.mat')).test1_feats;
                    axial_fat_feats_test2 = load(fullfile(results_root, P{axial_proxfat10_wilcoxon_rf_index}, 'top5_test2_feats.mat')).test2_feats;
                    axial_tumor_feats_train = load(fullfile(results_root, P{axial_tumor_wilcoxon_rf_index}, 'top5_training_feats.mat')).train_feats;
                    axial_tumor_feats_test1 = load(fullfile(results_root, P{axial_tumor_wilcoxon_rf_index}, 'top5_test1_feats.mat')).test1_feats;
                    axial_tumor_feats_test2 = load(fullfile(results_root, P{axial_tumor_wilcoxon_rf_index}, 'top5_test2_feats.mat')).test2_feats;
                    coronal_fat_feats_train = load(fullfile(results_root, P{coronal_proxfat10_wilcoxon_rf_index}, 'top5_training_feats.mat')).train_feats;
                    coronal_fat_feats_test1 = load(fullfile(results_root, P{coronal_proxfat10_wilcoxon_rf_index}, 'top5_test1_feats.mat')).test1_feats;
                    coronal_fat_feats_test2 = load(fullfile(results_root, P{coronal_proxfat10_wilcoxon_rf_index}, 'top5_test2_feats.mat')).test2_feats;
                    coronal_tumor_feats_train = load(fullfile(results_root, P{coronal_tumor_wilcoxon_rf_index}, 'top5_training_feats.mat')).train_feats;
                    coronal_tumor_feats_test1 = load(fullfile(results_root, P{coronal_tumor_wilcoxon_rf_index}, 'top5_test1_feats.mat')).test1_feats;
                    coronal_tumor_feats_test2 = load(fullfile(results_root, P{coronal_tumor_wilcoxon_rf_index}, 'top5_test2_feats.mat')).test2_feats;
                    axial_fat_feat_names = load(fullfile(results_root, P{axial_proxfat5_wilcoxon_rf_index}, 'Top_Features_Information.mat')).topResults.names';
                    coronal_fat_feat_names = load(fullfile(results_root, P{coronal_proxfat5_wilcoxon_rf_index}, 'Top_Features_Information.mat')).topResults.names';
                    axial_tumor_feat_names = load(fullfile(results_root, P{axial_tumor_wilcoxon_rf_index}, 'Top_Features_Information.mat')).topResults.names';
                    coronal_tumor_feat_names = load(fullfile(results_root, P{coronal_tumor_wilcoxon_rf_index}, 'Top_Features_Information.mat')).topResults.names';
                    feature_column_names = [axial_fat_feat_names, coronal_fat_feat_names, axial_tumor_feat_names, coronal_tumor_feat_names];
                    features_training = [axial_fat_feats_train, coronal_fat_feats_train, axial_tumor_feats_train, coronal_tumor_feats_train];
                    features_holdout1 = [axial_fat_feats_test1, coronal_fat_feats_test1, axial_tumor_feats_test1, coronal_tumor_feats_test1];
                    features_holdout2 = [axial_fat_feats_test2, coronal_fat_feats_test2, axial_tumor_feats_test2, coronal_tumor_feats_test2];
                elseif strcmp(best_fat_roi{1}, 'Fat')
                    axial_fat_feats_train = load(fullfile(results_root, P{axial_fat_wilcoxon_rf_index}, 'top5_training_feats.mat')).train_feats;
                    axial_fat_feats_test1 = load(fullfile(results_root, P{axial_fat_wilcoxon_rf_index}, 'top5_test1_feats.mat')).test1_feats;
                    axial_fat_feats_test2 = load(fullfile(results_root, P{axial_fat_wilcoxon_rf_index}, 'top5_test2_feats.mat')).test2_feats;
                    axial_tumor_feats_train = load(fullfile(results_root, P{axial_tumor_wilcoxon_rf_index}, 'top5_training_feats.mat')).train_feats;
                    axial_tumor_feats_test1 = load(fullfile(results_root, P{axial_tumor_wilcoxon_rf_index}, 'top5_test1_feats.mat')).test1_feats;
                    axial_tumor_feats_test2 = load(fullfile(results_root, P{axial_tumor_wilcoxon_rf_index}, 'top5_test2_feats.mat')).test2_feats;
                    coronal_fat_feats_train = load(fullfile(results_root, P{coronal_fat_wilcoxon_rf_index}, 'top5_training_feats.mat')).train_feats;
                    coronal_fat_feats_test1 = load(fullfile(results_root, P{coronal_fat_wilcoxon_rf_index}, 'top5_test1_feats.mat')).test1_feats;
                    coronal_fat_feats_test2 = load(fullfile(results_root, P{coronal_fat_wilcoxon_rf_index}, 'top5_test2_feats.mat')).test2_feats;
                    coronal_tumor_feats_train = load(fullfile(results_root, P{coronal_tumor_wilcoxon_rf_index}, 'top5_training_feats.mat')).train_feats;
                    coronal_tumor_feats_test1 = load(fullfile(results_root, P{coronal_tumor_wilcoxon_rf_index}, 'top5_test1_feats.mat')).test1_feats;
                    coronal_tumor_feats_test2 = load(fullfile(results_root, P{coronal_tumor_wilcoxon_rf_index}, 'top5_test2_feats.mat')).test2_feats;
                    axial_fat_feat_names = load(fullfile(results_root, P{axial_fat_wilcoxon_rf_index}, 'Top_Features_Information.mat')).topResults.names';
                    coronal_fat_feat_names = load(fullfile(results_root, P{coronal_fat_wilcoxon_rf_index}, 'Top_Features_Information.mat')).topResults.names';
                    axial_tumor_feat_names = load(fullfile(results_root, P{axial_tumor_wilcoxon_rf_index}, 'Top_Features_Information.mat')).topResults.names';
                    coronal_tumor_feat_names = load(fullfile(results_root, P{coronal_tumor_wilcoxon_rf_index}, 'Top_Features_Information.mat')).topResults.names';
                    feature_column_names = [axial_fat_feat_names, coronal_fat_feat_names, axial_tumor_feat_names, coronal_tumor_feat_names];
                    features_training = [axial_fat_feats_train, coronal_fat_feats_train, axial_tumor_feats_train, coronal_tumor_feats_train];
                    features_holdout1 = [axial_fat_feats_test1, coronal_fat_feats_test1, axial_tumor_feats_test1, coronal_tumor_feats_test1];
                    features_holdout2 = [axial_fat_feats_test2, coronal_fat_feats_test2, axial_tumor_feats_test2, coronal_tumor_feats_test2];
                end
            elseif length(rois) == 1
                if strcmp(rois{1}, 'Proximal_Fat5')
                    axial_feats_train = load(fullfile(results_root, P{axial_proxfat5_wilcoxon_rf_index}, 'top5_training_feats.mat')).train_feats;
                    axial_feats_test1 = load(fullfile(results_root, P{axial_proxfat5_wilcoxon_rf_index}, 'top5_test1_feats.mat')).test1_feats;
                    axial_feats_test2 = load(fullfile(results_root, P{axial_proxfat5_wilcoxon_rf_index}, 'top5_test2_feats.mat')).test2_feats;
                    coronal_feats_train = load(fullfile(results_root, P{coronal_proxfat5_wilcoxon_rf_index}, 'top5_training_feats.mat')).train_feats;
                    coronal_feats_test1 = load(fullfile(results_root, P{coronal_proxfat5_wilcoxon_rf_index}, 'top5_test1_feats.mat')).test1_feats;
                    coronal_feats_test2 = load(fullfile(results_root, P{coronal_proxfat5_wilcoxon_rf_index}, 'top5_test2_feats.mat')).test2_feats;
                    axial_feat_names = load(fullfile(results_root, P{axial_proxfat5_wilcoxon_rf_index}, 'Top_Features_Information.mat')).topResults.names';
                    coronal_feat_names = load(fullfile(results_root, P{coronal_proxfat5_wilcoxon_rf_index}, 'Top_Features_Information.mat')).topResults.names';
                    feature_column_names = [axial_feat_names, coronal_feat_names];
                    features_training = [axial_feats_train, coronal_feats_train];
                    features_holdout1 = [axial_feats_test1, coronal_feats_test1];
                    features_holdout2 = [axial_feats_test2, coronal_feats_test2];
                elseif strcmp(rois{1}, 'Proximal_Fat10')
                    axial_feats_train = load(fullfile(results_root, P{axial_proxfat10_wilcoxon_rf_index}, 'top5_training_feats.mat')).train_feats;
                    axial_feats_test1 = load(fullfile(results_root, P{axial_proxfat10_wilcoxon_rf_index}, 'top5_test1_feats.mat')).test1_feats;
                    axial_feats_test2 = load(fullfile(results_root, P{axial_proxfat10_wilcoxon_rf_index}, 'top5_test2_feats.mat')).test2_feats;
                    coronal_feats_train = load(fullfile(results_root, P{coronal_proxfat10_wilcoxon_rf_index}, 'top5_training_feats.mat')).train_feats;
                    coronal_feats_test1 = load(fullfile(results_root, P{coronal_proxfat10_wilcoxon_rf_index}, 'top5_test1_feats.mat')).test1_feats;
                    coronal_feats_test2 = load(fullfile(results_root, P{coronal_proxfat10_wilcoxon_rf_index}, 'top5_test2_feats.mat')).test2_feats;
                    axial_feat_names = load(fullfile(results_root, P{axial_proxfat10_wilcoxon_rf_index}, 'Top_Features_Information.mat')).topResults.names';
                    coronal_feat_names = load(fullfile(results_root, P{coronal_proxfat10_wilcoxon_rf_index}, 'Top_Features_Information.mat')).topResults.names';
                    feature_column_names = [axial_feat_names, coronal_feat_names];
                    features_training = [axial_feats_train, coronal_feats_train];
                    features_holdout1 = [axial_feats_test1, coronal_feats_test1];
                    features_holdout2 = [axial_feats_test2, coronal_feats_test2];
                elseif strcmp(rois{1}, 'Fat')
                    axial_feats_train = load(fullfile(results_root, P{axial_fat_wilcoxon_rf_index}, 'top5_training_feats.mat')).train_feats;
                    axial_feats_test1 = load(fullfile(results_root, P{axial_fat_wilcoxon_rf_index}, 'top5_test1_feats.mat')).test1_feats;
                    axial_feats_test2 = load(fullfile(results_root, P{axial_fat_wilcoxon_rf_index}, 'top5_test2_feats.mat')).test2_feats;
                    coronal_feats_train = load(fullfile(results_root, P{coronal_fat_wilcoxon_rf_index}, 'top5_training_feats.mat')).train_feats;
                    coronal_feats_test1 = load(fullfile(results_root, P{coronal_fat_wilcoxon_rf_index}, 'top5_test1_feats.mat')).test1_feats;
                    coronal_feats_test2 = load(fullfile(results_root, P{coronal_fat_wilcoxon_rf_index}, 'top5_test2_feats.mat')).test2_feats;
                    axial_feat_names = load(fullfile(results_root, P{axial_fat_wilcoxon_rf_index}, 'Top_Features_Information.mat')).topResults.names';
                    coronal_feat_names = load(fullfile(results_root, P{coronal_fat_wilcoxon_rf_index}, 'Top_Features_Information.mat')).topResults.names';
                    feature_column_names = [axial_feat_names, coronal_feat_names];
                    features_training = [axial_feats_train, coronal_feats_train];
                    features_holdout1 = [axial_feats_test1, coronal_feats_test1];
                    features_holdout2 = [axial_feats_test2, coronal_feats_test2];
                elseif strcmp(rois{1}, 'Tumor')
                    axial_feats_train = load(fullfile(results_root, P{axial_tumor_wilcoxon_rf_index}, 'top5_training_feats.mat')).train_feats;
                    axial_feats_test1 = load(fullfile(results_root, P{axial_tumor_wilcoxon_rf_index}, 'top5_test1_feats.mat')).test1_feats;
                    axial_feats_test2 = load(fullfile(results_root, P{axial_tumor_wilcoxon_rf_index}, 'top5_test2_feats.mat')).test2_feats;
                    coronal_feats_train = load(fullfile(results_root, P{coronal_tumor_wilcoxon_rf_index}, 'top5_training_feats.mat')).train_feats;
                    coronal_feats_test1 = load(fullfile(results_root, P{coronal_tumor_wilcoxon_rf_index}, 'top5_test1_feats.mat')).test1_feats;
                    coronal_feats_test2 = load(fullfile(results_root, P{coronal_tumor_wilcoxon_rf_index}, 'top5_test2_feats.mat')).test2_feats;
                    axial_feat_names = load(fullfile(results_root, P{axial_tumor_wilcoxon_rf_index}, 'Top_Features_Information.mat')).topResults.names';
                    coronal_feat_names = load(fullfile(results_root, P{coronal_tumor_wilcoxon_rf_index}, 'Top_Features_Information.mat')).topResults.names';
                    feature_column_names = [axial_feat_names, coronal_feat_names];
                    features_training = [axial_feats_train, coronal_feats_train];
                    features_holdout1 = [axial_feats_test1, coronal_feats_test1];
                    features_holdout2 = [axial_feats_test2, coronal_feats_test2];
                end
            end
        elseif strcmp(plane{1}, 'Axial')
            if strcmp(best_fat_roi{1}, 'Proximal_Fat5')
                axial_fat_feats_train = load(fullfile(results_root, P{axial_proxfat5_wilcoxon_rf_index}, 'top5_training_feats.mat')).train_feats;
                axial_fat_feats_test1 = load(fullfile(results_root, P{axial_proxfat5_wilcoxon_rf_index}, 'top5_test1_feats.mat')).test1_feats;
                axial_fat_feats_test2 = load(fullfile(results_root, P{axial_proxfat5_wilcoxon_rf_index}, 'top5_test2_feats.mat')).test2_feats;
                axial_tumor_feats_train = load(fullfile(results_root, P{axial_tumor_wilcoxon_rf_index}, 'top5_training_feats.mat')).train_feats;
                axial_tumor_feats_test1 = load(fullfile(results_root, P{axial_tumor_wilcoxon_rf_index}, 'top5_test1_feats.mat')).test1_feats;
                axial_tumor_feats_test2 = load(fullfile(results_root, P{axial_tumor_wilcoxon_rf_index}, 'top5_test2_feats.mat')).test2_feats;
                axial_fat_feat_names = load(fullfile(results_root, P{axial_proxfat5_wilcoxon_rf_index}, 'Top_Features_Information.mat')).topResults.names';
                axial_tumor_feat_names = load(fullfile(results_root, P{axial_tumor_wilcoxon_rf_index}, 'Top_Features_Information.mat')).topResults.names';
                feature_column_names = [axial_fat_feat_names, axial_tumor_feat_names];
                features_training = [axial_fat_feats_train, axial_tumor_feats_train];
                features_holdout1 = [axial_fat_feats_test1, axial_tumor_feats_test1];
                features_holdout2 = [axial_fat_feats_test2, axial_tumor_feats_test2];
            elseif strcmp(best_fat_roi{1}, 'Proximal_Fat10')
                axial_fat_feats_train = load(fullfile(results_root, P{axial_proxfat10_wilcoxon_rf_index}, 'top5_training_feats.mat')).train_feats;
                axial_fat_feats_test1 = load(fullfile(results_root, P{axial_proxfat10_wilcoxon_rf_index}, 'top5_test1_feats.mat')).test1_feats;
                axial_fat_feats_test2 = load(fullfile(results_root, P{axial_proxfat10_wilcoxon_rf_index}, 'top5_test2_feats.mat')).test2_feats;
                axial_tumor_feats_train = load(fullfile(results_root, P{axial_tumor_wilcoxon_rf_index}, 'top5_training_feats.mat')).train_feats;
                axial_tumor_feats_test1 = load(fullfile(results_root, P{axial_tumor_wilcoxon_rf_index}, 'top5_test1_feats.mat')).test1_feats;
                axial_tumor_feats_test2 = load(fullfile(results_root, P{axial_tumor_wilcoxon_rf_index}, 'top5_test2_feats.mat')).test2_feats;
                axial_fat_feat_names = load(fullfile(results_root, P{axial_proxfat5_wilcoxon_rf_index}, 'Top_Features_Information.mat')).topResults.names';
                axial_tumor_feat_names = load(fullfile(results_root, P{axial_tumor_wilcoxon_rf_index}, 'Top_Features_Information.mat')).topResults.names';
                feature_column_names = [axial_fat_feat_names, axial_tumor_feat_names];
                features_training = [axial_fat_feats_train, axial_tumor_feats_train];
                features_holdout1 = [axial_fat_feats_test1, axial_tumor_feats_test1];
                features_holdout2 = [axial_fat_feats_test2, axial_tumor_feats_test2];
            elseif strcmp(best_fat_roi{1}, 'Fat')
                axial_fat_feats_train = load(fullfile(results_root, P{axial_fat_wilcoxon_rf_index}, 'top5_training_feats.mat')).train_feats;
                axial_fat_feats_test1 = load(fullfile(results_root, P{axial_fat_wilcoxon_rf_index}, 'top5_test1_feats.mat')).test1_feats;
                axial_fat_feats_test2 = load(fullfile(results_root, P{axial_fat_wilcoxon_rf_index}, 'top5_test2_feats.mat')).test2_feats;
                axial_tumor_feats_train = load(fullfile(results_root, P{axial_tumor_wilcoxon_rf_index}, 'top5_training_feats.mat')).train_feats;
                axial_tumor_feats_test1 = load(fullfile(results_root, P{axial_tumor_wilcoxon_rf_index}, 'top5_test1_feats.mat')).test1_feats;
                axial_tumor_feats_test2 = load(fullfile(results_root, P{axial_tumor_wilcoxon_rf_index}, 'top5_test2_feats.mat')).test2_feats;
                axial_fat_feat_names = load(fullfile(results_root, P{axial_fat_wilcoxon_rf_index}, 'Top_Features_Information.mat')).topResults.names';
                axial_tumor_feat_names = load(fullfile(results_root, P{axial_tumor_wilcoxon_rf_index}, 'Top_Features_Information.mat')).topResults.names';
                feature_column_names = [axial_fat_feat_names, axial_tumor_feat_names];
                features_training = [axial_fat_feats_train, axial_tumor_feats_train];
                features_holdout1 = [axial_fat_feats_test1, axial_tumor_feats_test1];
                features_holdout2 = [axial_fat_feats_test2, axial_tumor_feats_test2];
            end
        elseif strcmp(plane{1}, 'Coronal')
            if strcmp(best_fat_roi{1}, 'Proximal_Fat5')
                coronal_fat_feats_train = load(fullfile(results_root, P{coronal_proxfat5_wilcoxon_rf_index}, 'top5_training_feats.mat')).train_feats;
                coronal_fat_feats_test1 = load(fullfile(results_root, P{coronal_proxfat5_wilcoxon_rf_index}, 'top5_test1_feats.mat')).test1_feats;
                coronal_fat_feats_test2 = load(fullfile(results_root, P{coronal_proxfat5_wilcoxon_rf_index}, 'top5_test2_feats.mat')).test2_feats;
                coronal_tumor_feats_train = load(fullfile(results_root, P{coronal_tumor_wilcoxon_rf_index}, 'top5_training_feats.mat')).train_feats;
                coronal_tumor_feats_test1 = load(fullfile(results_root, P{coronal_tumor_wilcoxon_rf_index}, 'top5_test1_feats.mat')).test1_feats;
                coronal_tumor_feats_test2 = load(fullfile(results_root, P{coronal_tumor_wilcoxon_rf_index}, 'top5_test2_feats.mat')).test2_feats;
                coronal_fat_feat_names = load(fullfile(results_root, P{coronal_proxfat5_wilcoxon_rf_index}, 'Top_Features_Information.mat')).topResults.names';
                coronal_tumor_feat_names = load(fullfile(results_root, P{coronal_tumor_wilcoxon_rf_index}, 'Top_Features_Information.mat')).topResults.names';
                feature_column_names = [coronal_fat_feat_names, coronal_tumor_feat_names];
                features_training = [coronal_fat_feats_train, coronal_tumor_feats_train];
                features_holdout1 = [coronal_fat_feats_test1, coronal_tumor_feats_test1];
                features_holdout2 = [coronal_fat_feats_test2, coronal_tumor_feats_test2];
            elseif strcmp(best_fat_roi{1}, 'Proximal_Fat10')
                coronal_fat_feats_train = load(fullfile(results_root, P{coronal_proxfat10_wilcoxon_rf_index}, 'top5_training_feats.mat')).train_feats;
                coronal_fat_feats_test1 = load(fullfile(results_root, P{coronal_proxfat10_wilcoxon_rf_index}, 'top5_test1_feats.mat')).test1_feats;
                coronal_fat_feats_test2 = load(fullfile(results_root, P{coronal_proxfat10_wilcoxon_rf_index}, 'top5_test2_feats.mat')).test2_feats;
                coronal_tumor_feats_train = load(fullfile(results_root, P{coronal_tumor_wilcoxon_rf_index}, 'top5_training_feats.mat')).train_feats;
                coronal_tumor_feats_test1 = load(fullfile(results_root, P{coronal_tumor_wilcoxon_rf_index}, 'top5_test1_feats.mat')).test1_feats;
                coronal_tumor_feats_test2 = load(fullfile(results_root, P{coronal_tumor_wilcoxon_rf_index}, 'top5_test2_feats.mat')).test2_feats;
                coronal_fat_feat_names = load(fullfile(results_root, P{coronal_proxfat5_wilcoxon_rf_index}, 'Top_Features_Information.mat')).topResults.names';
                coronal_tumor_feat_names = load(fullfile(results_root, P{coronal_tumor_wilcoxon_rf_index}, 'Top_Features_Information.mat')).topResults.names';
                feature_column_names = [coronal_fat_feat_names, coronal_tumor_feat_names];
                features_training = [coronal_fat_feats_train, coronal_tumor_feats_train];
                features_holdout1 = [coronal_fat_feats_test1, coronal_tumor_feats_test1];
                features_holdout2 = [coronal_fat_feats_test2, coronal_tumor_feats_test2];
            elseif strcmp(best_fat_roi{1}, 'Fat')
                coronal_fat_feats_train = load(fullfile(results_root, P{coronal_fat_wilcoxon_rf_index}, 'top5_training_feats.mat')).train_feats;
                coronal_fat_feats_test1 = load(fullfile(results_root, P{coronal_fat_wilcoxon_rf_index}, 'top5_test1_feats.mat')).test1_feats;
                coronal_fat_feats_test2 = load(fullfile(results_root, P{coronal_fat_wilcoxon_rf_index}, 'top5_test2_feats.mat')).test2_feats;
                coronal_tumor_feats_train = load(fullfile(results_root, P{coronal_tumor_wilcoxon_rf_index}, 'top5_training_feats.mat')).train_feats;
                coronal_tumor_feats_test1 = load(fullfile(results_root, P{coronal_tumor_wilcoxon_rf_index}, 'top5_test1_feats.mat')).test1_feats;
                coronal_tumor_feats_test2 = load(fullfile(results_root, P{coronal_tumor_wilcoxon_rf_index}, 'top5_test2_feats.mat')).test2_feats;
                coronal_fat_feat_names = load(fullfile(results_root, P{coronal_fat_wilcoxon_rf_index}, 'Top_Features_Information.mat')).topResults.names';
                coronal_tumor_feat_names = load(fullfile(results_root, P{coronal_tumor_wilcoxon_rf_index}, 'Top_Features_Information.mat')).topResults.names';
                feature_column_names = [coronal_fat_feat_names, coronal_tumor_feat_names];
                features_training = [coronal_fat_feats_train, coronal_tumor_feats_train];
                features_holdout1 = [coronal_fat_feats_test1, coronal_tumor_feats_test1];
                features_holdout2 = [coronal_fat_feats_test2, coronal_tumor_feats_test2];
            end
        end
    end
elseif strcmp(selector(1), 'mrmr')
    if strcmp(classifier{1}, 'qda')
        if strcmp(plane{1}, 'Multi-Plane')
            if length(rois) == 2
                if strcmp(best_fat_roi{1}, 'Proximal_Fat5')
                    axial_fat_feats_train = load(fullfile(results_root, P{axial_proxfat5_mrmr_qda_index}, 'top5_training_feats.mat')).train_feats;
                    axial_fat_feats_test1 = load(fullfile(results_root, P{axial_proxfat5_mrmr_qda_index}, 'top5_test1_feats.mat')).test1_feats;
                    axial_fat_feats_test2 = load(fullfile(results_root, P{axial_proxfat5_mrmr_qda_index}, 'top5_test2_feats.mat')).test2_feats;
                    axial_tumor_feats_train = load(fullfile(results_root, P{axial_tumor_mrmr_qda_index}, 'top5_training_feats.mat')).train_feats;
                    axial_tumor_feats_test1 = load(fullfile(results_root, P{axial_tumor_mrmr_qda_index}, 'top5_test1_feats.mat')).test1_feats;
                    axial_tumor_feats_test2 = load(fullfile(results_root, P{axial_tumor_mrmr_qda_index}, 'top5_test2_feats.mat')).test2_feats;
                    coronal_fat_feats_train = load(fullfile(results_root, P{coronal_proxfat5_mrmr_qda_index}, 'top5_training_feats.mat')).train_feats;
                    coronal_fat_feats_test1 = load(fullfile(results_root, P{coronal_proxfat5_mrmr_qda_index}, 'top5_test1_feats.mat')).test1_feats;
                    coronal_fat_feats_test2 = load(fullfile(results_root, P{coronal_proxfat5_mrmr_qda_index}, 'top5_test2_feats.mat')).test2_feats;
                    coronal_tumor_feats_train = load(fullfile(results_root, P{coronal_tumor_mrmr_qda_index}, 'top5_training_feats.mat')).train_feats;
                    coronal_tumor_feats_test1 = load(fullfile(results_root, P{coronal_tumor_mrmr_qda_index}, 'top5_test1_feats.mat')).test1_feats;
                    coronal_tumor_feats_test2 = load(fullfile(results_root, P{coronal_tumor_mrmr_qda_index}, 'top5_test2_feats.mat')).test2_feats;
                    axial_fat_feat_names = load(fullfile(results_root, P{axial_proxfat5_mrmr_qda_index}, 'Top_Features_Information.mat')).topResults.names';
                    coronal_fat_feat_names = load(fullfile(results_root, P{coronal_proxfat5_mrmr_qda_index}, 'Top_Features_Information.mat')).topResults.names';
                    axial_tumor_feat_names = load(fullfile(results_root, P{axial_tumor_mrmr_qda_index}, 'Top_Features_Information.mat')).topResults.names';
                    coronal_tumor_feat_names = load(fullfile(results_root, P{coronal_tumor_mrmr_qda_index}, 'Top_Features_Information.mat')).topResults.names';
                    feature_column_names = [axial_fat_feat_names, coronal_fat_feat_names, axial_tumor_feat_names, coronal_tumor_feat_names];
                    features_training = [axial_fat_feats_train, coronal_fat_feats_train, axial_tumor_feats_train, coronal_tumor_feats_train];
                    features_holdout1 = [axial_fat_feats_test1, coronal_fat_feats_test1, axial_tumor_feats_test1, coronal_tumor_feats_test1];
                    features_holdout2 = [axial_fat_feats_test2, coronal_fat_feats_test2, axial_tumor_feats_test2, coronal_tumor_feats_test2];
                elseif strcmp(best_fat_roi{1}, 'Proximal_Fat10')
                    axial_fat_feats_train = load(fullfile(results_root, P{axial_proxfat10_mrmr_qda_index}, 'top5_training_feats.mat')).train_feats;
                    axial_fat_feats_test1 = load(fullfile(results_root, P{axial_proxfat10_mrmr_qda_index}, 'top5_test1_feats.mat')).test1_feats;
                    axial_fat_feats_test2 = load(fullfile(results_root, P{axial_proxfat10_mrmr_qda_index}, 'top5_test2_feats.mat')).test2_feats;
                    axial_tumor_feats_train = load(fullfile(results_root, P{axial_tumor_mrmr_qda_index}, 'top5_training_feats.mat')).train_feats;
                    axial_tumor_feats_test1 = load(fullfile(results_root, P{axial_tumor_mrmr_qda_index}, 'top5_test1_feats.mat')).test1_feats;
                    axial_tumor_feats_test2 = load(fullfile(results_root, P{axial_tumor_mrmr_qda_index}, 'top5_test2_feats.mat')).test2_feats;
                    coronal_fat_feats_train = load(fullfile(results_root, P{coronal_proxfat10_mrmr_qda_index}, 'top5_training_feats.mat')).train_feats;
                    coronal_fat_feats_test1 = load(fullfile(results_root, P{coronal_proxfat10_mrmr_qda_index}, 'top5_test1_feats.mat')).test1_feats;
                    coronal_fat_feats_test2 = load(fullfile(results_root, P{coronal_proxfat10_mrmr_qda_index}, 'top5_test2_feats.mat')).test2_feats;
                    coronal_tumor_feats_train = load(fullfile(results_root, P{coronal_tumor_mrmr_qda_index}, 'top5_training_feats.mat')).train_feats;
                    coronal_tumor_feats_test1 = load(fullfile(results_root, P{coronal_tumor_mrmr_qda_index}, 'top5_test1_feats.mat')).test1_feats;
                    coronal_tumor_feats_test2 = load(fullfile(results_root, P{coronal_tumor_mrmr_qda_index}, 'top5_test2_feats.mat')).test2_feats;
                    axial_fat_feat_names = load(fullfile(results_root, P{axial_proxfat5_mrmr_qda_index}, 'Top_Features_Information.mat')).topResults.names';
                    coronal_fat_feat_names = load(fullfile(results_root, P{coronal_proxfat5_mrmr_qda_index}, 'Top_Features_Information.mat')).topResults.names';
                    axial_tumor_feat_names = load(fullfile(results_root, P{axial_tumor_mrmr_qda_index}, 'Top_Features_Information.mat')).topResults.names';
                    coronal_tumor_feat_names = load(fullfile(results_root, P{coronal_tumor_mrmr_qda_index}, 'Top_Features_Information.mat')).topResults.names';
                    feature_column_names = [axial_fat_feat_names, coronal_fat_feat_names, axial_tumor_feat_names, coronal_tumor_feat_names];
                    features_training = [axial_fat_feats_train, coronal_fat_feats_train, axial_tumor_feats_train, coronal_tumor_feats_train];
                    features_holdout1 = [axial_fat_feats_test1, coronal_fat_feats_test1, axial_tumor_feats_test1, coronal_tumor_feats_test1];
                    features_holdout2 = [axial_fat_feats_test2, coronal_fat_feats_test2, axial_tumor_feats_test2, coronal_tumor_feats_test2];
                elseif strcmp(best_fat_roi{1}, 'Fat')
                    axial_fat_feats_train = load(fullfile(results_root, P{axial_fat_mrmr_qda_index}, 'top5_training_feats.mat')).train_feats;
                    axial_fat_feats_test1 = load(fullfile(results_root, P{axial_fat_mrmr_qda_index}, 'top5_test1_feats.mat')).test1_feats;
                    axial_fat_feats_test2 = load(fullfile(results_root, P{axial_fat_mrmr_qda_index}, 'top5_test2_feats.mat')).test2_feats;
                    axial_tumor_feats_train = load(fullfile(results_root, P{axial_tumor_mrmr_qda_index}, 'top5_training_feats.mat')).train_feats;
                    axial_tumor_feats_test1 = load(fullfile(results_root, P{axial_tumor_mrmr_qda_index}, 'top5_test1_feats.mat')).test1_feats;
                    axial_tumor_feats_test2 = load(fullfile(results_root, P{axial_tumor_mrmr_qda_index}, 'top5_test2_feats.mat')).test2_feats;
                    coronal_fat_feats_train = load(fullfile(results_root, P{coronal_fat_mrmr_qda_index}, 'top5_training_feats.mat')).train_feats;
                    coronal_fat_feats_test1 = load(fullfile(results_root, P{coronal_fat_mrmr_qda_index}, 'top5_test1_feats.mat')).test1_feats;
                    coronal_fat_feats_test2 = load(fullfile(results_root, P{coronal_fat_mrmr_qda_index}, 'top5_test2_feats.mat')).test2_feats;
                    coronal_tumor_feats_train = load(fullfile(results_root, P{coronal_tumor_mrmr_qda_index}, 'top5_training_feats.mat')).train_feats;
                    coronal_tumor_feats_test1 = load(fullfile(results_root, P{coronal_tumor_mrmr_qda_index}, 'top5_test1_feats.mat')).test1_feats;
                    coronal_tumor_feats_test2 = load(fullfile(results_root, P{coronal_tumor_mrmr_qda_index}, 'top5_test2_feats.mat')).test2_feats;
                    axial_fat_feat_names = load(fullfile(results_root, P{axial_fat_mrmr_qda_index}, 'Top_Features_Information.mat')).topResults.names';
                    coronal_fat_feat_names = load(fullfile(results_root, P{coronal_fat_mrmr_qda_index}, 'Top_Features_Information.mat')).topResults.names';
                    axial_tumor_feat_names = load(fullfile(results_root, P{axial_tumor_mrmr_qda_index}, 'Top_Features_Information.mat')).topResults.names';
                    coronal_tumor_feat_names = load(fullfile(results_root, P{coronal_tumor_mrmr_qda_index}, 'Top_Features_Information.mat')).topResults.names';
                    feature_column_names = [axial_fat_feat_names, coronal_fat_feat_names, axial_tumor_feat_names, coronal_tumor_feat_names];
                    features_training = [axial_fat_feats_train, coronal_fat_feats_train, axial_tumor_feats_train, coronal_tumor_feats_train];
                    features_holdout1 = [axial_fat_feats_test1, coronal_fat_feats_test1, axial_tumor_feats_test1, coronal_tumor_feats_test1];
                    features_holdout2 = [axial_fat_feats_test2, coronal_fat_feats_test2, axial_tumor_feats_test2, coronal_tumor_feats_test2];
                end
            elseif length(rois) == 1
                if strcmp(rois{1}, 'Proximal_Fat5')
                    axial_feats_train = load(fullfile(results_root, P{axial_proxfat5_mrmr_qda_index}, 'top5_training_feats.mat')).train_feats;
                    axial_feats_test1 = load(fullfile(results_root, P{axial_proxfat5_mrmr_qda_index}, 'top5_test1_feats.mat')).test1_feats;
                    axial_feats_test2 = load(fullfile(results_root, P{axial_proxfat5_mrmr_qda_index}, 'top5_test2_feats.mat')).test2_feats;
                    coronal_feats_train = load(fullfile(results_root, P{coronal_proxfat5_mrmr_qda_index}, 'top5_training_feats.mat')).train_feats;
                    coronal_feats_test1 = load(fullfile(results_root, P{coronal_proxfat5_mrmr_qda_index}, 'top5_test1_feats.mat')).test1_feats;
                    coronal_feats_test2 = load(fullfile(results_root, P{coronal_proxfat5_mrmr_qda_index}, 'top5_test2_feats.mat')).test2_feats;
                    axial_feat_names = load(fullfile(results_root, P{axial_proxfat5_mrmr_qda_index}, 'Top_Features_Information.mat')).topResults.names';
                    coronal_feat_names = load(fullfile(results_root, P{coronal_proxfat5_mrmr_qda_index}, 'Top_Features_Information.mat')).topResults.names';
                    feature_column_names = [axial_feat_names, coronal_feat_names];
                    features_training = [axial_feats_train, coronal_feats_train];
                    features_holdout1 = [axial_feats_test1, coronal_feats_test1];
                    features_holdout2 = [axial_feats_test2, coronal_feats_test2];
                elseif strcmp(rois{1}, 'Proximal_Fat10')
                    axial_feats_train = load(fullfile(results_root, P{axial_proxfat10_mrmr_qda_index}, 'top5_training_feats.mat')).train_feats;
                    axial_feats_test1 = load(fullfile(results_root, P{axial_proxfat10_mrmr_qda_index}, 'top5_test1_feats.mat')).test1_feats;
                    axial_feats_test2 = load(fullfile(results_root, P{axial_proxfat10_mrmr_qda_index}, 'top5_test2_feats.mat')).test2_feats;
                    coronal_feats_train = load(fullfile(results_root, P{coronal_proxfat10_mrmr_qda_index}, 'top5_training_feats.mat')).train_feats;
                    coronal_feats_test1 = load(fullfile(results_root, P{coronal_proxfat10_mrmr_qda_index}, 'top5_test1_feats.mat')).test1_feats;
                    coronal_feats_test2 = load(fullfile(results_root, P{coronal_proxfat10_mrmr_qda_index}, 'top5_test2_feats.mat')).test2_feats;
                    axial_feat_names = load(fullfile(results_root, P{axial_proxfat10_mrmr_qda_index}, 'Top_Features_Information.mat')).topResults.names';
                    coronal_feat_names = load(fullfile(results_root, P{coronal_proxfat10_mrmr_qda_index}, 'Top_Features_Information.mat')).topResults.names';
                    feature_column_names = [axial_feat_names, coronal_feat_names];
                    features_training = [axial_feats_train, coronal_feats_train];
                    features_holdout1 = [axial_feats_test1, coronal_feats_test1];
                    features_holdout2 = [axial_feats_test2, coronal_feats_test2];
                elseif strcmp(rois{1}, 'Fat')
                    axial_feats_train = load(fullfile(results_root, P{axial_fat_mrmr_qda_index}, 'top5_training_feats.mat')).train_feats;
                    axial_feats_test1 = load(fullfile(results_root, P{axial_fat_mrmr_qda_index}, 'top5_test1_feats.mat')).test1_feats;
                    axial_feats_test2 = load(fullfile(results_root, P{axial_fat_mrmr_qda_index}, 'top5_test2_feats.mat')).test2_feats;
                    coronal_feats_train = load(fullfile(results_root, P{coronal_fat_mrmr_qda_index}, 'top5_training_feats.mat')).train_feats;
                    coronal_feats_test1 = load(fullfile(results_root, P{coronal_fat_mrmr_qda_index}, 'top5_test1_feats.mat')).test1_feats;
                    coronal_feats_test2 = load(fullfile(results_root, P{coronal_fat_mrmr_qda_index}, 'top5_test2_feats.mat')).test2_feats;
                    axial_feat_names = load(fullfile(results_root, P{axial_fat_mrmr_qda_index}, 'Top_Features_Information.mat')).topResults.names';
                    coronal_feat_names = load(fullfile(results_root, P{coronal_fat_mrmr_qda_index}, 'Top_Features_Information.mat')).topResults.names';
                    feature_column_names = [axial_feat_names, coronal_feat_names];
                    features_training = [axial_feats_train, coronal_feats_train];
                    features_holdout1 = [axial_feats_test1, coronal_feats_test1];
                    features_holdout2 = [axial_feats_test2, coronal_feats_test2];
                elseif strcmp(rois{1}, 'Tumor')
                    axial_feats_train = load(fullfile(results_root, P{axial_tumor_mrmr_qda_index}, 'top5_training_feats.mat')).train_feats;
                    axial_feats_test1 = load(fullfile(results_root, P{axial_tumor_mrmr_qda_index}, 'top5_test1_feats.mat')).test1_feats;
                    axial_feats_test2 = load(fullfile(results_root, P{axial_tumor_mrmr_qda_index}, 'top5_test2_feats.mat')).test2_feats;
                    coronal_feats_train = load(fullfile(results_root, P{coronal_tumor_mrmr_qda_index}, 'top5_training_feats.mat')).train_feats;
                    coronal_feats_test1 = load(fullfile(results_root, P{coronal_tumor_mrmr_qda_index}, 'top5_test1_feats.mat')).test1_feats;
                    coronal_feats_test2 = load(fullfile(results_root, P{coronal_tumor_mrmr_qda_index}, 'top5_test2_feats.mat')).test2_feats;
                    axial_feat_names = load(fullfile(results_root, P{axial_tumor_mrmr_qda_index}, 'Top_Features_Information.mat')).topResults.names';
                    coronal_feat_names = load(fullfile(results_root, P{coronal_tumor_mrmr_qda_index}, 'Top_Features_Information.mat')).topResults.names';
                    feature_column_names = [axial_feat_names, coronal_feat_names];
                    features_training = [axial_feats_train, coronal_feats_train];
                    features_holdout1 = [axial_feats_test1, coronal_feats_test1];
                    features_holdout2 = [axial_feats_test2, coronal_feats_test2];
                end
            end
        elseif strcmp(plane{1}, 'Axial')
            if strcmp(best_fat_roi{1}, 'Proximal_Fat5')
                axial_fat_feats_train = load(fullfile(results_root, P{axial_proxfat5_mrmr_qda_index}, 'top5_training_feats.mat')).train_feats;
                axial_fat_feats_test1 = load(fullfile(results_root, P{axial_proxfat5_mrmr_qda_index}, 'top5_test1_feats.mat')).test1_feats;
                axial_fat_feats_test2 = load(fullfile(results_root, P{axial_proxfat5_mrmr_qda_index}, 'top5_test2_feats.mat')).test2_feats;
                axial_tumor_feats_train = load(fullfile(results_root, P{axial_tumor_mrmr_qda_index}, 'top5_training_feats.mat')).train_feats;
                axial_tumor_feats_test1 = load(fullfile(results_root, P{axial_tumor_mrmr_qda_index}, 'top5_test1_feats.mat')).test1_feats;
                axial_tumor_feats_test2 = load(fullfile(results_root, P{axial_tumor_mrmr_qda_index}, 'top5_test2_feats.mat')).test2_feats;
                axial_fat_feat_names = load(fullfile(results_root, P{axial_proxfat5_mrmr_qda_index}, 'Top_Features_Information.mat')).topResults.names';
                axial_tumor_feat_names = load(fullfile(results_root, P{axial_tumor_mrmr_qda_index}, 'Top_Features_Information.mat')).topResults.names';
                feature_column_names = [axial_fat_feat_names, axial_tumor_feat_names];
                features_training = [axial_fat_feats_train, axial_tumor_feats_train];
                features_holdout1 = [axial_fat_feats_test1, axial_tumor_feats_test1];
                features_holdout2 = [axial_fat_feats_test2, axial_tumor_feats_test2];
            elseif strcmp(best_fat_roi{1}, 'Proximal_Fat10')
                axial_fat_feats_train = load(fullfile(results_root, P{axial_proxfat10_mrmr_qda_index}, 'top5_training_feats.mat')).train_feats;
                axial_fat_feats_test1 = load(fullfile(results_root, P{axial_proxfat10_mrmr_qda_index}, 'top5_test1_feats.mat')).test1_feats;
                axial_fat_feats_test2 = load(fullfile(results_root, P{axial_proxfat10_mrmr_qda_index}, 'top5_test2_feats.mat')).test2_feats;
                axial_tumor_feats_train = load(fullfile(results_root, P{axial_tumor_mrmr_qda_index}, 'top5_training_feats.mat')).train_feats;
                axial_tumor_feats_test1 = load(fullfile(results_root, P{axial_tumor_mrmr_qda_index}, 'top5_test1_feats.mat')).test1_feats;
                axial_tumor_feats_test2 = load(fullfile(results_root, P{axial_tumor_mrmr_qda_index}, 'top5_test2_feats.mat')).test2_feats;
                axial_fat_feat_names = load(fullfile(results_root, P{axial_proxfat5_mrmr_qda_index}, 'Top_Features_Information.mat')).topResults.names';
                axial_tumor_feat_names = load(fullfile(results_root, P{axial_tumor_mrmr_qda_index}, 'Top_Features_Information.mat')).topResults.names';
                feature_column_names = [axial_fat_feat_names, axial_tumor_feat_names];
                features_training = [axial_fat_feats_train, axial_tumor_feats_train];
                features_holdout1 = [axial_fat_feats_test1, axial_tumor_feats_test1];
                features_holdout2 = [axial_fat_feats_test2, axial_tumor_feats_test2];
            elseif strcmp(best_fat_roi{1}, 'Fat')
                axial_fat_feats_train = load(fullfile(results_root, P{axial_fat_mrmr_qda_index}, 'top5_training_feats.mat')).train_feats;
                axial_fat_feats_test1 = load(fullfile(results_root, P{axial_fat_mrmr_qda_index}, 'top5_test1_feats.mat')).test1_feats;
                axial_fat_feats_test2 = load(fullfile(results_root, P{axial_fat_mrmr_qda_index}, 'top5_test2_feats.mat')).test2_feats;
                axial_tumor_feats_train = load(fullfile(results_root, P{axial_tumor_mrmr_qda_index}, 'top5_training_feats.mat')).train_feats;
                axial_tumor_feats_test1 = load(fullfile(results_root, P{axial_tumor_mrmr_qda_index}, 'top5_test1_feats.mat')).test1_feats;
                axial_tumor_feats_test2 = load(fullfile(results_root, P{axial_tumor_mrmr_qda_index}, 'top5_test2_feats.mat')).test2_feats;
                axial_fat_feat_names = load(fullfile(results_root, P{axial_fat_mrmr_qda_index}, 'Top_Features_Information.mat')).topResults.names';
                axial_tumor_feat_names = load(fullfile(results_root, P{axial_tumor_mrmr_qda_index}, 'Top_Features_Information.mat')).topResults.names';
                feature_column_names = [axial_fat_feat_names, axial_tumor_feat_names];
                features_training = [axial_fat_feats_train, axial_tumor_feats_train];
                features_holdout1 = [axial_fat_feats_test1, axial_tumor_feats_test1];
                features_holdout2 = [axial_fat_feats_test2, axial_tumor_feats_test2];
            end
        elseif strcmp(plane{1}, 'Coronal')
            if strcmp(best_fat_roi{1}, 'Proximal_Fat5')
                coronal_fat_feats_train = load(fullfile(results_root, P{coronal_proxfat5_mrmr_qda_index}, 'top5_training_feats.mat')).train_feats;
                coronal_fat_feats_test1 = load(fullfile(results_root, P{coronal_proxfat5_mrmr_qda_index}, 'top5_test1_feats.mat')).test1_feats;
                coronal_fat_feats_test2 = load(fullfile(results_root, P{coronal_proxfat5_mrmr_qda_index}, 'top5_test2_feats.mat')).test2_feats;
                coronal_tumor_feats_train = load(fullfile(results_root, P{coronal_tumor_mrmr_qda_index}, 'top5_training_feats.mat')).train_feats;
                coronal_tumor_feats_test1 = load(fullfile(results_root, P{coronal_tumor_mrmr_qda_index}, 'top5_test1_feats.mat')).test1_feats;
                coronal_tumor_feats_test2 = load(fullfile(results_root, P{coronal_tumor_mrmr_qda_index}, 'top5_test2_feats.mat')).test2_feats;
                coronal_fat_feat_names = load(fullfile(results_root, P{coronal_proxfat5_mrmr_qda_index}, 'Top_Features_Information.mat')).topResults.names';
                coronal_tumor_feat_names = load(fullfile(results_root, P{coronal_tumor_mrmr_qda_index}, 'Top_Features_Information.mat')).topResults.names';
                feature_column_names = [coronal_fat_feat_names, coronal_tumor_feat_names];
                features_training = [coronal_fat_feats_train, coronal_tumor_feats_train];
                features_holdout1 = [coronal_fat_feats_test1, coronal_tumor_feats_test1];
                features_holdout2 = [coronal_fat_feats_test2, coronal_tumor_feats_test2];
            elseif strcmp(best_fat_roi{1}, 'Proximal_Fat10')
                coronal_fat_feats_train = load(fullfile(results_root, P{coronal_proxfat10_mrmr_qda_index}, 'top5_training_feats.mat')).train_feats;
                coronal_fat_feats_test1 = load(fullfile(results_root, P{coronal_proxfat10_mrmr_qda_index}, 'top5_test1_feats.mat')).test1_feats;
                coronal_fat_feats_test2 = load(fullfile(results_root, P{coronal_proxfat10_mrmr_qda_index}, 'top5_test2_feats.mat')).test2_feats;
                coronal_tumor_feats_train = load(fullfile(results_root, P{coronal_tumor_mrmr_qda_index}, 'top5_training_feats.mat')).train_feats;
                coronal_tumor_feats_test1 = load(fullfile(results_root, P{coronal_tumor_mrmr_qda_index}, 'top5_test1_feats.mat')).test1_feats;
                coronal_tumor_feats_test2 = load(fullfile(results_root, P{coronal_tumor_mrmr_qda_index}, 'top5_test2_feats.mat')).test2_feats;
                coronal_fat_feat_names = load(fullfile(results_root, P{coronal_proxfat5_mrmr_qda_index}, 'Top_Features_Information.mat')).topResults.names';
                coronal_tumor_feat_names = load(fullfile(results_root, P{coronal_tumor_mrmr_qda_index}, 'Top_Features_Information.mat')).topResults.names';
                feature_column_names = [coronal_fat_feat_names, coronal_tumor_feat_names];
                features_training = [coronal_fat_feats_train, coronal_tumor_feats_train];
                features_holdout1 = [coronal_fat_feats_test1, coronal_tumor_feats_test1];
                features_holdout2 = [coronal_fat_feats_test2, coronal_tumor_feats_test2];
            elseif strcmp(best_fat_roi{1}, 'Fat')
                coronal_fat_feats_train = load(fullfile(results_root, P{coronal_fat_mrmr_qda_index}, 'top5_training_feats.mat')).train_feats;
                coronal_fat_feats_test1 = load(fullfile(results_root, P{coronal_fat_mrmr_qda_index}, 'top5_test1_feats.mat')).test1_feats;
                coronal_fat_feats_test2 = load(fullfile(results_root, P{coronal_fat_mrmr_qda_index}, 'top5_test2_feats.mat')).test2_feats;
                coronal_tumor_feats_train = load(fullfile(results_root, P{coronal_tumor_mrmr_qda_index}, 'top5_training_feats.mat')).train_feats;
                coronal_tumor_feats_test1 = load(fullfile(results_root, P{coronal_tumor_mrmr_qda_index}, 'top5_test1_feats.mat')).test1_feats;
                coronal_tumor_feats_test2 = load(fullfile(results_root, P{coronal_tumor_mrmr_qda_index}, 'top5_test2_feats.mat')).test2_feats;
                coronal_fat_feat_names = load(fullfile(results_root, P{coronal_fat_mrmr_qda_index}, 'Top_Features_Information.mat')).topResults.names';
                coronal_tumor_feat_names = load(fullfile(results_root, P{coronal_tumor_mrmr_qda_index}, 'Top_Features_Information.mat')).topResults.names';
                feature_column_names = [coronal_fat_feat_names, coronal_tumor_feat_names];
                features_training = [coronal_fat_feats_train, coronal_tumor_feats_train];
                features_holdout1 = [coronal_fat_feats_test1, coronal_tumor_feats_test1];
                features_holdout2 = [coronal_fat_feats_test2, coronal_tumor_feats_test2];
            end
        end
    elseif strcmp(classifier(1), 'lda')
        if strcmp(plane{1}, 'Multi-Plane')
            if length(rois) == 2
                if strcmp(best_fat_roi{1}, 'Proximal_Fat5')
                    axial_fat_feats_train = load(fullfile(results_root, P{axial_proxfat5_mrmr_lda_index}, 'top5_training_feats.mat')).train_feats;
                    axial_fat_feats_test1 = load(fullfile(results_root, P{axial_proxfat5_mrmr_lda_index}, 'top5_test1_feats.mat')).test1_feats;
                    axial_fat_feats_test2 = load(fullfile(results_root, P{axial_proxfat5_mrmr_lda_index}, 'top5_test2_feats.mat')).test2_feats;
                    axial_tumor_feats_train = load(fullfile(results_root, P{axial_tumor_mrmr_lda_index}, 'top5_training_feats.mat')).train_feats;
                    axial_tumor_feats_test1 = load(fullfile(results_root, P{axial_tumor_mrmr_lda_index}, 'top5_test1_feats.mat')).test1_feats;
                    axial_tumor_feats_test2 = load(fullfile(results_root, P{axial_tumor_mrmr_lda_index}, 'top5_test2_feats.mat')).test2_feats;
                    coronal_fat_feats_train = load(fullfile(results_root, P{coronal_proxfat5_mrmr_lda_index}, 'top5_training_feats.mat')).train_feats;
                    coronal_fat_feats_test1 = load(fullfile(results_root, P{coronal_proxfat5_mrmr_lda_index}, 'top5_test1_feats.mat')).test1_feats;
                    coronal_fat_feats_test2 = load(fullfile(results_root, P{coronal_proxfat5_mrmr_lda_index}, 'top5_test2_feats.mat')).test2_feats;
                    coronal_tumor_feats_train = load(fullfile(results_root, P{coronal_tumor_mrmr_lda_index}, 'top5_training_feats.mat')).train_feats;
                    coronal_tumor_feats_test1 = load(fullfile(results_root, P{coronal_tumor_mrmr_lda_index}, 'top5_test1_feats.mat')).test1_feats;
                    coronal_tumor_feats_test2 = load(fullfile(results_root, P{coronal_tumor_mrmr_lda_index}, 'top5_test2_feats.mat')).test2_feats;
                    axial_fat_feat_names = load(fullfile(results_root, P{axial_proxfat5_mrmr_lda_index}, 'Top_Features_Information.mat')).topResults.names';
                    coronal_fat_feat_names = load(fullfile(results_root, P{coronal_proxfat5_mrmr_lda_index}, 'Top_Features_Information.mat')).topResults.names';
                    axial_tumor_feat_names = load(fullfile(results_root, P{axial_tumor_mrmr_lda_index}, 'Top_Features_Information.mat')).topResults.names';
                    coronal_tumor_feat_names = load(fullfile(results_root, P{coronal_tumor_mrmr_lda_index}, 'Top_Features_Information.mat')).topResults.names';
                    feature_column_names = [axial_fat_feat_names, coronal_fat_feat_names, axial_tumor_feat_names, coronal_tumor_feat_names];
                    features_training = [axial_fat_feats_train, coronal_fat_feats_train, axial_tumor_feats_train, coronal_tumor_feats_train];
                    features_holdout1 = [axial_fat_feats_test1, coronal_fat_feats_test1, axial_tumor_feats_test1, coronal_tumor_feats_test1];
                    features_holdout2 = [axial_fat_feats_test2, coronal_fat_feats_test2, axial_tumor_feats_test2, coronal_tumor_feats_test2];
                elseif strcmp(best_fat_roi{1}, 'Proximal_Fat10')
                    axial_fat_feats_train = load(fullfile(results_root, P{axial_proxfat10_mrmr_lda_index}, 'top5_training_feats.mat')).train_feats;
                    axial_fat_feats_test1 = load(fullfile(results_root, P{axial_proxfat10_mrmr_lda_index}, 'top5_test1_feats.mat')).test1_feats;
                    axial_fat_feats_test2 = load(fullfile(results_root, P{axial_proxfat10_mrmr_lda_index}, 'top5_test2_feats.mat')).test2_feats;
                    axial_tumor_feats_train = load(fullfile(results_root, P{axial_tumor_mrmr_lda_index}, 'top5_training_feats.mat')).train_feats;
                    axial_tumor_feats_test1 = load(fullfile(results_root, P{axial_tumor_mrmr_lda_index}, 'top5_test1_feats.mat')).test1_feats;
                    axial_tumor_feats_test2 = load(fullfile(results_root, P{axial_tumor_mrmr_lda_index}, 'top5_test2_feats.mat')).test2_feats;
                    coronal_fat_feats_train = load(fullfile(results_root, P{coronal_proxfat10_mrmr_lda_index}, 'top5_training_feats.mat')).train_feats;
                    coronal_fat_feats_test1 = load(fullfile(results_root, P{coronal_proxfat10_mrmr_lda_index}, 'top5_test1_feats.mat')).test1_feats;
                    coronal_fat_feats_test2 = load(fullfile(results_root, P{coronal_proxfat10_mrmr_lda_index}, 'top5_test2_feats.mat')).test2_feats;
                    coronal_tumor_feats_train = load(fullfile(results_root, P{coronal_tumor_mrmr_lda_index}, 'top5_training_feats.mat')).train_feats;
                    coronal_tumor_feats_test1 = load(fullfile(results_root, P{coronal_tumor_mrmr_lda_index}, 'top5_test1_feats.mat')).test1_feats;
                    coronal_tumor_feats_test2 = load(fullfile(results_root, P{coronal_tumor_mrmr_lda_index}, 'top5_test2_feats.mat')).test2_feats;
                    axial_fat_feat_names = load(fullfile(results_root, P{axial_proxfat5_mrmr_lda_index}, 'Top_Features_Information.mat')).topResults.names';
                    coronal_fat_feat_names = load(fullfile(results_root, P{coronal_proxfat5_mrmr_lda_index}, 'Top_Features_Information.mat')).topResults.names';
                    axial_tumor_feat_names = load(fullfile(results_root, P{axial_tumor_mrmr_lda_index}, 'Top_Features_Information.mat')).topResults.names';
                    coronal_tumor_feat_names = load(fullfile(results_root, P{coronal_tumor_mrmr_lda_index}, 'Top_Features_Information.mat')).topResults.names';
                    feature_column_names = [axial_fat_feat_names, coronal_fat_feat_names, axial_tumor_feat_names, coronal_tumor_feat_names];
                    features_training = [axial_fat_feats_train, coronal_fat_feats_train, axial_tumor_feats_train, coronal_tumor_feats_train];
                    features_holdout1 = [axial_fat_feats_test1, coronal_fat_feats_test1, axial_tumor_feats_test1, coronal_tumor_feats_test1];
                    features_holdout2 = [axial_fat_feats_test2, coronal_fat_feats_test2, axial_tumor_feats_test2, coronal_tumor_feats_test2];
                elseif strcmp(best_fat_roi{1}, 'Fat')
                    axial_fat_feats_train = load(fullfile(results_root, P{axial_fat_mrmr_lda_index}, 'top5_training_feats.mat')).train_feats;
                    axial_fat_feats_test1 = load(fullfile(results_root, P{axial_fat_mrmr_lda_index}, 'top5_test1_feats.mat')).test1_feats;
                    axial_fat_feats_test2 = load(fullfile(results_root, P{axial_fat_mrmr_lda_index}, 'top5_test2_feats.mat')).test2_feats;
                    axial_tumor_feats_train = load(fullfile(results_root, P{axial_tumor_mrmr_lda_index}, 'top5_training_feats.mat')).train_feats;
                    axial_tumor_feats_test1 = load(fullfile(results_root, P{axial_tumor_mrmr_lda_index}, 'top5_test1_feats.mat')).test1_feats;
                    axial_tumor_feats_test2 = load(fullfile(results_root, P{axial_tumor_mrmr_lda_index}, 'top5_test2_feats.mat')).test2_feats;
                    coronal_fat_feats_train = load(fullfile(results_root, P{coronal_fat_mrmr_lda_index}, 'top5_training_feats.mat')).train_feats;
                    coronal_fat_feats_test1 = load(fullfile(results_root, P{coronal_fat_mrmr_lda_index}, 'top5_test1_feats.mat')).test1_feats;
                    coronal_fat_feats_test2 = load(fullfile(results_root, P{coronal_fat_mrmr_lda_index}, 'top5_test2_feats.mat')).test2_feats;
                    coronal_tumor_feats_train = load(fullfile(results_root, P{coronal_tumor_mrmr_lda_index}, 'top5_training_feats.mat')).train_feats;
                    coronal_tumor_feats_test1 = load(fullfile(results_root, P{coronal_tumor_mrmr_lda_index}, 'top5_test1_feats.mat')).test1_feats;
                    coronal_tumor_feats_test2 = load(fullfile(results_root, P{coronal_tumor_mrmr_lda_index}, 'top5_test2_feats.mat')).test2_feats;
                    axial_fat_feat_names = load(fullfile(results_root, P{axial_fat_mrmr_lda_index}, 'Top_Features_Information.mat')).topResults.names';
                    coronal_fat_feat_names = load(fullfile(results_root, P{coronal_fat_mrmr_lda_index}, 'Top_Features_Information.mat')).topResults.names';
                    axial_tumor_feat_names = load(fullfile(results_root, P{axial_tumor_mrmr_lda_index}, 'Top_Features_Information.mat')).topResults.names';
                    coronal_tumor_feat_names = load(fullfile(results_root, P{coronal_tumor_mrmr_lda_index}, 'Top_Features_Information.mat')).topResults.names';
                    feature_column_names = [axial_fat_feat_names, coronal_fat_feat_names, axial_tumor_feat_names, coronal_tumor_feat_names];
                    features_training = [axial_fat_feats_train, coronal_fat_feats_train, axial_tumor_feats_train, coronal_tumor_feats_train];
                    features_holdout1 = [axial_fat_feats_test1, coronal_fat_feats_test1, axial_tumor_feats_test1, coronal_tumor_feats_test1];
                    features_holdout2 = [axial_fat_feats_test2, coronal_fat_feats_test2, axial_tumor_feats_test2, coronal_tumor_feats_test2];
                end
            elseif length(rois) == 1
                if strcmp(rois{1}, 'Proximal_Fat5')
                    axial_feats_train = load(fullfile(results_root, P{axial_proxfat5_mrmr_lda_index}, 'top5_training_feats.mat')).train_feats;
                    axial_feats_test1 = load(fullfile(results_root, P{axial_proxfat5_mrmr_lda_index}, 'top5_test1_feats.mat')).test1_feats;
                    axial_feats_test2 = load(fullfile(results_root, P{axial_proxfat5_mrmr_lda_index}, 'top5_test2_feats.mat')).test2_feats;
                    coronal_feats_train = load(fullfile(results_root, P{coronal_proxfat5_mrmr_lda_index}, 'top5_training_feats.mat')).train_feats;
                    coronal_feats_test1 = load(fullfile(results_root, P{coronal_proxfat5_mrmr_lda_index}, 'top5_test1_feats.mat')).test1_feats;
                    coronal_feats_test2 = load(fullfile(results_root, P{coronal_proxfat5_mrmr_lda_index}, 'top5_test2_feats.mat')).test2_feats;
                    axial_feat_names = load(fullfile(results_root, P{axial_proxfat5_mrmr_lda_index}, 'Top_Features_Information.mat')).topResults.names';
                    coronal_feat_names = load(fullfile(results_root, P{coronal_proxfat5_mrmr_lda_index}, 'Top_Features_Information.mat')).topResults.names';
                    feature_column_names = [axial_feat_names, coronal_feat_names];
                    features_training = [axial_feats_train, coronal_feats_train];
                    features_holdout1 = [axial_feats_test1, coronal_feats_test1];
                    features_holdout2 = [axial_feats_test2, coronal_feats_test2];
                elseif strcmp(rois{1}, 'Proximal_Fat10')
                    axial_feats_train = load(fullfile(results_root, P{axial_proxfat10_mrmr_lda_index}, 'top5_training_feats.mat')).train_feats;
                    axial_feats_test1 = load(fullfile(results_root, P{axial_proxfat10_mrmr_lda_index}, 'top5_test1_feats.mat')).test1_feats;
                    axial_feats_test2 = load(fullfile(results_root, P{axial_proxfat10_mrmr_lda_index}, 'top5_test2_feats.mat')).test2_feats;
                    coronal_feats_train = load(fullfile(results_root, P{coronal_proxfat10_mrmr_lda_index}, 'top5_training_feats.mat')).train_feats;
                    coronal_feats_test1 = load(fullfile(results_root, P{coronal_proxfat10_mrmr_lda_index}, 'top5_test1_feats.mat')).test1_feats;
                    coronal_feats_test2 = load(fullfile(results_root, P{coronal_proxfat10_mrmr_lda_index}, 'top5_test2_feats.mat')).test2_feats;
                    axial_feat_names = load(fullfile(results_root, P{axial_proxfat10_mrmr_lda_index}, 'Top_Features_Information.mat')).topResults.names';
                    coronal_feat_names = load(fullfile(results_root, P{coronal_proxfat10_mrmr_lda_index}, 'Top_Features_Information.mat')).topResults.names';
                    feature_column_names = [axial_feat_names, coronal_feat_names];
                    features_training = [axial_feats_train, coronal_feats_train];
                    features_holdout1 = [axial_feats_test1, coronal_feats_test1];
                    features_holdout2 = [axial_feats_test2, coronal_feats_test2];
                elseif strcmp(rois{1}, 'Fat')
                    axial_feats_train = load(fullfile(results_root, P{axial_fat_mrmr_lda_index}, 'top5_training_feats.mat')).train_feats;
                    axial_feats_test1 = load(fullfile(results_root, P{axial_fat_mrmr_lda_index}, 'top5_test1_feats.mat')).test1_feats;
                    axial_feats_test2 = load(fullfile(results_root, P{axial_fat_mrmr_lda_index}, 'top5_test2_feats.mat')).test2_feats;
                    coronal_feats_train = load(fullfile(results_root, P{coronal_fat_mrmr_lda_index}, 'top5_training_feats.mat')).train_feats;
                    coronal_feats_test1 = load(fullfile(results_root, P{coronal_fat_mrmr_lda_index}, 'top5_test1_feats.mat')).test1_feats;
                    coronal_feats_test2 = load(fullfile(results_root, P{coronal_fat_mrmr_lda_index}, 'top5_test2_feats.mat')).test2_feats;
                    axial_feat_names = load(fullfile(results_root, P{axial_fat_mrmr_lda_index}, 'Top_Features_Information.mat')).topResults.names';
                    coronal_feat_names = load(fullfile(results_root, P{coronal_fat_mrmr_lda_index}, 'Top_Features_Information.mat')).topResults.names';
                    feature_column_names = [axial_feat_names, coronal_feat_names];
                    features_training = [axial_feats_train, coronal_feats_train];
                    features_holdout1 = [axial_feats_test1, coronal_feats_test1];
                    features_holdout2 = [axial_feats_test2, coronal_feats_test2];
                elseif strcmp(rois{1}, 'Tumor')
                    axial_feats_train = load(fullfile(results_root, P{axial_tumor_mrmr_lda_index}, 'top5_training_feats.mat')).train_feats;
                    axial_feats_test1 = load(fullfile(results_root, P{axial_tumor_mrmr_lda_index}, 'top5_test1_feats.mat')).test1_feats;
                    axial_feats_test2 = load(fullfile(results_root, P{axial_tumor_mrmr_lda_index}, 'top5_test2_feats.mat')).test2_feats;
                    coronal_feats_train = load(fullfile(results_root, P{coronal_tumor_mrmr_lda_index}, 'top5_training_feats.mat')).train_feats;
                    coronal_feats_test1 = load(fullfile(results_root, P{coronal_tumor_mrmr_lda_index}, 'top5_test1_feats.mat')).test1_feats;
                    coronal_feats_test2 = load(fullfile(results_root, P{coronal_tumor_mrmr_lda_index}, 'top5_test2_feats.mat')).test2_feats;
                    axial_feat_names = load(fullfile(results_root, P{axial_tumor_mrmr_lda_index}, 'Top_Features_Information.mat')).topResults.names';
                    coronal_feat_names = load(fullfile(results_root, P{coronal_tumor_mrmr_lda_index}, 'Top_Features_Information.mat')).topResults.names';
                    feature_column_names = [axial_feat_names, coronal_feat_names];
                    features_training = [axial_feats_train, coronal_feats_train];
                    features_holdout1 = [axial_feats_test1, coronal_feats_test1];
                    features_holdout2 = [axial_feats_test2, coronal_feats_test2];
                end
            end
        elseif strcmp(plane{1}, 'Axial')
            if strcmp(best_fat_roi{1}, 'Proximal_Fat5')
                axial_fat_feats_train = load(fullfile(results_root, P{axial_proxfat5_mrmr_lda_index}, 'top5_training_feats.mat')).train_feats;
                axial_fat_feats_test1 = load(fullfile(results_root, P{axial_proxfat5_mrmr_lda_index}, 'top5_test1_feats.mat')).test1_feats;
                axial_fat_feats_test2 = load(fullfile(results_root, P{axial_proxfat5_mrmr_lda_index}, 'top5_test2_feats.mat')).test2_feats;
                axial_tumor_feats_train = load(fullfile(results_root, P{axial_tumor_mrmr_lda_index}, 'top5_training_feats.mat')).train_feats;
                axial_tumor_feats_test1 = load(fullfile(results_root, P{axial_tumor_mrmr_lda_index}, 'top5_test1_feats.mat')).test1_feats;
                axial_tumor_feats_test2 = load(fullfile(results_root, P{axial_tumor_mrmr_lda_index}, 'top5_test2_feats.mat')).test2_feats;
                axial_fat_feat_names = load(fullfile(results_root, P{axial_proxfat5_mrmr_lda_index}, 'Top_Features_Information.mat')).topResults.names';
                axial_tumor_feat_names = load(fullfile(results_root, P{axial_tumor_mrmr_lda_index}, 'Top_Features_Information.mat')).topResults.names';
                feature_column_names = [axial_fat_feat_names, axial_tumor_feat_names];
                features_training = [axial_fat_feats_train, axial_tumor_feats_train];
                features_holdout1 = [axial_fat_feats_test1, axial_tumor_feats_test1];
                features_holdout2 = [axial_fat_feats_test2, axial_tumor_feats_test2];
            elseif strcmp(best_fat_roi{1}, 'Proximal_Fat10')
                axial_fat_feats_train = load(fullfile(results_root, P{axial_proxfat10_mrmr_lda_index}, 'top5_training_feats.mat')).train_feats;
                axial_fat_feats_test1 = load(fullfile(results_root, P{axial_proxfat10_mrmr_lda_index}, 'top5_test1_feats.mat')).test1_feats;
                axial_fat_feats_test2 = load(fullfile(results_root, P{axial_proxfat10_mrmr_lda_index}, 'top5_test2_feats.mat')).test2_feats;
                axial_tumor_feats_train = load(fullfile(results_root, P{axial_tumor_mrmr_lda_index}, 'top5_training_feats.mat')).train_feats;
                axial_tumor_feats_test1 = load(fullfile(results_root, P{axial_tumor_mrmr_lda_index}, 'top5_test1_feats.mat')).test1_feats;
                axial_tumor_feats_test2 = load(fullfile(results_root, P{axial_tumor_mrmr_lda_index}, 'top5_test2_feats.mat')).test2_feats;
                axial_fat_feat_names = load(fullfile(results_root, P{axial_proxfat5_mrmr_lda_index}, 'Top_Features_Information.mat')).topResults.names';
                axial_tumor_feat_names = load(fullfile(results_root, P{axial_tumor_mrmr_lda_index}, 'Top_Features_Information.mat')).topResults.names';
                feature_column_names = [axial_fat_feat_names, axial_tumor_feat_names];
                features_training = [axial_fat_feats_train, axial_tumor_feats_train];
                features_holdout1 = [axial_fat_feats_test1, axial_tumor_feats_test1];
                features_holdout2 = [axial_fat_feats_test2, axial_tumor_feats_test2];
            elseif strcmp(best_fat_roi{1}, 'Fat')
                axial_fat_feats_train = load(fullfile(results_root, P{axial_fat_mrmr_lda_index}, 'top5_training_feats.mat')).train_feats;
                axial_fat_feats_test1 = load(fullfile(results_root, P{axial_fat_mrmr_lda_index}, 'top5_test1_feats.mat')).test1_feats;
                axial_fat_feats_test2 = load(fullfile(results_root, P{axial_fat_mrmr_lda_index}, 'top5_test2_feats.mat')).test2_feats;
                axial_tumor_feats_train = load(fullfile(results_root, P{axial_tumor_mrmr_lda_index}, 'top5_training_feats.mat')).train_feats;
                axial_tumor_feats_test1 = load(fullfile(results_root, P{axial_tumor_mrmr_lda_index}, 'top5_test1_feats.mat')).test1_feats;
                axial_tumor_feats_test2 = load(fullfile(results_root, P{axial_tumor_mrmr_lda_index}, 'top5_test2_feats.mat')).test2_feats;
                axial_fat_feat_names = load(fullfile(results_root, P{axial_fat_mrmr_lda_index}, 'Top_Features_Information.mat')).topResults.names';
                axial_tumor_feat_names = load(fullfile(results_root, P{axial_tumor_mrmr_lda_index}, 'Top_Features_Information.mat')).topResults.names';
                feature_column_names = [axial_fat_feat_names, axial_tumor_feat_names];
                features_training = [axial_fat_feats_train, axial_tumor_feats_train];
                features_holdout1 = [axial_fat_feats_test1, axial_tumor_feats_test1];
                features_holdout2 = [axial_fat_feats_test2, axial_tumor_feats_test2];
            end
        elseif strcmp(plane{1}, 'Coronal')
            if strcmp(best_fat_roi{1}, 'Proximal_Fat5')
                coronal_fat_feats_train = load(fullfile(results_root, P{coronal_proxfat5_mrmr_lda_index}, 'top5_training_feats.mat')).train_feats;
                coronal_fat_feats_test1 = load(fullfile(results_root, P{coronal_proxfat5_mrmr_lda_index}, 'top5_test1_feats.mat')).test1_feats;
                coronal_fat_feats_test2 = load(fullfile(results_root, P{coronal_proxfat5_mrmr_lda_index}, 'top5_test2_feats.mat')).test2_feats;
                coronal_tumor_feats_train = load(fullfile(results_root, P{coronal_tumor_mrmr_lda_index}, 'top5_training_feats.mat')).train_feats;
                coronal_tumor_feats_test1 = load(fullfile(results_root, P{coronal_tumor_mrmr_lda_index}, 'top5_test1_feats.mat')).test1_feats;
                coronal_tumor_feats_test2 = load(fullfile(results_root, P{coronal_tumor_mrmr_lda_index}, 'top5_test2_feats.mat')).test2_feats;
                coronal_fat_feat_names = load(fullfile(results_root, P{coronal_proxfat5_mrmr_lda_index}, 'Top_Features_Information.mat')).topResults.names';
                coronal_tumor_feat_names = load(fullfile(results_root, P{coronal_tumor_mrmr_lda_index}, 'Top_Features_Information.mat')).topResults.names';
                feature_column_names = [coronal_fat_feat_names, coronal_tumor_feat_names];
                features_training = [coronal_fat_feats_train, coronal_tumor_feats_train];
                features_holdout1 = [coronal_fat_feats_test1, coronal_tumor_feats_test1];
                features_holdout2 = [coronal_fat_feats_test2, coronal_tumor_feats_test2];
            elseif strcmp(best_fat_roi{1}, 'Proximal_Fat10')
                coronal_fat_feats_train = load(fullfile(results_root, P{coronal_proxfat10_mrmr_lda_index}, 'top5_training_feats.mat')).train_feats;
                coronal_fat_feats_test1 = load(fullfile(results_root, P{coronal_proxfat10_mrmr_lda_index}, 'top5_test1_feats.mat')).test1_feats;
                coronal_fat_feats_test2 = load(fullfile(results_root, P{coronal_proxfat10_mrmr_lda_index}, 'top5_test2_feats.mat')).test2_feats;
                coronal_tumor_feats_train = load(fullfile(results_root, P{coronal_tumor_mrmr_lda_index}, 'top5_training_feats.mat')).train_feats;
                coronal_tumor_feats_test1 = load(fullfile(results_root, P{coronal_tumor_mrmr_lda_index}, 'top5_test1_feats.mat')).test1_feats;
                coronal_tumor_feats_test2 = load(fullfile(results_root, P{coronal_tumor_mrmr_lda_index}, 'top5_test2_feats.mat')).test2_feats;
                coronal_fat_feat_names = load(fullfile(results_root, P{coronal_proxfat5_mrmr_lda_index}, 'Top_Features_Information.mat')).topResults.names';
                coronal_tumor_feat_names = load(fullfile(results_root, P{coronal_tumor_mrmr_lda_index}, 'Top_Features_Information.mat')).topResults.names';
                feature_column_names = [coronal_fat_feat_names, coronal_tumor_feat_names];
                features_training = [coronal_fat_feats_train, coronal_tumor_feats_train];
                features_holdout1 = [coronal_fat_feats_test1, coronal_tumor_feats_test1];
                features_holdout2 = [coronal_fat_feats_test2, coronal_tumor_feats_test2];
            elseif strcmp(best_fat_roi{1}, 'Fat')
                coronal_fat_feats_train = load(fullfile(results_root, P{coronal_fat_mrmr_lda_index}, 'top5_training_feats.mat')).train_feats;
                coronal_fat_feats_test1 = load(fullfile(results_root, P{coronal_fat_mrmr_lda_index}, 'top5_test1_feats.mat')).test1_feats;
                coronal_fat_feats_test2 = load(fullfile(results_root, P{coronal_fat_mrmr_lda_index}, 'top5_test2_feats.mat')).test2_feats;
                coronal_tumor_feats_train = load(fullfile(results_root, P{coronal_tumor_mrmr_lda_index}, 'top5_training_feats.mat')).train_feats;
                coronal_tumor_feats_test1 = load(fullfile(results_root, P{coronal_tumor_mrmr_lda_index}, 'top5_test1_feats.mat')).test1_feats;
                coronal_tumor_feats_test2 = load(fullfile(results_root, P{coronal_tumor_mrmr_lda_index}, 'top5_test2_feats.mat')).test2_feats;
                coronal_fat_feat_names = load(fullfile(results_root, P{coronal_fat_mrmr_lda_index}, 'Top_Features_Information.mat')).topResults.names';
                coronal_tumor_feat_names = load(fullfile(results_root, P{coronal_tumor_mrmr_lda_index}, 'Top_Features_Information.mat')).topResults.names';
                feature_column_names = [coronal_fat_feat_names, coronal_tumor_feat_names];
                features_training = [coronal_fat_feats_train, coronal_tumor_feats_train];
                features_holdout1 = [coronal_fat_feats_test1, coronal_tumor_feats_test1];
                features_holdout2 = [coronal_fat_feats_test2, coronal_tumor_feats_test2];
            end
        end
    elseif strcmp(classifier(1), 'rf')
        if strcmp(plane{1}, 'Multi-Plane')
            if length(rois) == 2
                if strcmp(best_fat_roi{1}, 'Proximal_Fat5')
                    axial_fat_feats_train = load(fullfile(results_root, P{axial_proxfat5_mrmr_rf_index}, 'top5_training_feats.mat')).train_feats;
                    axial_fat_feats_test1 = load(fullfile(results_root, P{axial_proxfat5_mrmr_rf_index}, 'top5_test1_feats.mat')).test1_feats;
                    axial_fat_feats_test2 = load(fullfile(results_root, P{axial_proxfat5_mrmr_rf_index}, 'top5_test2_feats.mat')).test2_feats;
                    axial_tumor_feats_train = load(fullfile(results_root, P{axial_tumor_mrmr_rf_index}, 'top5_training_feats.mat')).train_feats;
                    axial_tumor_feats_test1 = load(fullfile(results_root, P{axial_tumor_mrmr_rf_index}, 'top5_test1_feats.mat')).test1_feats;
                    axial_tumor_feats_test2 = load(fullfile(results_root, P{axial_tumor_mrmr_rf_index}, 'top5_test2_feats.mat')).test2_feats;
                    coronal_fat_feats_train = load(fullfile(results_root, P{coronal_proxfat5_mrmr_rf_index}, 'top5_training_feats.mat')).train_feats;
                    coronal_fat_feats_test1 = load(fullfile(results_root, P{coronal_proxfat5_mrmr_rf_index}, 'top5_test1_feats.mat')).test1_feats;
                    coronal_fat_feats_test2 = load(fullfile(results_root, P{coronal_proxfat5_mrmr_rf_index}, 'top5_test2_feats.mat')).test2_feats;
                    coronal_tumor_feats_train = load(fullfile(results_root, P{coronal_tumor_mrmr_rf_index}, 'top5_training_feats.mat')).train_feats;
                    coronal_tumor_feats_test1 = load(fullfile(results_root, P{coronal_tumor_mrmr_rf_index}, 'top5_test1_feats.mat')).test1_feats;
                    coronal_tumor_feats_test2 = load(fullfile(results_root, P{coronal_tumor_mrmr_rf_index}, 'top5_test2_feats.mat')).test2_feats;
                    axial_fat_feat_names = load(fullfile(results_root, P{axial_proxfat5_mrmr_rf_index}, 'Top_Features_Information.mat')).topResults.names';
                    coronal_fat_feat_names = load(fullfile(results_root, P{coronal_proxfat5_mrmr_rf_index}, 'Top_Features_Information.mat')).topResults.names';
                    axial_tumor_feat_names = load(fullfile(results_root, P{axial_tumor_mrmr_rf_index}, 'Top_Features_Information.mat')).topResults.names';
                    coronal_tumor_feat_names = load(fullfile(results_root, P{coronal_tumor_mrmr_rf_index}, 'Top_Features_Information.mat')).topResults.names';
                    feature_column_names = [axial_fat_feat_names, coronal_fat_feat_names, axial_tumor_feat_names, coronal_tumor_feat_names];
                    features_training = [axial_fat_feats_train, coronal_fat_feats_train, axial_tumor_feats_train, coronal_tumor_feats_train];
                    features_holdout1 = [axial_fat_feats_test1, coronal_fat_feats_test1, axial_tumor_feats_test1, coronal_tumor_feats_test1];
                    features_holdout2 = [axial_fat_feats_test2, coronal_fat_feats_test2, axial_tumor_feats_test2, coronal_tumor_feats_test2];
                elseif strcmp(best_fat_roi{1}, 'Proximal_Fat10')
                    axial_fat_feats_train = load(fullfile(results_root, P{axial_proxfat10_mrmr_rf_index}, 'top5_training_feats.mat')).train_feats;
                    axial_fat_feats_test1 = load(fullfile(results_root, P{axial_proxfat10_mrmr_rf_index}, 'top5_test1_feats.mat')).test1_feats;
                    axial_fat_feats_test2 = load(fullfile(results_root, P{axial_proxfat10_mrmr_rf_index}, 'top5_test2_feats.mat')).test2_feats;
                    axial_tumor_feats_train = load(fullfile(results_root, P{axial_tumor_mrmr_rf_index}, 'top5_training_feats.mat')).train_feats;
                    axial_tumor_feats_test1 = load(fullfile(results_root, P{axial_tumor_mrmr_rf_index}, 'top5_test1_feats.mat')).test1_feats;
                    axial_tumor_feats_test2 = load(fullfile(results_root, P{axial_tumor_mrmr_rf_index}, 'top5_test2_feats.mat')).test2_feats;
                    coronal_fat_feats_train = load(fullfile(results_root, P{coronal_proxfat10_mrmr_rf_index}, 'top5_training_feats.mat')).train_feats;
                    coronal_fat_feats_test1 = load(fullfile(results_root, P{coronal_proxfat10_mrmr_rf_index}, 'top5_test1_feats.mat')).test1_feats;
                    coronal_fat_feats_test2 = load(fullfile(results_root, P{coronal_proxfat10_mrmr_rf_index}, 'top5_test2_feats.mat')).test2_feats;
                    coronal_tumor_feats_train = load(fullfile(results_root, P{coronal_tumor_mrmr_rf_index}, 'top5_training_feats.mat')).train_feats;
                    coronal_tumor_feats_test1 = load(fullfile(results_root, P{coronal_tumor_mrmr_rf_index}, 'top5_test1_feats.mat')).test1_feats;
                    coronal_tumor_feats_test2 = load(fullfile(results_root, P{coronal_tumor_mrmr_rf_index}, 'top5_test2_feats.mat')).test2_feats;
                    axial_fat_feat_names = load(fullfile(results_root, P{axial_proxfat5_mrmr_rf_index}, 'Top_Features_Information.mat')).topResults.names';
                    coronal_fat_feat_names = load(fullfile(results_root, P{coronal_proxfat5_mrmr_rf_index}, 'Top_Features_Information.mat')).topResults.names';
                    axial_tumor_feat_names = load(fullfile(results_root, P{axial_tumor_mrmr_rf_index}, 'Top_Features_Information.mat')).topResults.names';
                    coronal_tumor_feat_names = load(fullfile(results_root, P{coronal_tumor_mrmr_rf_index}, 'Top_Features_Information.mat')).topResults.names';
                    feature_column_names = [axial_fat_feat_names, coronal_fat_feat_names, axial_tumor_feat_names, coronal_tumor_feat_names];
                    features_training = [axial_fat_feats_train, coronal_fat_feats_train, axial_tumor_feats_train, coronal_tumor_feats_train];
                    features_holdout1 = [axial_fat_feats_test1, coronal_fat_feats_test1, axial_tumor_feats_test1, coronal_tumor_feats_test1];
                    features_holdout2 = [axial_fat_feats_test2, coronal_fat_feats_test2, axial_tumor_feats_test2, coronal_tumor_feats_test2];
                elseif strcmp(best_fat_roi{1}, 'Fat')
                    axial_fat_feats_train = load(fullfile(results_root, P{axial_fat_mrmr_rf_index}, 'top5_training_feats.mat')).train_feats;
                    axial_fat_feats_test1 = load(fullfile(results_root, P{axial_fat_mrmr_rf_index}, 'top5_test1_feats.mat')).test1_feats;
                    axial_fat_feats_test2 = load(fullfile(results_root, P{axial_fat_mrmr_rf_index}, 'top5_test2_feats.mat')).test2_feats;
                    axial_tumor_feats_train = load(fullfile(results_root, P{axial_tumor_mrmr_rf_index}, 'top5_training_feats.mat')).train_feats;
                    axial_tumor_feats_test1 = load(fullfile(results_root, P{axial_tumor_mrmr_rf_index}, 'top5_test1_feats.mat')).test1_feats;
                    axial_tumor_feats_test2 = load(fullfile(results_root, P{axial_tumor_mrmr_rf_index}, 'top5_test2_feats.mat')).test2_feats;
                    coronal_fat_feats_train = load(fullfile(results_root, P{coronal_fat_mrmr_rf_index}, 'top5_training_feats.mat')).train_feats;
                    coronal_fat_feats_test1 = load(fullfile(results_root, P{coronal_fat_mrmr_rf_index}, 'top5_test1_feats.mat')).test1_feats;
                    coronal_fat_feats_test2 = load(fullfile(results_root, P{coronal_fat_mrmr_rf_index}, 'top5_test2_feats.mat')).test2_feats;
                    coronal_tumor_feats_train = load(fullfile(results_root, P{coronal_tumor_mrmr_rf_index}, 'top5_training_feats.mat')).train_feats;
                    coronal_tumor_feats_test1 = load(fullfile(results_root, P{coronal_tumor_mrmr_rf_index}, 'top5_test1_feats.mat')).test1_feats;
                    coronal_tumor_feats_test2 = load(fullfile(results_root, P{coronal_tumor_mrmr_rf_index}, 'top5_test2_feats.mat')).test2_feats;
                    axial_fat_feat_names = load(fullfile(results_root, P{axial_fat_mrmr_rf_index}, 'Top_Features_Information.mat')).topResults.names';
                    coronal_fat_feat_names = load(fullfile(results_root, P{coronal_fat_mrmr_rf_index}, 'Top_Features_Information.mat')).topResults.names';
                    axial_tumor_feat_names = load(fullfile(results_root, P{axial_tumor_mrmr_rf_index}, 'Top_Features_Information.mat')).topResults.names';
                    coronal_tumor_feat_names = load(fullfile(results_root, P{coronal_tumor_mrmr_rf_index}, 'Top_Features_Information.mat')).topResults.names';
                    feature_column_names = [axial_fat_feat_names, coronal_fat_feat_names, axial_tumor_feat_names, coronal_tumor_feat_names];
                    features_training = [axial_fat_feats_train, coronal_fat_feats_train, axial_tumor_feats_train, coronal_tumor_feats_train];
                    features_holdout1 = [axial_fat_feats_test1, coronal_fat_feats_test1, axial_tumor_feats_test1, coronal_tumor_feats_test1];
                    features_holdout2 = [axial_fat_feats_test2, coronal_fat_feats_test2, axial_tumor_feats_test2, coronal_tumor_feats_test2];
                end
            elseif length(rois) == 1
                if strcmp(rois{1}, 'Proximal_Fat5')
                    axial_feats_train = load(fullfile(results_root, P{axial_proxfat5_mrmr_rf_index}, 'top5_training_feats.mat')).train_feats;
                    axial_feats_test1 = load(fullfile(results_root, P{axial_proxfat5_mrmr_rf_index}, 'top5_test1_feats.mat')).test1_feats;
                    axial_feats_test2 = load(fullfile(results_root, P{axial_proxfat5_mrmr_rf_index}, 'top5_test2_feats.mat')).test2_feats;
                    coronal_feats_train = load(fullfile(results_root, P{coronal_proxfat5_mrmr_rf_index}, 'top5_training_feats.mat')).train_feats;
                    coronal_feats_test1 = load(fullfile(results_root, P{coronal_proxfat5_mrmr_rf_index}, 'top5_test1_feats.mat')).test1_feats;
                    coronal_feats_test2 = load(fullfile(results_root, P{coronal_proxfat5_mrmr_rf_index}, 'top5_test2_feats.mat')).test2_feats;
                    axial_feat_names = load(fullfile(results_root, P{axial_proxfat5_mrmr_rf_index}, 'Top_Features_Information.mat')).topResults.names';
                    coronal_feat_names = load(fullfile(results_root, P{coronal_proxfat5_mrmr_rf_index}, 'Top_Features_Information.mat')).topResults.names';
                    feature_column_names = [axial_feat_names, coronal_feat_names];
                    features_training = [axial_feats_train, coronal_feats_train];
                    features_holdout1 = [axial_feats_test1, coronal_feats_test1];
                    features_holdout2 = [axial_feats_test2, coronal_feats_test2];
                elseif strcmp(rois{1}, 'Proximal_Fat10')
                    axial_feats_train = load(fullfile(results_root, P{axial_proxfat10_mrmr_rf_index}, 'top5_training_feats.mat')).train_feats;
                    axial_feats_test1 = load(fullfile(results_root, P{axial_proxfat10_mrmr_rf_index}, 'top5_test1_feats.mat')).test1_feats;
                    axial_feats_test2 = load(fullfile(results_root, P{axial_proxfat10_mrmr_rf_index}, 'top5_test2_feats.mat')).test2_feats;
                    coronal_feats_train = load(fullfile(results_root, P{coronal_proxfat10_mrmr_rf_index}, 'top5_training_feats.mat')).train_feats;
                    coronal_feats_test1 = load(fullfile(results_root, P{coronal_proxfat10_mrmr_rf_index}, 'top5_test1_feats.mat')).test1_feats;
                    coronal_feats_test2 = load(fullfile(results_root, P{coronal_proxfat10_mrmr_rf_index}, 'top5_test2_feats.mat')).test2_feats;
                    axial_feat_names = load(fullfile(results_root, P{axial_proxfat10_mrmr_rf_index}, 'Top_Features_Information.mat')).topResults.names';
                    coronal_feat_names = load(fullfile(results_root, P{coronal_proxfat10_mrmr_rf_index}, 'Top_Features_Information.mat')).topResults.names';
                    feature_column_names = [axial_feat_names, coronal_feat_names];
                    features_training = [axial_feats_train, coronal_feats_train];
                    features_holdout1 = [axial_feats_test1, coronal_feats_test1];
                    features_holdout2 = [axial_feats_test2, coronal_feats_test2];
                elseif strcmp(rois{1}, 'Fat')
                    axial_feats_train = load(fullfile(results_root, P{axial_fat_mrmr_rf_index}, 'top5_training_feats.mat')).train_feats;
                    axial_feats_test1 = load(fullfile(results_root, P{axial_fat_mrmr_rf_index}, 'top5_test1_feats.mat')).test1_feats;
                    axial_feats_test2 = load(fullfile(results_root, P{axial_fat_mrmr_rf_index}, 'top5_test2_feats.mat')).test2_feats;
                    coronal_feats_train = load(fullfile(results_root, P{coronal_fat_mrmr_rf_index}, 'top5_training_feats.mat')).train_feats;
                    coronal_feats_test1 = load(fullfile(results_root, P{coronal_fat_mrmr_rf_index}, 'top5_test1_feats.mat')).test1_feats;
                    coronal_feats_test2 = load(fullfile(results_root, P{coronal_fat_mrmr_rf_index}, 'top5_test2_feats.mat')).test2_feats;
                    axial_feat_names = load(fullfile(results_root, P{axial_fat_mrmr_rf_index}, 'Top_Features_Information.mat')).topResults.names';
                    coronal_feat_names = load(fullfile(results_root, P{coronal_fat_mrmr_rf_index}, 'Top_Features_Information.mat')).topResults.names';
                    feature_column_names = [axial_feat_names, coronal_feat_names];
                    features_training = [axial_feats_train, coronal_feats_train];
                    features_holdout1 = [axial_feats_test1, coronal_feats_test1];
                    features_holdout2 = [axial_feats_test2, coronal_feats_test2];
                elseif strcmp(rois{1}, 'Tumor')
                    axial_feats_train = load(fullfile(results_root, P{axial_tumor_mrmr_rf_index}, 'top5_training_feats.mat')).train_feats;
                    axial_feats_test1 = load(fullfile(results_root, P{axial_tumor_mrmr_rf_index}, 'top5_test1_feats.mat')).test1_feats;
                    axial_feats_test2 = load(fullfile(results_root, P{axial_tumor_mrmr_rf_index}, 'top5_test2_feats.mat')).test2_feats;
                    coronal_feats_train = load(fullfile(results_root, P{coronal_tumor_mrmr_rf_index}, 'top5_training_feats.mat')).train_feats;
                    coronal_feats_test1 = load(fullfile(results_root, P{coronal_tumor_mrmr_rf_index}, 'top5_test1_feats.mat')).test1_feats;
                    coronal_feats_test2 = load(fullfile(results_root, P{coronal_tumor_mrmr_rf_index}, 'top5_test2_feats.mat')).test2_feats;
                    axial_feat_names = load(fullfile(results_root, P{axial_tumor_mrmr_rf_index}, 'Top_Features_Information.mat')).topResults.names';
                    coronal_feat_names = load(fullfile(results_root, P{coronal_tumor_mrmr_rf_index}, 'Top_Features_Information.mat')).topResults.names';
                    feature_column_names = [axial_feat_names, coronal_feat_names];
                    features_training = [axial_feats_train, coronal_feats_train];
                    features_holdout1 = [axial_feats_test1, coronal_feats_test1];
                    features_holdout2 = [axial_feats_test2, coronal_feats_test2];
                end
            end
        elseif strcmp(plane{1}, 'Axial')
            if strcmp(best_fat_roi{1}, 'Proximal_Fat5')
                axial_fat_feats_train = load(fullfile(results_root, P{axial_proxfat5_mrmr_rf_index}, 'top5_training_feats.mat')).train_feats;
                axial_fat_feats_test1 = load(fullfile(results_root, P{axial_proxfat5_mrmr_rf_index}, 'top5_test1_feats.mat')).test1_feats;
                axial_fat_feats_test2 = load(fullfile(results_root, P{axial_proxfat5_mrmr_rf_index}, 'top5_test2_feats.mat')).test2_feats;
                axial_tumor_feats_train = load(fullfile(results_root, P{axial_tumor_mrmr_rf_index}, 'top5_training_feats.mat')).train_feats;
                axial_tumor_feats_test1 = load(fullfile(results_root, P{axial_tumor_mrmr_rf_index}, 'top5_test1_feats.mat')).test1_feats;
                axial_tumor_feats_test2 = load(fullfile(results_root, P{axial_tumor_mrmr_rf_index}, 'top5_test2_feats.mat')).test2_feats;
                axial_fat_feat_names = load(fullfile(results_root, P{axial_proxfat5_mrmr_rf_index}, 'Top_Features_Information.mat')).topResults.names';
                axial_tumor_feat_names = load(fullfile(results_root, P{axial_tumor_mrmr_rf_index}, 'Top_Features_Information.mat')).topResults.names';
                feature_column_names = [axial_fat_feat_names, axial_tumor_feat_names];
                features_training = [axial_fat_feats_train, axial_tumor_feats_train];
                features_holdout1 = [axial_fat_feats_test1, axial_tumor_feats_test1];
                features_holdout2 = [axial_fat_feats_test2, axial_tumor_feats_test2];
            elseif strcmp(best_fat_roi{1}, 'Proximal_Fat10')
                axial_fat_feats_train = load(fullfile(results_root, P{axial_proxfat10_mrmr_rf_index}, 'top5_training_feats.mat')).train_feats;
                axial_fat_feats_test1 = load(fullfile(results_root, P{axial_proxfat10_mrmr_rf_index}, 'top5_test1_feats.mat')).test1_feats;
                axial_fat_feats_test2 = load(fullfile(results_root, P{axial_proxfat10_mrmr_rf_index}, 'top5_test2_feats.mat')).test2_feats;
                axial_tumor_feats_train = load(fullfile(results_root, P{axial_tumor_mrmr_rf_index}, 'top5_training_feats.mat')).train_feats;
                axial_tumor_feats_test1 = load(fullfile(results_root, P{axial_tumor_mrmr_rf_index}, 'top5_test1_feats.mat')).test1_feats;
                axial_tumor_feats_test2 = load(fullfile(results_root, P{axial_tumor_mrmr_rf_index}, 'top5_test2_feats.mat')).test2_feats;
                axial_fat_feat_names = load(fullfile(results_root, P{axial_proxfat5_mrmr_rf_index}, 'Top_Features_Information.mat')).topResults.names';
                axial_tumor_feat_names = load(fullfile(results_root, P{axial_tumor_mrmr_rf_index}, 'Top_Features_Information.mat')).topResults.names';
                feature_column_names = [axial_fat_feat_names, axial_tumor_feat_names];
                features_training = [axial_fat_feats_train, axial_tumor_feats_train];
                features_holdout1 = [axial_fat_feats_test1, axial_tumor_feats_test1];
                features_holdout2 = [axial_fat_feats_test2, axial_tumor_feats_test2];
            elseif strcmp(best_fat_roi{1}, 'Fat')
                axial_fat_feats_train = load(fullfile(results_root, P{axial_fat_mrmr_rf_index}, 'top5_training_feats.mat')).train_feats;
                axial_fat_feats_test1 = load(fullfile(results_root, P{axial_fat_mrmr_rf_index}, 'top5_test1_feats.mat')).test1_feats;
                axial_fat_feats_test2 = load(fullfile(results_root, P{axial_fat_mrmr_rf_index}, 'top5_test2_feats.mat')).test2_feats;
                axial_tumor_feats_train = load(fullfile(results_root, P{axial_tumor_mrmr_rf_index}, 'top5_training_feats.mat')).train_feats;
                axial_tumor_feats_test1 = load(fullfile(results_root, P{axial_tumor_mrmr_rf_index}, 'top5_test1_feats.mat')).test1_feats;
                axial_tumor_feats_test2 = load(fullfile(results_root, P{axial_tumor_mrmr_rf_index}, 'top5_test2_feats.mat')).test2_feats;
                axial_fat_feat_names = load(fullfile(results_root, P{axial_fat_mrmr_rf_index}, 'Top_Features_Information.mat')).topResults.names';
                axial_tumor_feat_names = load(fullfile(results_root, P{axial_tumor_mrmr_rf_index}, 'Top_Features_Information.mat')).topResults.names';
                feature_column_names = [axial_fat_feat_names, axial_tumor_feat_names];
                features_training = [axial_fat_feats_train, axial_tumor_feats_train];
                features_holdout1 = [axial_fat_feats_test1, axial_tumor_feats_test1];
                features_holdout2 = [axial_fat_feats_test2, axial_tumor_feats_test2];
            end
        elseif strcmp(plane{1}, 'Coronal')
            if strcmp(best_fat_roi{1}, 'Proximal_Fat5')
                coronal_fat_feats_train = load(fullfile(results_root, P{coronal_proxfat5_mrmr_rf_index}, 'top5_training_feats.mat')).train_feats;
                coronal_fat_feats_test1 = load(fullfile(results_root, P{coronal_proxfat5_mrmr_rf_index}, 'top5_test1_feats.mat')).test1_feats;
                coronal_fat_feats_test2 = load(fullfile(results_root, P{coronal_proxfat5_mrmr_rf_index}, 'top5_test2_feats.mat')).test2_feats;
                coronal_tumor_feats_train = load(fullfile(results_root, P{coronal_tumor_mrmr_rf_index}, 'top5_training_feats.mat')).train_feats;
                coronal_tumor_feats_test1 = load(fullfile(results_root, P{coronal_tumor_mrmr_rf_index}, 'top5_test1_feats.mat')).test1_feats;
                coronal_tumor_feats_test2 = load(fullfile(results_root, P{coronal_tumor_mrmr_rf_index}, 'top5_test2_feats.mat')).test2_feats;
                coronal_fat_feat_names = load(fullfile(results_root, P{coronal_proxfat5_mrmr_rf_index}, 'Top_Features_Information.mat')).topResults.names';
                coronal_tumor_feat_names = load(fullfile(results_root, P{coronal_tumor_mrmr_rf_index}, 'Top_Features_Information.mat')).topResults.names';
                feature_column_names = [coronal_fat_feat_names, coronal_tumor_feat_names];
                features_training = [coronal_fat_feats_train, coronal_tumor_feats_train];
                features_holdout1 = [coronal_fat_feats_test1, coronal_tumor_feats_test1];
                features_holdout2 = [coronal_fat_feats_test2, coronal_tumor_feats_test2];
            elseif strcmp(best_fat_roi{1}, 'Proximal_Fat10')
                coronal_fat_feats_train = load(fullfile(results_root, P{coronal_proxfat10_mrmr_rf_index}, 'top5_training_feats.mat')).train_feats;
                coronal_fat_feats_test1 = load(fullfile(results_root, P{coronal_proxfat10_mrmr_rf_index}, 'top5_test1_feats.mat')).test1_feats;
                coronal_fat_feats_test2 = load(fullfile(results_root, P{coronal_proxfat10_mrmr_rf_index}, 'top5_test2_feats.mat')).test2_feats;
                coronal_tumor_feats_train = load(fullfile(results_root, P{coronal_tumor_mrmr_rf_index}, 'top5_training_feats.mat')).train_feats;
                coronal_tumor_feats_test1 = load(fullfile(results_root, P{coronal_tumor_mrmr_rf_index}, 'top5_test1_feats.mat')).test1_feats;
                coronal_tumor_feats_test2 = load(fullfile(results_root, P{coronal_tumor_mrmr_rf_index}, 'top5_test2_feats.mat')).test2_feats;
                coronal_fat_feat_names = load(fullfile(results_root, P{coronal_proxfat5_mrmr_rf_index}, 'Top_Features_Information.mat')).topResults.names';
                coronal_tumor_feat_names = load(fullfile(results_root, P{coronal_tumor_mrmr_rf_index}, 'Top_Features_Information.mat')).topResults.names';
                feature_column_names = [coronal_fat_feat_names, coronal_tumor_feat_names];
                features_training = [coronal_fat_feats_train, coronal_tumor_feats_train];
                features_holdout1 = [coronal_fat_feats_test1, coronal_tumor_feats_test1];
                features_holdout2 = [coronal_fat_feats_test2, coronal_tumor_feats_test2];
            elseif strcmp(best_fat_roi{1}, 'Fat')
                coronal_fat_feats_train = load(fullfile(results_root, P{coronal_fat_mrmr_rf_index}, 'top5_training_feats.mat')).train_feats;
                coronal_fat_feats_test1 = load(fullfile(results_root, P{coronal_fat_mrmr_rf_index}, 'top5_test1_feats.mat')).test1_feats;
                coronal_fat_feats_test2 = load(fullfile(results_root, P{coronal_fat_mrmr_rf_index}, 'top5_test2_feats.mat')).test2_feats;
                coronal_tumor_feats_train = load(fullfile(results_root, P{coronal_tumor_mrmr_rf_index}, 'top5_training_feats.mat')).train_feats;
                coronal_tumor_feats_test1 = load(fullfile(results_root, P{coronal_tumor_mrmr_rf_index}, 'top5_test1_feats.mat')).test1_feats;
                coronal_tumor_feats_test2 = load(fullfile(results_root, P{coronal_tumor_mrmr_rf_index}, 'top5_test2_feats.mat')).test2_feats;
                coronal_fat_feat_names = load(fullfile(results_root, P{coronal_fat_mrmr_rf_index}, 'Top_Features_Information.mat')).topResults.names';
                coronal_tumor_feat_names = load(fullfile(results_root, P{coronal_tumor_mrmr_rf_index}, 'Top_Features_Information.mat')).topResults.names';
                feature_column_names = [coronal_fat_feat_names, coronal_tumor_feat_names];
                features_training = [coronal_fat_feats_train, coronal_tumor_feats_train];
                features_holdout1 = [coronal_fat_feats_test1, coronal_tumor_feats_test1];
                features_holdout2 = [coronal_fat_feats_test2, coronal_tumor_feats_test2];
            end
        end
    end
end

fprintf("Loaded in feature matrices! \n");

output_path = "/Users/leobao/Documents/MultiPlanePipeline/AACR2023/Results_Multi/";
experiment_date = strrep(string(datetime("now")), " ", "_");
experiment_date = strrep(experiment_date, ":", "_");
output_path = strcat(output_path, experiment_date, "_", region, "_", scheme, "/");
if(~exist(output_path, "dir"))
    mkdir(output_path);
end

fprintf("Created output directory for this experiment! \n");

label_path_root = '/Users/leobao/Documents/MultiPlanePipeline/AACR2023/Labels/';
matrix_path_root = '/Users/leobao/Documents/MultiPlanePipeline/AACR2023/';

training_label_path = string(strcat(label_path_root, 'Axial_Tumor', '/train_labels.csv'));
testing1_label_path = string(strcat(label_path_root, 'Axial_Tumor', '/test1_labels.csv'));
testing2_label_path = string(strcat(label_path_root, 'Axial_Tumor', '/test2_labels.csv'));

data_labels_training = readmatrix(training_label_path);
data_labels_holdout1 = readmatrix(testing1_label_path);
data_labels_holdout2 = readmatrix(testing2_label_path);

% feature_file_name = string(strcat(output_path, 'top_features.xlsx'));
% feature_table = cell2table(feature_column_names);
% writetable(feature_table, feature_file_name);

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