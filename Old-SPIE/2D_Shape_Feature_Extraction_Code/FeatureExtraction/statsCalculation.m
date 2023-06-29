% Code for constructing the feature matrix based statistical moments per shape feature
% Written by Charlems


% Load feature files
clear;clc;
dataset = 'UH-';
folders = {[dataset 'RectalCA-']};

for ss = 1:length(folders)
    folder = folders{ss};
    
    if strcmp(folder,'UH-RectalCA-')
        if isunix
            folder_path = ['/Volumes/GoogleDrive/Shared drives/INVent/tom/Conferences/2022_SPIE/Shape_Features/2D_Feature_Matrices/Axial/' folder]; % For Linux
        elseif ispc
            folder_path = ['Results\mat_' dataset '_train\2D_ShapeDescriptors\']; % For Windows
        end

    end

    % statistics from feature DIFFERENCES between consecutive slices
    %diff_statMom_lumen = getStatMom_Diff_2DFeats(folder_path, folder, '-Lumen', 1); % last parameter, the number of cases you need to analyse
    %diff_statMom_ERW = getStatMom_Diff_2DFeats(folder_path, folder, '-ERW',1); % last parameter, the number of cases you need to analyse
    
    % statistics from raw FEATURES along a volume
    raw_statMom_lumen = getStatMom_Raw_2DFeats(folder_path, folder, '-fat', 2); % last parameter, the number of cases you need to analyse
    %raw_statMom_ERW = getStatMom_Raw_2DFeats(folder_path, folder, '-ERW', 1); % last parameter, the number of cases you need to analyse

    if strcmp(folder,'UH-RectalCA-')
        if isunix
            save_path = ['/Volumes/GoogleDrive/Shared drives/INVent/tom/Conferences/2022_SPIE/Shape_Features/2D_Feature_Matrices/Axial_descriptors/' dataset '/']; % For Linux
        elseif ispc
            save_path = ['Results\mat_' dataset '_train\2D_ShapeStatistics\']; % For Windows
        end

    end
    
    %save([save_path 'Diff_StatMomFeatsLumen_train.mat'],'diff_statMom_lumen');
    %save([save_path 'Diff_StatMomFeatsERW_train.mat'],'diff_statMom_ERW');
    save([save_path 'Raw_StatMomFeatsLumen_train.mat'],'raw_statMom_lumen');
    %save([save_path 'Raw_StatMomFeatsERW_train.mat'],'raw_statMom_ERW');
end