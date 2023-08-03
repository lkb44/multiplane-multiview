% Code for extracting 2D shape feature
% Written by Charlems

clear;clc;
addpath(genpath('Libraries'));
addpath(genpath('/Users/leobao/Library/CloudStorage/GoogleDrive-lkb44@case.edu/Shared drives/INVent/tom/Conferences/2022-SPIE/Resampling'));
addpath(genpath('/Users/leobao/Library/CloudStorage/GoogleDrive-lkb44@case.edu/Shared drives/INVent/tom/Conferences/2022-SPIE/Feature_Selection_Classification_Code/Feature_Classifier/mha'));

%% Path to data
path_data = '/Users/leobao/Documents/MultiPlanePipeline/Data/AxialTrainingLargestSlice_TRG_T/Patient-*';

%dataset = 'UH-';

pathIn = '/Users/leobao/Documents/MultiPlanePipeline/Data/AxialTrainingLargestSlice_TRG_T/';
path = path_data; % For Linux
pathOut= '/Users/leobao/Documents/MultiPlanePipeline/Data/AxialTrainingShapeFeatures_TRG_T/';
%folders = {'RectalCA_'};
regions = {'proxfat10'};%, 'tumor', 'proxfat'}; %'Lumen' or 'ERW'

Masks = dir(fullfile(path));
theFilesvol = {Masks.name};

%% Extract features
for k = 1:length(theFilesvol)
    folder = theFilesvol(k);
    
    for rr = 1:length(regions)
        folderIn = dir(fullfile(pathIn, theFilesvol{k}, 'masks/*'));
        P = {folderIn.name};
        region = regions{rr};
        %for i=1:180 % Move case by case - 1:116
            %if i<=9 u_case = ['00' num2str(i)]; end
            %if i>=10 && i<=99 u_case = ['0' num2str(i)]; end
            %if i>=100 u_case = num2str(i); end

        fat_index = find(contains(P,'_fat'));
        proxfat5_index = find(contains(P,'proxfat5'));
        proxfat10_index = find(contains(P,'proxfat10'));
        proxfat15_index = find(contains(P,'proxfat15'));
        tumor_index = find(contains(P,'tumor'));
            
        case_mask_path_proxfat5 = fullfile(pathIn, theFilesvol{k}, 'masks/',P{proxfat5_index}); % For Linux
        case_mask_path_proxfat10 = fullfile(pathIn, theFilesvol{k}, 'masks/',P{proxfat10_index});
        case_mask_path_proxfat15 = fullfile(pathIn, theFilesvol{k}, 'masks/',P{proxfat15_index});
        case_mask_path_tumor = fullfile(pathIn, theFilesvol{k}, 'masks/',P{tumor_index}); % For Linux
        case_mask_path_fat = fullfile(pathIn, theFilesvol{k}, 'masks/', P{fat_index});

            
        if strcmp(region,'proxfat5')
            case_mask_path = case_mask_path_proxfat5;
        end
        if strcmp(region,'proxfat10')
            case_mask_path = case_mask_path_proxfat10;
        end
        if strcmp(region,'proxfat15')
            case_mask_path = case_mask_path_proxfat15;
        end
        if strcmp(region,'tumor')
            case_mask_path = case_mask_path_tumor;
        end
        if strcmp(region, 'fat')
            case_mask_path = case_mask_path_fat;
        end
            
        if exist(case_mask_path)
            vol_proxfat5 = vv(case_mask_path_proxfat5); %load_untouch_nii(case_mask_path_lumen); % use vv
            vol_proxfat10 = vv(case_mask_path_proxfat10); %load_untouch_nii(case_mask_path_lumen); % use vv
            vol_proxfat15 = vv(case_mask_path_proxfat15); %load_untouch_nii(case_mask_path_lumen); % use vv
            vol_tumor = vv(case_mask_path_tumor); %load_untouch_nii(case_mask_path_tumor); % use vv
            vol_fat = vv(case_mask_path_fat);
                %vol_rw = vv(case_mask_path_rw);
                
            if strcmp(region, 'tumor')
                feature_matrix=extract_2Dfeats(vol_tumor);
            end
            if strcmp(region,'proxfat5')
                feature_matrix=extract_2Dfeats(vol_proxfat5); %Function
            end
            if strcmp(region,'proxfat10')
                feature_matrix=extract_2Dfeats(vol_proxfat10); %Function
            end
            if strcmp(region,'proxfat15')
                feature_matrix=extract_2Dfeats(vol_proxfat15); %Function
            end
            if strcmp(region, 'fat')
                feature_matrix=extract_2Dfeats(vol_fat);
            end

            file_name = strcat(folder, '_', region, '.mat');

            %if strcmp(folder,[dataset 'RectalCA-'])
                
            file_path = string(strcat(pathOut, file_name)); % For Linux
                
            %end
            save(file_path, 'feature_matrix');

            disp(['Finished:' case_mask_path])
        else
            disp(['Path not found:' case_mask_path])
        end
        %end
    end
end


