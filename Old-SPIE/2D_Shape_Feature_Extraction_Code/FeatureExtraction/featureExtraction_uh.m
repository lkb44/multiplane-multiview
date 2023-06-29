% Code for extracting 2D shape feature
% Written by Charlems

clear;clc;
addpath(genpath('Libraries'));
addpath(genpath('/Volumes/GoogleDrive/Shared drives/INVent/tom/Conferences/2022_SPIE/Resampling/'));
addpath(genpath('/Volumes/GoogleDrive/Shared drives/INVent/tom/Conferences/2022_SPIE/2D_Shape_Feature_Extraction_Code/Libraries/images_lib/mha/'));

%% Path to data
path_data = '/Volumes/GoogleDrive/Shared drives/INVent/tom/Conferences/2022_SPIE/Data/Coronal_resampled/';

dataset = 'UH-';


if isunix
    path = [path_data]; % For Linux
elseif ispc
    path = [path_data dataset 'RectalCA\']; % For Windows
end

folders = {[dataset 'RectalCA-']};
regions = {'lumen'}; %'Lumen' or 'ERW'

%% Extract features
for ss = 1:length(folders)
    folder = folders{ss};
    for rr = 1:length(regions)
        region = regions{rr};
        for i=1:180 % Move case by case - 1:116
            if i<=9 u_case = ['00' num2str(i)]; end
            if i>=10 && i<=99 u_case = ['0' num2str(i)]; end
            if i>=100 u_case = num2str(i); end
            
            if isunix
                case_mask_path_lumen = [path folder u_case '/masks/mask_lumen_largest_slice.mha']; % For Linux
                case_mask_path_tumor = [path folder u_case '/masks/mask_tumor_largest_slice.mha']; % For Linux
                case_mask_path_fat = [path folder u_case '/masks/mask_fat_largest_slice.mha'];
                case_mask_path_rw = [path folder u_case '/masks/mask_rw_largest_slice.mha'];
                
            elseif ispc
                case_mask_path_lumen = [path folder u_case '\masks\mask-2.nii']; % For Windows
                case_mask_path_tumor = [path folder u_case '\masks\mask-ERW.nii']; % For Windows
            end
            
            if strcmp(region,'lumen')
                case_mask_path = case_mask_path_lumen;
            end
            if strcmp(region,'tumor')
                case_mask_path = case_mask_path_tumor;
            end
            if strcmp(region, 'fat')
                case_mask_path = case_mask_path_fat;
            end
            if strcmp(region, 'rw')
                case_mask_path = case_mask_path_rw;
            end
            
            if exist(case_mask_path)
                vol_lumen = vv(case_mask_path_lumen); %load_untouch_nii(case_mask_path_lumen); % use vv
                %vol_tumor = vv(case_mask_path_tumor); %load_untouch_nii(case_mask_path_tumor); % use vv
                vol_fat = vv(case_mask_path_fat);
                vol_rw = vv(case_mask_path_rw);
                
                if strcmp(region,'lumen')
                    feature_matrix=extract_2Dfeats(vol_lumen); %Function
                end
                if strcmp(region, 'tumor')
                    feature_matrix=extract_2Dfeats(vol_tumor);
                end
                if strcmp(region,'rw')
                    feature_matrix=extract_2Dfeats(vol_rw); %Function
                end
                if strcmp(region, 'fat')
                    feature_matrix=extract_2Dfeats(vol_fat);
                end
                
                file_name = [folder u_case '-' region '.mat'];
                if strcmp(folder,[dataset 'RectalCA-'])
                    if isunix
                        file_path = ['/Volumes/GoogleDrive/Shared drives/INVent/tom/Conferences/2022_SPIE/Shape_Features/2D_Feature_Matrices/Coronal/All/' file_name]; % For Linux
                    elseif ispc
                        file_path = ['Results\mat_' dataset '_train\2D_ShapeDescriptors\' file_name]; % For Windows
                    end
                end
                save(file_path, 'feature_matrix');
                disp(['Finished:' case_mask_path])
            else
                disp(['Path not found:' case_mask_path])
            end
        end
    end
end