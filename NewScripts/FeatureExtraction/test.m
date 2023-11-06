clc;
clear;
close all;
addpath(genpath('/Users/leobao/Library/CloudStorage/GoogleDrive-lkb44@case.edu/Shared drives/INVent/tom/Conferences/2022-SPIE/Resampling'));
addpath(genpath('/Users/leobao/Library/CloudStorage/GoogleDrive-lkb44@case.edu/Shared drives/INVent/tom/Conferences/2022-SPIE/Feature_Selection_Classification_Code/Feature_Classifier/mha'));
addpath(genpath('/Users/leobao/Library/CloudStorage/GoogleDrive-lkb44@case.edu/Shared drives/INVent/sxs2658/TextureFeatures/subfunctions'))
addpath(genpath('/Users/leobao/Library/CloudStorage/GoogleDrive-lkb44@case.edu/Shared drives/INVent/sxs2658/TextureFeatures/haralick'));
addpath(genpath('/Users/leobao/Library/CloudStorage/GoogleDrive-lkb44@case.edu/Shared drives/INVent/sxs2658/TextureFeatures/grey'));
addpath(genpath('/Users/leobao/Library/CloudStorage/GoogleDrive-lkb44@case.edu/Shared drives/INVent/sxs2658/TextureFeatures/gabor'));
addpath(genpath('/Users/leobao/Library/CloudStorage/GoogleDrive-lkb44@case.edu/Shared drives/INVent/sxs2658/TextureFeatures'));
addpath(genpath('/Users/leobao/Library/CloudStorage/GoogleDrive-lkb44@case.edu/Shared drives/INVent/sxs2658/TextureFeatures/haralick/cpp_files'));
addpath(genpath('/Users/leobao/Library/CloudStorage/GoogleDrive-lkb44@case.edu/Shared drives/INVent/sxs2658/TextureFeatures/haralick'));
% addpath('D:\Data\scripts-master\')
%% Dataset Path
%path ='G:\Shared drives\INVent\tom\Conferences\2022_SPIE\Data\Coronal_resampled';
path_data = '/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/Data/CoronalTrainingLargestSlice/UH-RectalCA-*';
% path ='G:\Shared drives\INVent\sxs2658\Rectal_Data\UH_Resamp\Post_Resamp';

pathIn = '/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/Data/CoronalTrainingLargestSlice/';
path = path_data;
pathOut= '/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/Data/CoronalTrainingTextureFeatures/';

regions = {'tumor'};

if strcmp(regions,'proxfat')
    Regions = 'Proximal_Fat';
end
if strcmp(regions,'tumor')
    Regions = 'Tumor';
end
if strcmp(regions, 'fat')
    Regions = 'Fat';
end

Masks = dir(fullfile(path_data));
theFilesvol = {Masks.name};

if exist('ws_options','var')~=1 || isempty(ws_options)
    ws_options = 3:2:11;
end

ws1 = ws_options(1);
ws2 = ws_options(end);

%% Feature name
class_options = {'gray','gradient','haralick','gabor','laws'}; %(optional) feature classes to extract

%% Extraction
feat = {};
featStats={};


for k = 1:length(theFilesvol)
    %% *---------- Reading image(s) with a specific format -----------*%
%     folder = dir(fullfile(path, theFilesvol{k},'*.mha'));
%     P = {folder.name};
    folder = theFilesvol(k);

    for rr = 1:length(regions)
        folderIn = dir(fullfile(pathIn, theFilesvol{k}, 'masks/*'));
        P = {folderIn.name};
        region = regions{rr};
        volIn = dir(fullfile(pathIn, theFilesvol{k}));
        V = {volIn.name};

        
        fat_index = find(contains(P,'_fat'));
        proxfat_index = find(contains(P,'proxfat'));
        tumor_index = find(contains(P,'tumor'));
        vol_index = find(contains(V,'largest'));

        case_mask_path_proxfat = fullfile(pathIn, theFilesvol{k}, 'masks/',P{proxfat_index}); % For Linux
        case_mask_path_tumor = fullfile(pathIn, theFilesvol{k}, 'masks/',P{tumor_index}); % For Linux
        case_mask_path_fat = fullfile(pathIn, theFilesvol{k}, 'masks/', P{fat_index});
        case_vol_path = fullfile(pathIn, theFilesvol{k}, V{vol_index});

        if strcmp(regions,'proxfat')
            case_mask_path = case_mask_path_proxfat;
        end
        if strcmp(regions,'tumor')
            case_mask_path = case_mask_path_tumor;
        end
        if strcmp(regions, 'fat')
            case_mask_path = case_mask_path_fat;
        end
        
        if exist(case_mask_path)
            mskFileName = char(case_mask_path);
            volFileName = char(case_vol_path);
            fprintf(1, 'Now reading %s\n', volFileName);

            [Jvol, Jmask] = vv(volFileName,mskFileName);
            vol = Jvol;
            mask = Jmask;

            %[vol,mask] = boundingbox2(vol,mask,max(ws_options),'cropz','off');
            file_name = strcat(folder, '_', region, '.mat');

            [featints, featnames, featstats, statnames] = extract2DFeatureInfo(vol,mask,class_options, ws_options);
            % feat{k} = {cell2mat(featints)};
            feature_matrix = reshape(featstats, 1, []);
            % featStats{k} = {feature_matrix};
            % PatientId{k} = theFilesvol{k};

            file_path = string(strcat(pathOut, Regions, '/', file_name));
            % save('feat_preUH_CoronalFat.mat', 'feat','-v7.3');
            save(file_path, 'feature_matrix');
            % save VAfeatnames featnames
            % save UHfeatStats_PreCoronalFat featStats
            disp(['Finished: ' case_mask_path])
            
        else
            disp(['Path not found: ' case_mask_path])
        end
    end
    
end