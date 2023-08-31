clc;
clear;
close all;
addpath(genpath('/Users/leobao/Library/CloudStorage/GoogleDrive-lkb44@case.edu/Shared drives/INVent/tom/Conferences/2022-SPIE/Resampling'));
addpath(genpath('/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/2D_Shape_Feature_Extraction_Code/scripts/mha'));
addpath(genpath('/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/2D_Shape_Feature_Extraction_Code/scripts/feature_extraction/subfunctions'))
addpath(genpath('/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/2D_Shape_Feature_Extraction_Code/scripts/feature_extraction/haralick'));
addpath(genpath('/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/2D_Shape_Feature_Extraction_Code/scripts/feature_extraction/grey'));
addpath(genpath('/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/2D_Shape_Feature_Extraction_Code/scripts/feature_extraction/gabor'));
addpath(genpath('/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/2D_Shape_Feature_Extraction_Code/scripts/feature_extraction/laws'));
addpath(genpath('/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/2D_Shape_Feature_Extraction_Code/scripts/feature_extraction/collage'));
addpath(genpath('/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/2D_Shape_Feature_Extraction_Code/scripts/feature_extraction'));

path_data = '/Volumes/Crucial X6/ReExtraction/Current/Patient-*';

pathIn = '/Volumes/Crucial X6/ReExtraction/Current/';
path = path_data;
pathOut= '/Volumes/Crucial X6/ReExtraction/CoronalTraining_TextureFeatures/';

regions = {'tumor'};

if strcmp(regions,'proxfat5')
    Regions = 'Proximal_Fat5';
end
if strcmp(regions,'proxfat10')
    Regions = 'Proximal_Fat10';
end
if strcmp(regions,'proxfat15')
    Regions = 'Proximal_Fat15';
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

class_options = {'gray','gradient','haralick','gabor','laws'};

feat = {};
featStats={};


for k = 1:length(theFilesvol)
    folder = theFilesvol(k);

    for rr = 1:length(regions)
        folderIn = dir(fullfile(pathIn, theFilesvol{k}, 'masks/*'));
        P = {folderIn.name};
        P = P(~startsWith(P, '.'));
        region = regions{rr};
        volIn = dir(fullfile(pathIn, theFilesvol{k}));
        V = {volIn.name};
        V = V(~startsWith(V, '.'));

        
        fat_index = find(contains(P,'_fat'));
        proxfat5_index = find(contains(P,'proxfat5'));
        tumor_index = find(contains(P,'tumor'));
        proxfat10_index = find(contains(P,'proxfat10'));
        proxfat15_index = find(contains(P,'proxfat15'));
        vol_index = find(contains(V,'ls'));

        case_mask_path_proxfat5 = fullfile(pathIn, theFilesvol{k}, 'masks/', P{proxfat5_index});
        case_mask_path_proxfat10 = fullfile(pathIn, theFilesvol{k}, 'masks/', P{proxfat10_index});
        case_mask_path_proxfat15 = fullfile(pathIn, theFilesvol{k}, 'masks/', P{proxfat15_index});
        case_mask_path_tumor = fullfile(pathIn, theFilesvol{k}, 'masks/', P{tumor_index});
        case_mask_path_fat = fullfile(pathIn, theFilesvol{k}, 'masks/', P{fat_index});
        case_vol_path = fullfile(pathIn, theFilesvol{k}, V{vol_index});

        if strcmp(regions,'proxfat5')
            case_mask_path = case_mask_path_proxfat5;
        end
        if strcmp(regions,'proxfat10')
            case_mask_path = case_mask_path_proxfat10;
        end
        if strcmp(regions,'proxfat15')
            case_mask_path = case_mask_path_proxfat15;
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

            file_name = strcat(folder, '_', region, '.mat');

            [featints, featnames, featstats, statnames] = extract2DFeatureInfo(vol,mask,class_options, ws_options);
            feature_matrix = reshape(featstats, 1, []);

            file_path = string(strcat(pathOut, Regions, '/', file_name));
            save(file_path, 'feature_matrix');
            disp(['Finished: ' case_mask_path])
            
        else
            disp(['Path not found: ' case_mask_path])
        end
    end
    
end



