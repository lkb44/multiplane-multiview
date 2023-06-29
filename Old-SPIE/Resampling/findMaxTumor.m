close all
clear
clc

addpath(genpath(pwd))
addpath(genpath("/Users/leobao/Library/CloudStorage/GoogleDrive-lkb44@case.edu/Shared drives/INVent/tom/Conferences/2022-SPIE/Feature_Selection_Classification_Code/Feature_Classifier/mha"));

pathIn = '/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/Data/AxialMasks/';

DirectoryIn= '/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/Data/AxialMasks/RectalCA_*';
pathOut= '/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/Data/MaxAxialTumorSlices/';

Imgs = dir(fullfile(DirectoryIn));
theFilesvol = {Imgs.name};

for k = 1:length(theFilesvol)
    folderIn = dir(fullfile(pathIn, theFilesvol{k}, 'RectalCA_*'));
    if size(folderIn,1)==0, continue; end
    P = {folderIn.name};
    
    tumor_index = find(contains(P, 'tumor'));

    tumor_file_in = fullfile(pathIn, theFilesvol{k}, P{tumor_index});

    tumor_info = mha_read_header(tumor_file_in);
    tumor_data = mha_read_volume(tumor_info);
    
    % Finds index of largest tumor slice
    max_tumor = findMaxNSlices(tumor_data, 1);

    % Selects individual slice from tumor mask
    slice = tumor_data(:, :, max_tumor);

    % Modify header information
    slice.Dimensions = [size(slice,1), size(slice,2), 1];
    slice.ElementSpacing = [slice.ElementSpacing(1), slice.ElementSpacing(2), 0];
    mkdir(fullfile(pathOut, theFilesvol{k}));
end