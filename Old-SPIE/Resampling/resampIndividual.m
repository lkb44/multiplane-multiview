clc;
clear;
close all;
setenv('PATH', getenv('PATH')+":/Applications/Slicer.app/Contents/lib/Slicer-5.2/cli-modules/")
addpath(genpath(pwd))
addpath(genpath("/Users/leobao/Library/CloudStorage/GoogleDrive-lkb44@case.edu/Shared drives/INVent/tom/Conferences/2022-SPIE/Feature_Selection_Classification_Code/Feature_Classifier/mha"));

problems = [];
pathIn = '/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/Data/AxialSingle';

DirectoryIn= '/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/Data/AxialSingle/RectalCA_*';
pathOut= '/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/Data/AxialProblems/';

Imgs = dir(fullfile(DirectoryIn));
theFilesvol = {Imgs.name};

for k= 1:length(theFilesvol)
    folderIn = dir(fullfile(pathIn, theFilesvol{k}, 'RectalCA_*'));
    if size(folderIn,1)==0, continue; end
    P = {folderIn.name};
    
    img_index = find(~contains(P,'label'));
    mask_index = find(contains(P,'label'));

    img_file_in = fullfile(pathIn, theFilesvol{k}, P{img_index});
    mask_file_in = fullfile(pathIn, theFilesvol{k},P{mask_index});

    mkdir(fullfile(pathOut, theFilesvol{k}));
    folderOut = dir(fullfile(pathOut, theFilesvol{k},'pre*'));
    img_path_out = fullfile(pathOut, theFilesvol{k});
    mask_path_out = fullfile(pathOut, theFilesvol{k});
    
    try
        resampleScalarVolume(1, 1, 1,img_file_in, img_path_out, mask_file_in,mask_path_out)
    catch
        warning(strcat("PROBLEM WITH ", img_file_in, ". MOVING ON..."));
        problems = [problems img_file_in];
    end
end