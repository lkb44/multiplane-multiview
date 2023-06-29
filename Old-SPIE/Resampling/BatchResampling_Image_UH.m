clc;
clear;
close all;
setenv('PATH', getenv('PATH')+":/Applications/Slicer.app/Contents/lib/Slicer-4.11/cli-modules/")
addpath(genpath(pwd))
addpath(genpath("/Volumes/GoogleDrive/My Drive/tom/Rectal Segmentation/Feature Classifier/mha"));

problems = [];
pathIn = '../Data/Original_Coronal/UH/';
% path ='G:\Shared drives\INVent\amir\Rectal\UH_preCRT_Resampled';
DirectoryIn= '../Data/Original_Coronal/UH/UH-RectalCA-*';
pathOut= '../Data/Coronal_resampled/';
Imgs = dir(fullfile(DirectoryIn));
theFilesvol = {Imgs.name};
for k= 1:length(theFilesvol)
    folderIn = dir(fullfile(pathIn, theFilesvol{k}, 'UH-RectalCA_*'));
    if size(folderIn,1)==0, continue; end
    P = {folderIn.name};
    
    img_index = find(~contains(P,'label'));
    mask_index = find(contains(P,'label'));

    img_file_in = fullfile(pathIn, theFilesvol{k}, P{img_index});
    mask_file_in = fullfile(pathIn, theFilesvol{k},P{mask_index});

    mkdir(fullfile(pathOut, theFilesvol{k}));
%     folderOut = dir(fullfile(pathOut, theFilesvol{k},'pre*'));
    img_path_out = fullfile(pathOut, theFilesvol{k});
    mask_path_out = fullfile(pathOut, theFilesvol{k});
    
    try
        resampleScalarVolume(1, 1, 1,img_file_in, img_path_out, mask_file_in,mask_path_out)
    catch
        warning(strcat("PROBLEM WITH ", img_file_in, ". MOVING ON..."));
        problems = [problems img_file_in];
    end
end

%https://www.mathworks.com/matlabcentral/answers/740062-use-command-line-in-matlab-using-and-zsh-1-command-not-found-is-displayed