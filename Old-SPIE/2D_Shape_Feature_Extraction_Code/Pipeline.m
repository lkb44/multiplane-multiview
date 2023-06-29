clear; close all; clc;

addpath(genpath("../FeatureExtraction/"));
addpath(genpath("../Libraries/"));

%% 2D shape feature extraction
featureExtraction;

%% 2D shape statistics calculation
statsCalculation;