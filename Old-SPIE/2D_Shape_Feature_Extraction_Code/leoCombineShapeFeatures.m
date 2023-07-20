clear all; close all; clc;

%% Specify parameters
view = "Coronal"; % change to 'Axial' or 'Coronal'
dataset = "TestingShapeFeatures"; % change to 'train', 'test', or 'va_test'
split = "_TRG_T/";
roi = "ProxFat10"; % for axail: change to 'lumen', 'tumor', or 'fat.' For coronal: change to 'lumen', 'rw', or 'fat'

%% Specify filepaths
root_path = "/Users/leobao/Documents/MultiPlanePipeline/Data/";

input_path = strcat(root_path, view, dataset, split, roi, '/'); 
output_path = strcat(root_path, view, dataset, split, roi, '/', view, '_',roi, "_2D_shape_features.mat");

id_matrix_output_path = strcat(root_path, view, dataset);

%% Get patient IDs
patients = dir(fullfile(input_path, strcat('*Patient-*')));
patient_ids = {patients.name};

%% Combine into one feature matrix

feature_matrix = [];

for i=1:length(patient_ids)
    
    pt_path = strcat(input_path, patient_ids{i}, '/');
    pt_features = load(pt_path);
    pt_features = pt_features.feature_matrix;
    
    feature_matrix = [feature_matrix; pt_features];
end

%% Save feature matrix

save(output_path, 'feature_matrix');

%patient_ids = transpose(patient_ids);
%patient_ids_table = cell2table(patient_ids);

%table_fname = strcat(view, '_', 'case_index_matrix.xlsx');
%writetable(patient_ids_table, fullfile(id_matrix_output_path, table_fname));
