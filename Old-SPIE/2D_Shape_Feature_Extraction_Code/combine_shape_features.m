clear all; close all; clc;

%% Specify parameters
view = "Coronal"; % change to 'Axial' or 'Coronal'
dataset = "va_test"; % change to 'train', 'test', or 'va_test'
roi = "fat"; % for axail: change to 'lumen', 'tumor', or 'fat.' For coronal: change to 'lumen', 'rw', or 'fat'

%% Specify filepaths
root_path = "/Volumes/GoogleDrive/Shared drives/INVent/tom/Conferences/2022_SPIE/Shape_Features/2D_Feature_Matrices/";

input_path = strcat(root_path, view, "/", dataset, "/"); 
output_path = strcat(root_path, view, '/', dataset, '/', view, "_", roi, "_", dataset, "_2D_shape_features.mat");

id_matrix_output_path = strcat(root_path, view, '/', dataset, '/');

%% Get patient IDs
patients = dir(fullfile(input_path, strcat('*RectalCA*')));
patient_ids = {patients.name};

%% Combine into one feature matrix

feature_matrix = [];

for i=1:length(patient_ids)
    
    pt_path = strcat(input_path, patient_ids{i}, '/', patient_ids{i}, '-', roi, '.mat');
    pt_features = load(pt_path);
    pt_features = pt_features.feature_matrix;
    
    feature_matrix = [feature_matrix; pt_features];
end

%% Save feature matrix

save(output_path, 'feature_matrix');

patient_ids = transpose(patient_ids);
patient_ids_table = cell2table(patient_ids);

table_fname = strcat(view, '_', 'case_index_matrix.xlsx');
writetable(patient_ids_table, fullfile(id_matrix_output_path, table_fname));
