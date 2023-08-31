clc;
clear;
close all;

view = "Coronal";
plane = 'cor';
roi = "ProxFat10";
region = 'proxfat10';

root_path = "/Volumes/Crucial X6/CombineTexture/";

matlab_input_path = strcat(root_path, view, '_', roi, '_Matlab/');
python_input_path = strcat(root_path, view, '_', roi, '_Python/');
output_path = strcat(root_path, view, '_', roi, '/');

% Get a list of .mat files in the directory
patient_list = dir(fullfile(python_input_path, 'Feats_Col_Patient-*.mat'));

% Initialize a cell array to store the patient IDs
patient_ids = cell(1, numel(patient_list));

% Extract patient IDs from the filenames
for i = 1:numel(patient_list)
    file_name = patient_list(i).name;
    [~, name_without_extension] = fileparts(file_name); % Remove the '.mat' extension
    patient_id = name_without_extension(19:21); % Extract the patient ID (assuming the format)
    patient_ids{i} = patient_id;
end

% Convert the cell array to a string array (optional)
patient_ids_string = string(patient_ids);


for i = 1 : length(patient_ids_string)
    matlab_patient_path = strcat(matlab_input_path, 'Patient-', patient_ids_string{i}, '_', region, '.mat');
    python_patient_path = strcat(python_input_path, 'Feats_Col_Patient-', patient_ids_string{i}, '_', plane, '.mat');
    output_patient_path = strcat(output_path, 'Patient-', patient_ids_string{i}, '.mat');

    matlab_features = load(matlab_patient_path);
    python_features = load(python_patient_path);

    feature_matrix = [];
    matlab = matlab_features.feature_matrix;
    python = python_features.array_data;

    feature_matrix = [matlab python];

    save(output_patient_path, "feature_matrix")

end