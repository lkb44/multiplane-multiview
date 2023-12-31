clear
clc

input_path = '/Users/leobao/Documents/MultiPlanePipeline/CombinedTexture/Axial_Tumor/Testing1/HighFrequency/';
output_path = '/Users/leobao/Documents/MultiPlanePipeline/CombinedTexture/Axial_Tumor/Testing1/HighFrequency/Axial_Tumor_Testing1.mat';

% Get a list of .mat files in the directory
file_list = dir(fullfile(input_path, 'Patient-*.mat'));

% Initialize a cell array to store the patient IDs
patient_ids = cell(1, numel(file_list));

% Extract patient IDs from the filenames
for i = 1:numel(file_list)
    file_name = file_list(i).name;
    [~, name_without_extension] = fileparts(file_name); % Remove the '.mat' extension
    patient_id = name_without_extension(9:end); % Extract the patient ID (assuming the format)
    patient_ids{i} = patient_id;
end

% Convert the cell array to a string array (optional)
patient_ids_string = string(patient_ids);
feature_matrix = [];

for i = 1 : length(patient_ids_string)

    patient_path = strcat(input_path, 'Patient-', patient_ids_string{i}, '.mat');
    features = load(patient_path);
    features = features.feature_matrix;

    feature_matrix = [feature_matrix; features];

end

save(output_path, "feature_matrix");