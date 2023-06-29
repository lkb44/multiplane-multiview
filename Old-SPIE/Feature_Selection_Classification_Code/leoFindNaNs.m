clear
clc
% Create a vector with NaN values
data = load('/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/Data/AxialTrainingTextureFeatures/Axial_Proximal_Fat_2D_shape_features.mat');
feature_matrices = data.feature_matrix;
features_training = whitenData(feature_matrices, 'modified');

feature_vector = features_training(1, :);
% Find the indices of NaN values
nan_indices = find(isnan(feature_vector));

% Display the NaN indices
disp(nan_indices);

% Create a matrix
my_matrix = [1 2 3; 4 5 6; 7 8 9; 10 11 12];

% Create a vector of indices
index_vector = [2 4];

% Select the rows of the matrix with the indices from the vector
selected_rows = my_matrix(index_vector, :);

% Display the selected rows
disp(selected_rows);
