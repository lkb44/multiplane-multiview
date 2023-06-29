clear
close all
clc

index = 202;

% Load the necessary data (example data)
training_data = load("/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/Data/AxialTrainingTextureFeatures/Axial_Proximal_Fat_2D_texture_features.mat");
testing_data = load("/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/Data/AxialTestingTextureFeatures/Axial_Proximal_Fat_2D_texture_features.mat");

training_feature_set = training_data.feature_matrix;
training_feature_set = training_feature_set(:, index);

testing_feature_set = testing_data.feature_matrix;
testing_feature_set = testing_feature_set(:, index);

combined_features = cat(1, training_feature_set, testing_feature_set);

training_responders = readtable("/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/Data/train_test_labels/train_labels.csv");
testing_responders = readtable("/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/Data/train_test_labels/test_labels.csv");

responders = cat(1, training_responders, testing_responders);
responders = table2array(responders);

% Reshape array into a column vector
responders = reshape(responders, [], 1);

responder_indices = find(responders == 0);
non_responder_indices = find(responders == 1);

data1 = combined_features(responder_indices);
data2 = combined_features(non_responder_indices);

% Define bin edges and number of bins for histogram
bin_edges = linspace(min([data1; data2]), max([data1; data2]), 20);
num_bins = length(bin_edges) - 1;

% Compute kernel density estimates for data
pdf1 = ksdensity(data1, bin_edges);
pdf2 = ksdensity(data2, bin_edges);

% Normalize kernel density estimates
pdf1 = pdf1 ./ sum(pdf1);
pdf2 = pdf2 ./ sum(pdf2);

% Compute histogram counts for data
counts1 = histcounts(data1, bin_edges);
counts2 = histcounts(data2, bin_edges);

% Compute widths of violin plots
width1 = 0.3;
width2 = 0.3;

% Create a new figure
figure;

% Plot violin plots
h1 = area(bin_edges(1:end-1), counts1./max(counts1).*width1, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'none');
hold on;
h2 = area(bin_edges(1:end-1), counts2./max(counts2).*width2, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'none');
alpha(h1,0.5);
alpha(h2,0.5);
plot(bin_edges(1:end-1) + width1/2, pdf1./max(pdf1).*width1, '-k', 'LineWidth', 1);
plot(bin_edges(1:end-1) + width2/2, pdf2./max(pdf2).*width2, '-k', 'LineWidth', 1);

% Add a title and axis labels
title('Violin Plots of Data 1 and Data 2');
xlabel('Data');
ylabel('Density');

% Add a legend
legend('Data 1', 'Data 2');
