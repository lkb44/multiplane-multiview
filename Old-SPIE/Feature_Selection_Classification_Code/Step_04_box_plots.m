close all; clear; clc;


%% Specify view and rois
view = {'Coronal'}; % {Can only be 'Axial'}, {'Coronal'} or {'Axial','Coronal'}
rois = {'rw'};

%% Get feature names
feature_column_names = readtable("../2D_Shape_Feature_Extraction_Code/2d_shape_feature_names.xlsx", 'ReadVariableNames',false);
feature_column_names = table2cell(feature_column_names);

%% Load datasets

% Set directory path to the feature matrices depending on the view
if(length(view) == 1 && view{1} == "Axial")
    train_matrix_path = "../Shape_Features/Feature_Matrices/Axial/train/";
    test_matrix_path = "../Shape_Features/Feature_Matrices/Axial/test/";
    va_test_matrix_path = "../Shape_Features/Feature_Matrices/Axial/va_test/";
elseif(length(view) == 1 && view{1} == "Coronal")
    train_matrix_path = "../Shape_Features/2D_Feature_Matrices/Coronal/train/"; 
    test_matrix_path = "../Shape_Features/2D_Feature_Matrices/Coronal/test/";
    va_test_matrix_path = "../Shape_Features/2D_Feature_Matrices/Coronal/va_test/";
elseif(length(view) == 2)
    % Assert Axial is first in list and Coronal is second. This is
    % important for data loading and naming output logs. If view is not in
    % this format, the results will be confusing!!
    assert(strcmp(view{1}, 'Axial') && strcmp(view{2}, 'Coronal'), 'View Cell array must contain Axial first and Coronal second');
    ax_train_matrix_path = "../Shape_Features/Feature_Matrices/Axial/train/";
    ax_test_matrix_path = "../Shape_Features/Feature_Matrices/Axial/test/";
    ax_va_test_matrix_path = "../Shape_Features/Feature_Matrices/Axial/va_test/";
    cor_train_matrix_path = "../Shape_Features/Feature_Matrices/Coronal/train/"; 
    cor_test_matrix_path = "../Shape_Features/Feature_Matrices/Coronal/test/";
    cor_va_test_matrix_path = "../Shape_Features/Feature_Matrices/Coronal/va_test/";
end

% If working with only 1 view
if(length(view) == 1)
    for i=1:length(rois)
        path_to_train = strcat(train_matrix_path, view{1}, "_", rois{i}, "_train_2D_shape_features.mat");
        datasets.(strcat("train_",rois{i})) = load(path_to_train);

        path_to_test = strcat(test_matrix_path, view{1}, "_", rois{i}, "_test_2D_shape_features.mat");
        datasets.(strcat("test_",rois{i})) = load(path_to_test);
        
        path_to_test = strcat(va_test_matrix_path, view{1}, "_", rois{i}, "_va_test_2D_shape_features.mat");
        datasets.(strcat("va_test_",rois{i})) = load(path_to_test);
    end
else % Working with both views
    for i=1:length(rois)
        
        % Axial datasets
        path_to_train = strcat(ax_train_matrix_path, view{1}, "_", rois{i}, "_train_3D_shape_features.mat");
        datasets.(strcat(view{1}, "_train_",rois{i})) = load(path_to_train);
        
        path_to_test = strcat(ax_test_matrix_path, view{1}, "_", rois{i}, "_test_3D_shape_features.mat");
        datasets.(strcat(view{1}, "_test_",rois{i})) = load(path_to_test);
        
        path_to_test = strcat(ax_va_test_matrix_path, view{1}, "_", rois{i}, "_va_test_3D_shape_features.mat");
        datasets.(strcat(view{1}, "_va_test_",rois{i})) = load(path_to_test);
        
        % Coronal Datasets
        path_to_train = strcat(cor_train_matrix_path, view{2}, "_", rois{i}, "_train_3D_shape_features.mat");
        datasets.(strcat(view{2}, "_train_",rois{i})) = load(path_to_train);
        
        path_to_test = strcat(cor_test_matrix_path, view{2}, "_", rois{i}, "_test_3D_shape_features.mat");
        datasets.(strcat(view{2}, "_test_",rois{i})) = load(path_to_test);
        
        path_to_test = strcat(cor_va_test_matrix_path, view{2}, "_", rois{i}, "_va_test_3D_shape_features.mat");
        datasets.(strcat(view{2}, "_va_test_",rois{i})) = load(path_to_test);
        
    end

end

%% Load in patient IDs and labels
train_labels_path = "../Data/train_test_labels/train_labels.csv";
test_labels_path = "../Data/train_test_labels/test_labels.csv";
va_test_labels = "../Data/train_test_labels/va_test_labels.csv";

data_labels_training = readtable(train_labels_path);
train_pts = data_labels_training.Patient;
data_labels_training = data_labels_training.ypT;

uh_data_labels_holdout = readtable(test_labels_path);
uh_test_pts = uh_data_labels_holdout.Patient;
uh_data_labels_holdout = uh_data_labels_holdout.ypT;

va_data_labels = readtable(va_test_labels);
va_test_pts = va_data_labels.Patient;
va_data_labels = va_data_labels.ypT;

data_labels_holdout = vertcat(uh_data_labels_holdout, va_data_labels);
test_pts = vertcat(uh_test_pts, va_test_pts);

%% Axial View, Multiple Regions

% Axial Training Metrics
axial_lumen_top_feats_table = readtable("../Shape_and_texture_Feature_Results/21-Aug-2022_22_34_30_shape_Coronal_fat_only/Top_Features_Information.xlsx");
axial_lumen_train_feat_array = datasets.train_rw.feature_matrix(:, axial_lumen_top_feats_table.indices);

axial_lumen_test_feat_array = datasets.test_rw.feature_matrix(:, axial_lumen_top_feats_table.indices);
axial_lumen_test_feat_array = vertcat(axial_lumen_test_feat_array, datasets.va_test_rw.feature_matrix(:, axial_lumen_top_feats_table.indices));

feat_array_names = axial_lumen_top_feats_table.names;
feat_array_names = strrep(feat_array_names, " ", "_");

[training_class_0_feats, training_class_1_feats] = get_classes(axial_lumen_train_feat_array, data_labels_training);
[testing_class_0_feats, testing_class_1_feats] = get_classes(axial_lumen_test_feat_array, data_labels_holdout);

training_class_0_boxes = create_class_struct(training_class_0_feats, feat_array_names, [1:5]);
training_class_1_boxes = create_class_struct(training_class_1_feats, feat_array_names, [1:5]);
testing_class_0_boxes = create_class_struct(testing_class_0_feats, feat_array_names, [1:5]);
testing_class_1_boxes = create_class_struct(testing_class_1_feats, feat_array_names, [1:5]);

for k=1:length(feat_array_names)
    
    [feature_cell,legend] = create_feature_classes(training_class_0_boxes, training_class_1_boxes,...
                        testing_class_0_boxes, testing_class_1_boxes,...
                        'class',feat_array_names{k}, "T0-T2 vs. T3-T4");

    color1 = 1/255 * [161 230 255];
    color2 = 1/255 * [224 161 255];
    colors = [color1; color2];
    
    if(k == 1)
        %ylim = [-95 65];
        ylim = [-100 80];
    elseif(k == 2)
        %ylim = [0.3 1.1];
        ylim = [0 1];
    elseif(k == 3)
        ylim = [0, 1.1];
    elseif(k == 4)
        %ylim = [-1000 10500]
        ylim = [0 1];
    elseif(k == 5)
        ylim = [0 1.1];
    end
    
    %ylim = [round(min(cell2mat(feature_cell(1))), 1) round(max(cell2mat(feature_cell(1))), 1)];

    test_feature_cell = cell(1,2);
    test_feature_cell{:, 1} = feature_cell{:, 3};
    test_feature_cell{:, 2} = feature_cell{:, 4};
    legend = [legend(:,3) legend(:,4)];


    featname = feat_array_names{k};
    featname_formatted = strrep(featname, "_", " ");
    %plotfeatureBoxPlot(feature_cell, featname_formatted, "Cohort", featname_formatted, legend, ylim, 20.00, 10.00, 0, colors, 'class')
    plotfeatureBoxPlot(test_feature_cell, featname_formatted, "Cohort", featname_formatted, legend, ylim, 20.00, 10.00, 0, colors, 'class')
    output_name = strcat('../Abstract/FIgures/Coronal_Shape_Boxplot/', featname, '_boxplot.png');
    
%     if(k == 4)
%         set(gca, 'YScale', 'log')
%     end
    saveas(gcf, output_name);

end



% ax_rw_train = readtable("../Shape_Feature_Results/09-Aug-2022_00_27_49_Axial_rw_only/CV_Training_Results.xlsx");
% ax_fat_train = readtable("../Shape_Feature_Results/09-Aug-2022_00_27_58_Axial_fat_only/CV_Training_Results.xlsx");
% ax_all_train = readtable("../Shape_Feature_Results/09-Aug-2022_00_29_08_Axial_all/CV_Training_Results.xlsx");
% 
% % Axial Holdout Testing UH + VA Metrics
% axial_lumen_uh_va = readtable("../Shape_Feature_Results/09-Aug-2022_00_27_39_Axial_lumen_only/Holdout_Testing_Dataset_Results.xlsx", "Sheet", "Entire_Test_Cohort_Results");
% axial_rw_uh_va = readtable("../Shape_Feature_Results/09-Aug-2022_00_27_49_Axial_rw_only/Holdout_Testing_Dataset_Results.xlsx", "Sheet", "Entire_Test_Cohort_Results");
% axial_fat_uh_va = readtable("../Shape_Feature_Results/09-Aug-2022_00_27_58_Axial_fat_only/Holdout_Testing_Dataset_Results.xlsx", "Sheet", "Entire_Test_Cohort_Results");
% axial_all_uh_va = readtable("../Shape_Feature_Results/09-Aug-2022_00_29_08_Axial_all/Holdout_Testing_Dataset_Results.xlsx", "Sheet", "Entire_Test_Cohort_Results");
% 
% % Axial Holdout Testing UH Only Metrics
% axial_lumen_uh = readtable("../Shape_Feature_Results/09-Aug-2022_00_27_39_Axial_lumen_only/Holdout_Testing_Dataset_Results.xlsx", "Sheet", "UH_Only");
% axial_rw_uh = readtable("../Shape_Feature_Results/09-Aug-2022_00_27_49_Axial_rw_only/Holdout_Testing_Dataset_Results.xlsx", "Sheet", "UH_Only");
% axial_fat_uh = readtable("../Shape_Feature_Results/09-Aug-2022_00_27_58_Axial_fat_only/Holdout_Testing_Dataset_Results.xlsx", "Sheet", "UH_Only");
% axial_all_uh = readtable("../Shape_Feature_Results/09-Aug-2022_00_29_08_Axial_all/Holdout_Testing_Dataset_Results.xlsx", "Sheet", "UH_Only");
% 
% % Axial Holdout Testing VA Only Metrics
% axial_lumen_va = readtable("../Shape_Feature_Results/09-Aug-2022_00_27_39_Axial_lumen_only/Holdout_Testing_Dataset_Results.xlsx", "Sheet", "VA_Only");
% axial_rw_va = readtable("../Shape_Feature_Results/09-Aug-2022_00_27_49_Axial_rw_only/Holdout_Testing_Dataset_Results.xlsx", "Sheet", "VA_Only");
% axial_fat_va = readtable("../Shape_Feature_Results/09-Aug-2022_00_27_58_Axial_fat_only/Holdout_Testing_Dataset_Results.xlsx", "Sheet", "VA_Only");
% axial_all_va = readtable("../Shape_Feature_Results/09-Aug-2022_00_29_08_Axial_all/Holdout_Testing_Dataset_Results.xlsx", "Sheet", "VA_Only");
% 
% % create barplot
% % The x-axes is the dataset
% x_axes = categorical({'Training';...
%                       'Holdout Testing UH + VA';...
%                       'Holdout Testing UH Only';...
%                       'Holdout Testing VA Only'});
% x_axes = reordercats(x_axes, {'Training';...
%                       'Holdout Testing UH + VA';...
%                       'Holdout Testing UH Only';...
%                       'Holdout Testing VA Only'});
%     
% % y-axes will be the mean metric per region grouped by dataset
% y_axes = [axial_lumen_train.AUC ax_rw_train.AUC ax_fat_train.AUC ax_all_train.AUC;...
%           axial_lumen_uh_va.AUC axial_rw_uh_va.AUC axial_fat_uh_va.AUC axial_all_uh_va.AUC;...
%           axial_lumen_uh.AUC axial_rw_uh.AUC axial_fat_uh.AUC axial_all_uh.AUC;...
%           axial_lumen_va.AUC axial_rw_va.AUC axial_fat_va.AUC axial_all_va.AUC];
%       
% figure();
% set(gca,'FontSize',18)
% b = bar(x_axes, y_axes);
% xlabel("Dataset", 'FontSize', 18);
% ylabel("Mean AUC", 'FontSize', 18);
% ylim([0,1.05]);
% title("AUC: Axial View", 'Fontsize', 20);
% b(1).FaceColor = [39/255, 174/255, 239/255];
% b(2).FaceColor = [135/255, 188/255, 69/255];
% b(3).FaceColor = [179/255, 61/255, 198/255];
% b(4).FaceColor = [239/255, 155/255, 32/255];
% 
% set(b, {'DisplayName'}, {'Lumen','Rectal Wall','Fat', 'All Regions'}')
% legend('Location', 'northeastoutside')
% 
% output_path = ("../Shape_Feature_Results/ax_auc_bar_plot.png");
% saveas(gcf, output_path);
% 
% %% Coronal View, Multiple Regions
% 
% % Coronal Training Metrics
% coronal_lumen_train = readtable("../Shape_Feature_Results/09-Aug-2022_10_23_25_Coronal_lumen_only/CV_Training_Results.xlsx");
% coronal_rw_train = readtable("../Shape_Feature_Results/09-Aug-2022_10_23_43_Coronal_rw_only/CV_Training_Results.xlsx");
% coronal_fat_train = readtable("../Shape_Feature_Results/09-Aug-2022_10_23_55_Coronal_fat_only/CV_Training_Results.xlsx");
% coronal_all_train = readtable("../Shape_Feature_Results/09-Aug-2022_10_46_38_Coronal_all/CV_Training_Results.xlsx");
%  
% % Coronal Holdout Testing UH + VA Metrics
% coronal_lumen_uh_va = readtable("../Shape_Feature_Results/09-Aug-2022_10_23_25_Coronal_lumen_only/Holdout_Testing_Dataset_Results.xlsx", "Sheet", "Entire_Test_Cohort_Results");
% coronal_rw_uh_va = readtable("../Shape_Feature_Results/09-Aug-2022_10_23_43_Coronal_rw_only/Holdout_Testing_Dataset_Results.xlsx", "Sheet", "Entire_Test_Cohort_Results");
% coronal_fat_uh_va = readtable("../Shape_Feature_Results/09-Aug-2022_10_23_55_Coronal_fat_only/Holdout_Testing_Dataset_Results.xlsx", "Sheet", "Entire_Test_Cohort_Results");
% coronal_all_uh_va = readtable("../Shape_Feature_Results/09-Aug-2022_10_46_38_Coronal_all/Holdout_Testing_Dataset_Results.xlsx", "Sheet", "Entire_Test_Cohort_Results");
% 
% % Coronal Holdout Testing UH Only Metrics
% coronal_lumen_uh = readtable("../Shape_Feature_Results/09-Aug-2022_10_23_25_Coronal_lumen_only/Holdout_Testing_Dataset_Results.xlsx", "Sheet", "UH_Only");
% coronal_rw_uh = readtable("../Shape_Feature_Results/09-Aug-2022_10_23_43_Coronal_rw_only/Holdout_Testing_Dataset_Results.xlsx", "Sheet", "UH_Only");
% coronal_fat_uh = readtable("../Shape_Feature_Results/09-Aug-2022_10_23_55_Coronal_fat_only/Holdout_Testing_Dataset_Results.xlsx", "Sheet", "UH_Only");
% coronal_all_uh = readtable("../Shape_Feature_Results/09-Aug-2022_10_46_38_Coronal_all/Holdout_Testing_Dataset_Results.xlsx", "Sheet", "UH_Only");
% 
% % Coronal Holdout Testing VA Only Metrics
% coronal_lumen_va = readtable("../Shape_Feature_Results/09-Aug-2022_10_23_25_Coronal_lumen_only/Holdout_Testing_Dataset_Results.xlsx", "Sheet", "VA_Only");
% coronal_rw_va = readtable("../Shape_Feature_Results/09-Aug-2022_10_23_43_Coronal_rw_only/Holdout_Testing_Dataset_Results.xlsx", "Sheet", "VA_Only");
% coronal_fat_va = readtable("../Shape_Feature_Results/09-Aug-2022_10_23_55_Coronal_fat_only/Holdout_Testing_Dataset_Results.xlsx", "Sheet", "VA_Only");
% coronal_all_va = readtable("../Shape_Feature_Results/09-Aug-2022_10_46_38_Coronal_all/Holdout_Testing_Dataset_Results.xlsx", "Sheet", "VA_Only");
% 
% % Create barplot
% % The x-axes is the dataset
% x_axes = categorical({'Training';...
%                       'Holdout Testing UH + VA';...
%                       'Holdout Testing UH Only';...
%                       'Holdout Testing VA Only'});
% x_axes = reordercats(x_axes, {'Training';...
%                       'Holdout Testing UH + VA';...
%                       'Holdout Testing UH Only';...
%                       'Holdout Testing VA Only'});
%                   
% % y-axes will be the mean metric per region grouped by dataset
% y_axes = [coronal_lumen_train.AUC coronal_rw_train.AUC coronal_fat_train.AUC coronal_all_train.AUC;...
%           coronal_lumen_uh_va.AUC coronal_rw_uh_va.AUC coronal_fat_uh_va.AUC coronal_all_uh_va.AUC;...
%           coronal_lumen_uh.AUC coronal_rw_uh.AUC coronal_fat_uh.AUC coronal_all_uh.AUC;...
%           coronal_lumen_va.AUC coronal_rw_va.AUC coronal_fat_va.AUC coronal_all_va.AUC];
% 
% figure();
% set(gca,'FontSize',18)
% b = bar(x_axes, y_axes);
% xlabel("Dataset", 'FontSize', 18);
% ylabel("Mean AUC", 'FontSize', 18);
% ylim([0,1.05]);
% title("AUC: Coronal View", 'Fontsize', 20);
% b(1).FaceColor = [39/255, 174/255, 239/255];
% b(2).FaceColor = [135/255, 188/255, 69/255];
% b(3).FaceColor = [179/255, 61/255, 198/255];
% b(4).FaceColor = [239/255, 155/255, 32/255];
% 
% set(b, {'DisplayName'}, {'Lumen','Rectal Wall','Fat', 'All Regions'}')
% legend('Location', 'northeastoutside');
% 
% output_path = ("../Shape_Feature_Results/cor_auc_bar_plot.png");
% saveas(gcf, output_path);
% 
% %% Axial and Coronal View, Multiple Regions
% % Coronal Training Metrics
% ax_cor_lumen_train = readtable("../Shape_Feature_Results/09-Aug-2022_11_30_08_Axial_Coronal_lumen_only/CV_Training_Results.xlsx");
% ax_cor_rw_train = readtable("../Shape_Feature_Results/09-Aug-2022_11_30_30_Axial_Coronal_rw_only/CV_Training_Results.xlsx");
% ax_cor_fat_train = readtable("../Shape_Feature_Results/09-Aug-2022_11_30_50_Axial_Coronal_fat_only/CV_Training_Results.xlsx");
% ax_cor_all_train = readtable("../Shape_Feature_Results/09-Aug-2022_11_32_29_Axial_Coronal_all/CV_Training_Results.xlsx");
%  
% % Coronal Holdout Testing UH + VA Metrics
% ax_cor_lumen_uh_va = readtable("../Shape_Feature_Results/09-Aug-2022_11_30_08_Axial_Coronal_lumen_only/Holdout_Testing_Dataset_Results.xlsx", "Sheet", "Entire_Test_Cohort_Results");
% ax_cor_rw_uh_va = readtable("../Shape_Feature_Results/09-Aug-2022_11_30_30_Axial_Coronal_rw_only/Holdout_Testing_Dataset_Results.xlsx", "Sheet", "Entire_Test_Cohort_Results");
% ax_cor_fat_uh_va = readtable("../Shape_Feature_Results/09-Aug-2022_11_30_50_Axial_Coronal_fat_only/Holdout_Testing_Dataset_Results.xlsx", "Sheet", "Entire_Test_Cohort_Results");
% ax_cor_all_uh_va = readtable("../Shape_Feature_Results/09-Aug-2022_11_32_29_Axial_Coronal_all/Holdout_Testing_Dataset_Results.xlsx", "Sheet", "Entire_Test_Cohort_Results");
%  
% % Coronal Holdout Testing UH Only Metrics
% ax_cor_lumen_uh = readtable("../Shape_Feature_Results/09-Aug-2022_11_30_08_Axial_Coronal_lumen_only/Holdout_Testing_Dataset_Results.xlsx", "Sheet", "UH_Only");
% ax_cor_rw_uh = readtable("../Shape_Feature_Results/09-Aug-2022_11_30_30_Axial_Coronal_rw_only/Holdout_Testing_Dataset_Results.xlsx", "Sheet", "UH_Only");
% ax_cor_fat_uh = readtable("../Shape_Feature_Results/09-Aug-2022_11_30_50_Axial_Coronal_fat_only/Holdout_Testing_Dataset_Results.xlsx", "Sheet", "UH_Only");
% ax_cor_all_uh = readtable("../Shape_Feature_Results/09-Aug-2022_11_32_29_Axial_Coronal_all/Holdout_Testing_Dataset_Results.xlsx", "Sheet", "UH_Only");
%  
% % Coronal Holdout Testing VA Only Metrics
% ax_cor_lumen_va = readtable("../Shape_Feature_Results/09-Aug-2022_11_30_08_Axial_Coronal_lumen_only/Holdout_Testing_Dataset_Results.xlsx", "Sheet", "VA_Only");
% ax_cor_rw_va = readtable("../Shape_Feature_Results/09-Aug-2022_11_30_30_Axial_Coronal_rw_only/Holdout_Testing_Dataset_Results.xlsx", "Sheet", "VA_Only");
% ax_cor_fat_va = readtable("../Shape_Feature_Results/09-Aug-2022_11_30_50_Axial_Coronal_fat_only/Holdout_Testing_Dataset_Results.xlsx", "Sheet", "VA_Only");
% ax_cor_all_va = readtable("../Shape_Feature_Results/09-Aug-2022_11_32_29_Axial_Coronal_all/Holdout_Testing_Dataset_Results.xlsx", "Sheet", "VA_Only");
% 
% % Create barplot
% % The x-axes is the dataset
% x_axes = categorical({'Training';...
%                       'Holdout Testing UH + VA';...
%                       'Holdout Testing UH Only';...
%                       'Holdout Testing VA Only'});
% x_axes = reordercats(x_axes, {'Training';...
%                       'Holdout Testing UH + VA';...
%                       'Holdout Testing UH Only';...
%                       'Holdout Testing VA Only'});
%                   
% % y-axes will be the mean metric per region grouped by dataset
% y_axes = [ax_cor_lumen_train.AUC ax_cor_rw_train.AUC ax_cor_fat_train.AUC ax_cor_all_train.AUC;...
%           ax_cor_lumen_uh_va.AUC ax_cor_rw_uh_va.AUC ax_cor_fat_uh_va.AUC ax_cor_all_uh_va.AUC;...
%           ax_cor_lumen_uh.AUC ax_cor_rw_uh.AUC ax_cor_fat_uh.AUC ax_cor_all_uh.AUC;...
%           ax_cor_lumen_va.AUC ax_cor_rw_va.AUC ax_cor_fat_va.AUC ax_cor_all_va.AUC];
%       
% figure();
% b = bar(x_axes, y_axes);
% xlabel("Dataset", 'FontSize', 18);
% ylabel("Mean AUC", 'FontSize', 18);
% ylim([0,1.05]);
% title("AUC: Axial AND Coronal Views", 'Fontsize', 20);
% b(1).FaceColor = [39/255, 174/255, 239/255];
% b(2).FaceColor = [135/255, 188/255, 69/255];
% b(3).FaceColor = [179/255, 61/255, 198/255];
% b(4).FaceColor = [239/255, 155/255, 32/255];
% 
% set(b, {'DisplayName'}, {'Lumen','Rectal Wall','Fat', 'All Regions'}')
% legend('Location', 'northeastoutside')
% 
% output_path = ("../Shape_Feature_Results/ax_cor_auc_bar_plot.png");
% saveas(gcf, output_path);

