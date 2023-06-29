close all; clear; clc;
%% Axial View, Multiple Regions
close all; clear; clc;

% Axial Training Metrics
axial_tumor_train = readtable("/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/Data/AxialTrainingShapeFeatures/15-Feb-2023_22_31_15_Axial_tumor_only_wilcoxon_qda/CV_Training_Results.xlsx");
ax_proxfat_train = readtable("/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/Data/AxialTrainingShapeFeatures/15-Feb-2023_22_36_59_Axial_proxfat_only_wilcoxon_qda/CV_Training_Results.xlsx");
ax_fat_train = readtable("/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/Data/AxialTrainingShapeFeatures/15-Feb-2023_22_29_19_Axial_fat_only_wilcoxon_qda/CV_Training_Results.xlsx");
ax_all_train = readtable("/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/Data/16-Feb-2023_12_29_05_all_wilcoxon_qda_withPruningAndOversampling/CV_Training_Results.xlsx");

% Axial Holdout Testing UH + VA Metrics
axial_tumor_test = readtable("/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/Data/AxialTrainingShapeFeatures/15-Feb-2023_22_31_15_Axial_tumor_only_wilcoxon_qda/Holdout_Testing_Dataset_Results.xlsx", "Sheet", "Entire_Test_Cohort_Results");
axial_proxfat_test = readtable("/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/Data/AxialTrainingShapeFeatures/15-Feb-2023_22_36_59_Axial_proxfat_only_wilcoxon_qda/Holdout_Testing_Dataset_Results.xlsx", "Sheet", "Entire_Test_Cohort_Results");
axial_fat_test = readtable("/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/Data/AxialTrainingShapeFeatures/15-Feb-2023_22_29_19_Axial_fat_only_wilcoxon_qda/Holdout_Testing_Dataset_Results.xlsx", "Sheet", "Entire_Test_Cohort_Results");
axial_all_test = readtable("/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/Data/16-Feb-2023_12_29_05_all_wilcoxon_qda_withPruningAndOversampling/Holdout_Testing_Dataset_Results.xlsx", "Sheet", "Entire_Test_Cohort_Results");

% create barplot
% The x-axes is the dataset
x_axes = categorical({'Training'; ...
                        'Holdout Testing'});
x_axes = reordercats(x_axes, {'Training';...
                      'Holdout Testing'});
% x_axes = reordercats(x_axes, {'Training', 'Holdout Testing'});
    
% y-axes will be the mean metric per region grouped by dataset
y_axes = [axial_tumor_train.AUC ax_proxfat_train.AUC ax_fat_train.AUC ax_all_train.AUC;...
        axial_tumor_test.AUC axial_proxfat_test.AUC axial_fat_test.AUC axial_all_test.AUC];
      
figure();
set(gca,'FontSize',18)
b = bar(x_axes, y_axes);
xlabel("Dataset", 'FontSize', 18);
ylabel("Mean AUC", 'FontSize', 18);
ylim([0,1.05]);
title("AUC: Axial View", 'Fontsize', 20);
b(1).FaceColor = [39/255, 174/255, 239/255];
b(2).FaceColor = [135/255, 188/255, 69/255];
b(3).FaceColor = [179/255, 61/255, 198/255];
b(4).FaceColor = [239/255, 155/255, 32/255];

set(b, {'DisplayName'}, {'Proximal Fat','Tumor','Fat', 'All Regions'}')
legend('Location', 'northeastoutside')

output_path = ("/Users/leobao/Documents/MultiPlanePipeline/2023-SPIE/Data/AxialTrainingShapeFeatures/ax_auc_bar_plot.png");
saveas(gcf, output_path);

%% Coronal View, Multiple Regions

% Coronal Training Metrics
coronal_lumen_train = readtable("../Shape_Feature_Results/09-Aug-2022_10_23_25_Coronal_lumen_only/CV_Training_Results.xlsx");
coronal_rw_train = readtable("../Shape_Feature_Results/09-Aug-2022_10_23_43_Coronal_rw_only/CV_Training_Results.xlsx");
coronal_fat_train = readtable("../Shape_Feature_Results/09-Aug-2022_10_23_55_Coronal_fat_only/CV_Training_Results.xlsx");
coronal_all_train = readtable("../Shape_Feature_Results/09-Aug-2022_10_46_38_Coronal_all/CV_Training_Results.xlsx");
 
% Coronal Holdout Testing UH + VA Metrics
coronal_lumen_uh_va = readtable("../Shape_Feature_Results/09-Aug-2022_10_23_25_Coronal_lumen_only/Holdout_Testing_Dataset_Results.xlsx", "Sheet", "Entire_Test_Cohort_Results");
coronal_rw_uh_va = readtable("../Shape_Feature_Results/09-Aug-2022_10_23_43_Coronal_rw_only/Holdout_Testing_Dataset_Results.xlsx", "Sheet", "Entire_Test_Cohort_Results");
coronal_fat_uh_va = readtable("../Shape_Feature_Results/09-Aug-2022_10_23_55_Coronal_fat_only/Holdout_Testing_Dataset_Results.xlsx", "Sheet", "Entire_Test_Cohort_Results");
coronal_all_uh_va = readtable("../Shape_Feature_Results/09-Aug-2022_10_46_38_Coronal_all/Holdout_Testing_Dataset_Results.xlsx", "Sheet", "Entire_Test_Cohort_Results");

% Coronal Holdout Testing UH Only Metrics
coronal_lumen_uh = readtable("../Shape_Feature_Results/09-Aug-2022_10_23_25_Coronal_lumen_only/Holdout_Testing_Dataset_Results.xlsx", "Sheet", "UH_Only");
coronal_rw_uh = readtable("../Shape_Feature_Results/09-Aug-2022_10_23_43_Coronal_rw_only/Holdout_Testing_Dataset_Results.xlsx", "Sheet", "UH_Only");
coronal_fat_uh = readtable("../Shape_Feature_Results/09-Aug-2022_10_23_55_Coronal_fat_only/Holdout_Testing_Dataset_Results.xlsx", "Sheet", "UH_Only");
coronal_all_uh = readtable("../Shape_Feature_Results/09-Aug-2022_10_46_38_Coronal_all/Holdout_Testing_Dataset_Results.xlsx", "Sheet", "UH_Only");

% Coronal Holdout Testing VA Only Metrics
coronal_lumen_va = readtable("../Shape_Feature_Results/09-Aug-2022_10_23_25_Coronal_lumen_only/Holdout_Testing_Dataset_Results.xlsx", "Sheet", "VA_Only");
coronal_rw_va = readtable("../Shape_Feature_Results/09-Aug-2022_10_23_43_Coronal_rw_only/Holdout_Testing_Dataset_Results.xlsx", "Sheet", "VA_Only");
coronal_fat_va = readtable("../Shape_Feature_Results/09-Aug-2022_10_23_55_Coronal_fat_only/Holdout_Testing_Dataset_Results.xlsx", "Sheet", "VA_Only");
coronal_all_va = readtable("../Shape_Feature_Results/09-Aug-2022_10_46_38_Coronal_all/Holdout_Testing_Dataset_Results.xlsx", "Sheet", "VA_Only");

% Create barplot
% The x-axes is the dataset
x_axes = categorical({'Training';...
                      'Holdout Testing UH + VA';...
                      'Holdout Testing UH Only';...
                      'Holdout Testing VA Only'});
x_axes = reordercats(x_axes, {'Training';...
                      'Holdout Testing UH + VA';...
                      'Holdout Testing UH Only';...
                      'Holdout Testing VA Only'});
                  
% y-axes will be the mean metric per region grouped by dataset
y_axes = [coronal_lumen_train.AUC coronal_rw_train.AUC coronal_fat_train.AUC coronal_all_train.AUC;...
          coronal_lumen_uh_va.AUC coronal_rw_uh_va.AUC coronal_fat_uh_va.AUC coronal_all_uh_va.AUC;...
          coronal_lumen_uh.AUC coronal_rw_uh.AUC coronal_fat_uh.AUC coronal_all_uh.AUC;...
          coronal_lumen_va.AUC coronal_rw_va.AUC coronal_fat_va.AUC coronal_all_va.AUC];

figure();
set(gca,'FontSize',18)
b = bar(x_axes, y_axes);
xlabel("Dataset", 'FontSize', 18);
ylabel("Mean AUC", 'FontSize', 18);
ylim([0,1.05]);
title("AUC: Coronal View", 'Fontsize', 20);
b(1).FaceColor = [39/255, 174/255, 239/255];
b(2).FaceColor = [135/255, 188/255, 69/255];
b(3).FaceColor = [179/255, 61/255, 198/255];
b(4).FaceColor = [239/255, 155/255, 32/255];

set(b, {'DisplayName'}, {'Lumen','Rectal Wall','Fat', 'All Regions'}')
legend('Location', 'northeastoutside');

output_path = ("../Shape_Feature_Results/cor_auc_bar_plot.png");
saveas(gcf, output_path);

%% Axial and Coronal View, Multiple Regions
% Coronal Training Metrics
ax_cor_lumen_train = readtable("../Shape_Feature_Results/09-Aug-2022_11_30_08_Axial_Coronal_lumen_only/CV_Training_Results.xlsx");
ax_cor_rw_train = readtable("../Shape_Feature_Results/09-Aug-2022_11_30_30_Axial_Coronal_rw_only/CV_Training_Results.xlsx");
ax_cor_fat_train = readtable("../Shape_Feature_Results/09-Aug-2022_11_30_50_Axial_Coronal_fat_only/CV_Training_Results.xlsx");
ax_cor_all_train = readtable("../Shape_Feature_Results/09-Aug-2022_11_32_29_Axial_Coronal_all/CV_Training_Results.xlsx");
 
% Coronal Holdout Testing UH + VA Metrics
ax_cor_lumen_uh_va = readtable("../Shape_Feature_Results/09-Aug-2022_11_30_08_Axial_Coronal_lumen_only/Holdout_Testing_Dataset_Results.xlsx", "Sheet", "Entire_Test_Cohort_Results");
ax_cor_rw_uh_va = readtable("../Shape_Feature_Results/09-Aug-2022_11_30_30_Axial_Coronal_rw_only/Holdout_Testing_Dataset_Results.xlsx", "Sheet", "Entire_Test_Cohort_Results");
ax_cor_fat_uh_va = readtable("../Shape_Feature_Results/09-Aug-2022_11_30_50_Axial_Coronal_fat_only/Holdout_Testing_Dataset_Results.xlsx", "Sheet", "Entire_Test_Cohort_Results");
ax_cor_all_uh_va = readtable("../Shape_Feature_Results/09-Aug-2022_11_32_29_Axial_Coronal_all/Holdout_Testing_Dataset_Results.xlsx", "Sheet", "Entire_Test_Cohort_Results");
 
% Coronal Holdout Testing UH Only Metrics
ax_cor_lumen_uh = readtable("../Shape_Feature_Results/09-Aug-2022_11_30_08_Axial_Coronal_lumen_only/Holdout_Testing_Dataset_Results.xlsx", "Sheet", "UH_Only");
ax_cor_rw_uh = readtable("../Shape_Feature_Results/09-Aug-2022_11_30_30_Axial_Coronal_rw_only/Holdout_Testing_Dataset_Results.xlsx", "Sheet", "UH_Only");
ax_cor_fat_uh = readtable("../Shape_Feature_Results/09-Aug-2022_11_30_50_Axial_Coronal_fat_only/Holdout_Testing_Dataset_Results.xlsx", "Sheet", "UH_Only");
ax_cor_all_uh = readtable("../Shape_Feature_Results/09-Aug-2022_11_32_29_Axial_Coronal_all/Holdout_Testing_Dataset_Results.xlsx", "Sheet", "UH_Only");
 
% Coronal Holdout Testing VA Only Metrics
ax_cor_lumen_va = readtable("../Shape_Feature_Results/09-Aug-2022_11_30_08_Axial_Coronal_lumen_only/Holdout_Testing_Dataset_Results.xlsx", "Sheet", "VA_Only");
ax_cor_rw_va = readtable("../Shape_Feature_Results/09-Aug-2022_11_30_30_Axial_Coronal_rw_only/Holdout_Testing_Dataset_Results.xlsx", "Sheet", "VA_Only");
ax_cor_fat_va = readtable("../Shape_Feature_Results/09-Aug-2022_11_30_50_Axial_Coronal_fat_only/Holdout_Testing_Dataset_Results.xlsx", "Sheet", "VA_Only");
ax_cor_all_va = readtable("../Shape_Feature_Results/09-Aug-2022_11_32_29_Axial_Coronal_all/Holdout_Testing_Dataset_Results.xlsx", "Sheet", "VA_Only");

% Create barplot
% The x-axes is the dataset
x_axes = categorical({'Training';...
                      'Holdout Testing UH + VA';...
                      'Holdout Testing UH Only';...
                      'Holdout Testing VA Only'});
x_axes = reordercats(x_axes, {'Training';...
                      'Holdout Testing UH + VA';...
                      'Holdout Testing UH Only';...
                      'Holdout Testing VA Only'});
                  
% y-axes will be the mean metric per region grouped by dataset
y_axes = [ax_cor_lumen_train.AUC ax_cor_rw_train.AUC ax_cor_fat_train.AUC ax_cor_all_train.AUC;...
          ax_cor_lumen_uh_va.AUC ax_cor_rw_uh_va.AUC ax_cor_fat_uh_va.AUC ax_cor_all_uh_va.AUC;...
          ax_cor_lumen_uh.AUC ax_cor_rw_uh.AUC ax_cor_fat_uh.AUC ax_cor_all_uh.AUC;...
          ax_cor_lumen_va.AUC ax_cor_rw_va.AUC ax_cor_fat_va.AUC ax_cor_all_va.AUC];
      
figure();
b = bar(x_axes, y_axes);
xlabel("Dataset", 'FontSize', 18);
ylabel("Mean AUC", 'FontSize', 18);
ylim([0,1.05]);
title("AUC: Axial AND Coronal Views", 'Fontsize', 20);
b(1).FaceColor = [39/255, 174/255, 239/255];
b(2).FaceColor = [135/255, 188/255, 69/255];
b(3).FaceColor = [179/255, 61/255, 198/255];
b(4).FaceColor = [239/255, 155/255, 32/255];

set(b, {'DisplayName'}, {'Lumen','Rectal Wall','Fat', 'All Regions'}')
legend('Location', 'northeastoutside')

output_path = ("../Shape_Feature_Results/ax_cor_auc_bar_plot.png");
saveas(gcf, output_path);

