close all; clear; clc;
%% Axial View, Multiple Regions
close all; clear; clc;

% Axial Training Metrics
axial_lumen_train = readtable("../Shape_Feature_Results/09-Aug-2022_00_27_39_Axial_lumen_only/CV_Training_Results.xlsx");
ax_rw_train = readtable("../Shape_Feature_Results/09-Aug-2022_00_27_49_Axial_rw_only/CV_Training_Results.xlsx");
ax_fat_train = readtable("../Shape_Feature_Results/09-Aug-2022_00_27_58_Axial_fat_only/CV_Training_Results.xlsx");
ax_all_train = readtable("../Shape_Feature_Results/09-Aug-2022_00_29_08_Axial_all/CV_Training_Results.xlsx");

% Axial Holdout Testing UH + VA Metrics
axial_lumen_uh_va = readtable("../Shape_Feature_Results/09-Aug-2022_00_27_39_Axial_lumen_only/Holdout_Testing_Dataset_Results.xlsx", "Sheet", "Entire_Test_Cohort_Results");
axial_rw_uh_va = readtable("../Shape_Feature_Results/09-Aug-2022_00_27_49_Axial_rw_only/Holdout_Testing_Dataset_Results.xlsx", "Sheet", "Entire_Test_Cohort_Results");
axial_fat_uh_va = readtable("../Shape_Feature_Results/09-Aug-2022_00_27_58_Axial_fat_only/Holdout_Testing_Dataset_Results.xlsx", "Sheet", "Entire_Test_Cohort_Results");
axial_all_uh_va = readtable("../Shape_Feature_Results/09-Aug-2022_00_29_08_Axial_all/Holdout_Testing_Dataset_Results.xlsx", "Sheet", "Entire_Test_Cohort_Results");

% Axial Holdout Testing UH Only Metrics
axial_lumen_uh = readtable("../Shape_Feature_Results/09-Aug-2022_00_27_39_Axial_lumen_only/Holdout_Testing_Dataset_Results.xlsx", "Sheet", "UH_Only");
axial_rw_uh = readtable("../Shape_Feature_Results/09-Aug-2022_00_27_49_Axial_rw_only/Holdout_Testing_Dataset_Results.xlsx", "Sheet", "UH_Only");
axial_fat_uh = readtable("../Shape_Feature_Results/09-Aug-2022_00_27_58_Axial_fat_only/Holdout_Testing_Dataset_Results.xlsx", "Sheet", "UH_Only");
axial_all_uh = readtable("../Shape_Feature_Results/09-Aug-2022_00_29_08_Axial_all/Holdout_Testing_Dataset_Results.xlsx", "Sheet", "UH_Only");

% Axial Holdout Testing VA Only Metrics
axial_lumen_va = readtable("../Shape_Feature_Results/09-Aug-2022_00_27_39_Axial_lumen_only/Holdout_Testing_Dataset_Results.xlsx", "Sheet", "VA_Only");
axial_rw_va = readtable("../Shape_Feature_Results/09-Aug-2022_00_27_49_Axial_rw_only/Holdout_Testing_Dataset_Results.xlsx", "Sheet", "VA_Only");
axial_fat_va = readtable("../Shape_Feature_Results/09-Aug-2022_00_27_58_Axial_fat_only/Holdout_Testing_Dataset_Results.xlsx", "Sheet", "VA_Only");
axial_all_va = readtable("../Shape_Feature_Results/09-Aug-2022_00_29_08_Axial_all/Holdout_Testing_Dataset_Results.xlsx", "Sheet", "VA_Only");

% create barplot
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
y_axes = [axial_lumen_train.AUC ax_rw_train.AUC ax_fat_train.AUC ax_all_train.AUC;...
          axial_lumen_uh_va.AUC axial_rw_uh_va.AUC axial_fat_uh_va.AUC axial_all_uh_va.AUC;...
          axial_lumen_uh.AUC axial_rw_uh.AUC axial_fat_uh.AUC axial_all_uh.AUC;...
          axial_lumen_va.AUC axial_rw_va.AUC axial_fat_va.AUC axial_all_va.AUC];
      
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

set(b, {'DisplayName'}, {'Lumen','Rectal Wall','Fat', 'All Regions'}')
legend('Location', 'northeastoutside')

output_path = ("../Shape_Feature_Results/ax_auc_bar_plot.png");
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