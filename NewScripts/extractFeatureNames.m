top_feature_path = '/Users/leobao/Documents/MultiPlanePipeline/Data/MissingCollageResults/13-Sep-2023_16_46_13_Axial_proxfat10_only_wilcoxon_qda/Top_Features_Information.mat';
feature_column_path = '/Users/leobao/Documents/MultiPlanePipeline/Feature_Names/Texture_Feature_Names.xlsx';
feature_column_names = readtable(feature_column_path, 'ReadVariableNames',false);
feature_column_names = table2cell(feature_column_names);

top_feats = load(top_feature_path).topResults;

top_indices = top_feats.indices';

top5_feature_names = cell(5, 1);

for i = 1:5
    index = top_indices(i);
    top5_feature_names(i) = feature_column_names(index);
end

disp(top5_feature_names)

% output_path = '';
% if(~exist(output_path, "dir"))
%     mkdir(output_path);
% end
% 
% feature_file_name = string(strcat(output_path, 'top5_features.xlsx'));
% feature_table = cell2table(feature_column_names);
% writetable(feature_table, feature_file_name);