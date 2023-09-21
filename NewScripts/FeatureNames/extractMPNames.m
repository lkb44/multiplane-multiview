% Prompt the user to enter the top_feature_path
top_feature_path = input('Enter the pathname for top features: ', 's');

% Check if the input pathname is valid
if ~exist(top_feature_path, 'file')
    error('The specified file does not exist.');
end

feature_column_path = '/Users/leobao/Documents/MultiPlanePipeline/Data/MissingCollageResults/15-Sep-2023_15_22_06_Multi-Region_mrmr_lda/top20_features.xlsx';
feature_column_names = readtable(feature_column_path, 'ReadVariableNames',false);
feature_column_names = table2cell(feature_column_names);

top_feats = load(top_feature_path).topResults;

top_indices = top_feats.indices';

top5_feature_names = cell(5, 1);

for i = 1:5
    index = top_indices(i);
    if index >= 1 && index <= 5
        new_name = string(strcat('Axial_ProxFat10_', feature_column_names(index)));
    elseif index >= 6 && index <= 10
        new_name = string(strcat('Coronal_ProxFat10_', feature_column_names(index)));
    elseif index >= 11 && index <= 15
        new_name = string(strcat('Axial_Tumor_', feature_column_names(index)));
    elseif index >= 16 && index <= 20
        new_name = string(strcat('Coronal_Tumor_', feature_column_names(index)));
    end
    top5_feature_names{i} = string(new_name);
end


% Convert the cell array to a string array
top5_feature_names_str = string(top5_feature_names);

% Remove brackets and single quotes
top5_feature_names_str = replace(top5_feature_names_str, '{', '');
top5_feature_names_str = replace(top5_feature_names_str, '}', '');
top5_feature_names_str = strtrim(top5_feature_names_str);

% Display the formatted output
disp(top5_feature_names_str)