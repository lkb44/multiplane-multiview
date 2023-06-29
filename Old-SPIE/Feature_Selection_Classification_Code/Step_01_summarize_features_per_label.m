% Code for importing shape features, built by Jhimli Mitra (c++) 
% In this code it was read a file per label
% Code written by Charlems

% The code moves case by case to obtain a feature matrix, it uses the addFeature function

%% Initialize workspace
clear;clc;

%% Set path variables
view = "Coronal";
dataset_type = "test";
dataset = 'UH-';
path_data = strcat("../Shape_Features/", view, "/", dataset_type); % location where features are saved
path_project = strcat("../Shape_Features/"); % location where combined features will be saved

if isunix
    path = [path_data dataset '-RectalCA/']; % For Linux
    dir_output = strcat(path_project,'Feature_Matrices/', view, "/", dataset_type, "/");
    %dir_output = strjoin(dir_output);
elseif ispc
    path = [path_data dataset '-RectalCA\']; % For Windows
    dir_output = strcat(path_project,'Feature_Matrices/', view, "\", dataset_type, "\");
    %dir_output = strjoin(dir_output);
end

%% Set labels to loop through
labels = {'lumen','rw', 'fat'}; %'label_ERW-3s','label_Lumen-3s'

%% Loop through cases and import shape features
cases = [];
for ii = 1:length(labels)

    label=labels{ii};
    feature_matrix=[];
    case_index_matrix=[];
    folder = [dataset 'RectalCA-']; 

    for i=1:176 %move case by case
        if i<=9 u_case = ['00' num2str(i)]; end
        if i>=10 && i<=99 u_case = ['0' num2str(i)]; end
        if i>=100 u_case = num2str(i); end

        if isunix
            case_path = strcat(path_data, '/', folder, u_case, '/shape_feats_label_', label, '.txt'); % For Linux
            %case_path = strrep(case_path, " ","");
        elseif ispc
            case_path = [path_data folder u_case '\shape_feats_label' label '.txt']; % For Windows
        end
        
        if exist(case_path)
            cases = [cases ; u_case];
            imported_data = importdata(case_path,'\n',100);
            feature_matrix_case=[];

            for j=2:size(imported_data,1) %move line by line into each case file
                feature_matrix_case = addFeature(imported_data{j}, feature_matrix_case);        
            end 

            case_index_matrix = [case_index_matrix;[folder u_case]];
            feature_matrix = [feature_matrix;feature_matrix_case];
            
            disp('Finished case');
        else
            disp('Path does not exist');
        end
    end
    filename = strcat(view, "_", label, "_", dataset_type, '_3D_shape_features.mat'); % Specify filename for concatenated features across all patients for one label
    filename = strrep(filename, " ","");
    
    if isunix
        save(fullfile(dir_output, filename),'feature_matrix'); % Save concatenated features across all patients for one label
        feature_matrix_table = table(feature_matrix);
        filename = strrep(filename, ".mat",".xlsx");
        writetable(feature_matrix_table, fullfile(dir_output, filename));
        
    elseif ispc
        save(fullfile(dir_output, filename),'feature_matrix');
    end
end

%% Save patient index
filename = strcat(view,'_case_index_matrix.mat');
save(fullfile(dir_output,filename),'case_index_matrix');
case_index_table = table(case_index_matrix);
filename = strrep(filename, ".mat",".xlsx");
writetable(case_index_table, fullfile(dir_output, filename));