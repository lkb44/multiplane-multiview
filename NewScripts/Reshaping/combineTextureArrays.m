clc;
clear;
close all;

view = "Coronal";
plane = 'cor';
roi = "ProxFat10";
region = 'proxfat10';

root_path = "/Users/leobao/Documents/MultiPlanePipeline/AACR2023/";
matlab_training_path = strcat(root_path, 'TrainingTextureFeatures/MATLAB/', view, '/', roi, '/');
% matlab_v1_path = strcat(root_path, 'Validation1TextureFeatures/MATLAB/', view, '/', roi, '/');
% matlab_v2_path = strcat(root_path, 'Validation2TextureFeatures/MATLAB/', view, '/', roi, '/');

collage_training_path = strcat(root_path, 'TrainingTextureFeatures/Collage/', view, '/', roi, '/');
% collage_v1_path = strcat(root_path, 'Validation1TextureFeatures/Collage/', view, '/', roi, '/');
% collage_v2_path = strcat(root_path, 'Validation2TextureFeatures/Collage/', view, '/', roi, '/');

training_output_path = strcat(root_path, 'TrainingTextureFeatures/Combined/', view, '/', roi, '/');
% v1_output_path = strcat(root_path, 'Validation1TextureFeatures/Combined/', view, '/', roi, '/');
% v2_output_path = strcat(root_path, 'Validation2TextureFeatures/Combined/', view, '/', roi, '/');

if(~exist(training_output_path, "dir"))
    mkdir(training_output_path);
end
% if(~exist(v1_output_path, "dir"))
%     mkdir(v1_output_path);
% end
% if(~exist(v2_output_path, "dir"))
%     mkdir(v2_output_path);
% end

training_axial_fat_ids = {006, 009, 011, 016, 019, 021, 022, 025, 027, 028, 029, 031, 032, 033, 034, 036, 037, 039, 040, 041, 042, 045, 047, 048, 053, 054, 055, 056, 057, 058, 059, 060, 061, 062, 063, 064, 065, 066, 067, 068, 069, 070, 071, 072, 073, 074, 075, 076, 077, 078, 079, 080, 081, 083, 084, 085, 086, 087, 088, 089, 090, 091, 092, 094};
% validation1_axial_fat_ids = {096, 097, 098, 099, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121};
% validation2_axial_fat_ids = {096, 097, 098, 099, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 127, 128, 129, 130, 131};

training_axial_proxfat5_ids = {006, 009, 011, 016, 019, 021, 022, 025, 027, 028, 029, 031, 032, 033, 034, 036, 037, 039, 040, 041, 042, 045, 047, 048, 053, 054, 055, 056, 057, 058, 059, 060, 061, 062, 063, 064, 065, 066, 067, 068, 069, 070, 071, 072, 073, 074, 075, 076, 077, 078, 079, 080, 081, 083, 084, 085, 086, 087, 088, 089, 090, 091, 092, 094};
% validation1_axial_proxfat5_ids = {095, 096, 097, 098, 099, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121};
% validation2_axial_proxfat5_ids = {095, 096, 097, 098, 099, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 127, 128, 129, 130, 131};

training_axial_proxfat10_ids = {006, 009, 011, 016, 019, 021, 022, 025, 027, 028, 029, 031, 032, 033, 034, 036, 037, 039, 040, 041, 042, 045, 047, 048, 053, 054, 055, 056, 057, 058, 059, 060, 061, 062, 063, 064, 065, 066, 067, 068, 069, 070, 071, 072, 073, 074, 075, 076, 077, 078, 079, 080, 081, 083, 084, 085, 086, 087, 088, 089, 090, 091, 092, 094};
% validation1_axial_proxfat10_ids = {095, 096, 097, 098, 099, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121};
% validation2_axial_proxfat10_ids = {095, 096, 097, 098, 099, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 127, 128, 129, 130, 131};

training_axial_proxfat15_ids = {006, 009, 011, 016, 019, 021, 022, 025, 027, 028, 029, 031, 032, 033, 034, 036, 037, 039, 040, 041, 042, 045, 047, 048, 053, 054, 055, 056, 057, 058, 059, 060, 061, 062, 063, 064, 065, 066, 067, 068, 069, 070, 071, 072, 073, 074, 075, 076, 077, 078, 079, 080, 081, 083, 084, 085, 086, 087, 088, 089, 090, 091, 092, 094};
% validation1_axial_proxfat15_ids = {095, 096, 097, 098, 099, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121};
% validation2_axial_proxfat15_ids = {095, 096, 097, 098, 099, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 127, 128, 129, 130, 131};

training_axial_tumor_ids = {006, 009, 011, 016, 019, 021, 022, 025, 027, 028, 029, 031, 032, 033, 034, 036, 037, 039, 040, 041, 042, 045, 047, 048, 053, 054, 055, 056, 057, 058, 059, 060, 061, 062, 063, 064, 065, 066, 067, 068, 069, 070, 071, 072, 073, 074, 075, 076, 077, 078, 079, 080, 081, 083, 084, 085, 086, 087, 088, 089, 090, 091, 092, 094};
% validation1_axial_tumor_ids = {095, 096, 097, 098, 099, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121};
% validation2_axial_tumor_ids = {095, 096, 097, 098, 099, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 127, 128, 129, 130, 131};

training_coronal_fat_ids = {006, 009, 016, 019, 021, 022, 025, 027, 028, 029, 031, 032, 033, 034, 036, 037, 039, 040, 041, 042, 045, 047, 048, 053, 054, 055, 056, 057, 058, 059, 060, 061, 062, 063, 064, 065, 066, 067, 068, 069, 073, 074, 075, 076, 077, 078, 079, 081, 083, 086, 087, 088, 089, 090, 091, 092, 094};
% validation1_coronal_fat_ids = {095, 096, 097, 098, 099, 100, 101, 102, 103, 104, 105, 107, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 121};
% validation2_coronal_fat_ids = {095, 096, 097, 098, 099, 100, 101, 102, 103, 104, 105, 107, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 121, 122, 123, 124, 125, 127, 128, 129, 130, 131};

training_coronal_proxfat5_ids = {006, 009, 011, 016, 019, 021, 022, 025, 027, 028, 029, 031, 032, 033, 034, 036, 037, 039, 040, 041, 042, 045, 047, 048, 053, 054, 055, 056, 057, 058, 059, 060, 061, 062, 063, 064, 065, 066, 067, 068, 069, 070, 071, 072, 073, 074, 075, 076, 077, 078, 079, 080, 081, 083, 084, 085, 086, 087, 088, 089, 090, 091, 092, 094};
% validation1_coronal_proxfat5_ids = {095, 096, 097, 098, 099, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121};
% validation2_coronal_proxfat5_ids = {095, 096, 097, 098, 099, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 127, 128, 129, 130, 131};

training_coronal_proxfat10_ids = {006, 009, 011, 016, 019, 021, 022, 025, 027, 028, 029, 031, 032, 033, 034, 036, 037, 039, 040, 041, 042, 045, 047, 048, 053, 054, 055, 056, 057, 058, 059, 060, 061, 062, 063, 064, 065, 066, 067, 068, 069, 070, 071, 072, 073, 074, 075, 076, 077, 078, 079, 080, 081, 083, 084, 085, 086, 087, 088, 089, 090, 091, 092, 094};
% validation1_coronal_proxfat10_ids = {095, 096, 097, 098, 099, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121};
% validation2_coronal_proxfat10_ids = {095, 096, 097, 098, 099, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 127, 128, 129, 130, 131};

training_coronal_proxfat15_ids = {006, 009, 011, 016, 019, 021, 022, 025, 027, 028, 029, 031, 032, 033, 034, 036, 037, 039, 040, 041, 042, 045, 047, 048, 053, 054, 055, 056, 057, 058, 059, 060, 061, 062, 063, 064, 065, 066, 067, 068, 069, 070, 071, 072, 073, 074, 075, 076, 077, 078, 079, 080, 081, 083, 084, 085, 086, 087, 088, 089, 090, 091, 092, 094};
% validation1_coronal_proxfat15_ids = {095, 096, 097, 098, 099, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121};
% validation2_coronal_proxfat15_ids = {095, 096, 097, 098, 099, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 127, 128, 129, 130, 131};

training_coronal_tumor_ids = {006, 009, 011, 016, 019, 021, 022, 025, 027, 028, 029, 031, 032, 033, 034, 036, 037, 039, 040, 041, 042, 045, 047, 048, 053, 054, 055, 056, 057, 058, 059, 060, 061, 062, 063, 064, 065, 066, 067, 068, 069, 070, 071, 072, 073, 074, 075, 076, 077, 078, 079, 080, 081, 083, 084, 085, 086, 087, 088, 089, 090, 091, 092, 094};
% validation1_coronal_tumor_ids = {095, 096, 097, 098, 099, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121};
% validation2_coronal_tumor_ids = {095, 096, 097, 098, 099, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 127, 128, 129, 130, 131};

if strcmp(plane, 'ax')
    if strcmp(region,'proxfat5')
        training_patient_ids = training_axial_proxfat5_ids;
        % v1_patient_ids = validation1_axial_proxfat5_ids;
        % v2_patient_ids = validation2_axial_proxfat5_ids;
    elseif strcmp(region,'proxfat10')
        training_patient_ids = training_axial_proxfat10_ids;
        % v1_patient_ids = validation1_axial_proxfat10_ids;
        % v2_patient_ids = validation2_axial_proxfat10_ids;
    elseif strcmp(region,'proxfat15')
        training_patient_ids = training_axial_proxfat15_ids;
        % v1_patient_ids = validation1_axial_proxfat15_ids;
        % v2_patient_ids = validation2_axial_proxfat15_ids;
    elseif strcmp(region,'tumor')
        training_patient_ids = training_axial_tumor_ids;
        % v1_patient_ids = validation1_axial_tumor_ids;
        % v2_patient_ids = validation2_axial_tumor_ids;
    elseif strcmp(region, 'fat')
        training_patient_ids = training_axial_fat_ids;
        % v1_patient_ids = validation1_axial_fat_ids;
        % v2_patient_ids = validation2_axial_fat_ids;
    end
elseif strcmp(plane, 'cor')
    if strcmp(region,'proxfat5')
        training_patient_ids = training_coronal_proxfat5_ids;
        % v1_patient_ids = validation1_coronal_proxfat5_ids;
        % v2_patient_ids = validation2_coronal_proxfat5_ids;
    elseif strcmp(region,'proxfat10')
        training_patient_ids = training_coronal_proxfat10_ids;
        % v1_patient_ids = validation1_coronal_proxfat10_ids;
        % v2_patient_ids = validation2_coronal_proxfat10_ids;
    elseif strcmp(region,'proxfat15')
        training_patient_ids = training_coronal_proxfat15_ids;
        % v1_patient_ids = validation1_coronal_proxfat15_ids;
        % v2_patient_ids = validation2_coronal_proxfat15_ids;
    elseif strcmp(region,'tumor')
        training_patient_ids = training_coronal_tumor_ids;
        % v1_patient_ids = validation1_coronal_tumor_ids;
        % v2_patient_ids = validation2_coronal_tumor_ids;
    elseif strcmp(region, 'fat')
        training_patient_ids = training_coronal_fat_ids;
        % v1_patient_ids = validation1_coronal_fat_ids;
        % v2_patient_ids = validation2_coronal_fat_ids;
    end
end


% Convert the cell array to a string array (optional)

training_patient_ids_string = cellfun(@(x) sprintf('%03d', x), training_patient_ids, 'UniformOutput', false);
% v1_patient_ids_string = cellfun(@(x) sprintf('%03d', x), v1_patient_ids, 'UniformOutput', false);
% v2_patient_ids_string = cellfun(@(x) sprintf('%03d', x), v2_patient_ids, 'UniformOutput', false);

for i = 1 : length(training_patient_ids_string)
    matlab_training_patient_path = strcat(matlab_training_path, 'Patient-', training_patient_ids_string{i}, '_', region, '.mat');
    python_training_patient_path = strcat(collage_training_path, 'Feats_Col_Patient-', training_patient_ids_string{i}, '_', plane, '.mat');
    output_patient_path = strcat(training_output_path, 'Patient-', training_patient_ids_string{i}, '.mat');

    matlab_features = load(matlab_training_patient_path);
    python_features = load(python_training_patient_path);

    feature_matrix = [];
    matlab = matlab_features.feature_matrix;
    python = python_features.array_data;

    feature_matrix = [matlab python];

    save(output_patient_path, "feature_matrix")
end
% 
% for i = 1 : length(v1_patient_ids_string)
%     fprintf('Combining Patient-%d\n', v1_patient_ids_string{i});
%     matlab_v1_patient_path = strcat(matlab_v1_path, 'Patient-', v1_patient_ids_string{i}, '_', region, '.mat');
%     python_v1_patient_path = strcat(collage_v1_path, 'Feats_Col_Patient-', v1_patient_ids_string{i}, '_', plane, '.mat');
% 
%     output_patient_path = strcat(v1_output_path, 'Patient-', v1_patient_ids_string{i}, '.mat');
% 
%     matlab_features = load(matlab_v1_patient_path);
%     python_features = load(python_v1_patient_path);
% 
%     feature_matrix = [];
%     matlab = matlab_features.feature_matrix;
%     python = python_features.array_data;
% 
%     feature_matrix = [matlab python];
% 
%     save(output_patient_path, "feature_matrix")
% end
% 
% for i = 1 : length(v2_patient_ids_string)
%     fprintf('Combining Patient-%d\n', v2_patient_ids_string{i});
%     matlab_v2_patient_path = strcat(matlab_v2_path, 'Patient-', v2_patient_ids_string{i}, '_', region, '.mat');
%     python_v2_patient_path = strcat(collage_v2_path, 'Feats_Col_Patient-', v2_patient_ids_string{i}, '_', plane, '.mat');
% 
%     output_patient_path = strcat(v2_output_path, 'Patient-', v2_patient_ids_string{i}, '.mat');
% 
%     matlab_features = load(matlab_v2_patient_path);
%     python_features = load(python_v2_patient_path);
% 
%     feature_matrix = [];
%     matlab = matlab_features.feature_matrix;
%     python = python_features.array_data;
% 
%     feature_matrix = [matlab python];
% 
%     save(output_patient_path, "feature_matrix")
% end
