clear
clc

view = "Coronal";
roi = "ProxFat10";
cohort = "Training";

root_path = "/Users/leobao/Documents/MultiPlanePipeline/AACR2023/";

input_path = strcat(root_path, cohort, 'TextureFeatures/Combined/', view, '/', roi, '/');

output_path = strcat(root_path, cohort, 'TextureFeatures/', view, '_', roi, '_', cohort, '.mat');

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

training_axial_fat_ids = {006, 009, 011, 016, 019, 021, 022, 025, 027, 028, 029, 031, 032, 033, 034, 036, 037, 039, 040, 041, 042, 045, 047, 048, 053, 054, 055, 056, 057, 058, 059, 060, 061, 062, 063, 064, 065, 066, 067, 068, 069, 070, 071, 072, 073, 074, 075, 076, 077, 078, 079, 080, 081, 083, 084, 085, 086, 087, 088, 089, 090, 091, 092, 094};
training_axial_proxfat5_ids = {006, 009, 011, 016, 019, 021, 022, 025, 027, 028, 029, 031, 032, 033, 034, 036, 037, 039, 040, 041, 042, 045, 047, 048, 053, 054, 055, 056, 057, 058, 059, 060, 061, 062, 063, 064, 065, 066, 067, 068, 069, 070, 071, 072, 073, 074, 075, 076, 077, 078, 079, 080, 081, 083, 084, 085, 086, 087, 088, 089, 090, 091, 092, 094};
training_axial_proxfat10_ids = {006, 009, 011, 016, 019, 021, 022, 025, 027, 028, 029, 031, 032, 033, 034, 036, 037, 039, 040, 041, 042, 045, 047, 048, 053, 054, 055, 056, 057, 058, 059, 060, 061, 062, 063, 064, 065, 066, 067, 068, 069, 070, 071, 072, 073, 074, 075, 076, 077, 078, 079, 080, 081, 083, 084, 085, 086, 087, 088, 089, 090, 091, 092, 094};
training_axial_proxfat15_ids = {006, 009, 011, 016, 019, 021, 022, 025, 027, 028, 029, 031, 032, 033, 034, 036, 037, 039, 040, 041, 042, 045, 047, 048, 053, 054, 055, 056, 057, 058, 059, 060, 061, 062, 063, 064, 065, 066, 067, 068, 069, 070, 071, 072, 073, 074, 075, 076, 077, 078, 079, 080, 081, 083, 084, 085, 086, 087, 088, 089, 090, 091, 092, 094};
training_axial_tumor_ids = {006, 009, 011, 016, 019, 021, 022, 025, 027, 028, 029, 031, 032, 033, 034, 036, 037, 039, 040, 041, 042, 045, 047, 048, 053, 054, 055, 056, 057, 058, 059, 060, 061, 062, 063, 064, 065, 066, 067, 068, 069, 070, 071, 072, 073, 074, 075, 076, 077, 078, 079, 080, 081, 083, 084, 085, 086, 087, 088, 089, 090, 091, 092, 094};
training_coronal_fat_ids = {006, 009, 016, 019, 021, 022, 025, 027, 028, 029, 031, 032, 033, 034, 036, 037, 039, 040, 041, 042, 045, 047, 048, 053, 054, 055, 056, 057, 058, 059, 060, 061, 062, 063, 064, 065, 066, 067, 068, 069, 073, 074, 075, 076, 077, 078, 079, 081, 083, 086, 087, 088, 089, 090, 091, 092, 094};
training_coronal_proxfat5_ids = {006, 009, 011, 016, 019, 021, 022, 025, 027, 028, 029, 031, 032, 033, 034, 036, 037, 039, 040, 041, 042, 045, 047, 048, 053, 054, 055, 056, 057, 058, 059, 060, 061, 062, 063, 064, 065, 066, 067, 068, 069, 070, 071, 072, 073, 074, 075, 076, 077, 078, 079, 080, 081, 083, 084, 085, 086, 087, 088, 089, 090, 091, 092, 094};
training_coronal_proxfat10_ids = {006, 009, 011, 016, 019, 021, 022, 025, 027, 028, 029, 031, 032, 033, 034, 036, 037, 039, 040, 041, 042, 045, 047, 048, 053, 054, 055, 056, 057, 058, 059, 060, 061, 062, 063, 064, 065, 066, 067, 068, 069, 070, 071, 072, 073, 074, 075, 076, 077, 078, 079, 080, 081, 083, 084, 085, 086, 087, 088, 089, 090, 091, 092, 094};
training_coronal_proxfat15_ids = {006, 009, 011, 016, 019, 021, 022, 025, 027, 028, 029, 031, 032, 033, 034, 036, 037, 039, 040, 041, 042, 045, 047, 048, 053, 054, 055, 056, 057, 058, 059, 060, 061, 062, 063, 064, 065, 066, 067, 068, 069, 070, 071, 072, 073, 074, 075, 076, 077, 078, 079, 080, 081, 083, 084, 085, 086, 087, 088, 089, 090, 091, 092, 094};
training_coronal_tumor_ids = {006, 009, 011, 016, 019, 021, 022, 025, 027, 028, 029, 031, 032, 033, 034, 036, 037, 039, 040, 041, 042, 045, 047, 048, 053, 054, 055, 056, 057, 058, 059, 060, 061, 062, 063, 064, 065, 066, 067, 068, 069, 070, 071, 072, 073, 074, 075, 076, 077, 078, 079, 080, 081, 083, 084, 085, 086, 087, 088, 089, 090, 091, 092, 094};

if strcmp(view, 'Axial')
    if strcmp(roi,'ProxFat5')
        patient_ids = training_axial_proxfat5_ids;
    elseif strcmp(roi,'ProxFat10')
        patient_ids = training_axial_proxfat10_ids;
    elseif strcmp(roi,'ProxFat15')
        patient_ids = training_axial_proxfat15_ids;
    elseif strcmp(roi,'Tumor')
        patient_ids = training_axial_tumor_ids;
    elseif strcmp(roi, 'Fat')
        patient_ids = training_axial_fat_ids;
    end
elseif strcmp(view, 'Coronal')
    if strcmp(roi,'ProxFat5')
        patient_ids = training_coronal_proxfat5_ids;
    elseif strcmp(roi,'ProxFat10')
        patient_ids = training_coronal_proxfat10_ids;
    elseif strcmp(roi,'ProxFat15')
        patient_ids = training_coronal_proxfat15_ids;
    elseif strcmp(roi,'Tumor')
        patient_ids = training_coronal_tumor_ids;
    elseif strcmp(roi, 'Fat')
        patient_ids = training_coronal_fat_ids;
    end
end

patient_ids_string = cellfun(@(x) sprintf('%03d', x), patient_ids, 'UniformOutput', false);

for i = 1 : length(patient_ids_string)

    patient_path = strcat(input_path, 'Patient-', patient_ids_string{i}, '.mat');
    features = load(patient_path);
    features = features.feature_matrix;

    feature_matrix = [feature_matrix; features];

end

save(output_path, "feature_matrix");