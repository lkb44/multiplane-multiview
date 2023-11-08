matlab = load('/Users/leobao/Documents/MultiPlanePipeline/AACR2023/TrainingTextureFeatures/MATLAB/Axial/Fat/Patient-021_fat.mat');
features = matlab.feature_matrix_fat;
output_patient_path = '/Users/leobao/Documents/MultiPlanePipeline/AACR2023/TrainingTextureFeatures/MATLAB/Axial/Fat/Patient-021_fat.mat';


feature_matrix = features;

save(output_patient_path, "feature_matrix")