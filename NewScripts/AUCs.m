clear
clc

MP_data = load("/Users/leobao/Documents/MultiPlanePipeline/AACR2023/Results_Multi/12-Nov-2023_03_18_22_multiregion_mrmr_lda/CV_Feature_Results.mat").feat_stats;
T_data = load("/Users/leobao/Documents/MultiPlanePipeline/AACR2023/Results/12-Nov-2023_02_47_58_Axial_tumor_only_mrmr_lda/CV_Feature_Results.mat").feat_stats;
PT_data = load("/Users/leobao/Documents/MultiPlanePipeline/AACR2023/Results/12-Nov-2023_02_46_49_Coronal_proxfat10_only_mrmr_lda/CV_Feature_Results.mat").feat_stats;

for i = 1:numel(T_data)
    aucValues = T_data(i).AUCs;
    meanAUC = mean(aucValues);
    T_data(i).AUCs = meanAUC;
end

for i = 1:numel(MP_data)
    aucValues = MP_data(i).AUCs;
    meanAUC = mean(aucValues);
    MP_data(i).AUCs = meanAUC;
end

for i = 1:numel(PT_data)
    aucValues = PT_data(i).AUCs;
    meanAUC = mean(aucValues);
    PT_data(i).AUCs = meanAUC;
end