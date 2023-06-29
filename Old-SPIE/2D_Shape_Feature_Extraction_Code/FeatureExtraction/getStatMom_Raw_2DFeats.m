function statistical_moments =  getStatMom_Raw_2DFeats(folder_path, folder, regionLabel, num_cases)

statistical_moments = [];
for i=1:num_cases %move case by case
    if i<=9 u_case = ['00' num2str(i)]; end
    if i>=10 && i<=99 u_case = ['0' num2str(i)]; end
    if i>=100 u_case = num2str(i); end
    
    file_path = [folder_path regionLabel '.mat'];
    
    %file_path_diff = [folder_path folder u_case regionLabel '-diff.mat'];
    
    if exist(file_path)
        disp(['Looking into: ' file_path]);
        load(file_path);
        statistical_moments_case = [];
        %diff_matrix = [];
        
        % Finding the first slide with annotation (up to bottom)
        for ud=1:size(feature_matrix,1)
           if sum(feature_matrix(ud,:))>0 break, end  
        end
        
        % Finding the first slide with annotation (bottom to up)
        for du=size(feature_matrix,1):-1:1
           if sum(feature_matrix(du,:))>0 break, end  
        end
        
        feature_matrix = feature_matrix(ud:du,:);
        
        for f=1:size(feature_matrix,2) %move feature by feature
            mean_feature = mean(feature_matrix(:,f));
            variance_feature = std(feature_matrix(:,f));
            skewness_feature = skewness(feature_matrix(:,f));
            kurtosis_feature = kurtosis(feature_matrix(:,f));
            statistical_moments_case = [statistical_moments_case mean_feature variance_feature skewness_feature kurtosis_feature];
        end
        statistical_moments = [statistical_moments; statistical_moments_case];
    else
        disp(['Path does not exist:' file_path]);
    end
end
