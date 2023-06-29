function statistical_moments =  getStatMom_Diff_2DFeats(folder_path, folder, regionLabel, num_cases)


statistical_moments = [];
for i=1:num_cases %move case by case - 1:130
    if i<=9 u_case = ['00' num2str(i)]; end
    if i>=10 && i<=99 u_case = ['0' num2str(i)]; end
    if i>=100 u_case = num2str(i); end
    
    file_path = [folder_path folder u_case regionLabel '.mat'];
    
    file_path_diff = [folder_path folder u_case regionLabel '-diff.mat'];
    
    if exist(file_path)
        disp(['Looking into: ' file_path]);
        load(file_path);
        statistical_moments_case = [];
        diff_matrix = [];
        
        % Finding the first slide with annotation (up to bottom)
        for ud=1:size(feature_matrix,1)
           if sum(feature_matrix(ud,:))>0 break, end  
        end
        
        % Finding the first slide with annotation (bottom to up)
        for du=size(feature_matrix,1):-1:1
           if sum(feature_matrix(du,:))>0 break, end  
        end
        
        feature_matrix = feature_matrix(ud:du,:);
        
        for s=2:size(feature_matrix,1)
            diff_slice = feature_matrix(s,:)-feature_matrix(s-1,:);
            %diff_slice = abs(feature_matrix(s,:)-feature_matrix(s-1,:));
            diff_matrix = [diff_matrix;diff_slice];
        end
        
        save(file_path_diff,'diff_matrix');
        
        %diff_matrix = feature_matrix;
        %diff_matrix = normc(diff_matrix);
        %diff_matrix = simplewhiten(diff_matrix);
        
        for f=1:size(diff_matrix,2) %move feature by feature
            mean_feature = mean(diff_matrix(:,f));
            variance_feature = std(diff_matrix(:,f));
            skewness_feature = skewness(diff_matrix(:,f));
            kurtosis_feature = kurtosis(diff_matrix(:,f));
            statistical_moments_case = [statistical_moments_case mean_feature variance_feature skewness_feature kurtosis_feature];
        end
        statistical_moments = [statistical_moments; statistical_moments_case];
    else
        disp(['Path does not exist:' file_path]);
    end
end
