% %%
% clear params;
% 
% load fisheriris;
% 
% data_set = meas;
% data_labels = (contains(species,'virginica'));
% params.classifier='LDA';
% params.fsname='wilcoxon';
% params.shuffle = 1;
% params.n = 3;
% params.nIter = 25;
% params.num_top_feats = 2;
% % params.remouts = 'quartile';
% stats = nFoldCV_withFS_v3(data_set,data_labels,params);

%%
% 

addpath(genpath('..'));

close all;
noise = 0.0; %a percentage, decimal form
% data = clusterincluster([],[],[],[],[],[],noise);
data = planes([],[],noise);
% data = halfkernel([],[],[],[],[],[],noise);
% data = outlier([],[],[],[],[],noise);


% data_set = [rand(size(data,1),33) data(:,1) rand(size(data,1),33) data(:,2) rand(size(data,1),34)]; %stick the best features in "semi random" spots: 34 and 68
data_set = data([1:50 451:500],1);
data_set = simplewhiten(data_set);
data_labels = data([1:50 451:500],3);
params.classifier='LDA';
params.svm_kernel='gaussian';
params.fsname='wilcoxon';
params.shuffle = 1;
params.n = length(data_labels);
params.nIter = 1;
params.num_top_feats = 1;
% params.remouts = 'quartile';

%
stats = nFoldCV_withFS_v3(data_set,data_labels,params);

%
medAUC = median(cat(1,stats.AUCs));    medAUC = medAUC(1);
stdAUC = std(cat(1,stats.AUCs));       stdAUC = stdAUC(1);

[ordered_count, ordered_featInds] = count_topfeatures(cat(1,stats.topfeatinds),1:size(data_set,2),'descend',1);
ordered_freq = ordered_count./(params.n*params.nIter); %proportion of times each feature was selected
f1 = ordered_freq(ordered_featInds==34);
f2 = ordered_freq(ordered_featInds==68);


%%
figure;
dotsize = 12;
scatter(data(:,1), data(:,2), dotsize, data(:,3)); axis equal; hold on;
colormap([1 0 .5;   % magenta
           0 0 .8;   % blue
           0 .6 0;   % dark green
           .3 1 0]); % bright green
title({'Synthetic Data (2 clear feats mixed w/ 100 noisy feats)',...
    sprintf('%s feature selection -- picking top %i features',params.fsname, params.num_top_feats),...
    sprintf('%s classifier 25 runs 3foldCV',params.classifier),...
    sprintf('medAUC = %.3f +/- %.3f',medAUC, stdAUC)});
if isempty(f1)
    xlabel(sprintf('Feature 1 - picked %i%% of time',0*100));
else
    xlabel(sprintf('Feature 1 - picked %i%% of time',f1*100));
end
if isempty(f2)
    ylabel(sprintf('Feature 2 - picked %i%% of time',0*100));
else
    ylabel(sprintf('Feature 2 - picked %i%% of time',f2*100));
end
