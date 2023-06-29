function topfeats = lasso_fs_binarized(X,y,n)
%INPUTS
%   X: m x f data of m observations and f features
%   y: corresponding class labels (binarized; i.e. 0 and 1)
%   n: number of top features to return
%OUTPUTS
%   topfeats: idxs (i.e. which columns of X) of top n predictors

%% Check inputs
%checking size
if size(X,1) ~= length(y)
    error('[lasso_fs_binarized]: num rows of data in "X" must equal num elements of labels in "y"');
end

if size(X,2) < n
    error('[lasso_fs_binarized]: num cols of data in "X" is less than num predictors, "n" requested');
end

%lasso can't handle any NAN features
idx_pool = 1:size(X,2);
idxs2rem = [];
[~,c] = ind2sub(size(X),find(isnan(X)));
idxs2rem = unique(c);
X(:,idxs2rem) = []; 
idx_pool(idxs2rem) = [];

%% LASSO regression
[B,FitInfo] = lassoglm(X,y,'binomial');

%optional: plot!
% figure;
% lassoPlot(B,FitInfo,'plottype','lambda','Xscale','log'); %not all that informative...
% legend('show') % Show legend

%Identify lambda at which minimum deviance occurs
% deviance = diff(FitInfo.Deviance); <-- doesnt work as expected
% [MinDeviance,idxLambdaMinDeviance] = min(deviance);

%Identify degrees of freedom at each lambda (i.e. #predictors with nonzero coefficients 
df = FitInfo.DF;
    
%identify the first lambda to contain the number of features requested
MinIdx = find(df<=n,1);

%Get feature indices at this lambda
if df(MinIdx) < n %in case we are actually less than n top features, go back 1 lambda and choose top n from sorted lit
   MinIdx = MinIdx - 1;
   candidateIdx = find(B(:,MinIdx)~=0);
   [~,sortIdx] = sort(B(candidateIdx,MinIdx),'descend');
   TopIdx = candidateIdx(sortIdx(1:n));
else
   TopIdx = find(B(:,MinIdx)~=0); 
end

%map back to original feature pool
topfeats = idx_pool(TopIdx);

%% LASSO regression - CV %is this more appropriate?
% [B_CV,FitInfo_CV] = lassoglm(X,y,'binomial','CV',3); 
% 
% %optional: plot!
% figure;
% lassoPlot(B_CV,FitInfo_CV,'plottype','CV'); %plot certainly looks much prettier
% legend('show') % Show legend
% 
% %Identify lambda at which minimum deviance occurs
% idxLambdaMinDeviance = FitInfo_CV.IndexMinDeviance;


end