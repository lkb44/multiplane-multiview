function visualize_feature_correlation(data)

% INPUTS 
% - data : M x N matrix of M observations and N features 
% OUTPUTS
% <none>
% just a heatmap visualization of correaltion between each feature across all observations

R = corrcoef(data(:,1:660));
R = abs(R);
% Ru = triu(R);
% Ru(tril(R)>0) = nan;
figure('Color','white');
imagesc(R); 
cmap = redblue; %myred? %myjet?
cmap(1,:) = [1 1 1];colormap(cmap);colorbar;set(gca,'Visible','off')