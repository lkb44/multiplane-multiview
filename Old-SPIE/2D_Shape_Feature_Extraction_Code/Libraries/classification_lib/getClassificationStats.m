function stats = getClassificationStats(predictions,truelabels,threshold,targetclass)

stats.truelabels= truelabels;
stats.threshold= threshold;
stats.prediction = predictions;
stats.decision = stats.prediction >= stats.threshold;
stats.tp = length(find(stats.decision(truelabels==targetclass)==targetclass)); 
stats.tn = length(find(stats.decision(truelabels~=targetclass)~=targetclass)); 
stats.fp = length(find(stats.decision(truelabels~=targetclass)==targetclass)); 
stats.fn = length(find(stats.decision(truelabels==targetclass)~=targetclass)); 
stats.acc = (stats.tp+stats.tn)/(stats.tp+stats.tn+stats.fp+stats.fn);
stats.ppv = stats.tp/(stats.tp+stats.fp);
stats.sens = stats.tp/(stats.tp+stats.fn);
stats.spec = stats.tn/(stats.fp+stats.tn);
Pre = ((stats.tp+stats.fp)*(stats.tp+stats.fn) + (stats.tn+stats.fn)*(stats.tn+stats.fp)) / (stats.tp+stats.tn+stats.fp+stats.fn)^2;
stats.kappa = (stats.acc - Pre) / (1 - Pre);        
stats.Fscore = 2*stats.tp/(2*stats.tp+stats.fp+stats.fn);

end