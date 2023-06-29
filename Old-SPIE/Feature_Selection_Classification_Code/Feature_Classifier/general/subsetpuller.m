subsets=struct;

for l=1:50
    subsets(l).training=stats_TI_full(l).subsets.training;
    subsets(l).testing=stats_TI_full(l).subsets.testing;
end