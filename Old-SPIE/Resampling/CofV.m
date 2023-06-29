function cv = CofV(data)
% Coefficient of Variation

s = std(data(:));
m = mean(data(:));

cv = s / m;

end