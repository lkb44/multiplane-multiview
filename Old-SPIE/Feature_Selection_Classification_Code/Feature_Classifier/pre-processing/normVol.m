%% Normalized Volume Intensities
%Jacob Antunes
%November 11, 2014

%% About
%Normalize intensities across entire volume

function v_hat = normVol(v, dMin, dMax)
%v = original volume
%dMin = desired minimum
%dMax = desired maximum

%v_hat = normalized volume
vMin = min(min(min(v)));
vMax = max(max(max(v)));

v = (dMax-dMin)*(v-vMin)/(vMax-vMin)+dMin;

v_hat = v;

end
