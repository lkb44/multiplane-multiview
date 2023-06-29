function [threshold, idx] = determine_threshold_euclidean_distance(varargin)
%% Determines best threshold using Euclidean Distance Formula
% Neil Sehgal @ 2018
% POSSIBLE INPUTS -- the outputs from perfcurve
%       1- X
%       2- Y
%       3- T

if nargin == 3
    x = varargin{1};
    y = varargin{2};
    y = y(:, 1);
    t = varargin{3};
    t = t(:, 1);
end



% x = stats.X;
% y = stats.Y(:,1);
% t = stats.T(:,1);
% figure;
% plot(x, y);
% hold on;

% synthetic data for testing
% x = 0:1/51:1;
% y = 0:1/51:1;
% y(10) = .99;
% figure;
% plot(x, y);

currDistance = 111;
bestDistance = 999;
bestCoordinate = 111;
top = [0, 1];

idx = -1;

for i = 1:size(x, 1)
    curr = [x(i), y(i)];
    currDistance = sqrt( (top(1)-curr(1))^2 + (top(2)-curr(2))^2 );
    if currDistance<bestDistance
        bestCoordinate = x(i);
        bestDistance = currDistance;
        idx = i;
    end
end

threshold = t(idx);

% plot(x(idx), y(idx), 'o');
