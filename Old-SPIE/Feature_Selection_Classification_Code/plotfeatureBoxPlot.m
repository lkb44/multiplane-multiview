function plotfeatureBoxPlot(feature_classes, featureName, x_label, y_label, legend, y_lim, x_width, y_width, xTickAngle, colors, format)
%PLOTFEATUREBOXPLOT Create a boxplot for a single feature.
%   The input cell array must contain the feature arrays for each class in
%   each cohort and follow a specific format. There are two options for the
%   format of the cell array:
%       1). class: features in the input cell array are arranged such that
%       [train class 0, train class 1, test class 0, test class 1, excluded
%       0, excluded class 1].
%       2). cohort: features in the input cell array are arranged such that
%       [train class 0, test class 0, excluded class 0, train class 1, test
%       class 1, excluded class 1].
%   Each class is given a specific color.
%
%   Parameters:
%       feature_classes: cell array
%           1x6 cell array containing the feature in each class of each
%           cohort. The cell array must follow either the "class" format or
%           the "cohort" format.
%       featureName: str
%           Name of feature being plotted. Used for naming the plot.
%       x_label: str
%           X-axis title for boxplot
%       y_label: str
%           Y-axis title for boxplot
%       legend: string array
%           Each string is a label for a box on the boxplot
%       y_lim: vector
%           Vector integers that define the range of the Y-axis
%       x_width: int
%           Length of box plot. Used to format the figure for saving.
%       y_width: int
%           Height of box plot. Used to format the figure for saving.
%       xTickAngle: int
%           Angle of X-axis labels. X-axis labels are derived from the
%           legend
%       colors: matrix
%           2x3 matrix of floats. Each row represents an RGB color in the
%           MATLAB RGB formatting (i.e. 1/255 * [R G B]). Columns 1, 2, and
%           3 represent the red, green, and blue values for each color.
%           One color represents the class 0 and the other color represent
%           class 1.
%       format: str
%           Determines how the features per class per cohort are arranged
%           in the cell array. Also used to color classes. 
%               If format = "class," boxes in odd number indices of the
%               boxplot will be one color, while the boxes in even number 
%               indices will be a different color.
%               If format = "cohort," boxes in the first 3 indices will be
%               one color, while the boxes in the last 3 indices will be a
%               different color.
%   Returns:
%       None. A boxplot figure is generated.

fprintf("Creating boxplot for %s\n", featureName);

% Create the boxplot
cell2boxplot(feature_classes, 1, featureName, x_label, y_label, legend, y_lim);
xtickangle(xTickAngle);
% set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 x_width y_width]);

% Set class indicies for coloring
h = findobj(gca,'Tag','Box');
if(isequal(format,'class'))
    h_class0 = h(1:2:end);
    h_class1 = h(2:2:end);
elseif(isequal(format,'cohort'))
    h_class0 = h(1:2);
    h_class1 = h(3:4);
else
    error("Format not recognized!");
end

% Color the class 0 boxes
for j=1:length(h_class0)
    patch(get(h_class0(j),'XData'),get(h_class0(j),'YData'),colors(1,:),'FaceAlpha',.3);
end

% Color the class 1 boxes
for j=1:length(h_class1)
    patch(get(h_class1(j),'XData'),get(h_class1(j),'YData'),colors(2,:),'FaceAlpha',.3);
end

fprintf("Boxplot created successfully!\n")

end