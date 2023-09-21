import openpyxl

# List of descriptors
descriptors = [
    "AngularSecondMoment",
    "Contrast",
    "Correlation",
    "SumOfSquareVariance",
    "SumAverage",
    "SumVariance",
    "SumEntropy",
    "Entropy",
    "DifferenceVariance",
    "DifferenceEntropy",
    "InformationMeasureOfCorrelation1",
    "InformationMeasureOfCorrelation2",
    "MaximalCorrelationCoefficient"
]

# List of statistics
statistics = ["median", "var", "kurtosis", "skewness"]

# List of window sizes
window_sizes = ["W3", "W5", "W7", "W9", "W11"]

# Create a new workbook and select the active worksheet
workbook = openpyxl.Workbook()
sheet = workbook.active

# Write headers
header_row = ["Feature Name"]
for window_size in window_sizes:
    for descriptor in descriptors:
        for stat in statistics:        
            header_row.append(f"{stat}-{descriptor}{window_size}")

sheet.append(header_row)

# Write feature names
for descriptor in descriptors:
    for stat in statistics:
        for window_size in window_sizes:
            feature_name = f"{stat}-{descriptor}{window_size}"
            sheet.append([feature_name])

# Save the workbook
workbook.save("/Users/leobao/Documents/MultiPlanePipeline/Feature_Names/Collage_Feature_Names.xlsx")