% Define parameters
inputDirectory = '/Volumes/Crucial X6/CollageFeatures/';
plane = 'Axial'; % Specify the plane: 'Axial' or 'Coronal'
planes = 'ax'; % Specify the plane: 'ax' or 'cor'
regions = {'Tumor'}; % Specify region: 'Fat', 'Tumor', 'ProxFat5', 'ProxFat10', or 'ProxFat15'
numPatients = 92;

% Initialize the combined array
combinedArray = zeros(numPatients, 260); % Assuming each patient has 260 values

% Loop through patients and regions
for patientID = 1 : numPatients
    for regionIdx = 1 : length(regions)
        region = regions{regionIdx};
        
        % Create the file path
        filename = sprintf('%s_%s/Feats_Col_Patient-%03d_%s.mat', plane, region, patientID, planes);
        filePath = fullfile(inputDirectory, filename);
        
        % Load the data from the file
        load(filePath); % Assuming the .mat file contains a variable named 'data'
        
        % Assuming 'data' is a row vector of size 1x260
        combinedArray(patientID, :) = array_data;
    end
end

% Save the combined array
output_filename = sprintf('%s_%s/Feats_Col.mat', plane, region);
output_filepath = fullfile(inputDirectory, output_filename);
save(output_filepath, 'combinedArray');
disp('Combining and saving complete.');
