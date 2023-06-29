function feature_matrix = extract_2Dfeats(vol)

feature_matrix = [];
%for j=1:size(vol)
    slice = uint8(vol(:,:));
    if max(max(slice))~=0
        BW = imbinarize(slice);
        BW = bwconncomp(BW);
        sf = regionprops(BW,'all');
        sf = sf(1);
        
        area = sf.Area;
        majorAxisLength = sf.MajorAxisLength;
        minorAxisLength = sf.MinorAxisLength;
        eccentricity = sf.Eccentricity;
        elongation = (minorAxisLength)/(majorAxisLength);
        orientation = sf.Orientation;
        perimeter = sf.Perimeter;
        convexPerimeter = regionprops(sf.ConvexImage,'perimeter');
        roundness = (4*pi*area)/(convexPerimeter.Perimeter*convexPerimeter.Perimeter);
        equivDiameter = sf.EquivDiameter;
        compactness = (4*pi*area)/(perimeter*perimeter);
        convexArea = sf.ConvexArea;
        solidity = sf.Solidity;
        convexity = convexPerimeter.Perimeter/perimeter;
        extent = sf.Extent;
        boundingBoxArea = sf.BoundingBox(1,3)*sf.BoundingBox(1,4);
        elongationBoundingBox = min(sf.BoundingBox(1,3),sf.BoundingBox(1,4))/max(sf.BoundingBox(1,3),sf.BoundingBox(1,4));
        filledArea = sf.FilledArea;
        
        feature_vector = [area majorAxisLength minorAxisLength eccentricity elongation ...
            orientation perimeter roundness equivDiameter compactness convexArea ...
            convexPerimeter.Perimeter solidity convexity extent boundingBoxArea ...
            elongationBoundingBox filledArea];
%         feature_vector = [area majorAxisLength minorAxisLength eccentricity elongation ...
%             orientation perimeter equivDiameter compactness ...
%             solidity extent];
        feature_matrix = [feature_matrix ; feature_vector];
    else
        feature_matrix = [feature_matrix ; zeros(1,18)]; %Esta es la linea de mas
        %feature_matrix = [feature_matrix ; zeros(1,11)]; %Esta es la linea de mas
    end
%end