% Se cambia fuente para cálculo de: EquivDiameter, ConvexArea,Solidity
function feature_matrix = extract_2Dfeats_erw_v2(volERW,volLumen)

feature_matrix = [];
for j=1:size(volERW.img,3)
    slice = uint8(volERW.img(:,:,j));
    slice_lumen = uint8(volLumen.img(:,:,j));
    slice_join = slice+slice_lumen;
    if max(max(slice))~=0
        BW = imbinarize(slice);
        BW = bwconncomp(BW);
        sf = regionprops(BW,'all');
        sf = sf(1);
        
        BW_join = imbinarize(slice_join);
        BW_join = bwconncomp(BW_join);
        sf_join = regionprops(BW_join,'all');
        sf_join = sf_join(1);
        
        area = sf.Area;
        majorAxisLength = sf_join.MajorAxisLength;%sf.MajorAxisLength; %sf_join.MajorAxisLength; % *******
        minorAxisLength = sf_join.MinorAxisLength;%sf.MinorAxisLength; %sf_join.MinorAxisLength; % *******
        eccentricity = sf_join.Eccentricity;%sf.Eccentricity; %sf_join.Eccentricity; % *******
        elongation = (minorAxisLength)/(majorAxisLength); % *******
        orientation = sf_join.Orientation;%sf.Orientation; %sf_join.Orientation; % *******
        perimeter = sf_join.Perimeter;%sf.Perimeter; %sf_join.Perimeter; % *******
        convexPerimeter = regionprops(sf_join.ConvexImage,'perimeter');%regionprops(sf.ConvexImage,'perimeter'); %regionprops(sf_join.ConvexImage,'perimeter'); % *******
        roundness = (4*pi*area)/(convexPerimeter.Perimeter*convexPerimeter.Perimeter);
        equivDiameter = sf_join.EquivDiameter;
        compactness = (4*pi*area)/(perimeter*perimeter);
        convexArea = sf_join.ConvexArea;
        solidity = area/convexArea; %sf.Solidity;
        convexity = convexPerimeter.Perimeter/perimeter;
        boundingBoxArea = sf_join.BoundingBox(1,3)*sf_join.BoundingBox(1,4);
        elongationBoundingBox = min(sf_join.BoundingBox(1,3),sf_join.BoundingBox(1,4))/max(sf_join.BoundingBox(1,3),sf_join.BoundingBox(1,4));
        extent = area/boundingBoxArea;%sf.Extent;
        filledArea = sf_join.FilledArea;
        
        centroid = sf.Centroid;
        convexImage = sf.ConvexImage;
        
        feature_vector = [area majorAxisLength minorAxisLength eccentricity elongation ...
            orientation perimeter roundness equivDiameter compactness convexArea ...
            convexPerimeter.Perimeter solidity convexity extent boundingBoxArea ...
            elongationBoundingBox filledArea];
%         feature_vector = [area majorAxisLength minorAxisLength eccentricity elongation ...
%             orientation perimeter equivDiameter compactness ...
%             solidity extent];
        feature_matrix = [feature_matrix ; feature_vector];
        
%         x = round(centroid(1,1)); % x, y:    Center of the circle
%         y = round(centroid(1,2));
%         r=1; % r:       Radius of the circle
%         imshow(BW);
%         hold on 
%         theta = 0 : (2 * pi / 10000) : (2 * pi);
%         pline_x = r * cos(theta) + x;
%         pline_y = r * sin(theta) + y;
%         %plot(pline_x, pline_y, '-');
    else
        feature_matrix = [feature_matrix ; zeros(1,18)]; %Esta es la linea de mas
        %feature_matrix = [feature_matrix ; zeros(1,11)]; %Esta es la linea de mas
    end
end