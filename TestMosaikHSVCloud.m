%A = imread('DJI_0486.JPG');
NewBinaryPumkinImage = imread('mosaik/KokogKilemark_transparent_mosaic_group1WithMaskOutputFiltered.tif');





%% Counts the number of pumpkins in the image. 

%GrayScalePumpkinImage = rgb2gray(A);
%sRegionProps = regionprops(NewBinaryPumkinImage, GrayScalePumpkinImage, {'Centroid','Area'});
sRegionProps = regionprops(NewBinaryPumkinImage, {'Centroid','Area'});
%%
%TablePumpkinPositions = AnalysisClass.ReadTableValues(TableName, TableNameHenrik,1); %TablePumpkinPositions;
%manual counted pumpkins
%%

%allCenters = Functionclass.ConvertFromStructToArray(sRegionProps);
%allCenters = MatrixPumpkins  % for first version filter values
%Functionclass.KNNTreeFindAverageDistance(TablePumpkinPositions,allCenters)

%%
AllArea = Functionclass.ConvertFromStructToArrayArea(sRegionProps);
%AreaProgramSum = sum(AllArea)

%%

figure(54)
imshow(BinaryPumkinImage2FilledErode)
title('Weighted (red) and Unweighted (blue) Centroids'); 
hold on
numObj = numel(sRegionProps);

for k = 1 : numObj
   % plot(sRegionProps(k).WeightedCentroid(1), sRegionProps(k).WeightedCentroid(2), 'r*');
    plot(sRegionProps(k).Centroid(1), sRegionProps(k).Centroid(2), 'bo');
end
hold off

%%
numObj = numel(sRegionProps);
%String = 'The manual counted pumpkins is counted to:'
%size(T)

String = 'The program counted pumpkins is counted to:'
numObj