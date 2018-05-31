%A = imread('DJI_0486.JPG');
%A = imread('pumpkinPictures2017/DJI_0255CropedInkvart.JPG');
A = imread('pumpkinPictures2017/DJI_0255Croped.JPG');
%DetectedPumpkins = imread('DJI_0486RedMarked.PNG');
%DetectedPumpkins = imread('pumpkinPictures2017/DJI_0255CropedInkvartRedMarked.png');
%A = imread('pumpkinPictures2017/DJI_0065.JPG');
%DetectedPumpkins = imread('pumpkinPictures2017/DJI_0065RedMarked.png');
%TableName = 'pumpkinPictures2017/DJI_0065_MarkedWithoutWhiteNoDiameter.csv';
LightImage = A + 0;
HSVImage=rgb2hsv(LightImage);
DetectedPumpkins = imread('pumpkinPictures2017/DJI_0255CropedMarkedPumpkins.png');
%TableName = 'DJI_0486ManualCountedDiameterFinal.csv';
%TableName = 'pumpkinPictures2017/DJI_0255CropedCountedWithDiameterupdated.csv';
TableName = 'pumpkinPictures2017/DJI_612_MarkedPumpkinsDiameterFinalNogreen.csv';

TableNameHenrik = 'DJI_0486_points.csv';

T = readtable(TableName);

%Channels = Functionclass.MarkedPumpkins(A,DetectedPumpkins);
%ChannelsRotated = rot90(Channels);
ChannelsHSV = Functionclass.MarkedPumpkinsHSV(LightImage,DetectedPumpkins);
ChannelsRotatedHSV = rot90(ChannelsHSV);
Shows = 'A_Original' ;
Folder = 'JustTesting';
%Functionclass.SaveFigures(A,Folder,Shows);
%%

BinaryPumkinImage2 =  Functionclass.FindPumpkinsInImageVersion2HSV(HSVImage,ChannelsRotatedHSV);
%BinaryImageInverted = imcomplement(BinaryPumkinImage2);
%figure(1)
imshow(BinaryPumkinImage2)
Shows = 'B_Binary' ;
Folder = 'JustTesting';
%Functionclass.SaveFigures(BinaryPumkinImage2,Folder,Shows);

%%
seDisk = strel('disk',2);
IM2 = imdilate(BinaryPumkinImage2,seDisk);   % my filter example
%IM2 = imdilate(BinaryImageKmean,seDisk); %BinaryImageKmean example
%IM2 = imdilate(BinImage,seDisk);   % superpixel eksample

%figure(77)
%imshow(IM2)
Shows = 'B_A_BinaryDialate' ;
Folder = 'JustTesting';
%Functionclass.SaveFigures(IM2,Folder,Shows);
%% Fill holes in pumpkins, caused by overbeslysning
BinaryPumkinImage2Filled = imfill(IM2,'holes');
%figure(2)
%imshow(BinaryPumkinImage2Filled), title('Filled binary image (BinaryPumkinImage2Filled)')
Shows = 'C_BinaryFilled' ;
Folder = 'JustTesting';
%Functionclass.SaveFigures(BinaryPumkinImage2Filled,Folder,Shows);
%%
seDisk = strel('disk',2);
BinaryPumkinImage2FilledErode = imerode(BinaryPumkinImage2Filled,seDisk);
%figure(78)
%imshow(BinaryPumkinImage2FilledErode)
Shows = 'C_C_BinaryDialate' ;
Folder = 'JustTesting';
%Functionclass.SaveFigures(BinaryPumkinImage2FilledErode,Folder,Shows);
%%  Filter erode, dialate (remove pumpkins with size on 3 pixels. )
se1 = strel('disk', 3);
CleanedBinaryImage = imopen(BinaryPumkinImage2FilledErode, se1);
%figure(3)
%imshow(CleanedBinaryImage), title('Opening (CleanedBinaryImage)')
Shows = 'D_DialateErode' ;
Folder = 'JustTesting';
%Functionclass.SaveFigures(CleanedBinaryImage,Folder,Shows);
%%
ColorBinaryCombinedImage = Functionclass.CombineBinaryAndColorImage(A,CleanedBinaryImage);
ColorBinaryCombinedGrayscaleImage = rgb2gray(ColorBinaryCombinedImage);
%figure(4)
%imshow(ColorBinaryCombinedGrayscaleImage), title('CombinedGrayscale (ColorBinaryCombinedGrayscaleImage)')
Shows = 'E_CombinedGrayBinary' ;
Folder = 'JustTesting';
%Functionclass.SaveFigures(ColorBinaryCombinedGrayscaleImage,Folder,Shows);
%% Filter the grayscaled pumpkins, in order to ignore noise before finding the gradient. 
se = strel('disk', 3);  %% tuned to 3
Ie = imerode(ColorBinaryCombinedGrayscaleImage, se);
Iobr2 = imreconstruct(Ie, ColorBinaryCombinedGrayscaleImage);
%figure(5)
%imshow(Iobr2), title('Opening-by-reconstruction (Iobr)')
Shows = 'F_CombinedGrayBinaryRemovedNoise' ;
Folder = 'JustTesting';
%Functionclass.SaveFigures(Iobr2,Folder,Shows);
%% Finding Gradient
hy = fspecial('sobel');
hx = hy';
Iy = imfilter(double(Iobr2), hy, 'replicate');
Ix = imfilter(double(Iobr2), hx, 'replicate');
gradmag = sqrt(Ix.^2 + Iy.^2);
%figure(5)
%imshow(gradmag,[])
Shows = 'G_Gradient' ;
Folder = 'JustTesting';
%Functionclass.SaveFigures(gradmag,Folder,Shows);
%% Trying to remove counting on pumpkin double, by removing single black pixels inside a pumpkin.
seDisk = strel('disk',1);
IM2 = imdilate(gradmag,seDisk);
%im2 = imopen(gradmag, seDisk);
%IM3 = imerode(BinaryPumkinImage2Filled,seDisk);

%figure(77)
%imshow(IM2)
Shows = 'G_H_BinaryDialate' ;
Folder = 'JustTesting';
%Functionclass.SaveFigures(IM2,Folder,Shows);
%% Trying to remove counting on pumpkin double, by removing single black pixels inside a pumpkin.
seDisk1 = strel('disk',1);
IM3 = imerode(IM2,seDisk1);

%figure(77)
%imshow(IM3)
Shows = 'G_I_BinaryDialate' ;
Folder = 'JustTesting';
%Functionclass.SaveFigures(IM3,Folder,Shows);
%%  Finding bwdistance 
[D1,IDX] = bwdist(~IM3,'euclidean');
%figure(5)
%imshow(D1)
Shows = 'H_BwDist' ;
Folder = 'JustTesting';

%Functionclass.SaveFigures(D1,Folder,Shows);
%% Watershed transform
L = watershed(D1);
 Lrgb = label2rgb(L);
 %figure(6), imshow(Lrgb), title('Watershed transform of gradient magnitude (Lrgb)')
 Shows = 'I_Watershed' ;
Folder = 'JustTesting';
%Functionclass.SaveFigures(Lrgb,Folder,Shows);
%% Creates a Binary Image after the watershed transform
NewBinaryPumkinImage =  L > 6; % not bigger than 12 
%figure(10), imshow(NewBinaryPumkinImage), title('Watershed transform of gradient magnitude (Lrgb)')
 Shows = 'J_BinaryWatershed' ;
Folder = 'JustTesting';
%Functionclass.SaveFigures(NewBinaryPumkinImage,Folder,Shows);

 %%
 %figure(7)
%p1=subplot(2, 1, 1);
%imshow(A);
%h = viscircles(MatrixPumpkins,RadiiRot);

%p2=subplot(2, 1, 2);
%imshow(BinaryPumkinImage2FilledErode)
%h = viscircles(MatrixPumpkins,RadiiRot);

%linkaxes([p1, p2]);

%%

% figure(8)
%p1=subplot(2, 1, 1);
%imshow(BinaryPumkinImage2);
%h = viscircles(MatrixPumpkins,RadiiRot);

%p2=subplot(2, 1, 2);
%imshow(BinaryPumkinImage2)
%h = viscircles(MatrixPumpkins,RadiiRot);

%linkaxes([p1, p2]);
%% Counts the number of pumpkins in the image. 

GrayScalePumpkinImage = rgb2gray(A);
sRegionProps = regionprops(NewBinaryPumkinImage, GrayScalePumpkinImage, {'Centroid','Area'});
%%
TablePumpkinPositions = AnalysisClass.ReadTableValues(TableName, TableNameHenrik,1); %TablePumpkinPositions;
%manual counted pumpkins 0 = ingen diamter    1 = med diameter 
%%

allCenters = Functionclass.ConvertFromStructToArray(sRegionProps);
%allCenters = MatrixPumpkins  % for first version filter values
Functionclass.KNNTreeFindAverageDistance(TablePumpkinPositions,allCenters)

%%
AllArea = Functionclass.ConvertFromStructToArrayArea(sRegionProps);
AreaProgramSum = sum(AllArea)
%%
'program found, i did not'
Name1 = '1';
Name2 = '2';

% TablePumpkinPositionsHenrik

Functionclass.FindOutlierPumpkins(TablePumpkinPositions,allCenters,A,7,Name1)
% TablePumpkinPositionsHenrik
'I found this, program did not'
Functionclass.FindOutlierPumpkins(allCenters,TablePumpkinPositions,A,7,Name2)
%Functionclass.FindOutlierPumpkins(TablePumpkinPositionsHenrik,TablePumpkinPositions,A,7)
%Functionclass.FindOutlierPumpkins(TablePumpkinPositions,TablePumpkinPositionsHenrik,A,7)
%%

ManualCountedArea = AnalysisClass.ManualCountedPumpkinArea(TableName,TableNameHenrik);
ManualCountedSumArea = sum(ManualCountedArea);
AveragePumpkinSize = AnalysisClass.OutliersAndArea(TablePumpkinPositions,allCenters,ManualCountedArea,AllArea)

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
String = 'The manual counted pumpkins is counted to:'
size(T)

String = 'The program counted pumpkins is counted to:'
numObj