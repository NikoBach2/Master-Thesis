A = imread('DJI_0486.JPG');
DetectedPumpkins = imread('DJI_0486RedMarked.PNG');

Channels = Functionclass.MarkedPumpkins(A,DetectedPumpkins);
ChannelsRotated = rot90(Channels);
Shows = 'A_Original' ;
Folder = 'imrodeDiskSize3';
Functionclass.SaveFigures(A,Folder,Shows);
%%
BinaryPumkinImage2 =  Functionclass.FindPumpkinsInImageVersion2(A,ChannelsRotated);
%BinaryImageInverted = imcomplement(BinaryPumkinImage2);
figure(1)
imshow(BinaryPumkinImage2)
Shows = 'B_Binary' ;
Folder = 'imrodeDiskSize3';
Functionclass.SaveFigures(BinaryPumkinImage2,Folder,Shows);
%% Fill holes in pumpkins, caused by overbeslysning
BinaryPumkinImage2Filled = imfill(BinaryPumkinImage2,'holes');
figure(2)
imshow(BinaryPumkinImage2Filled), title('Filled binary image (BinaryPumkinImage2Filled)')
Shows = 'C_BinaryFilled' ;
Folder = 'imrodeDiskSize3';
Functionclass.SaveFigures(BinaryPumkinImage2Filled,Folder,Shows);
%%  Filter erode, dialate (remove pumpkins with size on 3 pixels. )
se1 = strel('disk', 3);
CleanedBinaryImage = imopen(BinaryPumkinImage2Filled, se1);
figure(3)
imshow(CleanedBinaryImage), title('Opening (CleanedBinaryImage)')
Shows = 'D_DialateErode' ;
Folder = 'imrodeDiskSize3';
Functionclass.SaveFigures(CleanedBinaryImage,Folder,Shows);
%%
ColorBinaryCombinedImage = Functionclass.CombineBinaryAndColorImage(A,CleanedBinaryImage);
ColorBinaryCombinedGrayscaleImage = rgb2gray(ColorBinaryCombinedImage);
figure(4)
imshow(ColorBinaryCombinedGrayscaleImage), title('CombinedGrayscale (ColorBinaryCombinedGrayscaleImage)')
Shows = 'E_CombinedGrayBinary' ;
Folder = 'imrodeDiskSize3';
Functionclass.SaveFigures(ColorBinaryCombinedGrayscaleImage,Folder,Shows);
%% Filter the grayscaled pumpkins, in order to ignore noise before finding the gradient. 
se = strel('disk', 3);  %% tuned to 3
Ie = imerode(ColorBinaryCombinedGrayscaleImage, se);
Iobr2 = imreconstruct(Ie, ColorBinaryCombinedGrayscaleImage);
figure(5)
imshow(Iobr2), title('Opening-by-reconstruction (Iobr)')
Shows = 'F_CombinedGrayBinaryRemovedNoise' ;
Folder = 'imrodeDiskSize3';
Functionclass.SaveFigures(Iobr2,Folder,Shows);
%% Finding Gradient
hy = fspecial('sobel');
hx = hy';
Iy = imfilter(double(Iobr2), hy, 'replicate');
Ix = imfilter(double(Iobr2), hx, 'replicate');
gradmag = sqrt(Ix.^2 + Iy.^2);
figure(5)
imshow(gradmag,[])
Shows = 'G_Gradient' ;
Folder = 'imrodeDiskSize3';
Functionclass.SaveFigures(gradmag,Folder,Shows);
%%  Finding bwdistance 
[D1,IDX] = bwdist(~gradmag,'euclidean');
figure(5)
imshow(D1,[])
Shows = 'H_BwDist' ;
Folder = 'imrodeDiskSize3';

Functionclass.SaveFigures(D1,Folder,Shows);
%%
L = watershed(D1);
 Lrgb = label2rgb(L);
 figure(6), imshow(Lrgb), title('Watershed transform of gradient magnitude (Lrgb)')
 Shows = 'I_Watershed' ;
Folder = 'imrodeDiskSize3';
Functionclass.SaveFigures(Lrgb,Folder,Shows);
%%
NewBinaryPumkinImage =  L > 6; % not bigger than 12 
figure(10), imshow(NewBinaryPumkinImage), title('Watershed transform of gradient magnitude (Lrgb)')
 Shows = 'J_BinaryWatershed' ;
Folder = 'imrodeDiskSize3';
Functionclass.SaveFigures(NewBinaryPumkinImage,Folder,Shows);
 %%
 figure(7)
p1=subplot(2, 1, 1);
imshow(Lrgb);
%h = viscircles(MatrixPumpkins,RadiiRot);

p2=subplot(2, 1, 2);
imshow(A)
%h = viscircles(MatrixPumpkins,RadiiRot);

linkaxes([p1, p2]);

%%

 figure(8)
p1=subplot(2, 1, 1);
imshow(BinaryPumkinImage2Filled);
%h = viscircles(MatrixPumpkins,RadiiRot);

p2=subplot(2, 1, 2);
imshow(CleanedBinaryImage)
%h = viscircles(MatrixPumpkins,RadiiRot);

linkaxes([p1, p2]);
%%

GrayScalePumpkinImage = rgb2gray(A);
sRegionProps = regionprops(NewBinaryPumkinImage, GrayScalePumpkinImage, {'Centroid'});

figure(54)
imshow(GrayScalePumpkinImage)
title('Weighted (red) and Unweighted (blue) Centroids'); 
hold on
numObj = numel(sRegionProps);

for k = 1 : numObj
   % plot(sRegionProps(k).WeightedCentroid(1), sRegionProps(k).WeightedCentroid(2), 'r*');
    plot(sRegionProps(k).Centroid(1), sRegionProps(k).Centroid(2), 'bo');
end
hold off

String = 'The manual counted pumpkins is counted to:'
size(T)

String = 'The program counted pumpkins is counted to:'
numObj