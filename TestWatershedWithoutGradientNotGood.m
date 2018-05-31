A = imread('DJI_0486.JPG');
DetectedPumpkins = imread('DJI_0486RedMarked.PNG');

Channels = Functionclass.MarkedPumpkins(A,DetectedPumpkins);
ChannelsRotated = rot90(Channels);
%%
BinaryPumkinImage2 =  Functionclass.FindPumpkinsInImageVersion2(A,ChannelsRotated);
%BinaryImageInverted = imcomplement(BinaryPumkinImage2);
figure(1)
imshow(BinaryPumkinImage2)
%% Fill holes in pumpkins, caused by overbeslysning
BinaryPumkinImage2Filled = imfill(BinaryPumkinImage2,'holes');
figure(2)
imshow(BinaryPumkinImage2Filled), title('Filled binary image (BinaryPumkinImage2Filled)')

%%  Filter erode, dialate (remove pumpkins with size on 3 pixels. )
se1 = strel('disk', 3);
Io = imopen(BinaryPumkinImage2Filled, se1);
figure
imshow(Io), title('Opening (Io)')
%%
ColorBinaryCombinedImage = Functionclass.CombineBinaryAndColorImage(A,Io);
ColorBinaryCombinedGrayscaleImage = rgb2gray(ColorBinaryCombinedImage);
imshow(ColorBinaryCombinedGrayscaleImage), title('CombinedGrayscale (ColorBinaryCombinedGrayscaleImage)')
%% Filter the grayscaled pumpkins, in order to ignore noise before finding the gradient. 
se = strel('disk', 3);
Ie = imerode(ColorBinaryCombinedGrayscaleImage, se);
Iobr2 = imreconstruct(Ie, ColorBinaryCombinedGrayscaleImage);
figure(5)
imshow(Iobr2), title('Opening-by-reconstruction (Iobr)')
%% Finding Gradient
hy = fspecial('sobel');
hx = hy';
Iy = imfilter(double(Iobr2), hy, 'replicate');
Ix = imfilter(double(Iobr2), hx, 'replicate');
gradmag = sqrt(Ix.^2 + Iy.^2);
%%  Finding bwdistance 
[D1,IDX] = bwdist(~Iobr2,'euclidean');
figure(6)
imshow(D1,[]), title('Bw dist (D1)')
%%
L = watershed(D1);
 Lrgb = label2rgb(L);
 figure, imshow(Lrgb), title('Watershed transform of gradient magnitude (Lrgb)')
 %%
 figure(8)
p1=subplot(2, 1, 1);
imshow(Lrgb);
%h = viscircles(MatrixPumpkins,RadiiRot);

p2=subplot(2, 1, 2);
imshow(A)
%h = viscircles(MatrixPumpkins,RadiiRot);

linkaxes([p1, p2]);

%%

 figure(9)
p1=subplot(2, 1, 1);
imshow(A);
%h = viscircles(MatrixPumpkins,RadiiRot);

p2=subplot(2, 1, 2);
imshow(Io)
%h = viscircles(MatrixPumpkins,RadiiRot);

linkaxes([p1, p2]);