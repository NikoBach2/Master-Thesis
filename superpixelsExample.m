%[L,N] = superpixels(A,50000);
[L,N] = superpixels(LightImage,500);  % er valgt som 1/4 a 200000, da omkring dette antal virkede fint i test
%[L,N] = superpixels(BlurImage,10000);


figure
BW1 = boundarymask(L);

bws = bwmorph(BW1,'skel');
imshow(imoverlay(LightImage,bws,'cyan'),'InitialMagnification',67)


%%  
outputImage = zeros(size(A),'like',A);
idx = label2idx(L);
numRows = size(A,1);
numCols = size(A,2);
bwsCopy = bws;
for labelVal = 1:N
    redIdx = idx{labelVal};
    greenIdx = idx{labelVal}+numRows*numCols;
    blueIdx = idx{labelVal}+2*numRows*numCols;
    ColorsCombined = [LightImage(redIdx),LightImage(greenIdx),LightImage(blueIdx)];
    PumpkinFloating = ChannelsRotated / 255;
    ColorsCombinedDouble = im2double(ColorsCombined);
    MahalDistanceMatrice = mahal(ColorsCombinedDouble,PumpkinFloating);
    % =  Functionclass.FindPumpkinsInImageVersion2(ColorsCombined,ChannelsRotated)
    %if MahalDistanceMatrice < 52  %12
    outputImage(redIdx) = mean(LightImage(redIdx));
    outputImage(greenIdx) = mean(LightImage(greenIdx));
    outputImage(blueIdx) = mean(LightImage(blueIdx));
    
    %else
    %    outputImage(redIdx) = 0;
    %    outputImage(greenIdx) = 0;
    %    outputImage(blueIdx) = 0;
    %end
end    

figure(9)
imshow(outputImage,'InitialMagnification',67)
%imshow(imoverlay(outputImage,A),'InitialMagnification',67)


% Creates a Binary Image after the watershed transform
GrayScalePumpkinImage2 = rgb2gray(outputImage);
NewBinaryPumkinImage =  GrayScalePumpkinImage2 > 6; % not bigger than 12 
figure(10)
imshow(NewBinaryPumkinImage), title('Watershed transform of gradient magnitude (Lrgb)')

%% plot combined image superpixel pumpkins, and original image

BinImage = Functionclass.FindPumpkinsInImageVersion2(outputImage,ChannelsRotated);

BinImageMorph = BinImage & bwsCopy;

I2 = LightImage;
I2(BinImage) = 0;
figure
imshow(BinImage), title('Regional maxima superimposed on original image (I2)')

%%
%bws = bwmorph(BW1,'skel');
%Over = imoverlay(NewBinaryPumkinImage,bws,'black');
%figure(11)
%imshow(Over,'InitialMagnification',67)
%% Counts the number of pumpkins in the image. 

%GrayScalePumpkinImage = rgb2gray(A);
%GrayScalePumpkinImage3 = rgb2gray(Over);
%NewBinaryPumkinImage3 =  GrayScalePumpkinImage3 > 6; % not bigger than 12
%sRegionProps1 = regionprops(NewBinaryPumkinImage3, GrayScalePumpkinImage, {'Centroid'});


%figure(54)
%imshow(NewBinaryPumkinImage)
%title('Weighted (red) and Unweighted (blue) Centroids'); 
%hold on
%numObj = numel(sRegionProps1);

%for k = 1 : numObj
   % plot(sRegionProps(k).WeightedCentroid(1), sRegionProps(k).WeightedCentroid(2), 'r*');
 %   plot(sRegionProps1(k).Centroid(1), sRegionProps1(k).Centroid(2), 'bo');
%end
%hold off



%%
%Functionclass.SuperPixelsFunction(A,ChannelsRotated);
%%
%allCenters1 = Functionclass.ConvertFromStructToArray(sRegionProps1);
%Functionclass.KNNTreeFindAverageDistance(TablePumpkinPositions,allCenters1)


%%

%Functionclass.FindOutlierPumpkins(TablePumpkinPositions,allCenters1,A,10)
%Functionclass.FindOutlierPumpkins(allCenters1,TablePumpkinPositions,A,10)


%%
%redIdx = idx{1};
%A(redIdx)
%%
%ColorsCombined = [A(redIdx) A(redIdx) A(redIdx)];

%%
%figure(20)
%imshow(mod(L,200)/200)

% Lrgb1 = label2rgb(L,'jet',[1 1 1],'shuffle');
 
 %%
 %bwsInverted = imcomplement(bws);
 %bins = bwconncomp(bwsInverted,4)
 
 
 
 
 
 %%
%ArrayHej = bins.PixelIdxList(1,1)
%ArrayHej{1,1}
%bins(1).Centroid(2);
 %%
 %CC = bwconncomp(bws);
 %L1 = labelmatrix(bins);
 
 %figure(1)
 %imshow(L1);
%%

% figure(8)
%p1=subplot(2, 1, 1);
%imshow(Lrgb1);
%h = viscircles(MatrixPumpkins,RadiiRot);

%p2=subplot(2, 1, 2);
%imshow(imoverlay(A,bws,'cyan'),'InitialMagnification',67)
%h = viscircles(MatrixPumpkins,RadiiRot);

%linkaxes([p1, p2]);