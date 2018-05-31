% read the original photo into memory
%origp = imread('DJI_0486.tif');



%origp(:,:,[4]) = [];
%%
origp_filename = 'mosaik/Kok_og_Kilemark_transparent_mosaic_group1WithMask.tif';
%origp_filename = 'DJI_0486.tif';
derp_filename = 'mosaik/Kok_og_Kilemark_transparent_mosaic_group1WithMask_RGB_output.tif';
%derp_filename = 'testtt.tif';
% create a Gaussian low-pass filter
%h = fspecial('gaussian',5,2);

%BlockProcessFunc.BlockFunc()
A = imread('mosaik/Kok_og_Kilemark_transparent_mosaic_group1WithMask_udsnit.png');
DetectedPumpkins = imread('mosaik/Kok_og_Kilemark_transparent_mosaic_group1WithMask_udsnit_redMarked.png');
%%
%ChannelsHSV = Functionclass.MarkedPumpkinsHSV(A,DetectedPumpkins);
%ChannelsRotatedHSV = rot90(ChannelsHSV);
Channels = Functionclass.MarkedPumpkins(A,DetectedPumpkins);
ChannelsRotated = rot90(Channels);
% create a function handle that we will apply to each block
myFun = @(block_struct) BlockProcessFunc.BlockFunc(block_struct.data,ChannelsRotated);
%myFun = @(block_struct) BlockProcessFunc.RGB2gr(block_struct.data);
%myFun = @(block_struct) BlockProcessFunc.KmeanFunc(block_struct.data);
% setup block size
block_size = [800 800];
border_size = [50 50];
% specify a string filename as our input image
%origp_filename = 'cameraman.tif';
% specify a string filename as the 'Destination' of our output data
%derp_filename = 'output.tif';
% don't request an output argument this time!
blockproc(origp_filename,block_size,myFun,...
    'BorderSize',border_size,'Destination',derp_filename);

%%
imshow(derp_filename);


%%

%derp_filename1 = 'mosaik/kok_og_kile_markMaskEkstra_Working_Output.tif';
%derp_filenameGray = 'mosaik/kok_og_kile_markMaskEkstra_Working_Gray.tif';
derp_filename1 = imread('mosaik/Kok_og_Kilemark_transparent_mosaic_group1WithMask_RGB_output.tif');
derp_filenameGray = imread('mosaik/Kok_og_Kilemark_transparent_mosaic_group1WithMask_Gray.tif');

%GrayScalePumpkinImage = rgb2gray(A);
sRegionProps = regionprops(derp_filename1, derp_filenameGray, {'Centroid','Area'});
%%

figure(1)
imshow(derp_filenameGray)
title('Weighted (red) and Unweighted (blue) Centroids'); 
hold on
numObj = numel(sRegionProps);

for k = 1 : numObj
   % plot(sRegionProps(k).WeightedCentroid(1), sRegionProps(k).WeightedCentroid(2), 'r*');
    plot(sRegionProps(k).Centroid(1), sRegionProps(k).Centroid(2), 'bo');
end
hold off
