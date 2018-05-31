classdef BlockProcessFunc < handle
    %MAJORCLASS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods(Static)
        function f = BlockFunc(A,ChannelsRotatedHSV)

RGB = A(:,:,1:3);


BinaryPumkinImage2 =  Functionclass.FindPumpkinsInImageVersion2(RGB,ChannelsRotatedHSV);
%HSVImage=rgb2hsv(RGB);



%BinaryPumkinImage2 =  Functionclass.FindPumpkinsInImageVersion2HSV(HSVImage,ChannelsRotatedHSV);


seDisk = strel('disk',2);
IM2 = imdilate(BinaryPumkinImage2,seDisk);   % my filter example


BinaryPumkinImage2Filled = imfill(IM2,'holes');

seDisk = strel('disk',2);
BinaryPumkinImage2FilledErode = imerode(BinaryPumkinImage2Filled,seDisk);

se1 = strel('disk', 3);
CleanedBinaryImage = imopen(BinaryPumkinImage2FilledErode, se1);

ColorBinaryCombinedImage = Functionclass.CombineBinaryAndColorImage(RGB,CleanedBinaryImage);
ColorBinaryCombinedGrayscaleImage = rgb2gray(ColorBinaryCombinedImage);

se = strel('disk', 3);  %% tuned to 3
Ie = imerode(ColorBinaryCombinedGrayscaleImage, se);
Iobr2 = imreconstruct(Ie, ColorBinaryCombinedGrayscaleImage);

hy = fspecial('sobel');
hx = hy';
Iy = imfilter(double(Iobr2), hy, 'replicate');
Ix = imfilter(double(Iobr2), hx, 'replicate');
gradmag = sqrt(Ix.^2 + Iy.^2);

seDisk = strel('disk',1);
IM2 = imdilate(gradmag,seDisk);

seDisk1 = strel('disk',1);
IM3 = imerode(IM2,seDisk1);

[D1,IDX] = bwdist(~IM3,'euclidean');

L = watershed(D1);
 Lrgb = label2rgb(L);
 
 NewBinaryPumkinImage =  L > 6; % not bigger than 12

f = NewBinaryPumkinImage;

        end
        
         function f = KmeanFunc(A)
            RGB = A(:,:,1:3);
            
            lab_he = rgb2lab(RGB);

            ab = lab_he(:,:,2:3);
            nrows = size(ab,1);
            ncols = size(ab,2);
            ab = reshape(ab,nrows*ncols,2);

            nColors = 2;
            % repeat the clustering 3 times to avoid local minima
            [cluster_idx, cluster_center] = kmeans(ab,nColors,'distance','sqEuclidean', ...
                                      'Replicates',3);
            pixel_labels = reshape(cluster_idx,nrows,ncols);
            
            segmented_images = cell(1,3);
            rgb_label = repmat(pixel_labels,[1 1 3]);

            for k = 1:nColors
                color = RGB;
                color(rgb_label ~= k) = 0;
                segmented_images{k} = color;
            end
            
            PumpkinsGrayScaledKMeans = rgb2gray(segmented_images{2});
            BinaryImageKmean = PumpkinsGrayScaledKMeans > 6;
            f = BinaryImageKmean;
         end
        
        function f = RGB2gr(A)
            RGB = A(:,:,1:3);
            
            
            
            f = rgb2gray(RGB);
            
        end
    end
    
end


