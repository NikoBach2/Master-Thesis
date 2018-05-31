

%lab_he = rgb2lab(A);
lab_he = rgb2lab(LightImage);

ab = lab_he(:,:,2:3);
nrows = size(ab,1);
ncols = size(ab,2);
ab = reshape(ab,nrows*ncols,2);

nColors = 2;
% repeat the clustering 3 times to avoid local minima
[cluster_idx, cluster_center] = kmeans(ab,nColors,'distance','sqEuclidean', ...
                                      'Replicates',3);
%%
pixel_labels = reshape(cluster_idx,nrows,ncols);
imshow(pixel_labels,[]), title('image labeled by cluster index');


%%
segmented_images = cell(1,3);
rgb_label = repmat(pixel_labels,[1 1 3]);

for k = 1:nColors
    color = A;
    color(rgb_label ~= k) = 0;
    segmented_images{k} = color;
end

%%
figure(1)
imshow(segmented_images{1}), title('objects in cluster 1');
%%
figure(2)
imshow(segmented_images{2}), title('objects in cluster 2');
%%
%figure(3)
%imshow(segmented_images{3});

%%
%figure(4)
%imshow(segmented_images{4}), title('objects in cluster 4');
%%
%figure(5)
%imshow(segmented_images{5}), title('objects in cluster 5');
%%
PumpkinsGrayScaledKMeans = rgb2gray(segmented_images{2});
BinaryImageKmean = PumpkinsGrayScaledKMeans > 6;

figure(4)
imshow(BinaryImageKmean)

%figure(5)
%imshow(PumpkinsGrayScaledKMeans)
% MarkedPumpkinsBinaryImage = MarkedPumpkinsImage(:,:,1) > 240 & MarkedPumpkinsImage(:,:,2) < 7;
