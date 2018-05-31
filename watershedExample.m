I = rgb2gray(ColorBinaryCombinedImage);


[Gmag, Gdir] = imgradient(I,'prewitt');

figure(1)
imshow(ColorBinaryCombinedImage)
%%
figure
imshowpair(Gmag, Gdir, 'montage');
title('Gradient Magnitude, Gmag (left), and Gradient Direction, Gdir (right), using Prewitt method')

%%
[Gx, Gy] = imgradientxy(I);
figure
imshowpair(Gx, Gy, 'montage')
title('Directional Gradients, Gx and Gy, using Sobel method')
%%
%ColoredDetectedPumpkins = DetectedPumpkins(:,:) > 0 & A(:,:,:)

%MarkedPumpkinsBinaryImage = MarkedPumpkinsImage(:,:,1) > 250 & MarkedPumpkinsImage(:,:,2) < 5;


imshow(I)

text(732,501,'Image courtesy of Corel(R)',...
     'FontSize',7,'HorizontalAlignment','right')
 
 
 hy = fspecial('sobel');
hx = hy';
Iy = imfilter(double(I), hy, 'replicate');
Ix = imfilter(double(I), hx, 'replicate');
gradmag = sqrt(Ix.^2 + Iy.^2);
figure
imshow(gradmag,[]), title('Gradient magnitude (gradmag)')

L = watershed(gradmag);
Lrgb = label2rgb(L);
figure, imshow(L), title('Watershed transform of gradient magnitude (Lrgb)')

%%
se = strel('disk', 5);
Io = imopen(I, se);
figure
imagesc(Io), title('Opening (Io)') , colormap(gray)

%%
Ie = imerode(I, se);
Iobr = imreconstruct(Ie, I);
figure
imshow(Iobr), title('Opening-by-reconstruction (Iobr)')


Ioc = imclose(Io, se);
figure
imshow(Ioc), title('Opening-closing (Ioc)')

%%

Iobrd = imdilate(Iobr, se);
Iobrcbr = imreconstruct(imcomplement(Iobrd), imcomplement(Iobr));
Iobrcbr = imcomplement(Iobrcbr);
figure
imshow(Iobrcbr), title('Opening-closing by reconstruction (Iobrcbr)')

fgm = BinaryImageInverted;
figure
imshow(fgm), title('Regional maxima of opening-closing by reconstruction (fgm)')

%%
I2 = I;
I2(fgm) = 255;
figure
imshow(I2), title('Regional maxima superimposed on original image (I2)')


se2 = strel(ones(5,5));
fgm2 = imclose(fgm, se2);
fgm3 = imerode(fgm2, se2);

%%
fgm4 = bwareaopen(fgm3, 20);
I3 = I;
I3(fgm4) = 255;
figure
imshow(I3)
title('Modified regional maxima superimposed on original image (fgm4)')

%%
bw = imbinarize(Iobrcbr);%BinaryImageInverted
figure
imshow(bw), title('Thresholded opening-closing by reconstruction (bw)')


%%

D = bwdist(bw);
%DL = watershed(D); % Gy
DL = watershed(Gmag);
bgm = DL == 0;
figure
imshow(bgm), title('Watershed ridge lines (bgm)')


%%
gradmag2 = imimposemin(gradmag, bgm | fgm4);


L = watershed(gradmag2);  %Gy
%L = watershed(Gy);

I4 = I;
I4(imdilate(L == 0, ones(3, 3)) | bgm | fgm4) = 255;
figure
imshow(I4)
title('Markers and object boundaries superimposed on original image (I4)')


Lrgb = label2rgb(L, 'jet', 'w', 'shuffle');
figure
imshow(Lrgb)
title('Colored watershed label matrix (Lrgb)')
%%
figure
imshow(I)
hold on
himage = imshow(Lrgb);
himage.AlphaData = 0.3;
title('Lrgb superimposed transparently on original image')

%%

se3 = strel('cube',3);

erodedBW = imerode(BinaryImageInverted,se);

erodedBW1 = imerode(erodedBW,se);
erodedBW2 = imerode(erodedBW1,se);
erodedBW3 = imerode(erodedBW2,se);
erodedBW4 = imerode(erodedBW3,se);
%erodedBW5 = imerode(erodedBW4,se);
%erodedBW6 = imerode(erodedBW5,se);
%erodedBW7 = imerode(erodedBW6,se);
%erodedBW8 = imerode(erodedBW7,se);


%%

figure(40)
p1=subplot(2, 1, 1);
imshow(I4);
%h = viscircles(MatrixPumpkins,RadiiRot);

p2=subplot(2, 1, 2);
imshow(Lrgb)
%h = viscircles(MatrixPumpkins,RadiiRot);

linkaxes([p1, p2]);