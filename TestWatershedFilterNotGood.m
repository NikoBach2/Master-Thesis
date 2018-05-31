%%
I2 = A;
I2(Io) = 255;
figure
imshow(I2), title('Regional maxima superimposed on original image (I2)')

%%
se2 = strel(ones(5,5));
fgm2 = imclose(Io, se2);
fgm3 = imerode(fgm2, se2);
%%
fgm4 = bwareaopen(fgm3, 20);
I3 = A;
I3(fgm4) = 255;
figure
imshow(I3)
title('Modified regional maxima superimposed on original image (fgm4)')
%%
bw = imbinarize(Iobrcbr);
figure
imshow(bw), title('Thresholded opening-closing by reconstruction (bw)')