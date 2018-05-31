bw=imread('DJI_0486_croped.PNG'); % read image
bwMarker=imread('DJI_0486_croped_marked_pumpkins.PNG');% read marker

bw = bw(:,:,1);


bw(bwMarker(:,:,1)==255)=0; % impose the marker on the image
%%
L=watershed(bw);% do the watershed
%% plot
figure(1)
imshow(bw)

rgb = label2rgb(L,'jet',[.5 .5 .5]);
figure(2), imshow(0.5*rgb + bw*0.5,'InitialMagnification','fit')
title('Marker Controlled Watershed transform ')