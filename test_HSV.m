%RGB=reshape(ones(64,1)*reshape(jet(64),1,192),[64,64,3]);

HSV=rgb2hsv(ChannelsRotated);
%H=HSV(:,:,1);
%S=HSV(:,:,2);
%V=HSV(:,:,3);
%subplot(2,2,1), imshow(H)
%subplot(2,2,2), imshow(S)
%subplot(2,2,3), imshow(V)
%subplot(2,2,4), imshow(RGB)

% 360 er skaleret ned til mellem 0 og 1 

%%
%ChannelsHSV = Functionclass.MarkedPumpkinsHSV(A,DetectedPumpkins);
ChannelsHSV = Functionclass.MarkedPumpkinsHSV(A,B);
ChannelsRotatedHSV = rot90(ChannelsHSV);
Channels = Functionclass.MarkedPumpkins(A,B);
ChannelsRotated = rot90(Channels);
%%
figure(3)
scatter3(ChannelsRotatedHSV(:,1),ChannelsRotatedHSV(:,2),ChannelsRotatedHSV(:,3))
xlabel('H') % x-axis label
ylabel('S') % y-axis label
zlabel('V') % y-axis label

figure(4)
scatter3(ChannelsRotated(:,1),ChannelsRotated(:,2),ChannelsRotated(:,3))
xlabel('RED') % x-axis label
ylabel('GREEN') % y-axis label
zlabel('BLUE') % y-axis label

%%

LightImage = A + 200;

figure(89)
imshow(LightImage)

%%
H = fspecial('disk',15);
BlurImage = imfilter(A,H);

figure(90)
imshow(BlurImage)