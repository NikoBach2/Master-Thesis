Functionclass.add_two_numbers(1, 2)
Functionclass.add_three_numbers(1, 2, 3)

%%
A = imread('DJI_0084.JPG');
DetectedPumpkins = imread('DJI_0486SinglepumpkinsredDoublePumpkinsBlue.PNG');

Channels = Functionclass.MarkedPumpkins(A);
ChannelsRotated = rot90(Channels);

arrayTest=[77,106,91];
%%
%d1 = mahal(arrayTest,ChannelsRotated)
progresscomment('Starting filtering pumpkins. Mahanalobis');
BinaryPumkinImage1 = Functionclass.FindPumpkinsInImage(A,ChannelsRotated);
progresscomment('Finish filtering pumpkins. Mahanalobis');
%%
figure(3)
imshow(BinaryPumkinImage2);

%%

BinaryPumkinImage2 =  Functionclass.FindPumpkinsInImageVersion2(A,ChannelsRotated);

%%

%Functionclass.FindPumpkinsInImageVersion2(A,ChannelsRotated);

%D = bwdist(~BinaryPumkinImage2); % image B (above)

%D = -bwdist(~BinaryPumkinImage2); % image C (above)

%I = rgb2gray(A);
[D1,IDX] = bwdist(~Gmag,'euclidean');

%I2 = imtophat(BinaryPumkinImage2, strel('disk', 20));


%%
se = strel('disk', 2);
Ie = imerode(BinaryPumkinImage2, se);
Iobr = imreconstruct(Ie, BinaryPumkinImage2);
%%
figure(4)
imshow(Iobr)

%%
figure(5)
imshow(D1,[],'InitialMagnification','fit')
title('Distance transform of ~bw')

%%
Functionclass.ColorBinaryCombinedImage = CombineBinaryAndColorImage(A,BinaryPumkinImage2);