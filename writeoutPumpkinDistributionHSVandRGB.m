%A = imread('DJI_0486.JPG');
A = imread('pumpkinPictures2017/DJI_0255Croped.JPG');
%DetectedPumpkins = imread('DJI_0486RedMarked.PNG');
%LightImage = A + 25;
%HSVImage=rgb2hsv(LightImage);
DetectedPumpkins = imread('pumpkinPictures2017/DJI_0255CropedMarkedPumpkins.png');
%TableName = 'DJI_0486ManualCountedDiameterFinal.csv';
%TableName = 'pumpkinPictures2017/DJI_612_MarkedPumpkinsDiameterFinalNogreen.csv';

%TableNameHenrik = 'DJI_0486_points.csv';

%T = readtable(TableName);
%%

LightImage = A + 100;
Channels = Functionclass.MarkedPumpkins(LightImage,DetectedPumpkins);
ChannelsRotated = rot90(Channels);
%ChannelsHSV = Functionclass.MarkedPumpkinsHSV(A,DetectedPumpkins);
%ChannelsRotatedHSV = rot90(ChannelsHSV);


%%
figure(2)
scatter3(ChannelsRotated(:,1),ChannelsRotated(:,2),ChannelsRotated(:,3))
xlabel('R') % x-axis label
ylabel('G') % y-axis label
zlabel('B') % Z-axis label
view(90, 90);