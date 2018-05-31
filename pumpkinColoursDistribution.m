PumpkinImage = imread('DJI_0084.JPG');
MarkedPumpkinImage = imread('DJI_0084RedPumpkins.JPG');

[row,col] = size(PumpkinImage(:,:,1));
RedChannel = zeros();
GreenChannel = zeros;
BlueChannel = zeros;
for j = 1:row
   for i = 1:col   
       MarkedAreaIsRedChannel = MarkedPumpkinImage(j,i,1);
       MarkedAreaIsGreenChannel = MarkedPumpkinImage(j,i,2);
             if MarkedAreaIsRedChannel == 255 && MarkedAreaIsGreenChannel <5
                 RedChannel(end +1) = PumpkinImage(j,i,1);
                 GreenChannel(end +1) = PumpkinImage(j,i,2);
                 BlueChannel(end +1) = PumpkinImage(j,i,3);
             end
   end
end
Matrix = [RedChannel; GreenChannel; BlueChannel];
scatter3(RedChannel,GreenChannel,BlueChannel)
xlabel('Red') % x-axis label
ylabel('Green') % y-axis label
zlabel('Blue') % y-axis label

