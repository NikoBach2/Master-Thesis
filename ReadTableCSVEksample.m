%%
%filename = fullfile(matlabroot,'examples','matlab','DJI_0486ManualCountedDiameterFinal.csv');
T = readtable('pumpkinPictures2017/DJI_0079ManualCountedWithoutDiameter.csv');
THenrik = readtable('DJI_0486_points.csv');

arrayX1Henrik = table2array(THenrik(:,3));
arrayY1Henrik = table2array(THenrik(:,4));
arraySizeHenrik = table2array(THenrik(:,5));
arraySizeRadius = arraySizeHenrik.*0.5;
arrayAngleHenrik = table2array(THenrik(:,6));
TablePumpkinPositionsHenrik = [arrayX1Henrik arrayY1Henrik];

%%
new_x = arrayX1Henrik + (arraySizeRadius.*cos(arrayAngleHenrik)); % matlab kører baglæns 
new_y = arrayY1Henrik - (arraySizeRadius.*sin(arrayAngleHenrik));
TablePumpkinPositionsHenrik = [new_x new_y];
%%
%arrayX1 = T(:,3);
arrayX1 = table2array(T(:,3));
arrayX2 = table2array(T(:,5));
XPumpkinposition = arrayX1 + ((arrayX2 - arrayX1)/2);
XPumpkinRadius = (arrayX2 - arrayX1)/2;

arrayY1 = table2array(T(:,4));
arrayY2 = table2array(T(:,6));
YPumpkinposition = arrayY1 + ((arrayY2 - arrayY1)/2);
YPumpkinDiameter = arrayY2 - arrayY1;
%%
arrayX2(isnan(arrayX2)) = [0];
arrayY2(isnan(arrayY2)) = [0];
Radius = sqrt((arrayX2-arrayX1).^2 + (arrayY2-arrayY1).^2)/2;
ManualCountedArea = pi * (Radius.^2);  % array 
ManualCountedSumArea = sum(ManualCountedArea)

%%



%%
arrayX2(240)
arrayX1(240)
arrayY2(240)
arrayY1(240)
%%
Radius = sqrt((arrayX2-arrayX1).^2 + (arrayY2(240)-arrayY1(240))^2)/2
%%
 nnn = isnan(ManualCountedArea);
 
 if nnn == 1
     nnn
     
 end 
%%
TablePumpkinPositions = [XPumpkinposition YPumpkinposition];
%TablePumpkinPositions = [table2array(T(:,3)) table2array(T(:,4))];

%%
figure(11)
yourImage = A;
imshow(yourImage);
hold on;
plot(new_x, new_y, 'r*', 'LineWidth', 2, 'MarkerSize', 15);
title('Bottle wit coordinates on it', 'FontSize', 24);

%%
 figure(53)
imshow(A)
title('Manual (red) and Henrik (blue) Centroids'); 
hold on
numObj = numel(sRegionProps);
%XPumpkinposition YPumpkinposition
for k = 1 : 3693
    plot(new_x(k), new_y(k), 'r*');
   % plot(OutLierPumpkinsRotatedX(k), OutLierPumpkinsRotatedY(k), 'bo');
end

hold off
%%

yourImage = A;
imshow(yourImage);
hold on;
plot(arrayX2, arrayY2, 'r*', 'LineWidth', 2, 'MarkerSize', 15);
title('Bottle wit coordinates on it', 'FontSize', 24);