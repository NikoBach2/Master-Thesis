    [Idx, D] = knnsearch(TablePumpkinPositions,TablePumpkinPositionsHenrik);
    TablePumpkinPositions(Idx,:);
    [rowKTree colKTree] = size(D);
    OutLierPumpkinsX = []
    OutLierPumpkinsY = []
    %allCenters = zeros(row,2);
    for i = 1:rowKTree
        if D(i) > 7    % the maximum distance between pumpkins, before it is counted as outlier. 
        
        OutLierPumpkinsX(end + 1) = TablePumpkinPositionsHenrik(i,1);
        OutLierPumpkinsY(end + 1) = TablePumpkinPositionsHenrik(i,2);
        end    
    end
    OutLierPumpkinsRotatedX = rot90(OutLierPumpkinsX);
    OutLierPumpkinsRotatedY = rot90(OutLierPumpkinsY);
    %%
   % OutLierPumpkins = [OutLierPumpkinsRotatedX OutLierPumpkinsRotatedY]
 % TablePumpkinPositions(1,2)  
 
 figure(53)
imshow(A)
title('Manual (red) and Henrik (blue) Centroids'); 
hold on
%%
numObj = numel(OutLierPumpkinsRotatedX)
%%
%XPumpkinposition YPumpkinposition
for k = 1 : 3631
    plot(XPumpkinposition(k), YPumpkinposition(k), 'r*');
   % plot(OutLierPumpkinsRotatedX(k), OutLierPumpkinsRotatedY(k), 'bo');
end
for k = 1 : 1334
   % plot(XPumpkinposition(k), YPumpkinposition(k), 'r*');
    plot(OutLierPumpkinsRotatedX(k), OutLierPumpkinsRotatedY(k), 'bo');
end
hold off
  
