%load hospital;
X = TablePumpkinPositions; %manual counted
%Y = allCenters;   % program counted
Y = TablePumpkinPositionsHenrik;

[Idx, D] = knnsearch(X,Y);

X(Idx,:);

[rowKTree colKTree] = size(D)
S = sum(D)

AverageDistanceToPumpkinCenter = S/rowKTree  % pixels 


%%
figure(54)
imshow(A)
title('Weighted (red) and Unweighted (blue) Centroids'); 
hold on
numObj = numel(sRegionProps);

for k = 1 : numObj
   % XPumpkinposition, YPumpkinposition,
    plot(XPumpkinposition, YPumpkinposition, 'r*');
    plot(arrayX1Henrik, arrayY1Henrik, 'bo');
end
hold off
