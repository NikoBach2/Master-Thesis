TablePumpkinPositions;


Precision = 7; 
     [Idx, D] = knnsearch(TablePumpkinPositions,allCenters);
     
 %%
 
 
 
 
 %%
% Ag = [7 18 27 42 65 49 54 65 78 82 87 98];
[ng, bin] = histc(IdxInlierRotated, unique(IdxInlierRotated));
multiple = find(ng > 1);
index    = find(ismember(bin, multiple));   % her er arrayet sorteret 
     %%
     %A=[1;1;1;2;2;2;2;3;3;3];
     IdxInlierRotated(index)
%%
IdxClean = unique(IdxInlierRotated(index)); % which will give you the unique elements of A in array B

Ncount = histc(IdxInlierRotated(index), IdxClean); % this willgive the number of occurences of each unique element

[rowOriginal colOriginal] = size(IdxInlierRotated(index));
[rowCleaned colCleaned] = size(IdxClean);

NumberOfDoubleCountedPumpkins = rowOriginal - rowCleaned
%%
TablePumpkinPositions(IdxClean,:)

     %%
    TablePumpkinPositions(58,:)
    %%
    AreaDifferenceTotal = 0;
    [rowKTree colKTree] = size(D);
    OutLierPumpkinsX = [];
    OutLierPumpkinsY = [];
    %allCenters = zeros(row,2);
    IdxInlier = [];
    for i = 1:rowKTree
        if D(i) > Precision    % the maximum distance between pumpkins, before it is counted as outlier. 
        OutLierPumpkinsX(end + 1) = allCenters(i,1);
        OutLierPumpkinsY(end + 1) = allCenters(i,2);
        else
            RowNumber = Idx(i);
            AreaDifference = ManualCountedArea(RowNumber) - AllArea(i);
            IdxInlier(end +1) = Idx(i);
        end 
        AreaDifferenceTotal = AreaDifferenceTotal + AreaDifference; 
    end
    OutLierPumpkinsRotatedX = rot90(OutLierPumpkinsX);
    OutLierPumpkinsRotatedY = rot90(OutLierPumpkinsY);
    IdxInlierRotated = rot90(IdxInlier);
    AreaDifferenceTotal / rowKTree
    %%
 figure(56)
imshow(A)
title('Manual (red) and Henrik (blue) Centroids'); 
hold on
numObj = numel(TablePumpkinPositions(:,1))
numObj1 = numel(OutLierPumpkinsRotatedX);
%XPumpkinposition YPumpkinposition



for k = 1 : numObj
    plot(allCenters(k,1), allCenters(k,2), 'r*');
end
%for k = 1 : numObj1
    plot(TablePumpkinPositions(IdxClean,1), TablePumpkinPositions(IdxClean,2), 'bo');
%end
hold off
    f = numel(OutLierPumpkinsRotatedX);
    
    %%
    figure(100)
    cdfplot(D)
   % cdfplot(y)
    hold on
    xKmeans = -20:0.1:10;
    fKmeans = evcdf(D,0,3);
    plot(xKmeans,fKmeans,'m')
    %%
    wblsurv = 1-cdf('weibull',D,100,2);
    plot(D,wblsurv,'g-','LineWidth',2)
    legend('Empirical','LCB','UCB','Population',...
    'Location','NE')
