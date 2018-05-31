classdef AnalysisClass < handle
    %MAJORCLASS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods(Static)
function ReadTableAnd = ReadTableValues(TableName,TableNameHenrik,DiameterOn,A)
%T = readtable('pumpkinPictures2017/DJI_0079ManualCountedWithoutDiameter.csv');
%THenrik = readtable('DJI_0486_points.csv');
T = readtable(TableName);
THenrik = readtable(TableNameHenrik);

arrayX1Henrik = table2array(THenrik(:,3));
arrayY1Henrik = table2array(THenrik(:,4));
arraySizeHenrik = table2array(THenrik(:,5));
arraySizeRadius = arraySizeHenrik.*0.5;
arrayAngleHenrik = table2array(THenrik(:,6));

TablePumpkinPositionsHenrik = [arrayX1Henrik arrayY1Henrik];
new_x = arrayX1Henrik + (arraySizeRadius.*cos(arrayAngleHenrik)); % matlab kører baglæns 
new_y = arrayY1Henrik - (arraySizeRadius.*sin(arrayAngleHenrik));
TablePumpkinPositionsHenrik = [new_x new_y];

%arrayX1 = T(:,3);
arrayX1 = table2array(T(:,3));
arrayX2 = table2array(T(:,5));
XPumpkinposition = arrayX1 + ((arrayX2 - arrayX1)/2);
XPumpkinRadius = (arrayX2 - arrayX1)/2;

arrayY1 = table2array(T(:,4));
arrayY2 = table2array(T(:,6));
YPumpkinposition = arrayY1 + ((arrayY2 - arrayY1)/2);
YPumpkinDiameter = arrayY2 - arrayY1;

if DiameterOn == 1
TablePumpkinPositions = [XPumpkinposition YPumpkinposition];
else
TablePumpkinPositions = [table2array(T(:,3)) table2array(T(:,4))];
end

'False Positive'
Name1 = '1';
Name2 = '2';
%Functionclass.FindOutlierPumpkins(TablePumpkinPositions,allCenters,A,7,Name1)
Functionclass.FindOutlierPumpkins(TablePumpkinPositions,TablePumpkinPositionsHenrik,A,7,Name1)

% TablePumpkinPositionsHenrik
'False Negative'
%Functionclass.FindOutlierPumpkins(allCenters,TablePumpkinPositions,A,7,Name2)
Functionclass.FindOutlierPumpkins(TablePumpkinPositionsHenrik,TablePumpkinPositions,A,7,Name2)




ReadTableAnd = TablePumpkinPositions;
end

function f = ManualCountedPumpkinArea(TableName, TableNameHenrik)
      %  T = readtable('pumpkinPictures2017/DJI_0079ManualCountedWithoutDiameter.csv');
%THenrik = readtable('DJI_0486_points.csv');
 T = readtable(TableName);
THenrik = readtable(TableNameHenrik);


arrayX1Henrik = table2array(THenrik(:,3));
arrayY1Henrik = table2array(THenrik(:,4));
arraySizeHenrik = table2array(THenrik(:,5));
arraySizeRadius = arraySizeHenrik.*0.5;
arrayAngleHenrik = table2array(THenrik(:,6));
TablePumpkinPositionsHenrik = [arrayX1Henrik arrayY1Henrik];


new_x = arrayX1Henrik + (arraySizeRadius.*cos(arrayAngleHenrik)); % matlab kører baglæns 
new_y = arrayY1Henrik - (arraySizeRadius.*sin(arrayAngleHenrik));
TablePumpkinPositionsHenrik = [new_x new_y];

%arrayX1 = T(:,3);
arrayX1 = table2array(T(:,3));
arrayX2 = table2array(T(:,5));
XPumpkinposition = arrayX1 + ((arrayX2 - arrayX1)/2);
XPumpkinRadius = (arrayX2 - arrayX1)/2;

arrayY1 = table2array(T(:,4));
arrayY2 = table2array(T(:,6));
YPumpkinposition = arrayY1 + ((arrayY2 - arrayY1)/2);
YPumpkinDiameter = arrayY2 - arrayY1;

arrayX2(isnan(arrayX2)) = [0];
arrayY2(isnan(arrayY2)) = [0];
Radius = sqrt((arrayX2-arrayX1).^2 + (arrayY2-arrayY1).^2)/2;
ManualCountedArea = pi * (Radius.^2);

f = ManualCountedArea;
end
        
        function f = OutliersAndArea(TablePumpkinPositions,allCenters,ManualCountedArea,AllArea)
            [Idx, D] = knnsearch(TablePumpkinPositions,allCenters);
            Precision = 7; 
            
           
              AreaDifferenceTotal = 0;
              AreaDifference = 0;
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
     [ng, bin] = histc(IdxInlierRotated, unique(IdxInlierRotated));
multiple = find(ng > 1);
index    = find(ismember(bin, multiple));
            IdxClean = unique(IdxInlierRotated(index)); % which will give you the unique elements of A in array B

Ncount = histc(IdxInlierRotated(index), IdxClean); % this willgive the number of occurences of each unique element

[rowOriginal colOriginal] = size(IdxInlierRotated(index));
[rowCleaned colCleaned] = size(IdxClean);

NumberOfDoubleCountedPumpkins = rowOriginal - rowCleaned
    f = AreaDifferenceTotal / rowKTree;
        end
    end
    
end