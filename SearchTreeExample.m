Mdl = KDTreeSearcher(TablePumpkinPositions);

%%
%allCenters = sRegionProps(2).Centroid
%sRegionProps(1).Centroid
[row col] = size(sRegionProps)
allCenters = zeros(row,2);
%allCenters = sRegionProps(i).Centroid;
%allCenters1(1) = allCenters
sRegionProps(1).Centroid(1);
sRegionProps(1).Centroid(2);
for i = 1:size(sRegionProps)
   allCenters(i,1) = sRegionProps(i).Centroid(1);
   allCenters(i,2) = sRegionProps(i).Centroid(2);
end

%%
%Y = 
Idx = knnsearch(Mdl,allCenters);



