%allCenters;

%calculate distances of each point to each other point
dists = squareform( pdist(allCenters) );
%set diagonal to inf so point is not its own closest neighbor
dists(1:size(dists,1)+1:end) = inf;
%find index of closest point
[TY, minidx] = min(dists, [], 2);  % makes it n x 1 instead of n x 1. 
%[~, minidx] = min(dists);
%%
CloseObjects = allCenters(minidx,:);

DistanceToCloseObjects = sqrt((CloseObjects(:,1)-allCenters(:,1)).^2 + (CloseObjects(:,2)-allCenters(:,2)).^2);
 
TotalNumberOfPumpkinsX = zeros();
TotalNumberOfPumpkinsY = zeros();
TotalNumberOfPumpkinsXDouble = zeros();
TotalNumberOfPumpkinsYDouble = zeros();
for counter = 1:3737
if DistanceToCloseObjects(counter)>6
   TotalNumberOfPumpkinsX(end +1) = allCenters(counter,1);
   TotalNumberOfPumpkinsY(end +1) = allCenters(counter,2);
   
else
   allCenters(counter,1);
   allCenters(counter,2);
   TotalNumberOfPumpkinsXDouble(end +1) = allCenters(counter,1);
   TotalNumberOfPumpkinsYDouble(end +1) = allCenters(counter,2);
end   
end
