classdef Functionclass < handle
    %MAJORCLASS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods(Static)
        function res = add_two_numbers(a, b)
            res = a + b;
        end
        function res = add_three_numbers(a, b, c)
            res = Functionclass.add_two_numbers(a, b);
            res = Functionclass.add_two_numbers(res, c);
        end
        
        
        
        function f = FindPumpkin(n1,n2,x11)   % find out, if the pixel is in a pumpkin
coutner = 1;
%THis function checks if the chunk of pixels should be a pumpkin.
WindowSize = 20;
    for nq = 1:WindowSize
        for j = 1:WindowSize
          %  if x11(n1+nq,n2+j,1) >1 && x11(n1+nq,n2+j,3) > 1 % hvis koden melder fejl, 
                 %lad den så kører igennem uden funktionen, da den ikke har bygget matricen x11 færdig. 
                coutner = coutner + 1 ;
           %  end
        end
    end
    
    if coutner > 7    % if more than 1/8 is orange
        pumpkin = coutner;  % there is a pumpkin   ændrer tilbage til pumpkin = 1
        %Centerxy(end + 1) = n1,n2;
        %CentX(end + 1) = n1;
        %CentY(end + 1) = n2;
        %Weight(end + 1 )=pumpkin;
    else 
        pumpkin = 0;   % there is not a pumpkin 
    end
    f = pumpkin; %prod(1:n);  % f skal være antal sammenhængende pixels.
    % koerer et vindue henover billedet, og find ud af hvor stor procentdel,
    % der er orange. hvis stor procent = graeskar
end


function f = PumpkinFilter(imageHue, imageSaturation, imageValue)
[row,col] = size(imageHue);
Zeroes2 = [0;0;0 ];
White2 = [255;255;255 ];
for n = 1:row
    for ni = 1:col
        Hue = imageHue(n,ni);
        Saturation = imageSaturation(n,ni);
        Value = imageValue(n,ni);
       if 1
       if Hue < 59 && Saturation > 0.2 && Value >0.2 && (Saturation > 0.5 || Value >0.5)   % 64 
           
           x1bin(n,ni) = 1;
           array(n,ni) = n;
       else
        
           x1bin(n,ni) = 0;
           array(n,ni) = 0;
       end
       end
    end
   
end
f = x1bin;
end

function f = GaussianFunction(image)
    [row,col] = size(image);
    GaussianImage = size(image);
    %deviation = 1;
    for n = 1:row
       for ni = 1:col
           if n> 4 && n< row-4 && ni>4 && ni<col-4
           Mask = [image(n-2,ni-2) image(n-1,ni-2) image(n,ni-2) image(n+1,ni-2) image(n+2,ni-2); image(n-2,ni-1) image(n-1,ni-1) image(n,ni-1) image(n+1,ni-1) image(n+2,ni-1); image(n-2,ni) image(n-1,ni) image(n,ni) image(n+1,ni) image(n+2,ni); image(n-2,ni+1) image(n-1,ni+1) image(n,ni+1) image(n+1,ni+1) image(n+2,ni+1); image(n-2,ni+2) image(n-1,ni+2) image(n,ni+2) image(n+1,ni+2) image(n+2,ni+2)];
           %image(n,ni)
           
          % image(n-2,ni)
          % Mask = [image(n,ni-2)];
           GaussianImage(n,ni) = Functionclass.GaussianKernel(Mask);
           
           else 
               GaussianImage(n,ni) = 0;
           end
           
         % Gauss = 1/(2*pi*deviation) * exp(-((n^2 + ni^2 )/(2*deviation^2))); 
       end
    end
    
    f = GaussianImage ;
end

function f = GaussianKernel(mask)
          
    Kernel = [1 4 7 4 1; 4 16 26 16 4; 7 26 41 26 7; 4 16 26 16 4; 1 4 7 4 1];
           % 1  4  7  4 1
           % 4 16 26 16 4
           % 7 26 41 26 7
           % 4 16 26 16 4
           % 1  4  7  4 1
    NewValue = Kernel * mask;
    S = sum(NewValue);
    Q = sum(S);
    ReturnValue = Q / 25;
    f = ReturnValue;
end

function f = MarkedPumpkins(PumpkinImage,MarkedPumpkinImage)
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
%CombinedColorChannels = cat(3,RedChannel,GreenChannel,BlueChannel);
CombinedColorChannels = [RedChannel; GreenChannel; BlueChannel];
f = CombinedColorChannels;

end

function f = MarkedPumpkinsHSV(PumpkinImage,MarkedPumpkinImage)
[row,col] = size(PumpkinImage(:,:,1));
HChannel = zeros();
SChannel = zeros;
VChannel = zeros;
HSVPumpkininImage = rgb2hsv(PumpkinImage);
%HSVMarkedPumpinImage =rgb2hsv(MarkedPumpkinImage);
for j = 1:row
   for i = 1:col   
       MarkedAreaIsRedChannel = MarkedPumpkinImage(j,i,1);
       MarkedAreaIsGreenChannel = MarkedPumpkinImage(j,i,2);
             if MarkedAreaIsRedChannel == 255 && MarkedAreaIsGreenChannel <5
                 HChannel(end +1) = HSVPumpkininImage(j,i,1);
                 SChannel(end +1) = HSVPumpkininImage(j,i,2);
                 VChannel(end +1) = HSVPumpkininImage(j,i,3);
             end
   end
end
%CombinedColorChannels = cat(3,RedChannel,GreenChannel,BlueChannel);
CombinedHSVChannels = [HChannel; SChannel; VChannel];
f = CombinedHSVChannels;

end

function f = MahanalobisDistToPumpkinDistribution(Point,PumpkinDistribution)
    arrayTest=[91,106,77]; % is just made in order to cast the array Point. 
    PointCastet = cast(Point,'like',arrayTest);  % typecast til double
    MahanalobisDistance = mahal(PointCastet,PumpkinDistribution);
    f = MahanalobisDistance;
end

function f = FindPumpkinsInImage(PumpkinImage,PumpkinDistribution)
    [row,col] = size(PumpkinImage(:,:,1));
    BinaryPumkinImage = zeros(row,col);
    
    for i = 1:row
        for j = 1:col
             RGBPoint = [PumpkinImage(i,j,1),PumpkinImage(i,j,2),PumpkinImage(i,j,3)];
             DistanceToPumpkinDistribution = Functionclass.MahanalobisDistToPumpkinDistribution(RGBPoint,PumpkinDistribution);
             if DistanceToPumpkinDistribution < 6
                 BinaryPumkinImage(i,j) = 255;
             else
                 BinaryPumkinImage(i,j) = 0;
             end    
        end
    end
    
    f = BinaryPumkinImage;
end

function f = FindPumpkinsInImageVersion2(PumpkinImage,PumpkinDistribution)
    [row,col] = size(PumpkinImage(:,:,1));
    MatriceY = reshape(PumpkinImage,[],3);
    % keyboard
    Matrice2Double = im2double(MatriceY);
    PumpkinFloating = PumpkinDistribution / 255;
    %PumpkinFloating = PumpkinDistribution / 255.;
    MahalDistanceMatrice = mahal(Matrice2Double,PumpkinFloating);
    BinaryPumkinImage =  MahalDistanceMatrice < 6; % not bigger than 6  
    
    f = reshape(BinaryPumkinImage,row,col);            
end

function f = FindPumpkinsInImageVersion2HSV(PumpkinImage,PumpkinDistribution)
    [row,col] = size(PumpkinImage(:,:,1));
    MatriceY = reshape(PumpkinImage,[],3);
    % keyboard
    Matrice2Double = im2double(MatriceY);
    PumpkinFloating = PumpkinDistribution / 255;
    %MahalDistanceMatrice = mahal(Matrice2Double,PumpkinFloating);
    MahalDistanceMatrice = mahal(Matrice2Double,PumpkinDistribution);
    %keyboard
    BinaryPumkinImage =  MahalDistanceMatrice < 6; % not bigger than 6  
    
    f = reshape(BinaryPumkinImage,row,col);            
end

function f = FindDoubleCountedPumpkins(XManualCounted)
    
    [Idx, D] = knnsearch(XManualCounted,XManualCounted);
    XManualCounted(Idx,:);
    [rowKTree colKTree] = size(D);
    
    %SumOfDistances = sum(D);
    %AverageDistanceToPumpkinCenter = SumOfDistances/rowKTree;  % pixels 
    f = AverageDistanceToPumpkinCenter;
    
    
end 

function f = FindPumpkinsInImageVersion3(PumpkinImage,PumpkinDistribution,Morph) % is the same as version2, but with bigger mahal distance
    [row,col] = size(PumpkinImage(:,:,1));
    MatriceY = reshape(PumpkinImage,[],3);
    % keyboard
    Matrice2Double = im2double(MatriceY);
    PumpkinFloating = PumpkinDistribution / 255;
    MahalDistanceMatrice = mahal(Matrice2Double,PumpkinFloating);
    BinaryPumkinImage =  MahalDistanceMatrice < 6; % not bigger than 6  
    
    f = reshape(BinaryPumkinImage,row,col);            
end

function f = GoThroughImage(Image,PCADist) 
    
[row,col] = size(Image(:,:,1));   
    
    for i = 20:row-20
        for j = 20:col-20
            Window = PCAWindow(i,j,Image);
           % PCAWindow = pca(Window);
            
        end
    end
    % calculate PCA
end

function f = PCAWindow(x,y,Image) 
    
    WindowSize = 20;
    NewMatrix = zeroes(WindowSize,WindowSize);
    for nq = 1:WindowSize
        for j = 1:WindowSize
            NewMatrix(i,j) = Image(x+i,y+j); % matrix of the what values is in the window.
        end
    end
    f = NewMatrix;
    
end

function f = FinAverageAreaOfDoubleMarkedPumpkins(MarkedPumpkinsImage)
         [row,col] = size(MarkedPumpkinsImage)
         MarkedPumpkinsBinaryImage = MarkedPumpkinsImage(:,:,3) > 250 & MarkedPumpkinsImage(:,:,2) < 5;
         [rowBin,colBin] = size(MarkedPumpkinsBinaryImage)
         GrayScaleMarkedPumpkinImage = rgb2gray(MarkedPumpkinsImage);
         [rowGray,colGray] = size(GrayScaleMarkedPumpkinImage)
         AreaAndPerimeterOfMarkedPumpkins = regionprops(MarkedPumpkinsBinaryImage, GrayScaleMarkedPumpkinImage, {'Area', 'Perimeter','Centroid'});
         f = AreaAndPerimeterOfMarkedPumpkins;

end


function f = FindAverageAreaOfMarkedPumpkins(MarkedPumpkinsImage)
        [row,col] = size(MarkedPumpkinsImage)
         MarkedPumpkinsBinaryImage = MarkedPumpkinsImage(:,:,1) > 240 & MarkedPumpkinsImage(:,:,2) < 7;
         [rowBin,colBin] = size(MarkedPumpkinsBinaryImage)
         GrayScaleMarkedPumpkinImage = rgb2gray(MarkedPumpkinsImage);
         [rowGray,colGray] = size(GrayScaleMarkedPumpkinImage)
         AreaAndPerimeterOfMarkedPumpkins = regionprops(MarkedPumpkinsBinaryImage, GrayScaleMarkedPumpkinImage, {'Area', 'Perimeter','Centroid'});
         f = AreaAndPerimeterOfMarkedPumpkins;
end

function f = LoopThroughPumpkinsListAndFindClusterPumpkins(AreaPumpkinDouble,AreaPumpkinSingle,ListOfpumpkins)
        [row,col] = size(ListOfpumpkins);
       % NumberOfPumpkinsOfEachCoordinate = [col,1];
        for i = 1:col
           ListOfpumpkins(i);
            if ListOfpumpkins(i) > (AreaPumpkinSingle/2) 
            NumberOfPumpkinsOfEachCoordinate(i,1) = round(ListOfpumpkins(i)/AreaPumpkinSingle);
            NumberOfPumpkinsOfEachCoordinate(i,2) = ListOfpumpkins(i);
            else
            NumberOfPumpkinsOfEachCoordinate(i,1) = 1;
            NumberOfPumpkinsOfEachCoordinate(i,2) = ListOfpumpkins(i); %area   
                
            end
        end
         f = NumberOfPumpkinsOfEachCoordinate;
end

function f = CombineBinaryAndColorImage(PumpkinImage,BinaryImage)
    [row,col] = size(PumpkinImage(:,:,1));
    ColoredDetectedPumpkins = zeros(row,col,3);
    ZeroesForBackground = [0 0 0];
    PumpkinImage(5,1,:);
    for i = 1:row
        for j = 1:col
            if BinaryImage(i,j) == 0
           %    PumpkinImage(i,j,:) = PumpkinImage(i,j,:);
           % else
               PumpkinImage(i,j,1) = 0 ;
               PumpkinImage(i,j,2) = 0 ;
               PumpkinImage(i,j,3) = 0 ;
            end
        end
    end
    
    f = PumpkinImage;
    
end


function f = SaveFigure(Figure,xmin,ymin,width, height,filename)
    rect = [xmin ymin width height];
    I2 = imcrop(Figure,rect);
    imwrite(I2,filename);
end

function f = SaveFigures(Img,Folder,Showing)
   % Showing = 'Original' 
   % folder = 'TestImage'
   % filename=sprintf('%s/%sPumpkins3.PNG',Folder,k)
filename=sprintf('TestImage/%s/Pumpkins1%s.PNG',Folder,Showing);
Functionclass.SaveFigure(Img,290,320,100, 100,filename);
filename=sprintf('TestImage/%s/Pumpkins2%s.PNG',Folder,Showing);
Functionclass.SaveFigure(Img,450,480,100, 100,filename);
filename=sprintf('TestImage/%s/Pumpkins3%s.PNG',Folder,Showing);
Functionclass.SaveFigure(Img,170,30,100, 100,filename);
filename=sprintf('TestImage/%s/Pumpkins4%s.PNG',Folder,Showing);
Functionclass.SaveFigure(Img,1750,220,100, 100,filename);
filename=sprintf('TestImage/%s/Pumpkins5%s.PNG',Folder,Showing);
Functionclass.SaveFigure(Img,1810,350,100, 100,filename);
filename=sprintf('TestImage/%s/Pumpkins6%s.PNG',Folder,Showing);
Functionclass.SaveFigure(Img,3690,2550,150, 150,filename);
filename=sprintf('TestImage/%s/Pumpkins7%s.PNG',Folder,Showing);
Functionclass.SaveFigure(Img,2520,2850,150, 150,filename);
filename=sprintf('TestImage/%s/Pumpkins8%s.PNG',Folder,Showing);
Functionclass.SaveFigure(Img,1100,2740,150, 150,filename);
filename=sprintf('TestImage/%s/Pumpkins9%s.PNG',Folder,Showing);
Functionclass.SaveFigure(Img,80,2790,150, 150,filename);
filename=sprintf('TestImage/%s/Pumpkins10%s.PNG',Folder,Showing);
Functionclass.SaveFigure(Img,60,2490,150, 150,filename);
filename=sprintf('TestImage/%s/Pumpkins11%s.PNG',Folder,Showing);
Functionclass.SaveFigure(Img,2000,490,100, 100,filename);
filename=sprintf('TestImage/%s/Pumpkins12%s.PNG',Folder,Showing);
Functionclass.SaveFigure(Img,2800,20,100, 100,filename);

filename=sprintf('TestImage/%s/Pumpkins13%s.PNG',Folder,Showing);
Functionclass.SaveFigure(Img,1090,1042,400,400,filename);
filename=sprintf('TestImage/%s/Pumpkins14%s.PNG',Folder,Showing);
Functionclass.SaveFigure(Img,1980,2160,400, 400,filename);
filename=sprintf('TestImage/%s/Pumpkins15%s.PNG',Folder,Showing);
Functionclass.SaveFigure(Img,3460,1690,400, 400,filename);
filename=sprintf('TestImage/%s/Pumpkins16%s.PNG',Folder,Showing);
Functionclass.SaveFigure(Img,3350,1080,400, 400,filename);
filename=sprintf('TestImage/%s/Pumpkins17%s.PNG',Folder,Showing);
Functionclass.SaveFigure(Img,1150,2440,400, 400,filename);
    
end

function f = KNNTreeFindAverageDistance(XManualCounted,YProgramCounted)
    %X = TablePumpkinPositions; %manual counted
    %Y = allCenters;   % program counted
    %Y = TablePumpkinPositionsHenrik;

    [Idx, D] = knnsearch(XManualCounted,YProgramCounted);
    XManualCounted(Idx,:);
    [rowKTree colKTree] = size(D);
    SumOfDistances = sum(D);
    AverageDistanceToPumpkinCenter = SumOfDistances/rowKTree;  % pixels 
    f = AverageDistanceToPumpkinCenter;
end 

function f = KNNTreeFindDifferenceInManuelAndProgramCounted(XManualCounted,YProgramCounted)
    
    [Idx, D] = knnsearch(XManualCounted,YProgramCounted);
    XManualCounted(Idx,:);
    [rowKTree colKTree] = size(D);
    
    %SumOfDistances = sum(D);
    %AverageDistanceToPumpkinCenter = SumOfDistances/rowKTree;  % pixels 
    f = AverageDistanceToPumpkinCenter;
end 

function f = ConvertFromStructToArray(sRegionProps)
[row col] = size(sRegionProps)
allCenters = zeros(row,2);
sRegionProps(1).Centroid(1);
sRegionProps(1).Centroid(2);
for i = 1:size(sRegionProps)
   allCenters(i,1) = sRegionProps(i).Centroid(1);
   allCenters(i,2) = sRegionProps(i).Centroid(2);
end
f = allCenters;
end 
function f = ConvertFromStructToArrayArea(sRegionProps)
[row col] = size(sRegionProps)
AllArea = zeros(row,1);
%sRegionProps(1).Area(1);
%sRegionProps(1).Centroid(2);
for i = 1:size(sRegionProps)
   AllArea(i,1) = sRegionProps(i).Area;
  % AllArea(i,2) = sRegionProps(i).Centroid(2);
end
f = AllArea;
end 
function f = FindOutlierPumpkins(TablePumpkinPositions,TablePumpkinPositionsHenrik,A,Precision,name)
     [Idx, D] = knnsearch(TablePumpkinPositions,TablePumpkinPositionsHenrik);
    TablePumpkinPositions(Idx,:);
    [rowKTree colKTree] = size(D);
    OutLierPumpkinsX = [];
    OutLierPumpkinsY = [];
    %allCenters = zeros(row,2);
    for i = 1:rowKTree
        if D(i) > Precision    % the maximum distance between pumpkins, before it is counted as outlier. 
        OutLierPumpkinsX(end + 1) = TablePumpkinPositionsHenrik(i,1);
        OutLierPumpkinsY(end + 1) = TablePumpkinPositionsHenrik(i,2);
        end    
    end
    OutLierPumpkinsRotatedX = rot90(OutLierPumpkinsX);
    OutLierPumpkinsRotatedY = rot90(OutLierPumpkinsY);
figure()
imshow(A)
title('Manual (red) and Henrik (blue) Centroids'); 
hold on
numObj = numel(TablePumpkinPositions(:,1));
numObj1 = numel(OutLierPumpkinsRotatedX);
%XPumpkinposition YPumpkinposition
for k = 1 : numObj
    plot(TablePumpkinPositions(k,1), TablePumpkinPositions(k,2), 'r*');
end
for k = 1 : numObj1
    plot(OutLierPumpkinsRotatedX(k), OutLierPumpkinsRotatedY(k), 'bo');
end
hold off
    NumberOfOutLierPumpkins = numel(OutLierPumpkinsRotatedX)
    f = NumberOfOutLierPumpkins;
end 
  
function f = SuperPixelsFunction(A,ChannelsRotated)
   %[L,N] = superpixels(A,1000000);
   [L,N] = superpixels(A,100000);
BW1 = boundarymask(L);
bws = bwmorph(BW1,'skel');
figure(1)
imshow(imoverlay(A,bws,'cyan'),'InitialMagnification',67) 


outputImage = zeros(size(A),'like',A);
idx = label2idx(L);
numRows = size(A,1);
numCols = size(A,2);

for labelVal = 1:N
    redIdx = idx{labelVal};
    greenIdx = idx{labelVal}+numRows*numCols;
    blueIdx = idx{labelVal}+2*numRows*numCols;
    ColorsCombined = [A(redIdx),A(greenIdx),A(blueIdx)];
    PumpkinFloating = ChannelsRotated / 255;
    ColorsCombinedDouble = im2double(ColorsCombined);
    MahalDistanceMatrice = mahal(ColorsCombinedDouble,PumpkinFloating);
    % =  Functionclass.FindPumpkinsInImageVersion2(ColorsCombined,ChannelsRotated)
    if MahalDistanceMatrice < 6
    outputImage(redIdx) = mean(A(redIdx));
    outputImage(greenIdx) = mean(A(greenIdx));
    outputImage(blueIdx) = mean(A(blueIdx));
    
    else
        outputImage(redIdx) = 0;
        outputImage(greenIdx) = 0;
        outputImage(blueIdx) = 0;
    end
end    

figure
imshow(outputImage,'InitialMagnification',67)
%imshow(imoverlay(outputImage,A),'InitialMagnification',67)


% Creates a Binary Image after the watershed transform
GrayScalePumpkinImage2 = rgb2gray(outputImage);
NewBinaryPumkinImage =  GrayScalePumpkinImage2 > 6; % not bigger than 12 
figure(10)
imshow(NewBinaryPumkinImage), title('Watershed transform of gradient magnitude (Lrgb)')


bws = bwmorph(BW1,'skel');
Over = imoverlay(NewBinaryPumkinImage,bws,'black');
figure(11)
imshow(Over,'InitialMagnification',67)


%% Counts the number of pumpkins in the image. 

GrayScalePumpkinImage = rgb2gray(A);
GrayScalePumpkinImage3 = rgb2gray(Over);
NewBinaryPumkinImage3 =  GrayScalePumpkinImage3 > 6; % not bigger than 12
sRegionProps1 = regionprops(NewBinaryPumkinImage3, GrayScalePumpkinImage, {'Centroid'});


figure(54)
imshow(NewBinaryPumkinImage)
title('Weighted (red) and Unweighted (blue) Centroids'); 
hold on
numObj = numel(sRegionProps1);

for k = 1 : numObj
   % plot(sRegionProps(k).WeightedCentroid(1), sRegionProps(k).WeightedCentroid(2), 'r*');
    plot(sRegionProps1(k).Centroid(1), sRegionProps1(k).Centroid(2), 'bo');
end
hold off



end
        
    end
    
end