%%
A = imread('DJI_0486.JPG');
DetectedPumpkins = imread('DJI_0486SinglepumpkinsredDoublePumpkinsBlue.PNG');

hsv = rgb2hsv(A);
figure(1)
image(A)

hu = hsv(:,:,1);         % Extract the HSV values  % hu between 0 and 55 will show between red and yellow
s = hsv(:,:,2);            
v = hsv(:,:,3); 
%hu*360;
hu360 = hu.*360;

figure(2)
image(hsv)
%FilterNewPumpkins
hsv2 = cat(3,hu360,s,v);

hu2 = hsv2(:,:,1);         % Extract the HSV values  % hu between 0 and 55 will show between red and yellow
s2 = hsv2(:,:,2);            
v2 = hsv2(:,:,3);


figure(3)
image(hsv2)


rgb = hsv2rgb(hsv);

figure(4)
image(rgb)
% 2390   1368
% 2386   1368
R = A(:,:,1);
G = A(:,:,2);
B = A(:,:,3);
%centerx = [0,0];
WindowSize = 20;


a = A(:,:,1);

%%
progresscomment('Starting filtering pumpkins. Mahanalobis');
Channels = MarkedPumpkins(A);
scatter3(Channels(1,:),Channels(2,:),Channels(3,:))
xlabel('Red') % x-axis label
ylabel('Green') % y-axis label
zlabel('Blue') % y-axis label

ChannelsRotated = rot90(Channels);
%coeff = pca(Channels);
[coeff,score,latent] = pca(ChannelsRotated);
Xcentered = score*coeff';

arrayTest=[77,106,91];
d1 = mahal(arrayTest,ChannelsRotated)
BinaryPumkinImage2 = FindPumpkinsInImage(A,ChannelsRotated);


%%
figure(4)
image(BinaryPumkinImage2)

%MahanalobisDistToPumpkinDistribution(arrayTest,PumpkinDistribution)
%%  2 
BinaryImage = FindPumpkinsInImageVersion2(A,ChannelsRotated);

%figure(2)
%imshow(BinaryPumkinImage2)
BinaryImageInverted = imcomplement(BinaryImage); %% have no idea why it now has to be inverted.


figure(7)
p1=subplot(2, 1, 1);
imshow(BinaryImageInverted);
%h = viscircles(MatrixPumpkins,RadiiRot);

p2=subplot(2, 1, 2);
imshow(A)
%h = viscircles(MatrixPumpkins,RadiiRot);

linkaxes([p1, p2]);

progresscomment('Finished filtering pumpkins.');

%%



progresscomment('start Gaussian pumpkins.');

GaussianBluredImage = GaussianFunction(BinaryImageInverted);
%%
GrayScalePumpkinImage = rgb2gray(A);
GaussianBinaryImage =  GaussianBluredImage > 2;


sRegionProps = regionprops(BinaryImageInverted, GrayScalePumpkinImage, {'Centroid','WeightedCentroid','Area','EquivDiameter','MajorAxisLength','MinorAxisLength','Orientation','Perimeter'});

figure(54)
imshow(GrayScalePumpkinImage)
title('Weighted (red) and Unweighted (blue) Centroids'); 
hold on
numObj = numel(sRegionProps);
for k = 1 : numObj
    plot(sRegionProps(k).WeightedCentroid(1), sRegionProps(k).WeightedCentroid(2), 'r*');
    plot(sRegionProps(k).Centroid(1), sRegionProps(k).Centroid(2), 'bo');
end
hold off



figure(55)
imshow(GaussianBinaryImage)
sBin = regionprops(GaussianBinaryImage,'centroid');
centroids = cat(1, sBin.Centroid);
rmin = 1;
rmax = 50;
[centers,radii] = imfindcircles(GaussianBinaryImage,[rmin rmax]);
viscircles(centers,radii);
figure(56)
imshow(GaussianBinaryImage)
hold on
plot(centroids(:,1),centroids(:,2), 'b*')
hold off

progresscomment('Finished Gaussian pumpkins.');
%%
GaussianBinaryImageDistance = bwdist(~GaussianBinaryImage);

figure(59)
imshow(GaussianBinaryImageDistance,[],'InitialMagnification','fit')
title('Distance transform of ~bw');

%%
GaussianBinaryImageDistance = -GaussianBinaryImageDistance;
GaussianBinaryImageDistance(~GaussianBinaryImage) = Inf;

WatershedImage = watershed(GaussianBinaryImageDistance,8);

WatershedImage(~GaussianBinaryImage) = 0;
rgb = label2rgb(WatershedImage,'jet',[.5 .5 .5]);

figure(60)
imshow(rgb)

%%

 AreaAndPerimeterMarkedPumpkins = FindAverageAreaOfMarkedPumpkins(DetectedPumpkins);
allAreaMarked = [AreaAndPerimeterMarkedPumpkins.Area];
[rowPumpkinsMarked,colPumpkinsMarked] = size(allAreaMarked)
SumPumpkinsMarked = sum(allAreaMarked)
AveragePumpkinSize = SumPumpkinsMarked / colPumpkinsMarked
%%
 ColorBinaryCombinedImage = CombineBinaryAndColorImage(A,BinaryPumkinImage2);
 
 
 %%
 
 %value = getfield(AreaAndPerimeter, 'Area')
 AreaAndPerimeter(2).Area
 %%
allArea = [sRegionProps.Area];
Hej = LoopThroughPumpkinsListAndFindClusterPumpkins(500,AveragePumpkinSize,allArea);
 
 %%
 SumPumpkins = sum(Hej)
 %%
 

figure(8)
p1=subplot(2, 1, 1);
imshow(ColorBinaryCombinedImage);
%h = viscircles(MatrixPumpkins,RadiiRot);

p2=subplot(2, 1, 2);
imshow(DetectedPumpkins)
%h = viscircles(MatrixPumpkins,RadiiRot);

linkaxes([p1, p2]);

%%
progresscomment('Starting filtering pumpkins.');
%figure(5);
%image(B);
counter = 1;
%array(499,806);
%x1 = ones(150,150,3);
[row,col] = size(a);
%[row1,col1] = size(DetectedPumpkins);
LastRedValue = 255;
LastBlueValue = 0;
for ii = 1:row
    for jj = 1:col
        RGBRed = DetectedPumpkins(ii,jj,1);
        RGBGreen = DetectedPumpkins(ii,jj,2);
        RGBBlue = DetectedPumpkins(ii,jj,3);
        
        if RGBRed > 250 && RGBGreen == 0 && RGBBlue == 0
            RedInterval = A(ii,jj,1);
            BlueInterval = A(ii,jj,3);
            
            if LastRedValue > RedInterval 
               LastRedValue = RedInterval;
            end
            
            if LastBlueValue < BlueInterval 
               LastBlueValue = BlueInterval;
            end
        end
    end
end
LastBlueValue
LastRedValue

x1bin = PumpkinFilter(hu360,s,v);


figure(5)
imshow(x1bin)
progresscomment('Finished filtering pumpkins.');

%%

bw2 = bwareaopen(x1bin,50);

%figure(6)
%imshow(bw2);

figure(7)
p1=subplot(2, 1, 1);
imshow(bw2);
%h = viscircles(MatrixPumpkins,RadiiRot);

p2=subplot(2, 1, 2);
imshow(A)
%h = viscircles(MatrixPumpkins,RadiiRot);

linkaxes([p1, p2]);
%%
for n = 1:row
    for ni = 1:col
        Red = A(n,ni,1);
        Blue = A(n,ni,3);
        Green = A(n,ni,2);
        if 1
       if Red > LastRedValue && Blue < LastBlueValue %&& Green > 81 && Green < 222   
           %   186 Red   og 81 Blue    % 104 og 58
           % Green mellem 81 og 117
           %n;
           %ni;
           x1(n,ni,:) = A(n,ni,:);   % melder fejl 82 + 16 = 98 > 97 hmmmm får aldrig fuld størrelse
           array(n,ni) = n;
       else
           x1(n,ni,:) = Zeroes2 ;
           array(n,ni) = 0;
       end
        end
    end
end
%%
progresscomment('Pumpkins has been filtered.');
for n = 1:row    %  3078   % building a grayscale image, where white are pumpkins, and black are background
   for ni = 1:col    % 5472
       Red = x1(n,ni,1);
       Blue = x1(n,ni,3);
       Green = x1(n,ni,2);
       n;
       ni;
       if 1
       if Red > 1 && Blue > 1 %&& Green > 81 && Green < 222   %   186 Red   og 81 Blue    % 104 og 58
           % Green mellem 81 og 117
           n;
           ni;
          % x1(n,ni,:) = A(n,ni,:);
           %array(n,ni) = n; 
           pump = 0;
           if  n < row-WindowSize && ni < col-WindowSize
               pump = FindPumpkin(n,ni,x1);
               pump;
               
           end
           %pump
          % n
           %ni
           if 1
           if pump > 0 
               pump;
               %if LastPixelPump == 0 % Was last pixel a pumpkin?
               % a pumpkin has been detected
                    if counter < 2 
                       % secures that that not every pixel in a pumpkin is
                       % detected with a center point. 
                       counter = counter + 1
                       centerx(end + 1) = n;
                       centery(end + 1) = ni;
                       weight(end + 1 ) = pump;
                       Radii(end+1) = WindowSize/2;
                       col1 =3;
                    end
               
                   %[row1,col1] = size(centerx);
                    
                    for i = 1:col1
                    xiCent = centerx(i);
                    yiCent = centery(i);
                    centerx;
                    
                    weightiCent = weight(i);
                    
                   % if ni == 642 && n < 80
                   %        n
                   %        ni
                           
                   %     end
                    
   if  (n-WindowSize < xiCent) & (xiCent < n + WindowSize ) & (ni-WindowSize<yiCent) & (yiCent<ni + WindowSize )
                        ItsANewPumpkin = 0 ; 
                       
                        if pump >= weightiCent 
                            centerx(i) = n;
                            centery(i) = ni;
                            weight(i) = pump;
                            
                           if ni == 642 && n < 80
                           n
                           xiCent
                           ni
                           yiCent
                           weightiCent
                           pump
                           end
                           % weightiCent
                           % pump
                           % n
                           % ni
                           % xiCent
                           % yiCent
                        end    
                        break
                        
                        
                    
                    else
                        % It is a new pumpkin. 
                        ItsANewPumpkin = 1;  
                        %center(n,ni);
                    
                        LastPixelPump = 0;
                    end 
                    end
                    if ItsANewPumpkin == 1 
                   %har udkommenteret if her (skal tilbage, ved brug af programmet.)
                       counter = counter + 1;
                       centerx(end + 1) = n;
                       centery(end + 1) = ni;
                       weight(end + 1 ) = pump;
                       Radii(end+1) = WindowSize/2;
                       col1 = col1 +1;
                    end
                    
           end 
           end
   end
           %else 
            %   LastPixelPump = 0 ;
           %end
           
       end
       
    end
   
end
progresscomment('All pumpkins are found.');
%figure(6)
%imshow(x1);
Cen = vertcat(centerx,centery);
centerxRot = rot90(centerx);
centeryRot = rot90(centery);
weightRot = rot90(weight);
RadiiRot = rot90(Radii);
center = [centeryRot, centerxRot]; %% this shall be done in order to get the right orientation of the image. 
centerWeight = [centeryRot, centerxRot, weightRot];

[row,col] = size(centerWeight);
progresscomment('Start filtering double counted pumpkins.');
MatrixPumpkins = [row,col];
for n = 1:row
    xnCent = centerWeight(n,1);
    ynCent = centerWeight(n,2);
   for ni = n:row 
    
    xniCent = centerWeight(ni,1);
    yniCent = centerWeight(ni,2);
    xnCent;
        xniCent;
        ynCent;
        yniCent;
if  (xnCent-WindowSize < xniCent) & (xniCent < xnCent + WindowSize ) & (ynCent-WindowSize<yniCent) & (yniCent<ynCent + WindowSize )
        
        ItsANewPumpkin = 0;
        MatrixPumpkins(n,1) = 0;
        MatrixPumpkins(n,2) = 0;
    else
        ItsANewPumpkin = 1;
        FilterNewPumpkins = FilterNewPumpkins + 1;
        MatrixPumpkins(n,1) = xniCent + WindowSize/2;
        MatrixPumpkins(n,2) = yniCent + WindowSize/2;
        break
    end
    
   end
    
    
end

progresscomment('Finish.');
figure(7)
imshow(x1);
h = viscircles(center,RadiiRot);
size(Zeroes2);
figure(8);
imshow(array);
h = viscircles(MatrixPumpkins,RadiiRot);

figure(9)
imshow(A)
h = viscircles(MatrixPumpkins,RadiiRot);


%%

figure(18)
p1=subplot(2, 1, 1);
imshow(array);
h = viscircles(MatrixPumpkins,RadiiRot);

p2=subplot(2, 1, 2);
imshow(A)
h = viscircles(MatrixPumpkins,RadiiRot);

linkaxes([p1, p2]);
%%
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
           GaussianImage(n,ni) = GaussianKernel(Mask);
           
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

function f = MarkedPumpkins(PumpkinImage)
%PumpkinImage = imread('DJI_0084.JPG');
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
%CombinedColorChannels = cat(3,RedChannel,GreenChannel,BlueChannel);
CombinedColorChannels = [RedChannel; GreenChannel; BlueChannel];
f = CombinedColorChannels;

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
             DistanceToPumpkinDistribution = MahanalobisDistToPumpkinDistribution(RGBPoint,PumpkinDistribution);
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
    MahalDistanceMatrice = mahal(Matrice2Double,PumpkinFloating);
    BinaryPumkinImage =  MahalDistanceMatrice < 12; % not bigger than 12  
    
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

function pumpk = DefPump(p,o)  % finds out, if the pumpkin is counted earlier. 
WindowSize = 20;
newCenter = [0,0];
[rowCent,colCent] = size(center)
[row,col] = size(a);
for n = 1:row    %  3078
   for ni = 1:col    % 5472
       Red = A(n,ni,1);
       Blue = A(n,ni,3);
       Green = A(n,ni,2);
       if Red > 141 && Blue < 96 %&& Green > 81 && Green < 222   %   186 Red   og 81 Blue    % 104 og 58
           % Wellow mellem 81 og 117
           x1(n,ni,:) = A(n,ni,:);
           array(n,ni) = n; 
           pump = 0;
           if n > WindowSize && n < row-WindowSize && ni > WindowSize && ni < col-WindowSize
               pump = FindPumpkin(n,ni,x1);
               
           end
          
           if pump > 0 
               if LastPixelPump > 0 % Was last pixel a pumpkin?
                    if LastPixelPump > pump
                    counter = counter + 1;
                    %center(n,ni);
                    centerx(end + 1) = n;
                    centery(end + 1) = ni;
                    Radii(end + 1) = 30;
                    LastPixelPump = pump;
                    end
               
               
               end 
           else 
               LastPixelPump = 0 ;
           end
           
       end
       
    end
   
end
end
