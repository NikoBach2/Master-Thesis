A = imread('DJI_0084.JPG');

B = imread('DJI_0084RedBackground.JPG');
counter = 0;
bool = 0
%%
%[centers, radii, metric] = imfindcircles(A,[40,80]);figure(7)
%imshow(A);
%h = viscircles(centers,radii,'EdgeColor','b');
%fact(A,5);
%B= A(:,:,3); % your blue band!
%figure(2)
%imshow(B);
a = A(:,:,1);
[row,col] = size(a);
C = zeros;
Reds = zeros;
Greens = zeros;
Blue = zeros;
%countTable;
for j = 1:row
   for i = 1:col   
       Red = B(j,i,1);
       Green = B(j,i,2);
             if 1 | Red == 255 & Green < 5
                 [tabel] = size(Reds);
                 counter = counter + 1;
                 %for k = 1:tabel
                  %   if A(j,i,1) == Reds(k) && A(j,i,2)==Green(k) && A(j,i,3)==Blue
                   %      bool = 1; % den kommer aldrig herind
                    %     countTable(k)=A(j,i,:);
                  %   end
                     
                 %end
                 if bool == 0
                 %C(end + 1) = A(j,i,:);
                 Reds(end +1) = A(j,i,1);
                 Greens(end +1) = A(j,i,2);
                 Blue(end +1) = A(j,i,3);
                 end
             end
   end
end

Reds
%surf(Reds,Green,Blue)
scatter3(Reds,Greens,Blue,1,[Reds',Greens',Blue']/255)
xlabel('Red') % x-axis label
ylabel('Green') % y-axis label
zlabel('Blue') % y-axis label
%dagbog den. 15 november:
% bliver nød til at lave et array, fo alle tre farver. 
