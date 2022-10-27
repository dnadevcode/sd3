function [ nw,centers ] = get_remove_circles( imDenoised, imAverage,rMin, rMax )
% removes circles
if nargin < 4
  rMin = 12;
  rMax = 23;
end

% expect the circles to have higher intensity than all else
lim = multithresh(imDenoised,2);

circles = imDenoised > lim(2);
se = strel('disk',5);
dilatedBW = imdilate(circles,se);
%     dilatedBW = circles;
[centers] = imfindcircles(dilatedBW,[rMin rMax]);
 
% figure,imagesc(imDenoised)
% hold on
% for i=1:size(centers,1)
%     plot(centers(i,1),centers(i,2),'redx')
% end

%    compare2(circles,imAverage);

% circC = bwpropfilt(circles,'Area',[0 100]);
% compare2(circC,avg);

% dilate circles
%    compare2(circles,dilatedBW);

%
% circlesT = bwmorph(circles,'majority');
%
% circlesT = bwmorph(circlesT,'thicken',40);
% circlesT = bwmorph(circlesT,'majority');
% circlesT = bwmorph(circlesT,'fill');

% without circles
%withoutCircles = (not(dilatedBW)).*(imAverage);
%figure,imshow(withoutCircles,[])

%
ww = (not(circles)).*(imAverage);
%figure,imshow(ww,[])

nw = ww;
nw(nw==0) = nan;

end

