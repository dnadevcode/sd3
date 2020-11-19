function [img,ecc,aRat,length,width] = cont_draw(B)

rows = max(B(:,1))-min(B(:,1))+1;
cols = max(B(:,2))-min(B(:,2))+1;
img = zeros(rows,cols);

for i = 1:size(B,1)
  img(B(i,1),B(i,2)) = 1;
end

stats = regionprops(img,'Eccentricity','FilledArea','ConvexArea','MajorAxisLength','MinorAxisLength');
ecc = stats.Eccentricity;
aRat = stats.FilledArea/stats.ConvexArea;
length = stats.MajorAxisLength;
width = stats.MinorAxisLength;
