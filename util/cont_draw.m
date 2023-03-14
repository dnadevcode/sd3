function [img, stats] = cont_draw(B)
    % Draw contour and calculate statistics for potential molecule B
    %
    %       Args: 
    %           B - molecule edge mask
    %
    %       Returns:
    %           img - output image
    %           stats - statistics of image
    %
    minB = min(B);
    B=B-minB+1;
    maxB = max(B);

%     CC =[];
%     CC.Connectivity = 8;
%     CC.ImageSize = maxB;
%     CC.NumObjects = 1;
%     CC.PixelIdxList{1} = sub2ind(maxB,B(1:end-1,1),B(1:end-1,2));
%     stats = regionprops(CC,'Eccentricity','FilledArea','ConvexArea','MajorAxisLength','MinorAxisLength');

    
    img = zeros(maxB(1),maxB(2));
    
    for i = 1:size(B,1)
      img(B(i,1),B(i,2)) = 1;
    end
    
    % CC = bwconncomp(img);
    stats = regionprops(img,'Eccentricity','FilledArea','ConvexArea','MajorAxisLength','MinorAxisLength');


end
