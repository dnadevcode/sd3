function [curve, xF, yF, distance] = get_curve_parameters_spline( bw,mat )

    % skeletonize
    out = bwskel(bw==1,'MinBranchLength',20);

    % smoothen the skeleton
    [coor1 coor2]=find(out==1); % Finding all the co-ordinates for the corresponding component
    
    if length(coor1) > 3 % what is minimum number of points required
        yy1 = smooth(coor2,coor1,0.5,'loess');
        [curve, goodness, output] = fit(coor2,yy1,'smoothingspline');% ,'SmoothingParam',0.5


    % finally fit a spline//either x or y coord, depending which gives more
    % accurate!
%     [curve, goodness, output] = fit([coor2,yy1], ones(length(coor2),1),'lowess','Span',0.5);% ,'SmoothingParam',0.5

%         figure,imagesc(mat)
%         hold on
%         plot(curve)
    % find distances between points
    
        xd = 0.01;
        xvals = 1:xd:size(mat,2);
        yvals = curve(1:xd:size(mat,2));
        allD = zeros(1,length(xvals)-1);
        for i=1:length(xvals)-1
           allD(i) = sqrt((xvals(i+1)-xvals(i))^2+(yvals(i+1)-yvals(i))^2); % straight line distance
        end
        distAlong = cumsum(allD);
        % 
        pxPos = 0:floor(sum(allD));
        idxPx = zeros(1,length(pxPos));
        idxPx(1)=1;
        for idpxPos = 1:pxPos(end)
            idxPx(idpxPos+1) = find(distAlong >= idpxPos ,1,'first');
        end
%         find

        xF = xvals(idxPx);
        yF = yvals(idxPx);
        
%          figure,imagesc(mat)
%         hold on
%         plot(xF,yF)
        
        distance = zeros(1,length(xF)-1);
        for i=1:length(xF)-1
           distance(i) = sqrt((xF(i+1)-xF(i))^2+(yF(i+1)-yF(i))^2); % straight line distance
        end
        
    else
        curve = nan;
        xF = nan;
        yF = nan;
        distance = nan;
    end

        
        
%     % gets line parameters
%     nonnans = find(~isnan(mat));
%     [row, col] = ind2sub(size(bw), nonnans);
%     weight = mat(nonnans);
    
    
    
    % possible change: find skeletonized version of bw, that should go
    % through the center
    

    % use PCA to determine the flow of the molecule (this will determine  which
    % direction should we fit a line to.

    % RUN PCA: If coefficient very sharp, then need to switch 
%     [coeff,~,~,~,~,mu]= pca(vals);
%     % Calculate the slope of the line
%     %	(since the y axis is reversed in image analysis, a minus sign is introduced
%     %  in order for the number k to represent what is normally thought of as the slope)
%     k = -coeff(1,1)/coeff(2,1); 
% 
%     % intercept of first column (1 is top row)
%     b = mu(1) + (mu(2) - 1) * k;

%     mat2 = ones(size(bw));
%     mat2(nonnans) = mat(nonnans);
%     % should row/col be switched for this?
%     f=fit(row,col,'smoothingspline','Weights',weight);% make sure that the line is calculated in y direction (possible issue otherwise, need some tests..)
% %     f=fit(row,col,'smoothingspline');% make sure that the line is calculated in y direction (possible issue otherwise, need some tests..)
%     
%     % todo : more complicated structures 
%     
% %     k1 = f.p1;
% %     b1 = f.p2;
% %     c1 = f.p3;
% %     k = -1/k1;
% %     b = -b1/k1;
% %         c = -b1/k1;
% %         x = 1:max(row);
%         
%         xd = 0.1;
%         xvals = 1:xd:size(mat,1);
%         yvals = f(1:xd:size(mat,1));
%         allD = zeros(1,length(xvals)-1);
%         for i=1:length(xvals)-1
%            allD(i) = sqrt((xvals(i+1)-xvals(i))^2+(yvals(i+1)-yvals(i))^2); % straight line distance
%         end
%         distAlong = cumsum(allD);
%         % 
%         pxPos = 0:floor(sum(allD));
%         idxPx = zeros(1,length(pxPos));
%         idxPx(1)=1;
%         for idpxPos = 1:pxPos(end)
%             idxPx(idpxPos+1) = find(distAlong >= idpxPos ,1,'first');
%         end
% %         find
% 
%         xF = xvals(idxPx);
%         yF = yvals(idxPx);
%         
%         distance = zeros(1,length(xF)-1);
%         for i=1:length(xF)-1
%            distance(i) = sqrt((xF(i+1)-xF(i))^2+(yF(i+1)-yF(i))^2); % straight line distance
%         end
% 

%    figure
%         imagesc(mat')
%         hold on
%                 plot(f(x))
% % 
%   

%     plot1= 0;
%     if plot1
%         x = 1:max(row);
%         
%         figure
%         imagesc(mat')
%         hold on
%                 plot(f(x))
% % 
% %         plot(f.p1*x.^2+f.p2*x+f.p3)
%         hold on
% %         plot(row,col)
%         imagesc(row,col)
%         plot(-f.p1(
%         
%     end

end

