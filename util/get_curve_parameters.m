function [f, xF, yF] = get_curve_parameters( bw,mat )

    % gets line parameters
    nonnans = find(~isnan(mat));
    [row, col] = ind2sub(size(bw), nonnans);
    weight = mat(nonnans);

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
    % should row/col be switched for this?
    f=fit(row,col,'smoothingspline','Weights',weight);% make sure that the line is calculated in y direction (possible issue otherwise, need some tests..)
%     f=fit(row,col,'smoothingspline');% make sure that the line is calculated in y direction (possible issue otherwise, need some tests..)
    
    % todo : more complicated structures 
    
%     k1 = f.p1;
%     b1 = f.p2;
%     c1 = f.p3;
%     k = -1/k1;
%     b = -b1/k1;
%         c = -b1/k1;
%         x = 1:max(row);
        
        xd = 0.1;
        xvals = 1:xd:size(mat,1);
        yvals = f(1:xd:size(mat,1));
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

