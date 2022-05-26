function [kymos, lineParams, xyPars] = get_kymos_from_movies( molM, bwM, sPer, method)
% this function gets kymos from movies

if nargin < 4
    method = 1;
end
kymos = cell(1, length(molM));
lineParams = cell(1, length(molM));
xyPars = cell(1, length(molM));

%import PD.Core.Extraction.get_line_parameters;
%import PD.Core.Extraction.get_kymo;
%%
for i=1:length(molM)
  % get line parameters
  if method ==1
      [k,b] = get_line_parameters(bwM{i});
      kymos{i} = get_kymo(molM{i}, k , b, sPer);
      lineParams{i} = [k b];
  else
      [f,xF,yF] = get_curve_parameters(bwM{i},molM{i});
      tempKm = [];
        for lIdx=1:length(xF)
            Vq = interp2(molM{i}',xF(lIdx),yF(lIdx),'linear'); % Could change interpolation method
           	tempKm(1,lIdx) = Vq; % could be nansum
        end
        % also find the longest part without nan's, and bitmask the rest out 
%         bwconncomp(isnan( kymos{i)
        connCompoments = bwconncomp(~isnan(tempKm));
        conCompLengths = cellfun(@(x) length(x),connCompoments.PixelIdxList);
        [maxL, maxId] = max(conCompLengths);
        kymos{i} = nan(1,length(tempKm));
        kymos{i}(connCompoments.PixelIdxList{maxId}) = tempKm(connCompoments.PixelIdxList{maxId});
        lineParams{i} = f;
        xyPars{i} ={xF, yF};
  end
  %    import PD.Core.Extraction.save_kymos;
  %    save_kymos( kymos{i}, fold, i-1  )
end

    %   % plot
    plotExample = 0;
    if plotExample ==1
        %%
        idx=1
%            [k,b] = get_line_parameters(bwM{i});
%       kymos{i} = get_kymo(molM{i}, k , b, sPer);
        figure,
        tiledlayout(1,3);nexttile
        imagesc(molM{idx})
        nexttile
%         imagesc(dotM{idx})
        imagesc(molM{idx})

        hold on
        plot(-lineParams{idx}(1)*(1:size(molM{idx},2))+lineParams{idx}(2),'redx')
          nexttile
%         imagesc(dotM{idx})
        imagesc(molM{idx}')
            hold on
            
            plot(xF,yF,'redx')
%             
        colormap(gray)
        %%
%         xd = 0.1;
%         xvals = 1:xd:size(molM{idx},1);
%         yvals = f(1:xd:size(molM{idx},1));
%         allD = zeros(1,length(xvals)-1);
%         for i=1:length(xvals)-1
%            allD(i) = sqrt((xvals(i+1)-xvals(i))^2+(yvals(i+1)-yvals(i))^2); % straight line distance
%         end
%         distAlong = cumsum(allD);
%         % 
%         pxPos = 0:round(sum(allD));
%         idxPx = zeros(1,length(pxPos));
%         idxPx(1)=1;
%         for idpxPos = 1:pxPos(end)
%             idxPx(idpxPos+1) = find(distAlong >= idpxPos ,1,'first');
%         end
% %         find
% 
%         xF = xvals(idxPx);
%         yF = yvals(idxPx);
%         figure,plot(xF,yF)
        
        % then at each point, we can draw a diagonal (line) between two
        % points, and perpendicular to that will give us the perpendicular
        % points for the kymo

%         norm(
    end
end

