function [kymos,lineParams,xyPars, distance] = fit_spline_method(molM, bwM, sPer)
   
    %   Args:
    %       molM - molecule mini movie
    %       bwM - bitmask
    %       sPer - how many pixels perpendicularly to take
    %
    %   Returns:
    %       kymos - kymographs
    %       lineParams
    %       xyPars
    %       distance - distance between neighbouring pixels

    kymos = cell(1,length(molM));
    lineParams = cell(1,length(molM));
    xyPars = cell(1,length(molM));
    distance = cell(1,length(molM));

    for i=1:length(molM)
    
        [f,xF,yF, distance{i}] = get_curve_parameters_spline(bwM{i},molM{i}(:,:,1));
        tempKm =  zeros(size(molM{i},3),length(xF));
        for lIdx=1:length(xF)
            for tIdx =1:size(molM{i},3)
                Vq = interp2(molM{i}(:,:,tIdx)',yF(lIdx),xF(lIdx),'linear'); % Could change interpolation method
                tempKm(tIdx,lIdx) = Vq; % could be nansum
            end
        end
        % also find the longest part without nan's, and bitmask the rest out 
        connCompoments = bwconncomp(~isnan(tempKm));
        conCompLengths = cellfun(@(x) length(x),connCompoments.PixelIdxList);
        [maxL, maxId] = max(conCompLengths);
        kymos{i} = nan(size(tempKm));
        try
            kymos{i}(connCompoments.PixelIdxList{maxId}) = tempKm(connCompoments.PixelIdxList{maxId});
        catch
        end
        lineParams{i} = f;
        xyPars{i} ={yF,xF};
    end

end

