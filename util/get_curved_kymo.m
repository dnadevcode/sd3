function [kymo, boundaries] = get_curved_kymo(mol, xF, yF, sPer)
    % get_kymo - extract kymo from a movie
    % based on line parameters k,b and perpendicular length sPer
    % the line equations is  y = kx+b, and the normal equation is
    % todo: more complicated features, using spline.
    

    tempKm = zeros(1,length(xF));
    for lIdx=1:length(xF)
        Vq = interp2(mol',xF(lIdx),yF(lIdx),'linear'); % Could change interpolation method
        tempKm(1,lIdx) = Vq; % could be nansum
    end
    % also find the longest part without nan's, and bitmask the rest out 
%         bwconncomp(isnan( kymos{i)
    connCompoments = bwconncomp(~isnan(tempKm));
    conCompLengths = cellfun(@(x) length(x),connCompoments.PixelIdxList);
    [maxL, maxId] = max(conCompLengths);
    kymo = nan(1,length(tempKm));
    kymo(connCompoments.PixelIdxList{maxId}) = tempKm(connCompoments.PixelIdxList{maxId});
    
    boundaries = [size(mol,1) 1  0 size(mol,2)-1];

   
end

