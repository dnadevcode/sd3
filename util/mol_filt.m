function [acc,stats,img] = mol_filt(B, score, lowLim, highLim, elim, ratlim, lengthLims, widthLims)
    %   Args:
    %
    %   Returns:
    %

    % estimate length upper bound (length in x + length in y)
    lenX =  abs(max(B(:,2))-min(B(:,2)));
    lenY =  abs(max(B(:,1))-min(B(:,1)));
    l = lenX+lenY;

    % estimate width (initial upper bound 
    w = min(lenX,lenY);

    % length and width limits
    lOk = (l > lengthLims(1) && l < lengthLims(2));
    % this should be min of (max repeating elements in X, Y)

    %wOk = (w < widthLims(2));

 
    if score > lowLim && score < highLim && lOk %&& wOk
        
        [img,stats] = cont_draw(B);
        testofboundary = (stats.Eccentricity > elim && stats.FilledArea/stats.ConvexArea > ratlim && stats.MinorAxisLength < widthLims(2)) ;
        if testofboundary
            acc = true;
        else
            acc = false;
        end
    else
        acc = false;
        img =[];
    end
    stats.l = l;
    stats.w = w;
end