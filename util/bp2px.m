function [thyCurve_pxRes] = bp2px(thyCurve_bpRes, meanBpExt_pixels)
    % Convert from bp resolution to pixel resolution
    % with moving average window.
    if nargin < 2
        pxSize = 541.28;
    else
        pxSize = 1/meanBpExt_pixels;
    end
    thyCurve_pxRes(floor(length(thyCurve_bpRes)/pxSize)) = 0;
%	 if length(thyCurve_bpRes)  < round(pxSize)*2
%		 fprintf('thyCurve length = %i.\n',length(thyCurve_bpRes))
%		 fprintf('pxSize  = %.2f\n',pxSize)
%		 thyCurve_pxRes = mean(thyCurve_bpRes);
%	 else
	    xtraseq = cat(find(size(thyCurve_bpRes) - 1), ...
                  thyCurve_bpRes(end-round(pxSize):end), ...
                  thyCurve_bpRes, ...
                  thyCurve_bpRes(1:round(2*pxSize)));
   	 for i = 1:floor(length(thyCurve_bpRes)/pxSize)
      	  thyCurve_pxRes(i) = mean(xtraseq(floor((i*pxSize)-pxSize+1):floor((i*pxSize)+pxSize)));
	    end
%	 end
end
