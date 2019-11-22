function [ xcorrs ] = compute_pcc( shortVec,longVec,w1,w2,th )

if nargin < 5
	th = 2;
end

    % compute_pcc
    % Computes Pearson correlation coefficient using fft's
    %     Args:
    %         shortVec: short linear barcode
    %         longVec: long circular barcode
    %         w1: bitmask of first barcode
    %         w2: bitmask of second barcode
    % 
    %     Returns:
    %         xcorrs: PCC values
    %
    %     Example:
    %       
    % 
    %     % compute X*Y
    % shortLength = 100;
    % longLength = 1000;
    % 
    % shortVec = rand(1,shortLength);
    % longVec = rand(1,longLength);
    %     % bitmasks first bar
    %     w1 = ones(1,shortLength);
    %     w1(1)=0;
    %     w1(end)=0;
    % 
    %     % bitmasks second bar
    %     w2 = ones(1,longLength);
    %     w2(1)=0;
    %     w2(end)=0;
    %   [ xcorrs ] = compute_pcc( shortVec,longVec,w1,w2 )


	%%%%%%%%%%%%%%%%%%%%%%%%%%	
   % Ensure shortVec is the shortest
%	  if length(shortVec)>length(longVec)
%		  lv_temp = shortVec;
%   	  w2_temp = w1;
%		  shortVec = longVec;
%	     w1 = w2;
%		  w2 = w2_temp;
%        longVec = lv_temp;
%	 end
    longLength = length(longVec);
    
    % 0 where bitmask is 0
    shortVec(~w1) = 0;
    longVec(~w2) = 0;

    w2fft = fft(w2);
    % Compute the number of nonzero elements in each comparison
    conVec = conj(fft(w1,longLength));
    numForward =  (ifft(w2fft.*conVec));
	%numForward2 = ifft(w2fft.*conj(fft(fliplr(w1),longLength)));

    % ybar
    ybar =  (ifft(fft(longVec).*conVec))./numForward;
    ySq =  (ifft(fft(longVec.^2).*conVec));

    % main part. Compute xiyi
    conVec = conj(fft(shortVec,longLength));
    xiyi = (ifft(fft(longVec).*conVec));
    xbar =  (ifft(w2fft.*conVec))./numForward;

    % now, compute mean of the first vector sq
    conVec = conj(fft(shortVec.^2,longLength));
    xSq =  (ifft(w2fft.*conVec));

    denominator = sqrt(xSq-numForward.*xbar.^2).*sqrt(ySq-numForward.*ybar.^2);
    numerator = xiyi-numForward.*xbar.*ybar;
    % now, compute mean of the second vector sq
    ccF = numerator./denominator;
    ccF(numForward < th) = NaN;   % Added by Jens to avoid weird numbers

    % First part, just compute XY
    shortVecFlip = fliplr(shortVec);
    w1Flip = fliplr(w1);     

    % Compute the number of nonzero elements in each comparison
    conVec = conj(fft(w1Flip,longLength));
    numForward =  (ifft(w2fft.*conVec));
    % ybar
    ybar =  (ifft(fft(longVec).*conVec))./numForward;
    ySq =  (ifft(fft(longVec.^2).*conVec));

    % main part. Compute xiyi
    conVec = conj(fft(shortVecFlip,longLength));
    xiyi = (ifft(fft(longVec).*conVec));
    xbar =  (ifft(w2fft.*conVec))./numForward;

    % now, compute mean of the first vector sq
    conVec = conj(fft(shortVecFlip.^2,longLength));
    xSq =  (ifft(w2fft.*conVec));

    denominator = sqrt(xSq-numForward.*xbar.^2).*sqrt(ySq-numForward.*ybar.^2);
    numerator = xiyi-numForward.*xbar.*ybar;
    % now, compute mean of the second vector sq
    ccB = numerator./denominator;
   ccB(numForward < th) = NaN;   % Added by Jens to avoid weird numbers

    %xx2 = circshift(ccB,[0,length(shortVec)]);
    %xx1 = ccF;
    xcorrs = [ccF;ccB];
end

