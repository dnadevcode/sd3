function [kymo] = get_kymo(mol,k,b,sPer)

   [n,m] = size(mol(:,:,1));
   sz = -k*(m-1)+b;	
	if k > 0
		st = min(n,b);
		en = max(1,sz);
	elseif k < 0
		st = max(1,b);
		en = min(n,sz);
	else
		st = b;
		en = b;
	end

	yl = en - st;
	if k  == 0
		xl = 0;
	else
	xl =  yl/k;
	end
	l = ceil(sqrt(xl^2 + yl^2))+1;

	% Only scan over rows where the molecule is present
    kymo = zeros(size(mol,3),l);
	angle = atan(k);
    for i=1:size(mol,3)
        molC = mol(:,:,i);
		lSpc = linspace(-sPer,sPer,2*sPer+1);
        for lIdx=1:l
		  	xVals = (b-st)/k + 1 + cos(angle) * (lIdx - 1)  - sin(angle) * lSpc;
		  	yVals = st - sin(angle) * (lIdx - 1) -  cos(angle) * lSpc;
            Vq = interp2(molC,xVals,yVals,'linear'); % Could change interpolation method
            kymo(i,lIdx) = nanmean(Vq);
        end
    end
end

