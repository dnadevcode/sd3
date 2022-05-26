function [kymo, boundaries] = get_kymo(mol,k,b,sPer)
    % get_kymo - extract kymo from a movie
    % based on line parameters k,b and perpendicular length sPer
    % the line equations is  y = kx+b, and the normal equation is
    % todo: more complicated features, using spline.
    
    [n,m] = size(mol(:,:,1)); % y and x dimensions

    sz = -k*(m-1)+b; % for the second-last x coordinate, do we go out of bounds?

    % start and end bounds on the molecule (for y) to detect dots so that
    % we do not go out of the field of view
    if k > 0 
      stY = min(n,b);
      enY = max(1,sz);
    elseif k < 0
      stY = max(1,b);
      enY = min(n,sz);
    else
      stY = b; %constant
      enY = b;
      stX = 1;
      enX = m;
    end
    
    stX = -(stY-b)/k;
    enX = -(enY-b)/k;

    
    % boundaries for the molecule within the cut-out window
    boundaries = [stY enY stX enX];


    % y,x px between start and stop
    yl = enY - stY;
    xl = enX - stX;
% if k  == 0
%   xl = 0;
% else
%   xl =  yl/k;
% end
lenPx = ceil(sqrt(xl^2 + yl^2))+1; % total number of pixels

% Only scan over rows where the molecule is present
kymo = zeros(size(mol,3),lenPx);
angle = atan(k);
for i=1:size(mol,3)
  molC = mol(:,:,i);
  lSpc = linspace(-sPer,sPer,2*sPer+1); % take perpendicularly
  for lIdx=1:lenPx % if at least one point outside, produces a nan
    xVals = stX + cos(angle) * (lIdx - 1)  - sin(angle) * lSpc; % stX + cos(angle) * (lIdx - 1) is current x coordinate
    yVals = stY - sin(angle) * (lIdx - 1) -  cos(angle) * lSpc; % same for Y
    Vq = interp2(molC,xVals,yVals,'linear'); % Could change interpolation method
    kymo(i,lIdx) = nanmean(Vq); % could be nansum
  end
end
end

