function [fluorophorePosPx, fluorophorePosPxYOYO] = simulate_stretched_contour(molLen,stX, numPts, w)
% simulate stretched contour of DNA on glass. Allow to be slightly curve
% 
%       Args:
%   molLen - molecule length
%   stX - starting position
%   numPts - 
%   w - how much deviation from line

% Very simple random walk to get things going
% 

posPts = round(rand(1,numPts)*(molLen-50))+1; % for now will be integers;
x(1)=stX;
if nargin < 4
    w = randn(); % random walk parameter
end

for i=1:molLen-1
    dx(i)=w*randn();
    x(i+1)=x(i)+dx(i);
end

% DNA contour
fluorophorePosPxYOYO = [50+[1:molLen]' x'];
fluorophorePosPx =  [50+posPts' x(posPts)'];

end

