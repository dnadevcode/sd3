function [ avg ] = averageImages( registeredIm,n )
if nargin < 2
  n = length(registeredIm);
end
t = zeros(size(registeredIm{1},1),size(registeredIm{1},2),min(length(registeredIm),n));

for i=1:min(length(registeredIm),n)
  t(:,:,i) = registeredIm{i};
end


avg = nanmean(t,3);



end

