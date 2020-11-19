function [ kymo ] = get_kymo( mol, bw, k , b, sPer )

[n,m] = size(mol(:,:,1));
sz = ceil(-k*(m-1)+b);	% Could add extra space here to show that molecule ends.
st = max(1,min(sz,floor(b))); 	% Row where the molecule enters the frame
en = min(n,max(sz,floor(b)));	   % Row where the molecule exits the frame

% Only scan over rows where the molecule is present

kymo = zeros(size(mol,3),en-st+1);
angle = atan(k);
for i=1:size(mol,3)
  molC = mol(:,:,i);
  % what happens when the shape is different - maybe rotate to always
  for l=st:en
    j = (b-l)/k + 1;
    % only in the allowed columns
    newXvals = j + sin(angle) * linspace(-sPer,sPer,2*sPer+1);
    newYvals = l + cos(angle) * linspace(-sPer,sPer,2*sPer+1);
    Vq = interp2(molC,newXvals,newYvals,'nearest'); % Could change interpolation method
    kymo(i,l-st+1) = nanmean(Vq);
  end
  %for l=st:en
  %		j = (b-l)/k + 1;
  %      % only in the allowed columns
  %      newXvals = j + sin(angle) * linspace(-sPer,sPer,2*sPer+1);
  %      newYvals = l + cos(angle) * linspace(-sPer,sPer,2*sPer+1);
  %      Vq = interp2(molC,newXvals,newYvals,'nearest'); % Could change interpolation method
  %      kymo(i,l-st+1) = nanmean(Vq);
  %  end
end
end

