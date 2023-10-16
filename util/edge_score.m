function score = edge_score(B,im,Gdir,dist,stat)

% tic
szIM = size(im,[1 2]);
distVals = -dist:1:dist;
% bound = B;
h = zeros(1,size(B,1));
for point = 1:size(B,1)
    dir = Gdir(B(point,1),B(point,2)); % test to see if dirs were switched
    dx = cosd(dir); 
    dy = -sind(dir); 
    xser = round((distVals)*dx) + B(point,2);
    yser = round((distVals)*dy) + B(point,1);
    prof = zeros(1,2*dist+1);
    if min(min(yser),min(xser))>0 && max(xser) < szIM(2) && max(yser) < szIM(1)
        prof = arrayfun(@(x,y) im(y,x),xser,yser);
    end
    %h(point) = abs(sum(prof(1:dist))-sum(prof(dist+2:end)));
    h(point) = -sum(prof(1:dist))+sum(prof(dist+2:end));
end
% toc
score = stat(h);

end
