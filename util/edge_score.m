function score = edge_score(B,im,Gdir,dist,stat)

% tic
% bound = B;
h = zeros(1,size(B,1));
for point = 1:size(B,1)
    dir = Gdir(B(point,1),B(point,2)); % test to see if dirs were switched
    dx = cosd(dir); %magnus' way
    dy = -sind(dir); %magnus' way
    xser = round((-dist:1:dist)*dx) + B(point,2);
    yser = round((-dist:1:dist)*dy) + B(point,1);
    prof = zeros(1,2*dist+1);
    if min(min(yser),min(xser))>0 && max(xser) < size(im,2) && max(yser) < size(im,1)
        for j = 1:2*dist+1
            prof(j) = im(yser(j),xser(j));
        end
    end
    %h(point) = abs(sum(prof(1:dist))-sum(prof(dist+2:end)));
    h(point) = -sum(prof(1:dist))+sum(prof(dist+2:end));
end
% toc
score = stat(h);

end
