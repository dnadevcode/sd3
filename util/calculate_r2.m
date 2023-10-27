function [r2] = calculate_r2(bwM,lineparams)
    %   bwM - mask
    %   xy - xy coordinates of a line

    r2 = zeros(1,length(bwM));
    for i=1:length(bwM)
        px = find(bwM{i}(:)==1);
        [ypos,xpos] = ind2sub(size(bwM{i}),px);
        ymean = mean(ypos);

%         lineparams
        SSres = sum((ypos- (-lineparams{i}(1)*xpos+lineparams{i}(2))).^2);
        SStot = sum((ypos- ymean).^2);

        r2(i) = 1- SSres/SStot;


    end
%     r2 = 0;
end

