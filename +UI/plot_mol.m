function [] = plot_mol(output,idxMov,nameMol)


hold on



extractionMethod = 2;
if extractionMethod == 1
    angle = atan(output{idxMov}.lineParams{nameMol}(1));
end

     for j = 1:numel(output{idxMov}.dots{nameMol}.locations)
            if extractionMethod == 1
%                 dy = -sin(angle)*(barcodes.dots{molIdx}.locations(j)-1+barcodes.dots{molIdx}.leftOffset);
%                 dx = cos(angle)*(barcodes.dots{molIdx}.locations(j)-1+barcodes.dots{molIdx}.leftOffset);
%                 y = movies.pos{curIdx}(1)+vOff+dy-1;%barcodes.lineParams{molIdx}(2)+dy-1;
%                 x = movies.pos{curIdx}(2)+hOff+dx;
            else
                y =  output{idxMov}.xy{nameMol}{1}(output{idxMov}.dots{nameMol}.locations(j)+output{idxMov}.nanid(nameMol));
                x =  output{idxMov}.xy{nameMol}{2}(output{idxMov}.dots{nameMol}.locations(j)+output{idxMov}.nanid(nameMol));
            end
            plot(x,y,'rx'); % maybe too much to plot also this info
            str = sprintf('I = %.1f', output{idxMov}.dots{nameMol}.val(j)); %,barcodes.dots{molIdx}.depth(j)
            text(x-5,y-5,str,'Color','white');
     end

     
end

