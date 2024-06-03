function [dotLocsGlobal] = sdd_mark_bars(movies,barcodes,tiles,extractionMethod)
    % sdd_mark_bars
    %
    % marks result visually of detected molecules.
    %
    %   Args:
    %       movies,barcodes,tiles
    %
    %
    
    if nargin < 4
        extractionMethod = 1;
    end
 
    % center to the bottom of the molecule.
    voffList = zeros(1,length(barcodes.expBars));
    dotLocsGlobal = cell(1,length(barcodes.expBars));

if isfield(movies,'dotFigNum')
    
   
    axes(tiles.dotDet);
    hold on
    for molIdx = 1:length(barcodes.dots)
        curIdx = barcodes.idx(molIdx);
        str = sprintf('Mol. %i',curIdx);
        angle = atan(barcodes.lineParams{curIdx}(1));
        % 
        if extractionMethod == 1
            vOff = barcodes.boundaries{curIdx}(1);
            hOff = barcodes.boundaries{curIdx}(3);
        else
            vOff = barcodes.boundaries{molIdx}(1);
            hOff = barcodes.boundaries{molIdx}(3);
        end
        voffList(molIdx) = vOff;
        text(movies.pos{curIdx}(2),movies.pos{curIdx}(1)+vOff,str,'Color','white','Clipping','on');%
        for j = 1:numel(barcodes.dots{molIdx}.locations)
            try
                if extractionMethod == 1
                    dy = -sin(angle)*(barcodes.dots{molIdx}.locations(j)-1+barcodes.dots{molIdx}.leftOffset);
                    dx = cos(angle)*(barcodes.dots{molIdx}.locations(j)-1+barcodes.dots{molIdx}.leftOffset);
                    y = movies.pos{curIdx}(1)+vOff+dy-1;%barcodes.lineParams{molIdx}(2)+dy-1;
                    x = movies.pos{curIdx}(2)+hOff+dx;
                else
                    y =  movies.pos{curIdx}(1)+barcodes.xy{molIdx}{1}(barcodes.dots{molIdx}.locations(j)+barcodes.dots{molIdx}.leftOffset)-1;
                    x =  movies.pos{curIdx}(2)+barcodes.xy{molIdx}{2}(barcodes.dots{molIdx}.locations(j)+barcodes.dots{molIdx}.leftOffset)-1;
                end
                dotLocsGlobal{molIdx}{j} = [x,y];
    %             y = movies.pos{curIdx}(1)+vOff+dy-1;%barcodes.lineParams{molIdx}(2)+dy-1;
    %             x = movies.pos{curIdx}(2)+hOff+dx;
    
                plot(x,y,'rx'); % maybe too much to plot also this info
                str = sprintf('I = %.1f', barcodes.dots{molIdx}.val(j)); %,barcodes.dots{molIdx}.depth(j)
                text(x-5,y-5,str,'Color','white','Clipping','on');
            catch
                warning(['Mol_', num2str(molIdx), ' dot_', num2str(j), ' not plotted correctly']);
            end
        end
    end
  
    hold off
end
      

  
      %                 % Plot the estimated beginning of each molecule (the
      %                 % leftmost end)
      %                 x0 = movies.pos{molIdx}(2)+cos(angle)*(-1+barcodes.dots{molIdx}.leftOffset);
      %                 y0 = movies.pos{molIdx}(1)+barcodes.lineParams{molIdx}(2)-sin(angle)*(-1+barcodes.dots{molIdx}.leftOffset);
      %                 plot(x0,y0,'gd');
      %                 % Plot the estimated beginning of end molecule (the
      %                 % rightmost end)
      %                 kymLength = numel(barcodes.dotBars{molIdx});
      %                 xend2 = movies.pos{molIdx}(2) + cos(angle)*(kymLength-barcodes.dots{molIdx}.rightOffset-1);
      %                 yend2 = movies.pos{molIdx}(1)+barcodes.lineParams{molIdx}(2) - sin(angle)*(kymLength-barcodes.dots{molIdx}.rightOffset-1);
      %                 plot(xend2,yend2,'bd');
  

axes(tiles.molDet);
hold on
 for molIdx = 1:length(barcodes.expBars)
    curIdx = barcodes.idx(molIdx);
    str = sprintf('Mol. %i',curIdx);
    text(movies.pos{curIdx}(2),movies.pos{curIdx}(1)+voffList(molIdx),str,'Color','white','Clipping','on');
  end

hold off

%%
if isfield(movies,'dotFigNum2')
    
   
    axes(tiles.dotDet2);
    hold on
    for molIdx = 1:length(barcodes.dots2)
        curIdx = barcodes.idx(molIdx);
        str = sprintf('Mol. %i',curIdx);
        angle = atan(barcodes.lineParams{curIdx}(1));
        % 
        if extractionMethod == 1
            vOff = barcodes.boundaries{curIdx}(1);
            hOff = barcodes.boundaries{curIdx}(3);
        else
            vOff = barcodes.boundaries{molIdx}(1);
            hOff = barcodes.boundaries{molIdx}(3);
        end
        voffList(molIdx) = vOff;
        text(movies.pos{curIdx}(2),movies.pos{curIdx}(1)+vOff,str,'Color','white','Clipping','on');%
        for j = 1:numel(barcodes.dots{molIdx}.locations)
            try
                if extractionMethod == 1
                    dy = -sin(angle)*(barcodes.dots2{molIdx}.locations(j)-1+barcodes.dots2{molIdx}.leftOffset);
                    dx = cos(angle)*(barcodes.dots2{molIdx}.locations(j)-1+barcodes.dots2{molIdx}.leftOffset);
                    y = movies.pos{curIdx}(1)+vOff+dy-1;%barcodes.lineParams{molIdx}(2)+dy-1;
                    x = movies.pos{curIdx}(2)+hOff+dx;
                else
                    y =  movies.pos{curIdx}(1)+barcodes.xy{molIdx}{1}(barcodes.dots2{molIdx}.locations(j)+barcodes.dots2{molIdx}.leftOffset)-1;
                    x =  movies.pos{curIdx}(2)+barcodes.xy{molIdx}{2}(barcodes.dots2{molIdx}.locations(j)+barcodes.dots2{molIdx}.leftOffset)-1;
                end
                dotLocsGlobal{molIdx}{j} = [x,y];
    %             y = movies.pos{curIdx}(1)+vOff+dy-1;%barcodes.lineParams{molIdx}(2)+dy-1;
    %             x = movies.pos{curIdx}(2)+hOff+dx;
    
                plot(x,y,'rx'); % maybe too much to plot also this info
                str = sprintf('I = %.1f', barcodes.dots2{molIdx}.val(j)); %,barcodes.dots{molIdx}.depth(j)
                text(x-5,y-5,str,'Color','white','Clipping','on');
            catch
                warning(['Mol_', num2str(molIdx), ' dot_', num2str(j), ' not plotted correctly']);
            end
        end
    end
  
    hold off
end
      
%     axes(tiles.dotDet2);
% hold on
%  for molIdx = 1:length(barcodes.expBars)
%     curIdx = barcodes.idx(molIdx);
%     str = sprintf('Mol. %i',curIdx);
%     text(movies.pos{curIdx}(2),movies.pos{curIdx}(1)+voffList(molIdx),str,'Color','white','Clipping','on');
%   end
% 
% hold off