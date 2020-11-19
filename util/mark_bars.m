function mark_bars(movies,barcodes)

if isfield(movies,'dotFigNum')
  figure(movies.dotFigNum)
  molIdx = 0;
  hold on
  while molIdx < length(barcodes.dots) +1
    molIdx = molIdx + 1;
    while sum(barcodes.delid == molIdx) == 1
      molIdx = molIdx + 1;
    end
    if molIdx < length(barcodes.dots)+1
      str = sprintf('Mol. %i',molIdx);
      angle = atan(barcodes.lineParams{molIdx}(1));
      if barcodes.lineParams{molIdx}(1) < 0
        vOff = 0;
        hOff = (-barcodes.lineParams{molIdx}(1)) * cos(angle)-sin(angle)*barcodes.dots{molIdx}.leftOffset;
      elseif barcodes.lineParams{molIdx}(1) > 1
        vOff = size(movies.bwM{molIdx},1);
        hOff = cos(angle)*(barcodes.lineParams{molIdx}(2) - vOff);
      else
        vOff = barcodes.lineParams{molIdx}(2);
        hOff = 0;
      end
      text(movies.pos{molIdx}(2),movies.pos{molIdx}(1)+vOff,str,'Color','white');
      for j = 1:numel(barcodes.dots{molIdx}.locations)
        dy = -sin(angle)*(barcodes.dots{molIdx}.locations(j)-1+barcodes.dots{molIdx}.leftOffset);
        dx = cos(angle)*(barcodes.dots{molIdx}.locations(j)-1+barcodes.dots{molIdx}.leftOffset);
        y = movies.pos{molIdx}(1)+vOff+dy-1;%barcodes.lineParams{molIdx}(2)+dy-1;
        x = movies.pos{molIdx}(2)+hOff+dx;
        plot(x,y,'rx');
        str = sprintf('I = %.1f, depth = %.1f',barcodes.dots{molIdx}.val(j),barcodes.dots{molIdx}.depth(j));
        text(x-5,y-5,str,'Color','white');
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
    end
  end
  hold off
end

figure(movies.molFigNum)
molIdx = 0;
hold on
while molIdx < length(barcodes.expBars)+1
  molIdx = molIdx + 1;
  while sum(barcodes.delid == molIdx) == 1
    molIdx = molIdx + 1;
  end
  if molIdx < length(barcodes.expBars)+1
    str = sprintf('Mol. %i',molIdx);
    text(movies.pos{molIdx}(2),movies.pos{molIdx}(1),str,'Color','white');
  end
end
hold off
